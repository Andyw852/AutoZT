#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lineage_check.py —— 上游同源核对（V140）：各步结构是否都来自 S1 的当前 CONTCAR、派生产物是否旧于来源。

背景：S1 重新弛豫后整条 rerun 时，
  · autozt rerun 先 rm -rf 步骤目录再 gen，V122/V136 的"gen 见到旧产物 -> 下游失效"不触发
    （V140 起各 gen 改为无条件失效；在此之前已经 rerun 过的，用 --invalidate-from 补）；
  · S2.3 HSE 的 gen 从 S2.2 拷 POSCAR/WAVECAR，S2.2 没先重跑就 rerun S2.3，HSE 算的是旧结构；
  · S4 的 h5、S7.1 的形变势、S2 画图的 band_summary.json 仍是旧的"完成"。
S8/S8.4 的 gen 做同一核对（STRUCTURE_GUARD），不通过就不生成。本工具只读（除非 --invalidate-from）。
V148 起本工具另查（S8 的 gen 不查）：S8.1/S8.3 的产物是否旧于 S8.2/S8/S8.4，S8/S8.4 的带隙是否与 S2 画图一致。
V149 起结构比对扣除整体平移（S3/S3b 的 align_origin 是刚性平移，不算不同）。
V153：--invalidate-from 也接受 S1 / S2 链（下游 VASP 步骤整目录归档；在跑的不动）。
V175：另查 k 网格 —— S3 不合现行规则、S7 和 S3 不是同一张网格且不合规。gen 的闸门只在重新生成时起作用，
      已经算完的旧结果只能靠这里查出来（CrS₂ 的 15×15×1 当初就是这样"收口"的）。
V185：材料目录是 autozt fetch 拉回的本地副本（步骤目录里有 .tf_fetched）时逐个标出，末尾提示到集群上复核。
      fetch 只刷新 fetch_files 里的文件（扇出步骤如 S7 也只拉子目录里的同名文件），KPOINTS / INCAR / POSCAR
      往往是更早一次拉回的旧版；9p / drvfs 上修改时间也是拉回时刻。实测：本地 S7 的 ionrelax/KPOINTS 还是
      15×15×1，集群上早已是 48×48×3，本地跑出了一条假的"S7 网格不合规"。
V183：递归找材料时跳过 project_setting/、templates/、隐藏目录和 *.stale-* 归档目录（模板不是材料）。
V182：末尾按泛函分组列出本次查到的材料（step1 的 FUNC=，读不到从 INCAR 反推）；不止一种时提示。
      只是提示，不算上游问题、不改退出码 —— 出厂是 PBEsol，材料级 step.conf 可以改成 PBE，
      同一项目里混着两种很正常，但放进同一张表比较前要知道。

用法：
    python lineage_check.py <材料目录或项目根> [--glob "*/ke-dft-cpu"]
        列出每个材料的问题；退出码 1 = 有问题。
    python lineage_check.py <材料目录> --invalidate-from step3_uniform step7_deform ...
        已经 rerun 过这些上游：把它们下游的完成标记改名归档（*.stale-upstream-*，不删），
        之后 auto-advance 会在上游跑完后按顺序重做 S4 / S7.1 / S8 …。
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(HERE.parent.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


# [V183] 不是材料的目录：autozt init 把技能整套模板（含 step1_opt/ 等步骤目录）复制进
#   project_setting/templates/，递归找 step1_opt 时会把它当成材料（实测项目根多出 6 个"泛函未知"的"材料"）。
#   归档目录（*.stale-*）、隐藏目录同理跳过。
_NOT_MATERIAL_PARTS = ("project_setting", "templates")


def _is_material_path(p, root):
    try:
        parts = Path(p).relative_to(root).parts
    except ValueError:
        parts = Path(p).parts
    return not any(x in _NOT_MATERIAL_PARTS or x.startswith(".") or ".stale-" in x for x in parts)


FETCH_STAMP = ".tf_fetched"        # autozt fetch 在 result_dir/<步骤>/ 下写的回执（autozt.bootstrap.FETCH_STAMP）


def is_local_copy(mat):
    """[V185] 材料目录是不是 autozt fetch 拉回的本地副本。"""
    return any(Path(mat).glob("*/" + FETCH_STAMP))


def find_materials(root, pattern=None):
    root = Path(root)
    if any((root / d).is_dir() for d in kc.S1_DIRS):
        return [root]
    cands = root.glob(pattern) if pattern else (p.parent for d in kc.S1_DIRS for p in root.rglob(d))
    return sorted({Path(p) for p in cands
                   if any((Path(p) / d).is_dir() for d in kc.S1_DIRS) and _is_material_path(p, root)})


def _short(mat):
    """材料的短名：<材料>/ke-dft-cpu 这种布局取上一级目录名。"""
    mat = Path(mat)
    return mat.parent.name if mat.name.endswith("-cpu") and mat.parent.name else mat.name


def functional_report(mats):
    """[V182] 按泛函分组的几行文字；mats 为空返回 []。"""
    groups = {}
    for m in mats:
        groups.setdefault(kc.material_functional(m)["label"], []).append(_short(m))
    if not groups:
        return []
    lines = ["\n泛函（step1 的 FUNC=；GGA=PS 是 PBEsol，GGA=PE 是 PBE）："]
    for lab in sorted(groups):
        lines.append("  %s × %d：%s" % (lab, len(groups[lab]), "、".join(sorted(groups[lab]))))
    if len(groups) > 1:
        lines.append("  [提示] 这些材料不是同一个泛函：放进同一张表或互相比较前，先统一泛函，或在表里注明各自的泛函。")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("path")
    ap.add_argument("--glob", default=None, help="材料目录的 glob（相对 path），默认递归找含 step1_opt 的目录")
    ap.add_argument("--invalidate-from", nargs="+", default=None, metavar="STEP",
                    help="只对单个材料目录：这些上游已经重算 -> 归档它们下游的完成标记")
    a = ap.parse_args(argv)

    if a.invalidate_from:
        mat = Path(a.path)
        if not any((mat / d).is_dir() for d in kc.S1_DIRS):
            sys.exit("[ERROR] --invalidate-from 只接受单个材料目录（里面要有 step1_opt）：%s" % mat)
        unknown = [s for s in a.invalidate_from if s not in kc.DOWNSTREAM]
        if unknown:
            sys.exit("[ERROR] 不认识的步骤 %s；可用：%s" % (unknown, ", ".join(sorted(kc.DOWNSTREAM))))
        done = []
        for step in a.invalidate_from:
            # V148：S2 画图的下游按带隙变没变决定（gap 取它当前的 band_summary.json）
            done += kc.invalidate_downstream(mat, step, "lineage_check --invalidate-from",
                                             gap=kc.band_summary_gap(mat / step))
        print("[OK] 归档了 %d 个完成标记%s" % (len(done), "" if done else "（下游本来就没有完成标记）"))
        for step, m in done:
            print("     %s/%s" % (step, m))

    bad = 0
    mats = find_materials(a.path, a.glob)
    if not mats:
        print("没找到材料目录（含 step1_opt / step1_std_opt）：%s" % a.path)
        return 0
    local = [m for m in mats if is_local_copy(m)]               # V185
    for mat in mats:
        tag = "（本地拉回的副本）" if mat in local else ""
        ref, probs = kc.structure_lineage(mat)
        probs = probs + kc.post_lineage(mat)                 # V148：S8 之后的派生产物 + 带隙一致性
        probs = probs + kc.grid_lineage(mat)                 # V175：S3 / S7 的 k 网格
        if not probs:
            print("[OK]  %s%s%s" % (mat, tag, "" if ref else "（没有 S1 CONTCAR，只核对了派生产物）"))
            continue
        bad += 1
        print("[★]   %s%s：%d 处" % (mat, tag, len(probs)))
        for step, why in probs:
            print("        - %s：%s" % (step, why))
    if bad:
        print("\n%d/%d 个材料的上游不同源。结构不同的 VASP 步骤 rerun（S2.3 要等 S2.2、S2.2 要等 S2.1）；"
              "派生步骤在来源跑完后 rerun，或用 --invalidate-from 让 auto-advance 接手。" % (bad, len(mats)))
        print("k 网格问题：S3 不合规 -> retry S3 -> S4 -> S7 -> S7.1 -> S8/S8.4/S8.2/S8.1；"
              "只有 S7 不合规 -> retry S7 -> S7.1 -> S8/S8.4。")
    for ln in functional_report(mats):                       # V182：只提示，不改退出码
        print(ln)
    if local:                                                # V185
        print("\n[提示] %d 个材料目录是 autozt fetch 拉回的本地副本：fetch 只刷新 fetch_files 里的文件，"
              "KPOINTS / INCAR / POSCAR 等可能还是更早拉回的旧版，修改时间在 9p / drvfs 上也不可信。"
              "上面的结论（尤其网格、新旧）请在集群上的材料目录里再跑一次确认；"
              "要刷新本地副本，用 autozt fetch --all。" % len(local))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
