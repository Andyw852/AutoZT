#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""validate_dielectric.py —— step5_dielect 的**数值**校验（分析步，不提交 HPC）

为什么必须有这一步
--------------------------------------------------------------------------
step5_dielect 的完成标记是
    check: marker, marker: "OUTCAR:MACROSCOPIC STATIC DIELECTRIC TENSOR"
—— 只要 OUTCAR 里有这个**字符串**就算完成，**完全不管张量内容**。

实测（2026-09-14，jzzn 全工作树扫描）：
    qTPC20-r / Mg2C60 / Mo2S3(P1) / Pb2Sb2Te5 / Sn2Sb2Te5  五个工程的张量**全是 NaN**
    （HEAD 块还是单位矩阵 = DFPT 迭代初值签名），却全部被判为"完成"，一路走到 step8。
    NaN 随后被兜底成看似合理的标量（如 0.8614），或使 eps_static ≡ eps_inf
    -> Frohlich 耦合恒为 0 -> **POP 散射静默丢失**（实测签名 max|mobility/POP|=2.3e48）。

另外 skill 原有的 patch_dielec_guard（gen_step10_amset.py）写的是
    if any(eps_0[i][i] < 0 ...): raise
而 Python 里 float('nan') < 0 是 **False** —— 那个守卫对 NaN 恒不触发。

判据（任一项不过 -> 写 ok=false 并退出码 40，让 tf 判 error）
--------------------------------------------------------------------------
  1) 有 'MACROSCOPIC ... (including local field effects in DFT)' 块
  2) 无 NaN / Inf
  3) 对角 >= 1                      （静态电子响应不可能 < 1）
  4) 不是单位矩阵                    （DFPT 迭代初值签名 —— 零响应）
  5) 有 'IONIC CONTRIBUTION' 块       （缺它则 eps_static 只能等于 eps_inf）
  6) 非极性体系除外：eps_static != eps_inf （否则 Frohlich 耦合恒 0、POP 丢失）

用法：在工程根目录跑 python3 validate_dielectric.py [--step-dir step5_dielect]
"""
import argparse, json, math, os, re, sys
from pathlib import Path

TAG_INF = "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"
TAG_ION = "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"


def grab(txt, tag):
    i = txt.rfind(tag)
    if i < 0:
        return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        n = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|NaN|nan|[+-]?Inf(?:inity)?", ln)
        if len(n) >= 3:
            out = []
            for v in n[:3]:
                lv = v.lower()
                out.append(float("nan") if lv == "nan" else
                           (float("inf") if "inf" in lv else float(v)))
            rows.append(out)
        if len(rows) == 3:
            break
    return rows if len(rows) == 3 else None


def any_nan(M):
    return bool(M) and any(math.isnan(x) or math.isinf(x) for r in M for x in r)


def soc_isym_guard(step_dir):
    """SOC + DFPT 的参数组合哨兵（2026-09-15 新增）。

    VASP 6.4.3 / 6.5.0 在 LSORBIT=.TRUE. 时做 DFPT 有 spinor 旋转缺失 bug，
    官方绕法是 **ISYM=-1**。INCAR 的 incar_dfpt_*.tpl 默认写 ISYM=2，若某体系
    从 step1 继承了 LSORBIT=.TRUE. 而 ISYM 仍是 2，就会踩这个 bug ——
    表现为 DFPT 崩溃，或张量非对角项异常增大（实测 2.3e-3 vs 正常 1e-6）。

    gen_step8_dielect.py 的 _fix_soc_isym() 会在生成阶段自动改成 -1；本函数是
    第二道防线：即使输入被手工改过，也在校验阶段把危险组合报出来。
    """
    inc = step_dir / "INCAR"
    if not inc.is_file():
        return None, []
    try:
        t = inc.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None, []
    soc = bool(re.search(r"^\s*LSORBIT\s*=\s*\.TRUE\.", t, re.M | re.I))
    mi = re.search(r"^\s*ISYM\s*=\s*(-?\d+)", t, re.M | re.I)
    isym = int(mi.group(1)) if mi else None
    info = {"lsorbit": soc, "isym": isym}
    w = []
    if soc and isym is not None and isym != -1:
        w.append("LSORBIT=.TRUE. 但 ISYM=%d —— DFPT+SOC 的 spinor 旋转 bug 组合，"
                 "结果可能错误或崩溃；应改为 ISYM=-1" % isym)
    return info, w


def born_asr_and_crosscheck(txt, ei):
    """DFPT 产物的物理自洽检查（2026-09-15 新增，来源：A2B2Te5 四材料复核）。

    这三项都是**零额外算力**的强校验，用来回答"DFPT 到底算对没有"：

      1) 声学求和规则 Sigma Z* = 0。Born 有效电荷是 DFPT 的直接产物，
         求和为 0 要求整条线性响应链都算对，极难靠巧合满足。
         实测四材料（PBE+SOC, ISYM=-1, 9x9x1）全部精确到 1e-5。
         构造 Z* 是位移导数、求和规则被破坏到 1e-2 量级 -> DFPT 可疑。
      2) DFPT 与「独立粒子」两条独立算法交叉。VASP 在同一次运行里还会用
         能带求和的独立粒子法算一个 eps（TAG_IND），算法与 DFPT 完全不同。
         实测面内一致到 1.0~1.4%，面外差 8.6~10.3%（那是真实的局域场效应，
         层状材料 zz 方向局域场修正更大，不是错误）。
         若两者差 > 25%，多半是 NBANDS 太少或 k 网格太稀。
      3) 对角性。高对称晶胞的 eps 张量应当对角，非对角项应 ~1e-5。
         实测 ISYM=2 + SOC 的坏配置给出 0.1~0.4 量级的非对角项 ——
         这是 spinor 旋转 bug 的指纹，正是它逼我们把 SOC 锁到 ISYM=-1。

    返回值：{"physics": {...}}，并把可疑项追加进 reasons（由调用方合并）。
    """
    out = {"nkpts": None, "born_sum_max": None, "born_asr_ok": None,
           "cross_ratio": None, "offdiag_ratio": None}
    warns = []
    m = re.search(r"NKPTS\s*=\s*(\d+)", txt)
    if m:
        out["nkpts"] = int(m.group(1))
    # ---- 1) Born 有效电荷求和规则 ----
    try:
        i = txt.rfind("BORN EFFECTIVE CHARGES (including local field effects)")
        if i >= 0:
            Z = []
            for mo in re.finditer(
                    r"ion\s+(\d+)\s*\n((?:\s+\d+\s+[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n){3})",
                    txt[i:]):
                Z.append([[float(v) for v in ln.split()[1:4]]
                          for ln in mo.group(2).strip().split("\n")])
            if Z:
                n = len(Z)
                s = [sum(Z[k][r][c] for k in range(n)) for r in range(3) for c in range(3)]
                mx = max(abs(v) for v in s)
                out["born_sum_max"] = round(mx, 6)
                out["born_ions"] = n
                out["born_asr_ok"] = mx < 1e-2
                if mx >= 1e-1:
                    warns.append("Born 有效电荷求和规则 Sigma Z* = %.4f 被破坏（>0.1）"
                                 "—— DFPT 线性响应不可信" % mx)
                elif mx >= 1e-2:
                    warns.append("Born 有效电荷 Sigma Z* = %.4f 偏大（正常应 ~1e-5）" % mx)
    except Exception as e:
        warns.append("Born 求和规则解析失败：%s" % e)
    # ---- 2) DFPT vs 独立粒子 ----
    try:
        ind = None
        j = txt.rfind("HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE")
        if j >= 0:
            rows = []
            for ln in txt[j:].splitlines()[1:]:
                n2 = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", ln)
                if len(n2) >= 3:
                    rows.append([float(v) for v in n2[:3]])
                if len(rows) == 3:
                    break
            if len(rows) == 3:
                ind = rows
        if ind and ei and not any_nan(ei):
            d_d = [ei[k][k] for k in range(3)]
            d_i = [ind[k][k] for k in range(3)]
            rat = max(abs(d_d[k] - d_i[k]) / max(abs(d_d[k]), 1e-9) for k in range(3))
            out["cross_ratio"] = round(rat, 4)
            out["eps_inf_independent"] = [round(v, 4) for v in d_i]
            if rat > 0.25:
                warns.append("DFPT 与独立粒子 eps 差 %.1f%%（>25%%）—— 检查 NBANDS 是否够多、k 网格是否太稀" % (rat * 100))
        else:
            out["cross_ratio"] = "no-independent-particle-block"
    except Exception as e:
        warns.append("独立粒子交叉校验失败：%s" % e)
    # ---- 3) 对角性 ----
    try:
        if ei and not any_nan(ei):
            dg = max(abs(ei[k][k]) for k in range(3))
            od = max(abs(ei[r][c]) for r in range(3) for c in range(3) if r != c)
            out["offdiag_ratio"] = round(od / max(dg, 1e-9), 6)
            if od / max(dg, 1e-9) > 1e-4:
                warns.append("eps 非对角项相对大小 %.3g（高对称晶胞应 ~1e-6）"
                             "—— 典型病因是 SOC 下用了 ISYM != -1（spinor 旋转 bug）" % (od / max(dg, 1e-9)))
    except Exception:
        pass
    return out, warns

def _is_nonpolar(step_dir):
    """单元素体系 = 非极性：没有 IR 活性声子，eps_static ≡ eps_inf 是【物理正确】结果。

    与 step8_amset 的 _is_nonpolar() 同判据（读 POSCAR 第 6 行元素名，单元素即非极性）。
    这类体系里 Frohlich 耦合恒 0 不是 DFPT 失效的签名，不该判 FAIL。
    """
    for cand in (Path(step_dir) / "POSCAR", Path(step_dir).parent / "POSCAR"):
        try:
            syms = cand.read_text(errors="ignore").splitlines()[5].split()
            if syms and not any(ch.isdigit() for ch in syms[0]):
                return len(set(syms)) == 1
        except (OSError, IndexError):
            continue
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step-dir", default="step5_dielect")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    d = Path(a.step_dir)
    out = Path(a.json) if a.json else d / "dielectric_check.json"
    # run:gen 步骤的"步骤目录"没有任何人创建（gen 在技能根目录里跑）—— 显式补一下，
    # 否则 --json 指到步骤目录时会 FileNotFoundError，autozt 也就永远找不到 done_marker。
    out.parent.mkdir(parents=True, exist_ok=True)
    oc = d / "OUTCAR"
    res = {"step": "step5_dielect.validate", "outcar": str(oc), "ok": False, "reasons": []}
    if not oc.is_file():
        res["reasons"].append("缺 OUTCAR")
        out.write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n")
        print("[ERROR] %s" % res["reasons"][-1]); return 40
    txt = oc.read_text(errors="ignore")
    ei, io = grab(txt, TAG_INF), grab(txt, TAG_ION)
    res["eps_inf"] = ei
    res["eps_ionic"] = io
    if ei is None:
        res["reasons"].append("没有 '%s' 块（DFPT 没跑完 / 崩了 / 停在初值）" % TAG_INF[:52])
    else:
        if any_nan(ei):
            res["reasons"].append("eps_inf 含 NaN/Inf —— DFPT 发散或未收敛（HEAD 块为单位矩阵时即为迭代初值）")
        else:
            diag = [ei[i][i] for i in range(3)]
            if min(diag) < 1.0:
                res["reasons"].append("eps_inf 对角 %s 有 <1 的分量 —— 静态电子响应不可能小于 1" % diag)
            if max(abs(ei[i][j] - (1.0 if i == j else 0.0)) for i in range(3) for j in range(3)) < 1e-6:
                res["reasons"].append("eps_inf 恰为单位矩阵 —— DFPT 迭代初值签名（零响应）")
    if io is None:
        res["reasons"].append("没有 '%s' 块 —— 离子介电未算出，eps_static 只能退化成 eps_inf，POP 会静默丢失" % TAG_ION[:52])
    else:
        if any_nan(io):
            res["reasons"].append("IONIC 块含 NaN/Inf —— 离子介电段发散（注意：原 patch_dielec_guard 的 '<0' 判断对 NaN 恒为假）")
        else:
            if any(io[i][i] < 0 for i in range(3)):
                res["reasons"].append("离子介电对角为负 —— 非物理（常见于 Gamma 近零声学模污染）")
            elif ei is not None and not any_nan(ei):
                coup = max(abs(1.0 / ei[i][i] - 1.0 / (ei[i][i] + io[i][i])) for i in range(3))
                res["frohlich_coupling_max"] = coup
                if coup < 1e-9:
                    if _is_nonpolar(d):
                        # 非极性：POP 本就不适用，不是失败。记进 notes（warnings 稍后会被
                        # 物理自洽检查整体覆盖，所以单独放 notes）。
                        res.setdefault("notes", []).append(
                            "非极性体系（单元素）：eps_static ≡ eps_inf 是物理正确结果，"
                            "Frohlich 耦合恒 0 -> POP 不适用（不判 FAIL；"
                            "与 step8_amset 的非极性处理一致）")
                    else:
                        res["reasons"].append("eps_static == eps_inf，Frohlich 耦合恒 0 -> POP 静默丢失"
                                              "（非极性体系此项可忽略，需人工确认）")
    # ---- 物理自洽检查（求和规则 / 双算法交叉 / 对角性）----
    # 只有"求和规则被破坏到 1e-1"才判致命；其余进 warnings，不阻断流程 ——
    # 低对称(P1)晶胞的非对角项、宽隙体系的交叉偏差都可能是正常的。
    _phys, _warns = born_asr_and_crosscheck(txt, ei)
    _soc, _sw = soc_isym_guard(d)
    res["physics"] = _phys
    if _soc:
        res["physics"]["soc_isym"] = _soc
    res["warnings"] = _warns + _sw
    if _phys.get("born_sum_max") is not None and _phys["born_sum_max"] >= 0.1:
        res["reasons"].append(
            "Born 有效电荷求和规则 Sigma Z* = %.4f 被破坏（>0.1）—— DFPT 线性响应不可信"
            % _phys["born_sum_max"])
    res["ok"] = not res["reasons"]
    out.write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    if res["ok"]:
        print("[OK] 介电张量通过数值校验：eps_inf=%s ionic=%s Frohlich=%.6f"
              % ([round(ei[i][i], 3) for i in range(3)],
                 [round(io[i][i], 3) for i in range(3)], res.get("frohlich_coupling_max", 0)))
        _p = res.get("physics", {})
        print("     物理自洽：NKPTS=%s  Sigma Z* max=%s  双算法偏差=%s  非对角比=%s"
              % (_p.get("nkpts"), _p.get("born_sum_max"),
                 _p.get("cross_ratio"), _p.get("offdiag_ratio")))
        for w in res.get("warnings", []):
            print("     [WARN] %s" % w)
        return 0
    for w in res.get("warnings", []):
        print("[WARN] %s" % w)
    print("[ERROR] step5_dielect 的介电产物不可用：")
    for r in res["reasons"]:
        print("   - %s" % r)
    print("   产物: %s" % out)
    return 40


if __name__ == "__main__":
    sys.exit(main())
