#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stepconf.py —— 步骤配置文件 step.conf 的解析 / 合并 / 类型化取值
================================================================
一个步骤一个文件夹，文件夹里一个 step.conf，文件里分节。全文只有一种语法：
INCAR 风格的 `KEY = VALUE`，行尾 # 或 ! 之后是注释。没有 JSON，没有特殊符号。

    [params]        脚本行为参数（不写进 INCAR）
    [submit]        覆盖 submit.sh 的 Slurm 行；留空 = 沿用模板原值
    [incar]         覆盖继承来的 INCAR 标签
    [incar.delete]  删除继承来的标签（每行一个标签名，不写等号）
    [incar.final]   在脚本自动计算(NBANDS/KPAR/MAGMOM…)之后再覆盖，你说了算
    [<自定义>]      结构化数据，每行按空白切分，如 [kpath.extra] 的 "X 0.5 0.0 0.0"

合并：tf 在本地按 skill 默认 → project templates/ → templates/<步骤>/ 逐级叠加，
把结果连同来源注释写成一份 step.conf 推到超算。gen 脚本只读这一份，不做回落。
"""

import re
import sys
from pathlib import Path

CONF_NAME = "step.conf"
_COMMENT = re.compile(r"\s+[#!].*$")
_SECTION = re.compile(r"^\[([A-Za-z0-9_.\-]+)\]$")

# 驱动层保留键：写在 [params] 里、供 tf 决定步骤图（如 BANDGAP=pbe|hse 增删
# 整段 HSE），gen 脚本本身不消费。校验白名单时无条件放行，避免各 gen 脚本
# 都误报"不认识的键"。新增工作流级开关往这里加即可。
RESERVED_PARAMS = frozenset({"BANDGAP", "CONDA_SH", "CONDA_ENV", "MACE_MODEL_DIR", "AMSET_ENV", "POTCAR_DIR", "REFERENCES_DIR"})


def _strip(line):
    s = line.rstrip()
    if s.lstrip().startswith(("#", "!")):
        return ""
    return _COMMENT.sub("", s).strip()


def parse(text, src="<text>"):
    """-> {节名: [(key, value, 原始行号), ...]}；[incar.delete] 这类无等号的行 value=None。"""
    out, sec = {}, "params"
    for i, line in enumerate(text.splitlines(), 1):
        s = _strip(line)
        if not s:
            continue
        m = _SECTION.match(s)
        if m:
            sec = m.group(1).lower()
            out.setdefault(sec, [])
            continue
        if "=" in s:
            k, v = s.split("=", 1)
            out.setdefault(sec, []).append((k.strip(), v.strip(), i))
        else:
            out.setdefault(sec, []).append((s, None, i))
    return out


def read_submit(path, cwd=".", used_incar=False):
    """只读 step.conf 的 [submit] 节，返回 {key(小写, 连字符转下划线): value}。

    供那些不加载完整 step.conf（没有 [params] spec、只关心提交参数）的
    gen 脚本复用——绕开 StepConf 对 [params] 未知键的严格校验，只取 [submit]。
    文件不存在时返回 {}。value 过滤掉 None / 空。

    ★ used_incar：本函数的调用方**是否也会消费** [incar] / [incar.final] /
      [incar.delete]（走 StepConf.apply_incar 或等价逻辑）。默认 False。
      为 False 且 step.conf 里确实写了这几节时，打印告警 —— 因为这几个节
      只被 apply_incar 读取，只调 read_submit 的 gen 脚本会把用户写的覆盖
      **静默忽略**（"写了以为生效、其实没生效"是这条链上反复出现的失效模式：
      2026-09-15 之前 ke-dft-cpu 全部 gen 脚本都踩了这个洞）。
      真会消费它们的脚本请显式传 used_incar=True，避免误报。
    """
    p = Path(path)
    if not p.is_absolute():
        p = Path(cwd) / p
    if not p.is_file():
        return {}
    merged = parse(p.read_text(encoding="utf-8-sig"), str(p))
    if not used_incar:
        _ic = [s for s in ("incar", "incar.final", "incar.delete") if merged.get(s)]
        if _ic:
            print("[WARN] %s 里有 [%s] 段，但本步骤的 gen 脚本只读了 [submit] —— "
                  "这些 INCAR 覆盖会被**静默忽略**。要让它们生效，gen 脚本必须调 "
                  "StepConf.apply_incar()（或等价逻辑）；确实不需要就请删掉这几节。"
                  % (p, "]/[".join(_ic)), file=sys.stderr)
    return {k.lower().replace("-", "_"): v
            for k, v, _ in merged.get("submit", []) if v not in (None, "")}


_SUBMIT_FLAGS = {
    "nodes": "nodes", "ntasks_per_node": "ntasks-per-node",
    "ntasks": "ntasks", "cpus_per_task": "cpus-per-task",
    "qos": "qos", "partition": "partition", "time": "time",
    "job_name": "job-name", "gres": "gres", "mem": "mem",
}


def apply_submit(submit_path, sub_dict):
    """把 sub_dict 覆盖到 submit.sh 的 #SBATCH 行；无该行则在首个 #SBATCH 后补一行。

    统一的核数/队列覆盖出口：gen 脚本写 submit.sh 后调
        stepconf.apply_submit(submit_path, stepconf.read_submit(CONF_NAME))
    sub_dict 的 key 用下划线小写（read_submit 已把连字符转下划线），
    值为 None/空 = 跳过。sub_dict 空则不动。
    """
    changed = []
    if not sub_dict:
        return changed
    p = Path(submit_path)
    if not p.is_file():
        return changed
    text = p.read_text(encoding="utf-8")
    for k, v in sub_dict.items():
        if v in (None, ""):
            continue
        fl = _SUBMIT_FLAGS.get(k, k.replace("_", "-"))
        pat = re.compile(r"^(#SBATCH\s+--%s=)\S+.*$" % re.escape(fl), re.MULTILINE)
        if pat.search(text):
            text = pat.sub(r"\g<1>%s" % v, text)
        else:
            lines = text.splitlines()
            last = max([i for i, ln in enumerate(lines)
                        if ln.startswith("#SBATCH")] or [0])
            if last:
                lines.insert(last + 1, "#SBATCH --%s=%s" % (fl, v))
                text = "\n".join(lines) + "\n"
        changed.append("--%s=%s" % (fl, v))
    p.write_text(text, encoding="utf-8", newline="\n")
    return changed


def merge(sources):
    """sources = [(标签, 文本), ...]，越靠后优先级越高。
    -> (merged, prov)：merged 同 parse 的结构；prov[(节, key)] = 标签。"""
    merged, prov = {}, {}
    for tag, text in sources:
        for sec, items in parse(text, tag).items():
            cur = merged.setdefault(sec, [])
            for k, v, _ in items:
                if v is None:                       # 无等号的行（如 [incar.delete]）：追加去重
                    if all(x[0] != k for x in cur):
                        cur.append((k, None, 0))
                        prov[(sec, k)] = tag
                    continue
                for j, (ek, _, _) in enumerate(cur):
                    if ek.upper() == k.upper():
                        cur[j] = (ek, v, 0)
                        break
                else:
                    cur.append((k, v, 0))
                prov[(sec, k)] = tag
    return merged, prov


def dumps(merged, prov=None, header_lines=()):
    """把合并结果写回 step.conf 文本，每行标注来源。"""
    out = list(header_lines)
    order = (["params", "submit", "incar", "incar.delete", "incar.final"]
             + sorted(s for s in merged
                      if s not in ("params", "submit", "incar",
                                   "incar.delete", "incar.final")))
    for sec in order:
        items = merged.get(sec)
        if not items:
            continue
        out.append("")
        out.append("[%s]" % sec)
        w = max((len(k) for k, _, _ in items), default=1)
        wv = max((len(v) for _, v, _ in items if v is not None), default=1)
        for k, v, _ in items:
            src = (prov or {}).get((sec, k))
            if v is None:
                out.append("%-*s%s" % (w + wv + 3, k,
                                       ("  # <- %s" % src) if src else ""))
            else:
                out.append("%-*s = %-*s%s" % (w, k, wv, v,
                                              ("  # <- %s" % src) if src else ""))
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ 类型化取值
def _as_bool(s):
    low = str(s).strip().lower()
    if low in ("true", ".true.", "yes", "on", "1"):
        return True
    if low in ("false", ".false.", "no", "off", "0"):
        return False
    raise ValueError("不是布尔值: %r" % s)


def _as_elemmap(s):
    """'Mn:5.0 In:0.0' -> {'Mn': 5.0, 'In': 0.0}"""
    out = {}
    for tok in str(s).replace(",", " ").split():
        if ":" not in tok:
            raise ValueError("元素表要写成 元素:数值，收到 %r" % tok)
        el, val = tok.split(":", 1)
        out[el.strip()] = float(val)
    return out


_CAST = {
    "str": lambda s: str(s),
    "int": lambda s: int(str(s), 0),
    "float": float,
    "bool": _as_bool,
    "elemmap": _as_elemmap,
    "words": lambda s: str(s).split(),
}


class StepConf(object):
    def __init__(self, merged, spec, path=None, strict=True):
        self._m, self._spec, self.path = merged, spec, path
        self.params = {}
        raw = {k.upper(): v for k, v, _ in merged.get("params", [])}
        unknown = sorted(set(raw) - {k.upper() for k in spec} - RESERVED_PARAMS)
        # strict=False：只取自己 SPEC 里声明的键，其余键（别的步骤的参数，
        # 如 FUNC）忽略。材料级 step.conf 是**全技能共用**的一份（"每层只写
        # 跟上一层不一样的键"），只认自己那一两个键的脚本用严格模式必然被
        # 别的步骤的键打死（2026-09-16 MoS2 S3_uniform 实测：FUNC 直接让 gen 退出）。
        if strict and unknown:
            raise SystemExit("[ERROR] %s 的 [params] 里有本脚本不认识的键：%s\n"
                             "        可用键：%s"
                             % (path or CONF_NAME, ", ".join(unknown),
                                ", ".join(sorted(spec))))
        for key, (default, typ) in spec.items():
            s = raw.get(key.upper())
            if s is None or s == "":
                self.params[key] = None if (s == "" and default is not None
                                            and typ != "str") else default
                if s == "":
                    self.params[key] = None
                continue
            try:
                self.params[key] = _CAST[typ](s)
            except (ValueError, KeyError) as e:
                raise SystemExit("[ERROR] %s 的 %s=%r 解析失败（应为 %s）：%s"
                                 % (path or CONF_NAME, key, s, typ, e))

    def __getitem__(self, k):
        return self.params[k]

    # ---- dict 只读协议（2026-09-22）--------------------------------------
    # StepConf 原来只实现了 __getitem__。下游那些"通用比对/留档"代码习惯性
    # 按 dict 用（conf.get("ALM_CUT2") 之类），于是 gen_step4_disp.py 在
    # S4 生成位移数据集时直接
    #     AttributeError: 'StepConf' object has no attribute 'get'
    # —— 崩在写完 POSCAR-*/phono3py_disp.yaml 之后、写 disp_plan.json 之前，
    # 留下一套没有留档的半成品数据集（幂等检查随后也认不出它）。
    # 只补一个调用点不解决问题（同类写法会继续踩），所以这里把 StepConf
    # 补成**只读 dict**：get/keys/values/items/copy/in/iter/len 全部可用，
    # 写入仍需先构造 dict（本类没有 __setitem__）。
    # 值语义不变：params 的键就是 SPEC 里的键（大写）。
    def get(self, k, default=None):
        return self.params.get(k, default)

    def keys(self):
        return self.params.keys()

    def values(self):
        return self.params.values()

    def items(self):
        return self.params.items()

    def copy(self):
        return dict(self.params)

    def __contains__(self, k):
        return k in self.params

    def __iter__(self):
        return iter(self.params)

    def __len__(self):
        return len(self.params)

    def section(self, name):
        """任意节的原始行：[(key, value|None), ...]"""
        return [(k, v) for k, v, _ in self._m.get(name.lower(), [])]

    @property
    def incar(self):
        return {k.upper(): v for k, v, _ in self._m.get("incar", []) if v is not None}

    @property
    def incar_final(self):
        return {k.upper(): v for k, v, _ in self._m.get("incar.final", []) if v is not None}

    @property
    def incar_delete(self):
        return {k.upper() for k, _, _ in self._m.get("incar.delete", [])}

    @property
    def submit(self):
        return {k.lower(): v for k, v, _ in self._m.get("submit", [])
                if v not in (None, "")}

    def apply_incar(self, inherited, computed=None):
        """继承 → [incar] → 脚本计算 → [incar.final] → [incar.delete]"""
        out = dict(inherited)
        out.update(self.incar)
        out.update(computed or {})
        out.update(self.incar_final)
        for k in self.incar_delete:
            out.pop(k, None)
        return out


def load(spec, step_name=None, cwd=".", strict=True):
    """gen 脚本入口：读材料目录里 tf 推来的那一份 step.conf。

    strict=False 时忽略 SPEC 之外的键（共用 step.conf 里别的步骤的参数），
    只认自己声明的那几个 —— 给"只用一两个小开关"的 gen 脚本用。"""
    p = Path(cwd) / CONF_NAME
    if not p.is_file():
        raise SystemExit("[ERROR] 缺少 %s —— 该步骤的 gen_need 里漏了它？" % CONF_NAME)
    merged = parse(p.read_text(encoding="utf-8-sig"), str(p))
    got = {k.upper(): v for k, v, _ in merged.get("params", [])}.get("STEP")
    if step_name and got and got != step_name:
        raise SystemExit("[ERROR] %s 属于步骤 %r，本脚本是 %r —— gen_need 串了。"
                         % (p, got, step_name))
    spec = dict(spec)
    spec.setdefault("STEP", (step_name, "str"))
    return StepConf(merged, spec, str(p), strict=strict)


# =============================================================================
# step.conf 的 [incar] / [incar.final] / [incar.delete] —— 面向"已经写出 INCAR 文件"
# 的 gen 脚本的共用入口。
#
# 为什么要有这一段：ke-dft-cpu 的 gen 脚本分两派 ——
#   ① 自己用 build_incar(items, remove, incar_set) 拼字符串（step2.1_static、
#      step6_elastic、step2.3_hse）；
#   ② 用 ke_common.render_tpl + inherit_scf_tags + apply_parallel_tags 渲染文件
#      （step2.2_pbe、step3_uniform、step4_wave、step5_dielect、step7_deform、
#       step8_amset）。
# 两派此前都只调 read_submit()，于是 step.conf 里写的 [incar*] 三节**全是死的**：
# 用户以为改了 INCAR，实际被静默忽略。2026-09-15 在 step2.3_hse 上暴露并修复。
# 为避免"每个脚本各写一套"再走回老路，这里提供唯一实现，两派都调它。
#
# 语义与 StepConf.apply_incar 完全一致：
#     继承(=文件现有内容) → [incar] → [incar.final] → [incar.delete]
# =============================================================================

_INCAR_KV = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)(.*?)(\s*)$")


def incar_sections(path=CONF_NAME, cwd="."):
    """读 step.conf，返回 ([incar] 字典, [incar.final] 字典, [incar.delete] 集合)。

    键统一大写；大小写不合并不在此处理（INCAR 标签本身大小写不敏感，
    调用方按 .upper() 比对）。文件不存在返回三个空容器。
    """
    p = Path(path)
    if not p.is_absolute():
        p = Path(cwd) / p
    if not p.is_file():
        return {}, {}, set()
    merged = parse(p.read_text(encoding="utf-8-sig"), str(p))
    _up = lambda k: k.strip().upper()
    incar = {_up(k): v for k, v, _ in merged.get("incar", []) if v is not None}
    final = {_up(k): v for k, v, _ in merged.get("incar.final", []) if v is not None}
    dele = {_up(k) for k, _, _ in merged.get("incar.delete", [])}
    return incar, final, dele


def _split_incar_value(raw):
    """把 INCAR 的值与行尾注释分开：'1E-7  # 说明' -> ('1E-7', '  # 说明')。"""
    for i, ch in enumerate(raw):
        if ch in "#!":
            return raw[:i].rstrip(), raw[i:]
    return raw.rstrip(), ""


def apply_incar_file(incar_path, path=CONF_NAME, cwd=".", log=None):
    """把 step.conf 的 [incar] / [incar.final] / [incar.delete] 应用到已写出的 INCAR。

    保留原有行序、缩进、对齐与行尾注释；[incar.final] 里文件原本没有的标签
    追加到文件末尾（并标注来源）。返回改动清单
        [(节名, 键, 旧值-or-None, 新值-or-None), ...]
    供调用方写进自己的 log。step.conf 不存在或三节全空时**不碰文件**。
    """
    changed = []
    incar, final, dele = incar_sections(path, cwd)
    if not (incar or final or dele):
        return changed
    ip = Path(incar_path)
    if not ip.is_file():
        return changed
    lines = ip.read_text(encoding="utf-8").splitlines()

    # 文件里已有的键 -> 行号（后出现的覆盖先出现的，与 VASP 的"后者生效"一致）
    pos = {}
    for i, ln in enumerate(lines):
        m = _INCAR_KV.match(ln)
        if m:
            pos[m.group(2).upper()] = i

    def _put(section, key, val):
        if key in pos:
            i = pos[key]
            m = _INCAR_KV.match(lines[i])
            indent, _, eq, old_raw, trail = m.groups()
            old_val, comment = _split_incar_value(old_raw)
            keep = comment or trail.strip()
            lines[i] = "%s%s%s%s%s" % (indent, m.group(2), eq, val,
                                       ("  " + keep) if keep else "")
            if old_val != val:
                changed.append((section, key, old_val, val))
        else:
            lines.append("%-9s= %s" % (key, val))
            pos[key] = len(lines) - 1
            changed.append((section, key, None, val))

    for k, v in incar.items():
        _put("incar", k, v)
    for k, v in final.items():
        _put("incar.final", k, v)
    for k in sorted(dele):
        if k in pos:
            old_val, _ = _split_incar_value(_INCAR_KV.match(lines[pos[k]]).group(4))
            del lines[pos[k]]
            pos = {}
            for i, ln in enumerate(lines):
                m = _INCAR_KV.match(ln)
                if m:
                    pos[m.group(2).upper()] = i
            changed.append(("incar.delete", k, old_val, None))

    if changed:
        ip.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    if log is not None:
        for sec, k, o, n in changed:
            if n is None:
                log.append("step.conf [%s] 删除 %s（原 %s）" % (sec, k, o))
            elif o is None:
                log.append("step.conf [%s] 新增 %s = %s" % (sec, k, n))
            else:
                log.append("step.conf [%s] 覆盖 %s: %s -> %s" % (sec, k, o, n))
    return changed



def _spans(lines):
    """-> [(节名, 起始行, 结束行)]；第一个 [xxx] 之前的内容算 params。"""
    out, cur, start = [], "params", 0
    for i, line in enumerate(lines):
        s = _strip(line)
        if s and _SECTION.match(s):
            out.append((cur, start, i))
            cur, start = _SECTION.match(s).group(1).lower(), i + 1
    out.append((cur, start, len(lines)))
    return out


def set_value(path, section, key, value):
    """就地改一个键：文件/节/键不存在则新建，注释与其它行原样保留。
    value=None 表示删除该键。返回 (旧值, 新值)。"""
    p = Path(path)
    lines = (p.read_text(encoding="utf-8-sig").splitlines()
             if p.is_file() else ["[params]"])
    section, ku = section.lower(), key.upper()
    hit, tail = None, None
    for name, a, b in _spans(lines):
        if name != section:
            continue
        for i in range(a, b):
            s = _strip(lines[i])
            if s and "=" in s and s.split("=", 1)[0].strip().upper() == ku:
                hit = i
        tail = b
        while tail > a and not lines[tail - 1].strip():
            tail -= 1
    old = lines[hit].split("=", 1)[1].strip() if hit is not None else None
    if value is None:
        if hit is not None:
            lines.pop(hit)
    elif hit is not None:
        keep = _COMMENT.search(lines[hit])
        lines[hit] = "%s = %s%s" % (lines[hit].split("=", 1)[0].rstrip(),
                                    value, keep.group(0) if keep else "")
    elif tail is not None:
        lines.insert(tail, "%s = %s" % (key, value))
    else:
        if lines and lines[-1].strip():
            lines.append("")
        lines += ["[%s]" % section, "%s = %s" % (key, value)]
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return old, value
