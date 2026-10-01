# -*- coding: utf-8 -*-
"""autozt check —— 技能运行依赖探测（哪台机器能跑哪个技能、缺什么、怎么补）。

依赖声明在各技能 skill.yaml 的 ``requires:``：

    requires:
      python: [numpy, phonopy]          # 必需 python 模块（pip 名，自动换成 import 名）
      optional_python: [pheasy]         # 按参数选用（缺了只影响对应引擎/选项）
      amset_python: [amset]             # 在集群 yaml 的 amset_env 环境里检查（ke/zt）
      exe: [vasp_std, vaspkit]          # 必需可执行文件
      optional_exe: [vasp_ncl]
      conda: atomate2_p_a               # 集群 yaml 没写 conda_env 时用的环境名
      gpu: true                         # 需要 NVIDIA GPU（技能名以 -gpu 结尾时自动视为 true）

隐含依赖：exe 里有 vasp_std → 还要 POTCAR 库（集群 yaml 的 potcar_dir）；
python 里有 mace-torch → 还要 MACE 模型目录（mlff_model_dir）。
vasp_std/vasp_gam/vasp_ncl 也可由集群 yaml 的 vasp: 配置（可执行文件名或绝对路径）满足。

探测在目标机器上做（本机直接跑，集群经 ssh 一次），环境激活方式与提交模板一致：
CONDA_ENV 是 venv 目录就激活它，否则 source conda_sh 再 conda activate，都没有就用 PATH。
"""
import base64
import json
import os
import shlex

# pip 名 → import 名（不在表里的按 "-"→"_" 处理）
IMPORT_NAME = {
    "mace-torch": "mace", "scikit-fem": "skfem", "torch-geometric": "torch_geometric",
    "pyyaml": "yaml", "scikit-learn": "sklearn", "BoltzTraP2": "BoltzTraP2",
}
VASP_EXES = ("vasp_std", "vasp_gam", "vasp_ncl")

_PROBE_PY = r'''
import importlib.util, json, os, shutil, sys
cfg = json.loads(sys.argv[1])
def _mod(m):
    try:
        return importlib.util.find_spec(m) is not None
    except Exception:
        return False
def _exe(name, alts):
    if shutil.which(name):
        return True
    for a in alts:
        a = os.path.expandvars(os.path.expanduser(a))
        if (os.path.sep in a and os.path.isfile(a) and os.access(a, os.X_OK)) or shutil.which(a):
            return True
    return False
def _vaspkit_potcar():
    """~/.vaspkit 里 POTCAR_TYPE 对应的 <TYPE>_PATH 是否指向存在的目录。"""
    p = os.path.expanduser("~/.vaspkit")
    if not os.path.isfile(p):
        return False
    kv = {}
    for ln in open(p, errors="ignore"):
        ln = ln.split("#", 1)[0]
        if "=" in ln:
            k, v = ln.split("=", 1)
            kv[k.strip().upper()] = v.strip().strip("'\"")
    t = (kv.get("POTCAR_TYPE") or "PBE").upper()
    path = kv.get(t + "_PATH") or kv.get("PBE_PATH")
    return bool(path) and os.path.isdir(os.path.expanduser(path))
out = {"python": sys.executable, "version": sys.version.split()[0],
       "mods": {m: _mod(m) for m in cfg["mods"]},
       "exes": {e: _exe(e, cfg["alts"].get(e, [])) for e in cfg["exes"]},
       "paths": {p: os.path.isdir(os.path.expandvars(os.path.expanduser(p)))
                 for p in cfg["paths"]},
       "gpu": bool(shutil.which("nvidia-smi")),
       "vaspkit_potcar": _vaspkit_potcar(),
       "active_env": os.environ.get("VIRTUAL_ENV") or os.environ.get("CONDA_DEFAULT_ENV") or ""}
print("@@JSON " + json.dumps(out))
'''


def _import_name(pkg):
    return IMPORT_NAME.get(pkg, str(pkg).replace("-", "_"))


def _as_list(v):
    if not v:
        return []
    if isinstance(v, (list, tuple)):
        return [str(x) for x in v if str(x).strip()]
    return [str(v)]


def skill_requires(cfg, key):
    """合并后的技能骨架里的 requires（apply_skills 已把 skill.yaml 的 requires 挂到 _skill_requires）。"""
    t = (cfg.get("task_types") or {}).get(key) or {}
    req = dict(t.get("_skill_requires") or {})
    if not req:
        from autozt import _load_yaml_file, _PKG_ROOT
        p = os.path.join(_PKG_ROOT, "skill", key, "skill.yaml")
        req = dict((_load_yaml_file(p) or {}).get("requires") or {})
    return req


def load_cluster(name):
    from autozt import _load_yaml_file, pkg_setting_path
    p = pkg_setting_path(str(name) + ".yaml")
    return (_load_yaml_file(p) or {}) if p else None


def list_clusters():
    from autozt import _PKG_ROOT
    d = os.path.join(_PKG_ROOT, "setting")
    out = []
    for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        if f.endswith(".yaml") and not f.startswith(("tf", "template-")):
            c = load_cluster(f[:-5]) or {}
            if "ssh_host" in c:
                out.append(f[:-5])
    return out


def _env_for(cluster, req, is_local):
    """(conda_sh, conda_env) 与运行时一致：集群值优先；本机显式留空 = PATH。"""
    sh, env = cluster.get("conda_sh"), cluster.get("conda_env")
    if env:
        return str(sh or ""), str(env)
    if is_local and "conda_env" in cluster:
        return str(sh or ""), ""
    return str(sh or ""), str(req.get("conda") or "")


def _plan(cfg, keys, cluster_name, cluster):
    is_local = "ssh_host" in cluster and not cluster.get("ssh_host")
    vasp = cluster.get("vasp") or {}
    alts = {}
    for prof in (vasp.values() if isinstance(vasp, dict) else []):
        if isinstance(prof, dict):
            for v in VASP_EXES:
                k = v.split("_")[1]
                if prof.get(k):
                    alts.setdefault(v, []).append(str(prof[k]))
    groups, items = {}, []
    for key in keys:
        req = skill_requires(cfg, key)
        if not req:
            items.append({"skill": key, "req": None})
            continue
        sh, env = _env_for(cluster, req, is_local)
        g = groups.setdefault((sh, env), {"mods": set(), "exes": set(), "paths": set()})
        mods = [_import_name(x) for x in _as_list(req.get("python"))]
        omods = [_import_name(x) for x in _as_list(req.get("optional_python"))]
        exes = _as_list(req.get("exe"))
        oexes = _as_list(req.get("optional_exe"))
        g["mods"].update(mods + omods)
        g["exes"].update(exes + oexes)
        paths = {}
        # POTCAR 来源：用 vaspkit 的技能走 ~/.vaspkit 的 <TYPE>_PATH（vaspkit 103 生成）；
        # 不用 vaspkit 的（defect-dft-cpu）读集群 yaml 的 potcar_dir（注入成 POTCAR_DIR）。
        _vk_potcar = False
        if any(e in VASP_EXES for e in exes):
            if (req.get("potcar") or ("vaspkit" if "vaspkit" in exes else "potcar_dir")) \
                    == "vaspkit":
                _vk_potcar = True
            else:
                paths["potcar_dir"] = str(cluster.get("potcar_dir") or "")
        if "mace" in mods:
            paths["mlff_model_dir"] = str(cluster.get("mlff_model_dir") or "")
        g["paths"].update(p for p in paths.values() if p)
        amods = [_import_name(x) for x in _as_list(req.get("amset_python"))]
        aenv = str(cluster.get("amset_env") or "")
        if amods:
            ga = groups.setdefault((sh, aenv or env), {"mods": set(), "exes": set(),
                                                      "paths": set()})
            ga["mods"].update(amods)
        items.append({"skill": key, "req": req, "sh": sh, "env": env, "aenv": aenv or env,
                      "aenv_raw": aenv,
                      "mods": mods, "omods": omods, "amods": amods, "exes": exes,
                      "oexes": oexes, "paths": paths, "vk_potcar": _vk_potcar,
                      "gpu": bool(req.get("gpu")) or key.endswith("-gpu")})
    return groups, items, alts, is_local


def _probe_script(groups, alts):
    lines = ["_act() { local SH=\"$1\" ENV=\"$2\"; SH=\"${SH/#\\~/$HOME}\"; ENV=\"${ENV/#\\~/$HOME}\";",
             "  if [ -n \"$ENV\" ] && [ -x \"$ENV/bin/python\" ]; then source \"$ENV/bin/activate\";",
             "  elif [ -n \"$SH\" ] && [ -f \"$SH\" ]; then source \"$SH\";",
             "    if [ -n \"$ENV\" ]; then conda activate \"$ENV\"; fi;",
             "  elif [ -n \"$ENV\" ] && command -v conda >/dev/null 2>&1; then",
             "    eval \"$(conda shell.bash hook 2>/dev/null)\"; conda activate \"$ENV\"; fi; }",
             "_PYB64=%s" % shlex.quote(base64.b64encode(_PROBE_PY.encode()).decode())]
    for i, ((sh, env), g) in enumerate(sorted(groups.items())):
        payload = json.dumps({"mods": sorted(g["mods"]), "exes": sorted(g["exes"]),
                              "paths": sorted(g["paths"]), "alts": alts})
        lines.append("echo '@@GROUP %d'" % i)
        lines.append("( _act %s %s >/dev/null 2>&1; _py=$(command -v python3 || command -v python); "
                     "if [ -z \"$_py\" ]; then echo '@@NOPY'; else "
                     "echo \"$_PYB64\" | base64 -d > \"${TMPDIR:-/tmp}/.autozt_probe_$$.py\"; "
                     "\"$_py\" \"${TMPDIR:-/tmp}/.autozt_probe_$$.py\" %s; "
                     "rm -f \"${TMPDIR:-/tmp}/.autozt_probe_$$.py\"; fi ) 2>/dev/null"
                     % (shlex.quote(sh), shlex.quote(env), shlex.quote(payload)))
    return "\n".join(lines) + "\n"


def _parse(out, n):
    res, cur = {}, None
    for ln in (out or "").splitlines():
        if ln.startswith("@@GROUP "):
            cur = int(ln.split()[1])
        elif ln.startswith("@@NOPY") and cur is not None:
            res[cur] = {"python": None}
        elif ln.startswith("@@JSON ") and cur is not None:
            try:
                res[cur] = json.loads(ln[7:])
            except ValueError:
                pass
    return [res.get(i) for i in range(n)]


def _env_active(want, active):
    """请求的环境（名字或 venv 路径）是否真的激活了。"""
    if not want:
        return True
    if not active:
        return False
    w = os.path.basename(os.path.normpath(os.path.expanduser(str(want))))
    a = os.path.basename(os.path.normpath(str(active)))
    return w == a or str(want) == str(active)


def _hint(kind, name, cluster_name, env):
    where = ("（环境 %s）" % env) if env else "（PATH 里的 python）"
    if kind == "python":
        return "pip install %s%s" % (name, where)
    if kind == "exe" and name in VASP_EXES:
        return ("在 setting/%s.yaml 的 vasp.standard.%s 写 VASP 可执行文件（名字或绝对路径；"
                "VASP 需自行取得授权）" % (cluster_name, name.split("_")[1]))
    if kind == "exe":
        return "安装 %s 并放进 PATH%s" % (name, where)
    if kind == "potcar_dir":
        return "在 setting/%s.yaml 写 potcar_dir（VASP 赝势库根目录）" % cluster_name
    if kind == "mlff_model_dir":
        return "在 setting/%s.yaml 写 mlff_model_dir（MACE 模型目录）" % cluster_name
    if kind == "gpu":
        return "需要 NVIDIA GPU（nvidia-smi）；本机没有就把材料切到 GPU 集群（autozt hpc）"
    if kind == "python3":
        return "目标环境里找不到 python3/python"
    return name


def check(cfg, keys, cluster_name="local"):
    """返回 {"cluster", "host", "ok", "error", "skills": [...]}；每个技能
    status = ready | missing | unknown（未声明依赖）。"""
    from autozt.collect import LOCAL_HOST, run_remote
    cluster = load_cluster(cluster_name)
    if cluster is None:
        return {"cluster": cluster_name, "ok": False,
                "error": "没有集群配置 setting/%s.yaml" % cluster_name, "skills": []}
    groups, items, alts, is_local = _plan(cfg, keys, cluster_name, cluster)
    order = sorted(groups)
    host = LOCAL_HOST if is_local else (cluster.get("ssh_host") or "__default__")
    probed = []
    if order:
        rc, out = run_remote(dict(cfg, remote_path_prefix=cluster.get("remote_path_prefix")
                                  or cfg.get("remote_path_prefix")),
                             _probe_script({k: groups[k] for k in order}, alts),
                             host=host, use_stdin=True)
        probed = _parse(out, len(order))
        if not any(probed):
            return {"cluster": cluster_name, "host": cluster.get("ssh_host") or "本机",
                    "ok": False, "error": "探测失败（rc=%s）：%s" % (rc, (out or "")[-400:]),
                    "skills": []}
    by_env = dict(zip(order, probed))
    report = []
    for it in items:
        if it["req"] is None:
            report.append({"skill": it["skill"], "status": "unknown",
                           "missing": [], "optional_missing": [],
                           "hints": ["skill.yaml 没有 requires 声明，无法判断"]})
            continue
        r = by_env.get((it["sh"], it["env"])) or {}
        ra = by_env.get((it["sh"], it["aenv"])) or {}
        miss, omiss, hints = [], [], []
        if r.get("python") and it["env"] and not _env_active(it["env"], r.get("active_env")):
            hints.append("环境 %s 没能激活（conda_sh=%s），以上结果按 %s 判断；"
                         "检查 setting/%s.yaml 的 conda_sh/conda_env，或都留空用 PATH"
                         % (it["env"], it["sh"] or "未写", r.get("python"), cluster_name))
        if not r.get("python"):
            miss.append("python3")
            hints.append(_hint("python3", "", cluster_name, it["env"]))
        else:
            for m, pkg in zip(it["mods"], _as_list(it["req"].get("python"))):
                if not r["mods"].get(m):
                    miss.append("py:" + pkg)
                    hints.append(_hint("python", pkg, cluster_name, it["env"]))
            for m, pkg in zip(it["omods"], _as_list(it["req"].get("optional_python"))):
                if not r["mods"].get(m):
                    omiss.append("py:" + pkg)
            for m, pkg in zip(it["amods"], _as_list(it["req"].get("amset_python"))):
                if not (ra.get("mods") or {}).get(m):
                    miss.append("py:%s@amset_env" % pkg)
                    hints.append(("在 setting/%s.yaml 写 amset_env，并在该环境里 pip install %s"
                                  % (cluster_name, pkg)) if not it["aenv_raw"]
                                 else _hint("python", pkg, cluster_name, it["aenv"]))
            for e in it["exes"]:
                if not r["exes"].get(e):
                    miss.append("exe:" + e)
                    hints.append(_hint("exe", e, cluster_name, it["env"]))
            for e in it["oexes"]:
                if not r["exes"].get(e):
                    omiss.append("exe:" + e)
            for kind, p in it["paths"].items():
                if not p or not r["paths"].get(p):
                    miss.append(kind)
                    hints.append(_hint(kind, p, cluster_name, it["env"]))
            if it["vk_potcar"] and r["exes"].get("vaspkit") and not r.get("vaspkit_potcar"):
                miss.append("vaspkit_potcar")
                hints.append("在 ~/.vaspkit 配置 POTCAR 库路径（如 PBE_PATH = ~/POTCAR/PBE），"
                             "vaspkit 103 要用它生成 POTCAR")
            if it["gpu"] and not r.get("gpu"):
                miss.append("gpu")
                hints.append(_hint("gpu", "", cluster_name, it["env"]))
        for kind in ("potcar_dir", "mlff_model_dir"):
            if kind in it["paths"] and not it["paths"][kind] and kind not in miss:
                miss.append(kind)
                hints.append(_hint(kind, "", cluster_name, it["env"]))
        report.append({"skill": it["skill"], "status": "missing" if miss else "ready",
                       "missing": miss, "optional_missing": omiss,
                       "python": r.get("python"), "env": it["env"] or "(PATH)",
                       "env_active": bool(not it["env"] or _env_active(it["env"],
                                                                      r.get("active_env"))),
                       "hints": sorted(set(hints))})
    return {"cluster": cluster_name, "host": cluster.get("ssh_host") or "本机",
            "ok": True, "skills": report}


_CACHE = {}


def missing_for(cfg, key, cluster_name):
    """本进程内缓存的单技能探测：返回缺失的必需项列表（探测失败返回 None）。"""
    ck = (cluster_name, key)
    if ck not in _CACHE:
        rep = check(cfg, [key], cluster_name)
        sk = (rep.get("skills") or [None])[0] if rep.get("ok") else None
        _CACHE[ck] = sk
    sk = _CACHE[ck]
    if sk is None or sk.get("status") == "unknown":
        return None
    return sk


def render(rep):
    lines = ["集群 %s（%s）" % (rep.get("cluster"), rep.get("host") or "未配置")]
    if not rep.get("ok"):
        return "\n".join(lines + ["  " + str(rep.get("error"))])
    mark = {"ready": "✅", "missing": "❌", "unknown": "？"}
    for s in rep["skills"]:
        line = "  %s %-18s" % (mark.get(s["status"], "?"), s["skill"])
        if s["status"] == "ready":
            line += " 可运行"
            if s.get("optional_missing"):
                line += "（可选缺：%s）" % ", ".join(s["optional_missing"])
            if not s.get("env_active", True):
                line += "  ⚠ 环境 %s 未激活，按 %s 判断" % (s.get("env"), s.get("python"))
        elif s["status"] == "missing":
            line += " 缺：%s" % ", ".join(s["missing"])
        else:
            line += " 未声明依赖"
        lines.append(line)
        if s["status"] == "missing":
            lines += ["       → %s" % h for h in s.get("hints") or []]
    n = sum(1 for s in rep["skills"] if s["status"] == "ready")
    lines.append("  合计：%d/%d 个技能在 %s 上可运行" % (n, len(rep["skills"]), rep["cluster"]))
    return "\n".join(lines)
