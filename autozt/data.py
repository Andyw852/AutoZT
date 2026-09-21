"""data —— 数据采集 / 过滤 / 缓存（真实深模块）。

数据层：collect_data、过滤、状态缓存和快照。
被 15_hpc / 16_watch / 17_cli 共同依赖——抽成真模块后环2 即破
（15/16 依赖 data 而非 17）。跨分片依赖用函数内延迟 import 从包命名空间取。
"""
import os
import sys
import time
import hashlib
import json

def _dbg_t(label, t0):
    """AUTOZT_DEBUG_TIME=1 时向 stderr 打印各阶段耗时。"""
    if os.environ.get("AUTOZT_DEBUG_TIME"):
        import time as _t
        print("[计时] %s: %.1fs" % (label, _t.time() - t0), file=sys.stderr)

def _state_cache_path(cfg):
    return os.path.join(cfg.get("_config_dir") or os.getcwd(),
                        ".tf_state_cache.json")

def _state_cache_sig(cfg, types, tt, root):
    """缓存键：主配置 mtime+size + 采集范围指纹。"""
    cst = None
    cp = cfg.get("_config_path")
    if cp and os.path.isfile(cp):
        try:
            _s = os.stat(cp)
            cst = (os.path.realpath(cp), _s.st_mtime_ns, _s.st_size)
        except OSError:
            pass
    import hashlib
    h = hashlib.sha1()
    for t in types:
        h.update(("%s|%s|%s\n" % (t.get("key", ""), t.get("root", ""),
                                  t.get("local_root", ""))).encode("utf-8"))
    return (cst, tt, root, h.hexdigest())

def _state_cache_save(cfg, data, types, tt, root):
    try:
        import tempfile as _tf
        p = _state_cache_path(cfg)
        d = os.path.dirname(p) or "."
        os.makedirs(d, exist_ok=True)
        payload = {"ts": time.time(),
                   "sig": _state_cache_sig(cfg, types, tt, root),
                   "data": data}
        _fd, _tmp = _tf.mkstemp(dir=d, prefix=".tf_state.", suffix=".tmp")
        try:
            with os.fdopen(_fd, "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=False)
            os.replace(_tmp, p)
        finally:
            if os.path.exists(_tmp):
                try:
                    os.remove(_tmp)
                except OSError:
                    pass
    except Exception:
        pass

def _state_cache_load(cfg, types, tt, root, ttl):
    if ttl <= 0:
        return None
    try:
        with open(_state_cache_path(cfg), encoding="utf-8") as f:
            payload = json.load(f)
    except Exception:
        return None
    if payload.get("sig") != _state_cache_sig(cfg, types, tt, root):
        return None
    if time.time() - float(payload.get("ts") or 0) > ttl:
        return None
    from autozt.bootstrap import filter_config_conflicts
    return filter_config_conflicts(payload.get("data") or {})

def collect_data(cfg, types):
    """按类型采集全部材料状态（远端/本地两段路径），供各命令及 watch 循环复用。"""
    from autozt import (collect, collect_v3_batch, _dedup_segments, annotate,
                       _queue_total, check_duplicates)
    data_types = []
    queue_by_host = {}
    v2 = [t for t in types if not t.get("local_root")]
    if v2:
        _d = collect(cfg, v2)
        data_types.extend(_d["types"])
        queue_by_host[cfg.get("host") or "local"] = _d.get("queue") or {}
    segs = [t for t in types if t.get("local_root")]
    if segs:   # v3.20：跨段批量采集（一次 ssh 完成全部材料）
        _te_list, _qbh = collect_v3_batch(cfg, segs)
        queue_by_host.update(_qbh)
        for te in _te_list:
            exist = next((x for x in data_types
                          if x["key"] == te["key"] and x.get("local")), None)
            if exist:
                exist["materials"] = _dedup_segments(exist["materials"]
                                                     + te["materials"])
            else:
                data_types.append(te)
    from autozt.bootstrap import filter_config_conflicts
    data = annotate(filter_config_conflicts({"host": cfg.get("host") or "local", "types": data_types}))
    if queue_by_host:
        data["queue"] = _queue_total(queue_by_host)
    local_by_key = {t["key"]: t for t in types}
    for t in data["types"]:
        lc = local_by_key.get(t["key"], {})
        t["steps_cfg"] = lc.get("steps") or []
        t["gen_dir"] = lc.get("gen_dir")
        t["gen_need"] = lc.get("gen_need")
        t["skill_dir"] = lc.get("skill_dir")
    check_duplicates(data)
    return data

def apply_exclude(data, exclude):
    """-x：跳过指定项目（全名 / basename / <项目名>/<完整名>，逗号分隔）。"""
    from autozt import _name_matches
    if not exclude:
        return
    ex = [x.strip() for x in exclude.split(",") if x.strip()]
    for t in data["types"]:
        t["materials"] = [m for m in t["materials"]
                          if not any(_name_matches(m, x) for x in ex)]

def filter_projs(data, projs):
    """只保留指定材料（全名 / basename / <项目名>/<完整名>）；空列表 = 不过滤。"""
    from autozt import _name_matches
    if not projs:
        return
    want = [x for x in projs if x]
    for t in data["types"]:
        t["materials"] = [m for m in t["materials"]
                          if any(_name_matches(m, x) for x in want)]

def filter_status(data, spec):
    """v1.4 -status：只保留"有步骤处于指定状态"的材料（对任意命令生效：
    status 只看它们，start/retry/rerun/stop 只操作它们）。
    状态词大小写不限、支持别名，逗号分隔多个。"""
    from autozt import STATUS_ALIAS
    kinds = set()
    for x in str(spec or "").split(","):
        x = x.strip()
        # patch_cell_word：ready 展开成 TODO + PREP 两个 kind
        if x and STATUS_ALIAS.get(x.lower()) == "TODO+PREP":
            kinds.update(("TODO", "PREP"))
            continue
        if x:
            kinds.add(STATUS_ALIAS.get(x.lower(), x.upper()))
    if not kinds:
        return
    for t in data["types"]:
        t["materials"] = [m for m in t["materials"]
                          if any(s["kind"] in kinds for s in m["steps"])]

def status_spec_has_scancel(spec):
    """-status 里显式含 scancel → start/retry 放行 SCANCEL 步骤。"""
    return any(x.strip().lower() in ("scancel", "scancelled", "cancel",
                                     "cancelled", "canceled")
               for x in str(spec or "").split(","))

def _snapshot(data):
    """状态指纹：材料 → 各步骤 (label, kind, 作业状态)，watch 据此判断有无变化。"""
    return json.dumps({m["name"]: [(s["label"], s["kind"],
                                    (s.get("job") or {}).get("state"))
                                   for s in m["steps"]]
                       for t in data["types"] for m in t["materials"]},
                      sort_keys=True)
