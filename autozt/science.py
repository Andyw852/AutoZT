# -*- coding: utf-8 -*-
"""面向模型的科学工作流层。

这里只读整理已有技能契约和已回收结果，不创建项目、不改配置、不提交作业。
"""
from __future__ import annotations

import json
import math
import os
import hashlib
from typing import Any, Dict, Iterable, List, Optional, Tuple

RESULT_SCHEMA = "autozt/science-result/1"
CONVERSATION_SCHEMA = "autozt/conversation/1"
_DIMENSIONS = (None, "0D", "1D", "2D", "3D")
_PROPERTY_ALIASES = {
    "zt": {"zt", "zt_value", "z_t"},
    "seebeck": {"seebeck", "seebeck_uv/k", "s"},
    "sigma": {"sigma", "conductivity", "sigma_s/m"},
    "kappa_e": {"kappa_e", "kappa_e_w/mk", "electronic_thermal_conductivity"},
    "kappa_l": {"kappa_l", "kappa_l_w/mk", "kappa_lattice", "kappa_l_w/m/k"},
    "pf": {"pf", "pf_w/mk2", "power_factor"},
}
_PROPERTY_UNITS = {"zt": "1", "seebeck": "uV/K", "sigma": "S/m",
                   "kappa_e": "W/m/K", "kappa_l": "W/m/K", "pf": "W/m/K^2"}
_ROW_KEYS = {"zt": "ZT", "seebeck": "seebeck_uV/K", "sigma": "sigma_S/m",
             "kappa_e": "kappa_e_W/mK", "kappa_l": "kappa_L_W/mK", "pf": "PF_W/mK2"}


def _plan_id(goal: str, *, dimension: Optional[str], temperature: Iterable[Any],
             carrier: Iterable[Any], material: Optional[str], skills: Iterable[str]) -> str:
    """Create a stable identifier so a plan can be reviewed across transports."""
    payload = {"goal": goal, "dimension": dimension, "temperature": list(temperature),
               "carrier": list(carrier), "material": material, "skills": list(skills)}
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True,
                          separators=(",", ":")).encode("utf-8")
    return "plan-" + hashlib.sha256(encoded).hexdigest()[:16]


def _conversation(*, state: str, message: str, next_action: str,
                  requires_confirmation: bool, actions: Optional[Iterable[Any]] = None,
                  risk: str = "read", will_submit_jobs: bool = False,
                  plan_id: Optional[str] = None) -> Dict[str, Any]:
    """Return the small common envelope shared by science-facing interfaces."""
    scope = list(actions or [])
    return {
        "schema_version": CONVERSATION_SCHEMA,
        "conversation_state": state,
        "message": message,
        "next_action": next_action,
        "requires_user_confirmation": bool(requires_confirmation),
        "confirmation_payload": {
            "plan_id": plan_id,
            "actions": scope,
            "risk": risk,
            "will_submit_jobs": bool(will_submit_jobs),
        },
    }


def _as_list(value: Optional[Iterable[Any]]) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [value]
    return list(value)


def _float(value: Any) -> Optional[float]:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _close(a: Any, b: Any, tol: float = 1e-6) -> bool:
    x, y = _float(a), _float(b)
    if x is None or y is None:
        return str(a) == str(b)
    return abs(x - y) <= tol * max(1.0, abs(x), abs(y))


def _contains(values: Iterable[Any], wanted: Any, tol: float = 1e-6) -> bool:
    return any(_close(value, wanted, tol) for value in values)


def _canonical_property(name: str) -> str:
    key = str(name or "").strip().lower()
    for canonical, aliases in _PROPERTY_ALIASES.items():
        if key in aliases:
            return canonical
    return key


def _skill_name(item: Any) -> str:
    if not isinstance(item, dict):
        return str(item)
    return str(item.get("skill") or item.get("name") or item.get("key") or "")


def research_plan(goal: str, skills: Iterable[Dict[str, Any]], *, dimension: Optional[str] = None,
                  temperature: Optional[Iterable[Any]] = None,
                  carrier: Optional[Iterable[Any]] = None,
                  material: Optional[str] = None) -> Dict[str, Any]:
    """生成可审查、可交给动作层执行的计划，但绝不隐式提交作业。"""
    text = (goal or "").strip()
    available = [_skill_name(x) for x in skills if _skill_name(x)]
    low = text.lower()
    wants_zt = any(k in low for k in ("zt", "热电", "thermoelectric", "thermoelectricity"))
    selected: List[str] = []
    for candidate in ("zt-dft-cpu", "ke-dft-cpu", "kl-dft-cpu"):
        if candidate in available and (candidate == "zt-dft-cpu" or wants_zt):
            selected.append(candidate)
    if not selected and wants_zt:
        selected = [x for x in ("ke-dft-cpu", "kl-dft-cpu") if x in available]

    temps, carriers = _as_list(temperature), _as_list(carrier)
    gaps: List[Dict[str, str]] = []
    if not text:
        gaps.append({"code": "missing_goal", "message": "需要研究目标，例如计算 zT(T,n)。"})
    if wants_zt and not selected:
        gaps.append({"code": "no_zt_skill", "message": "技能目录中没有可用的 zT/电子输运/晶格热导技能。"})
    if dimension not in _DIMENSIONS:
        gaps.append({"code": "invalid_dimension", "message": "dimension 必须是 0D/1D/2D/3D。"})
    if wants_zt and not temps:
        gaps.append({"code": "missing_temperature_grid", "message": "请提供目标温度网格，避免跨步骤温度不一致。"})
    if wants_zt and not carriers:
        gaps.append({"code": "missing_carrier_grid", "message": "请提供载流子浓度/化学势网格，才能定义 zT(n,T)。"})
    if dimension == "2D" and not material:
        gaps.append({"code": "missing_material", "message": "二维工作流需要材料名，便于绑定 POSCAR 和 thickness_2d.json。"})

    stages: List[Dict[str, Any]] = []
    if selected:
        stages = [
            {"id": "structure", "stage": "structure", "skills": selected[:1],
             "requires": ["POSCAR"], "produces": ["CONTCAR", "workflow_method.txt"],
             "validator": "结构优化完成且维度已确定"},
            {"id": "electronic_transport", "stage": "electronic_transport",
             "skills": [x for x in selected if x in ("zt-dft-cpu", "ke-dft-cpu")],
             "requires": ["CONTCAR", "band structure", "dielectric", "deformation potential"],
             "produces": ["transport.json"], "validator": "S/σ/κe 通过技能判据"},
            {"id": "lattice_transport", "stage": "lattice_transport",
             "skills": [x for x in selected if x in ("zt-dft-cpu", "kl-dft-cpu")],
             "requires": ["CONTCAR", "force constants", "phonon stability"],
             "produces": ["kappa_summary.json"], "validator": "KAPPA_DONE 且声子/口径检查通过"},
            {"id": "combine", "stage": "combine", "skills": ["zt-dft-cpu"] if "zt-dft-cpu" in selected else [],
             "requires": ["transport.json", "kappa_summary.json"],
             "produces": ["zt_summary.json", "zt_summary.txt"],
             "validator": "ZT_DONE、温度/载流子/单位和二维厚度口径一致"},
        ]

    primary = "zt-dft-cpu" if "zt-dft-cpu" in selected else (selected[0] if selected else None)
    actions = []
    if primary and material:
        actions = [
            {"action": "init_project", "command": ["autozt", "-tt", primary, "-p", material, "init"],
             "risk": "local_project_setup", "execute": "explicit_only"},
            {"action": "start_workflow", "command": ["autozt", "-tt", primary, "-p", material, "start"],
             "risk": "may_submit_hpc_jobs", "execute": "explicit_only"},
            {"action": "advance_ready", "command": ["autozt", "-tt", primary, "-p", material, "auto", "on"],
             "risk": "may_submit_ready_jobs", "execute": "explicit_only"},
        ]
    plan_id = _plan_id(text, dimension=dimension, temperature=temps, carrier=carriers,
                       material=material, skills=selected)
    conversation = _conversation(
        state="needs_user_input" if gaps else "awaiting_confirmation",
        message=("研究计划缺少必要输入：" + "；".join(item["message"] for item in gaps)
                 if gaps else "研究计划已生成，当前只读预览，不会提交超算作业。"),
        next_action="provide_inputs" if gaps else "confirm_plan",
        requires_confirmation=not gaps,
        actions=actions,
        risk="mutate" if actions else "read",
        will_submit_jobs=bool(actions), plan_id=plan_id)
    return {
        "schema_version": "autozt/research-plan/1", "status": "needs_input" if gaps else "ready",
        "plan_id": plan_id, **conversation,
        "plan_summary": ("二维材料热电 zT 工作流" if dimension == "2D" else "热电 zT 工作流"),
        "plan_steps": stages,
        "goal": text, "selected_skills": selected, "available_skills": available, "material": material,
        "assumptions": {"dimension": dimension, "temperature": temps, "carrier": carriers},
        "required_inputs": [
            {"name": "POSCAR", "source": "user", "required": True, "units": None},
            {"name": "temperature", "source": "user", "required": wants_zt, "units": "K", "values": temps},
            {"name": "carrier", "source": "user", "required": wants_zt, "units": "cm-3", "values": carriers},
            {"name": "thickness_2d.json", "source": "workflow", "required": dimension == "2D", "units": "A"},
        ],
        "stages": stages, "actions": actions,
        "execution": {"mode": "dry_run", "submits_jobs": False,
                       "next": "review gaps, then call inspect/cycle or the explicit actions"},
        "gaps": gaps,
        "review": ["确认上游结构和单位", "确认二维厚度约定", "确认温度/载流子网格",
                   "确认提交集群和资源后再执行 actions"],
    }


def _json_files(root: str) -> Iterable[str]:
    if not root or not os.path.isdir(root):
        return
    for base, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d not in (".git", "__pycache__")]
        for name in files:
            path = os.path.join(base, name)
            try:
                if name.endswith(".json") and os.path.getsize(path) <= 8 * 1024 * 1024:
                    yield path
            except OSError:
                continue


def _load_jsons(root: str) -> Tuple[List[str], Dict[str, Any], List[Dict[str, Any]]]:
    paths = list(_json_files(root))
    loaded: Dict[str, Any] = {}
    errors: List[Dict[str, Any]] = []
    for path in paths:
        try:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = json.load(fh)
        except (OSError, ValueError, UnicodeError) as exc:
            errors.append({"path": path, "error": str(exc)})
    return paths, loaded, errors


def _by_basename(paths: Iterable[str], name: str) -> List[str]:
    return [p for p in paths if os.path.basename(p) == name]


def _first_loaded(paths: Iterable[str], loaded: Dict[str, Any], name: str) -> Tuple[Optional[str], Any]:
    for path in _by_basename(paths, name):
        if path in loaded:
            return path, loaded[path]
    return None, None


def _check(checks: List[Dict[str, Any]], code: str, ok: bool, message: str,
           evidence: Optional[Iterable[Any]] = None, actual: Any = None,
           expected: Any = None, severity: str = "error") -> None:
    checks.append({"code": code, "ok": bool(ok), "severity": severity,
                   "message": message, "evidence": list(evidence or []),
                   "actual": actual, "expected": expected})


def _number_series(obj: Any, *keys: str) -> List[Any]:
    if not isinstance(obj, dict):
        return []
    for key in keys:
        value = obj.get(key)
        if isinstance(value, list):
            return value
    return []


def preflight(root: str, *, dimension: Optional[str] = None, thickness: Optional[float] = None,
              temperature: Optional[Iterable[Any]] = None,
              carrier: Optional[Iterable[Any]] = None) -> Dict[str, Any]:
    """检查结果内容、口径和 validator，而不把文件存在误报成科学正确。"""
    paths, loaded, load_errors = _load_jsons(root)
    checks: List[Dict[str, Any]] = []
    _check(checks, "result_dir", bool(root and os.path.isdir(root)), "结果目录存在", [root])
    _check(checks, "json_parse", not load_errors, "JSON 产物均可解析", load_errors,
           actual=load_errors, expected="valid JSON")
    zt_path, zt = _first_loaded(paths, loaded, "zt_summary.json")
    tr_path, transport = _first_loaded(paths, loaded, "transport.json")
    kl_path, kappa = _first_loaded(paths, loaded, "kappa_summary.json")
    th_path, thick = _first_loaded(paths, loaded, "thickness_2d.json")
    _check(checks, "zt_file", zt is not None, "已找到 zT 汇总产物", [zt_path] if zt_path else [])
    _check(checks, "electronic_file", transport is not None, "已找到电子输运产物 transport.json", [tr_path] if tr_path else [])
    _check(checks, "lattice_file", kappa is not None, "已找到晶格热导产物 kappa_summary.json", [kl_path] if kl_path else [])

    if isinstance(zt, dict):
        _check(checks, "zt_validator", zt.get("ZT_DONE") is True,
               "zT 汇总 validator 已通过", [zt_path], zt.get("ZT_DONE"), True)
        expected_units = {"S": "uV/K", "sigma": "S/m", "kappa_e": "W/m/K",
                          "kappa_L": "W/m/K", "T": "K", "PF": "W/m/K^2"}
        units = zt.get("units") if isinstance(zt.get("units"), dict) else {}
        _check(checks, "units", all(units.get(k) == v for k, v in expected_units.items()),
               "zT 单位契约完整且匹配", [zt_path], units, expected_units)
        zt_dim = str(zt.get("dim") or "").upper()
        if dimension:
            _check(checks, "dimension", zt_dim == dimension, "结果维度与请求一致", [zt_path], zt_dim, dimension)
        zt_temps = _number_series(zt, "temperatures")
        wanted_temps = _as_list(temperature)
        if wanted_temps:
            _check(checks, "temperature_grid", all(_contains(zt_temps, x) for x in wanted_temps),
                   "zT 温度网格覆盖请求", [zt_path], zt_temps, wanted_temps)
        wanted_carrier = _as_list(carrier)
        zt_carrier = _number_series(zt, "doping_cm-3", "carrier_cm-3")
        if wanted_carrier:
            _check(checks, "carrier_grid", all(_contains(zt_carrier, x, 1e-5) for x in wanted_carrier),
                   "zT 载流子网格覆盖请求", [zt_path], zt_carrier, wanted_carrier)
        sources = zt.get("sources") if isinstance(zt.get("sources"), dict) else {}
        _check(checks, "provenance", bool(sources) and bool(zt.get("formula")),
               "zT 结果包含公式和来源信息", [zt_path], sources, "formula + sources")
        caliber = zt.get("cell_caliber_check")
        if isinstance(caliber, dict) and "ok" in caliber:
            _check(checks, "cell_caliber", caliber.get("ok") is True,
                   "电子/晶格元胞口径一致", [zt_path], caliber, {"ok": True})

    if isinstance(transport, dict):
        has_transport = any(k in transport for k in
                            ("temperatures", "seebeck", "sigma", "conductivity",
                             "electronic_thermal_conductivity", "doping"))
        _check(checks, "transport_validator", has_transport,
               "电子输运 JSON 含温度/输运字段", [tr_path], list(transport), "transport fields")
    if isinstance(kappa, dict):
        _check(checks, "kappa_validator", kappa.get("KAPPA_DONE") is True,
               "晶格热导 validator 已通过", [kl_path], kappa.get("KAPPA_DONE"), True)
        ktemps = _number_series(kappa, "temperatures")
        if isinstance(zt, dict) and zt.get("temperatures") and ktemps:
            ztemps = _number_series(zt, "temperatures")
            outside = [x for x in ztemps if not _contains(ktemps, x)]
            rows = zt.get("rows") if isinstance(zt.get("rows"), list) else []
            null_declared = any(isinstance(row, dict) and any(v is None for v in row.get("kappa_L_W/mK", []))
                                for row in rows)
            _check(checks, "kappa_temperature_coverage", not outside or null_declared,
                   "κL 温区不足时已显式标记不可计算点", [kl_path, zt_path],
                   {"outside": outside, "null_declared": null_declared}, "null for uncovered")

    if dimension == "2D":
        _check(checks, "2d_thickness_file", thick is not None,
               "二维结果包含共享厚度契约 thickness_2d.json", [th_path] if th_path else [])
        actual = thick.get("thickness_d_A", thick.get("thickness_d_ang")) if isinstance(thick, dict) else None
        _check(checks, "2d_thickness_value", _float(actual) is not None and _float(actual) > 0,
               "二维厚度是正数并带有数值口径", [th_path] if th_path else [], actual, "> 0 A")
        if thickness is not None:
            _check(checks, "2d_thickness_match", _close(actual, thickness, 1e-3),
                   "二维厚度与请求一致", [th_path] if th_path else [], actual, thickness)

    passed = all(item["ok"] for item in checks)
    conversation = _conversation(
        state="ready_to_execute" if passed else "preflight_review",
        message=("预检查通过，可进入执行前确认。" if passed else
                 "预检查发现需要复核的输入、单位或结果证据。"),
        next_action="confirm_execution" if passed else "review_preflight",
        requires_confirmation=passed, risk="mutate", will_submit_jobs=False)
    return {"schema_version": "autozt/preflight/1",
            "status": "pass" if passed else "review", **conversation,
            "root": root, "checks": checks, "files_scanned": len(paths),
            "loaded_files": len(loaded), "parse_errors": load_errors,
            "requested": {"dimension": dimension, "temperature": _as_list(temperature),
                          "carrier": _as_list(carrier), "thickness_A": thickness},
            "note": "检查包含文件内容、单位、validator、网格和二维厚度；最终物理正确性仍由技能 validator 决定。"}


def _result_record(canonical: str, value: Any, *, source: str, path: str,
                   context: Optional[Dict[str, Any]] = None, unit: Optional[str] = None,
                   method: Optional[str] = None, validator: Optional[Dict[str, Any]] = None,
                   provenance: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    ctx = dict(context or {})
    return {"property": canonical, "value": value, "unit": unit or _PROPERTY_UNITS.get(canonical),
            "temperature": ctx.get("temperature"), "carrier": ctx.get("carrier"),
            "direction": ctx.get("direction"), "method": method or "unknown",
            "status": "complete" if validator and validator.get("ok") else "observed",
            "validator": validator or {"ok": None, "status": "not_declared"},
            "provenance": provenance or {"source": source, "field_path": path},
            "source": source, "path": path, "context": ctx}


def _zt_results(path: str, data: Dict[str, Any], wanted: str,
                temperature: Any, carrier: Any, direction: Optional[str]) -> List[Dict[str, Any]]:
    rows = data.get("rows") if isinstance(data.get("rows"), list) else []
    temps = data.get("temperatures") if isinstance(data.get("temperatures"), list) else []
    dopings = data.get("doping_cm-3") if isinstance(data.get("doping_cm-3"), list) else []
    dim = str(data.get("dim") or "").lower()
    key = _ROW_KEYS.get(wanted)
    if not key:
        return []
    validator = {"ok": data.get("ZT_DONE") is True,
                 "status": "passed" if data.get("ZT_DONE") is True else "failed",
                 "marker": "ZT_DONE", "source": path}
    sources = data.get("sources") if isinstance(data.get("sources"), dict) else {}
    notes = data.get("notes") if isinstance(data.get("notes"), list) else []
    out: List[Dict[str, Any]] = []
    for i, row in enumerate(rows):
        if not isinstance(row, dict) or not isinstance(row.get(key), list):
            continue
        row_carrier = row.get("doping_cm-3", dopings[i] if i < len(dopings) else None)
        if carrier is not None and not _close(row_carrier, carrier, 1e-5):
            continue
        for j, value in enumerate(row[key]):
            temp = temps[j] if j < len(temps) else None
            if temperature is not None and not _close(temp, temperature):
                continue
            ctx = {"temperature": temp, "carrier": row_carrier,
                   "direction": direction or ("in_plane" if dim == "2d" else "cell_average")}
            field_path = "/rows/%d/%s/%d" % (i, key, j)
            out.append(_result_record(wanted, value, source=path, path=field_path, context=ctx,
                                      method=str(data.get("formula") or data.get("method") or "zt_summary"),
                                      validator=validator,
                                      provenance={"source": path, "field_path": field_path,
                                                  "sources": sources, "notes": notes}))
    return out


def query_results(root: str, property_name: str, *, temperature: Optional[float] = None,
                  carrier: Optional[float] = None, direction: Optional[str] = None) -> Dict[str, Any]:
    """按标准结果对象查询；没有标准产物时保留通用 JSON 兜底。"""
    paths, loaded, errors = _load_jsons(root)
    wanted = _canonical_property(property_name)
    hits: List[Dict[str, Any]] = []
    for path in _by_basename(paths, "zt_summary.json"):
        if isinstance(loaded.get(path), dict):
            hits.extend(_zt_results(path, loaded[path], wanted, temperature, carrier, direction))

    if not hits:
        aliases = _PROPERTY_ALIASES.get(wanted, {str(property_name).lower()})
        seen = set()
        for path, data in loaded.items():
            if not isinstance(data, (dict, list)):
                continue
            units = data.get("units", {}) if isinstance(data, dict) else {}
            validator = {"ok": data.get("KAPPA_DONE") is True if isinstance(data, dict) and "KAPPA_DONE" in data else None,
                         "status": "passed" if isinstance(data, dict) and data.get("KAPPA_DONE") is True else "not_declared",
                         "source": path}

            def walk(obj: Any, trail: str = "", context: Optional[Dict[str, Any]] = None) -> None:
                ctx = dict(context or {})
                if isinstance(obj, dict):
                    for ck in ("temperature", "Temperature", "T", "carrier", "Carrier", "doping", "Doping", "direction"):
                        if ck in obj and isinstance(obj[ck], (int, float, str)):
                            ctx[ck.lower()] = obj[ck]
                    for key, value in obj.items():
                        low = str(key).lower()
                        if (low in aliases or any(alias in low for alias in aliases)) and isinstance(value, (int, float, str)):
                            temp = ctx.get("temperature", ctx.get("t"))
                            car = ctx.get("carrier", ctx.get("doping"))
                            if (temperature is None or _close(temp, temperature)) and (carrier is None or _close(car, carrier, 1e-5)):
                                fp = trail + "/" + str(key)
                                ident = (path, fp, repr(value))
                                if ident not in seen:
                                    seen.add(ident)
                                    unit = units.get(key) if isinstance(units, dict) else None
                                    hits.append(_result_record(wanted, value, source=path, path=fp,
                                                                context={"temperature": temp, "carrier": car,
                                                                         "direction": ctx.get("direction") or direction},
                                                                unit=unit, method=str(data.get("method") or os.path.basename(path)),
                                                                validator=validator,
                                                                provenance={"source": path, "field_path": fp}))
                        walk(value, trail + "/" + str(key), ctx)
                elif isinstance(obj, list):
                    for i, value in enumerate(obj[:10000]):
                        walk(value, trail + "/" + str(i), ctx)
            walk(data)

    found = bool(hits)
    conversation = _conversation(
        state="completed" if found else "blocked",
        message=("已找到结果，并保留单位、validator 和 provenance。" if found else
                 "没有找到符合条件的结果，需检查产物、性质名称或筛选条件。"),
        next_action="review_result" if found else "inspect_results",
        requires_confirmation=False, risk="read", will_submit_jobs=False)
    return {"schema_version": RESULT_SCHEMA, "status": "found" if found else "not_found",
            **conversation,
            "property": property_name,
            "filters": {"temperature": temperature, "carrier": carrier, "direction": direction},
            "results": hits[:200], "count": min(len(hits), 200), "parse_errors": errors,
            "note": "结果包含数值、单位、方法、validator 状态和来源；科学解释仍需结合技能说明。"}
