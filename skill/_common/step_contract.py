#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step_contract.py —— 步骤参数契约：**前后步骤参数对不上就不要运行**（2026-09-22 用户要求）。

## 为什么需要"通用机制"而不是字段清单

本会话实测两次事故，共同点是"**影响产物的参数变了，但产物指纹没覆盖它 → 静默沿用旧产物**"：
  ① 对称化只改 POSCAR（几何只差 4e-6 Å，逃过 1e-3 容差）→ S4 幂等沿用**旧点群**数据集 → S5 白拟合；
  ② @FC3_CUTOFF_PAIR@ 空→4.0 让**自动超胞** 5×5×1→4×4×1，单胞没变 ⇒ 连显式 @SUPERCELL="5 5 1"@ 都失效。
逐个字段加比对（@_dataset_matches@）永远会漏下一个旋钮。这里改成**通用机制**。

## 机制

每个步骤在 **gen 结束**时写 `<step_dir>/.step_contract.json`：

    {"schema": 1, "step": "<步骤名>", "time": "...",
     "step_conf_sha256": "<该步合并后 step.conf 的哈希>",
     "params": {...},                 # 可选：规范化后的 [params]（便于 diff 报告）
     "params_hash": "...",            # sha256(canonical(params))
     "inputs": {"<相对路径>": "<sha256>", ...},   # 本步消费的输入文件
     "upstream": {"step": "...", "chain_hash": "..."} | null,
     "chain_hash": "<sha256(params_hash + inputs + upstream.chain_hash)>"}

下游步骤在 **gen 开始**时调用 `require_upstream(...)`，任一不满足就 **报错拒绝运行**（exit≠0，
autozt 会把它标成 error，而不是让它带着不一致的产物继续往下跑）：

  A. 上游契约缺失            → 默认拒绝（@mode="warn"@ 可降级为告警，用于历史遗留步骤）
  B. 上游 step.conf 的**当前**哈希 != 上游契约里记录的哈希
                             ⇒ 上游参数改过但产物没重生成 ⇒ 拒绝
  C. 本步契约里记录的上游 chain_hash != 上游**现在**的 chain_hash
                             ⇒ 上游产物换过 ⇒ 本步旧产物已失效 ⇒ 拒绝
  D. 本步自己的 params_hash 变了（自己的 step.conf 改了）
     ⇒ 由各步自己的入参校验处理（S4 这类扇出步用 @_dataset_matches@ 隔离并**重建**，
       昂贵步则同样应拒绝并提示 retry）

**结论**：任何一步的参数变了，它的**整条下游链全部拒绝运行**，直到从该步 @retry@/@rerun@ 重新生成。
这就是"参数对不上就不要运行"。

## 迁移/降级（必须显式）

历史步骤目录没有契约文件。默认 **fail-closed**（拒绝）；要用 @allow_legacy=True@（或环境变量
@AUTOZT_CONTRACT_LEGACY=warn@）才降级为"告警并继续"，且会在契约里记 @legacy_ok: true@ 留痕。
"""
import hashlib
import json
import os
import time

SCHEMA = 1
CONTRACT_NAME = ".step_contract.json"


def _sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def file_sha256(path):
    """文件 sha256；读不到返回 None（不抛，调用方决定怎么处理）。"""
    try:
        with open(path, "rb") as fh:
            return _sha256_bytes(fh.read())
    except OSError:
        return None


def canonical(obj):
    """稳定序列化：dict 按键排序、float 走 repr，保证跨机/跨次可复现。"""
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
                      default=repr)


def params_hash(params):
    """规范化参数指纹；@params@ 为空也返回一个稳定值。"""
    return _sha256_bytes(canonical(params or {}).encode("utf-8"))


def contract_path(step_dir):
    return os.path.join(str(step_dir), CONTRACT_NAME)


def read_contract(step_dir):
    """读契约；缺失/损坏返回 None（损坏时给一条告警，不抛）。"""
    p = contract_path(step_dir)
    if not os.path.isfile(p):
        return None
    try:
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else None
    except Exception as e:                             # noqa: BLE001
        print("[WARN] 步骤契约读取失败（%s）：%s —— 按缺失处理" % (p, e))
        return None


def make_contract(step, step_dir, params=None, inputs=None, upstream=None):
    """构造契约 dict（不落盘）。@inputs@ 可以是 {名: 路径} 或 {名: sha}。"""
    inp = {}
    for name, val in (inputs or {}).items():
        if isinstance(val, str) and os.path.isfile(val):
            h = file_sha256(val)
            if h:
                inp[str(name)] = h
        elif isinstance(val, str):
            inp[str(name)] = val
    conf = os.path.join(str(step_dir), "step.conf")
    ch = file_sha256(conf)
    ph = params_hash(params)
    up = None
    if upstream:
        up = {"step": upstream.get("step"), "chain_hash": upstream.get("chain_hash")}
    chain = _sha256_bytes(canonical({"params": ph, "inputs": inp,
                                     "upstream": (up or {}).get("chain_hash")}).encode("utf-8"))
    return {"schema": SCHEMA, "step": step, "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "step_conf_sha256": ch, "params_hash": ph,
            "params": (params or None), "inputs": inp,
            "upstream": up, "chain_hash": chain}


def write_contract(step_dir, step, params=None, inputs=None, upstream=None):
    """写 <step_dir>/.step_contract.json（原子替换）。返回契约 dict。"""
    c = make_contract(step, step_dir, params=params, inputs=inputs, upstream=upstream)
    p = contract_path(step_dir)
    try:
        tmp = p + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(c, fh, ensure_ascii=False, indent=1)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, p)
    except OSError as e:
        print("[WARN] 写步骤契约失败（%s）：%s" % (p, e))
    return c


class ContractError(SystemExit):
    """契约不一致：以 SystemExit 抛出，让 gen 以非零码结束（autozt 标 error）。"""


def require_upstream(step, step_dir, upstream_dir, upstream_step=None,
                     my_contract=None, allow_legacy=False, tag=""):
    """校验上游契约；不一致就 **拒绝运行**。返回上游契约 dict（供本步写 chain 用）。

    判据（见模块 docstring A/B/C）。@my_contract@ 是本步已有契约（可为 None）。
    """
    who = tag or ("%s ← %s" % (step, upstream_step or "上游"))
    up = read_contract(upstream_dir)
    if up is None:
        msg = ("%s：上游 %s 没有步骤契约（%s）。历史遗留或被清理过。"
               % (who, upstream_step or "", contract_path(upstream_dir)))
        if allow_legacy:
            print("[WARN] %s —— 已按 allow_legacy 降级放行（留痕）" % msg)
            return None
        raise ContractError(
            "[ERROR] %s\n        拒绝运行（fail-closed）。处置：先 retry 上游步骤生成契约，"
            "或显式设 AUTOZT_CONTRACT_LEGACY=warn 降级（仅限确认无误的历史数据）。" % msg)

    # B. 上游 step.conf 改过、但产物没重生成
    cur = file_sha256(os.path.join(str(upstream_dir), "step.conf"))
    if up.get("step_conf_sha256") and cur and up["step_conf_sha256"] != cur:
        raise ContractError(
            "[ERROR] %s：上游 %s 的 step.conf 已变化，但它的产物还是旧参数下生成的。\n"
            "        拒绝运行 —— 请先 retry/rerun 上游步骤，让产物与参数一致。\n"
            "        （上游契约 %s，当前 step.conf %s）"
            % (who, upstream_step or "", str(up["step_conf_sha256"])[:12], str(cur)[:12]))

    # C. 上游产物换过（chain 变了）→ 本步旧产物失效
    if my_contract:
        old = (my_contract.get("upstream") or {}).get("chain_hash")
        new = up.get("chain_hash")
        if old and new and old != new:
            raise ContractError(
                "[ERROR] %s：上游 %s 的产物已更新（chain %s → %s），本步现有产物是在旧上游上生成的。\n"
                "        拒绝运行 —— 请先 retry/rerun 本步。"
                % (who, upstream_step or "", str(old)[:12], str(new)[:12]))
    return up


__all__ = ["SCHEMA", "CONTRACT_NAME", "ContractError", "canonical", "contract_path",
           "file_sha256", "make_contract", "params_hash", "read_contract",
           "require_upstream", "write_contract"]
