#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step_contract.py（步骤参数契约闸）回归。

用户要求："这是一个保护措施，前后步骤参数对不上就不要运行"。
本套件把该机制的**全部判据**钉成硬闸（纯本地、不碰集群）：
  A. 上游契约缺失 → fail-closed 拒绝；allow_legacy=True 才降级放行；
  B. 上游 step.conf 改过、产物没重生成 → 拒绝；
  C. 本步记录的上游 chain_hash != 上游现在的 chain_hash → 拒绝；
  + 契约构造/落盘/回读/损坏容错、chain_hash 的决定因素（params/inputs/upstream）。
"""
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skill" / "_common"))
import step_contract as sc  # noqa: E402

FAIL = []


def ck(cond, msg):
    if cond:
        print("  PASS  %s" % msg)
    else:
        print("  FAIL  %s" % msg)
        FAIL.append(msg)


def ck_raises(fn, exc, msg):
    try:
        fn()
    except exc as e:
        ck(True, msg + "（%s）" % type(e).__name__)
        return e
    except Exception as e:                                  # noqa: BLE001
        ck(False, msg + "：抛了 %s 而不是 %s（%s）" % (type(e).__name__, exc.__name__, e))
        return None
    ck(False, msg + "：没有抛异常")
    return None


ROOT.joinpath("tmp").mkdir(exist_ok=True)
D = Path(tempfile.mkdtemp(dir=str(ROOT / "tmp"), prefix="stepcontract-"))
(D / "step.conf").write_text("ALM_CUT3 = 6.0\n", encoding="utf-8")
IN1 = D / "in1.txt"
IN1.write_text("hello", encoding="utf-8")

print("== A. canonical / params_hash / file_sha256 ==")
ck(sc.canonical({"b": 1, "a": 2}) == sc.canonical({"a": 2, "b": 1}), "canonical 键序无关")
ck(sc.canonical({"a": 1}) != sc.canonical({"a": 2}), "canonical 区分值")
ck(sc.params_hash({"A": 1, "B": 2}) == sc.params_hash({"B": 2, "A": 1}), "params_hash 键序无关")
ck(sc.params_hash({"A": 1}) != sc.params_hash({"A": 2}), "params_hash 值变则变")
ck(sc.params_hash(None) == sc.params_hash({}), "params_hash(None) == params_hash({})")
ck(len(sc.params_hash({})) == 64, "params_hash 是 sha256 十六进制")
ck(sc.file_sha256(IN1) == hashlib.sha256(b"hello").hexdigest(), "file_sha256 正确")
ck(sc.file_sha256(D / "nope.txt") is None, "file_sha256 读不到返回 None（不抛）")

print("== B. make_contract 的字段与 chain_hash 决定因素 ==")
args = dict(params={"A": 1}, inputs={"poscar": str(IN1), "known": "deadbeef"})
c = sc.make_contract("S4_disp", D, **args)
ck(c["schema"] == sc.SCHEMA == 1, "schema 正确")
ck(c["step"] == "S4_disp", "step 名写入")
ck(len(c["step_conf_sha256"]) == 64 and c["step_conf_sha256"] == sc.file_sha256(D / "step.conf"),
   "记录本步 step.conf 的 sha256")
ck(c["inputs"]["poscar"] == hashlib.sha256(b"hello").hexdigest(), "路径型输入自动取哈希")
ck(c["inputs"]["known"] == "deadbeef", "非路径输入按已给哈希原样保留")
ck(c["upstream"] is None and c["params"] == {"A": 1}, "upstream=None / params 留档")
ck(sc.make_contract("S4_disp", D, **args)["chain_hash"] == c["chain_hash"], "同 params/inputs/upstream → chain 相同")
ck(sc.make_contract("S9_other", D, **args)["chain_hash"] == c["chain_hash"], "chain 不含步骤名（只由 params/inputs/upstream 决定）")
ck(sc.make_contract("S4_disp", D, params={"A": 2},
                    inputs=args["inputs"])["chain_hash"] != c["chain_hash"], "params 变 → chain 变")
ck(sc.make_contract("S4_disp", D, params={"A": 1},
                    inputs={"poscar": str(IN1), "known": "cafe"})["chain_hash"] != c["chain_hash"],
   "inputs 变 → chain 变")
c_up = sc.make_contract("S4_disp", D, upstream={"step": "S3_nac", "chain_hash": "UP1"}, **args)
ck(c_up["upstream"] == {"step": "S3_nac", "chain_hash": "UP1"}, "upstream 归一化为 step+chain_hash")
ck(c_up["chain_hash"] != c["chain_hash"], "upstream 变 → chain 变")

print("== C. write_contract / read_contract ==")
w = sc.write_contract(D, "S4_disp", **args)
ck(os.path.isfile(sc.contract_path(D)), "契约落盘到 <step_dir>/.step_contract.json")
ck(not os.path.exists(sc.contract_path(D) + ".tmp"), "原子替换，无 .tmp 残留")
r = sc.read_contract(D)
ck(r["chain_hash"] == w["chain_hash"] and r["step"] == "S4_disp", "回读一致")
open(sc.contract_path(D), "w", encoding="utf-8").write("{ this is not json")
ck(sc.read_contract(D) is None, "损坏契约按缺失处理（不抛）")

print("== D. require_upstream 判据 A/B/C ==")
ck(issubclass(sc.ContractError, SystemExit), "ContractError 是 SystemExit（gen 非零退出 → autozt 标 error）")
up_missing = D / "up_missing"
up_missing.mkdir()
e = ck_raises(lambda: sc.require_upstream("S5_fc", D, up_missing, "S4_disp"),
              sc.ContractError, "A 上游无契约 → 拒绝运行（fail-closed）")
ck(bool(getattr(e, "code", None)), "A 拒绝时带非空 code")
ck(sc.require_upstream("S5_fc", D, up_missing, "S4_disp", allow_legacy=True) is None,
   "A allow_legacy=True → 降级放行并返回 None")

up = D / "up"
up.mkdir()
(up / "step.conf").write_text("U = 1\n", encoding="utf-8")
sc.write_contract(up, "S4_disp", params={"U": 1}, inputs={})
got = sc.require_upstream("S5_fc", D, up, "S4_disp")
ck(got is not None and got["step"] == "S4_disp", "B/C 未改 → 放行并返回上游契约")

(up / "step.conf").write_text("U = 2\n", encoding="utf-8")
ck_raises(lambda: sc.require_upstream("S5_fc", D, up, "S4_disp"),
          sc.ContractError, "B 上游 step.conf 改了但产物没重生成 → 拒绝")

up_before = sc.write_contract(up, "S4_disp", params={"U": 2}, inputs={})
mine = {"upstream": {"step": "S4_disp", "chain_hash": up_before["chain_hash"]}}
ck(sc.require_upstream("S5_fc", D, up, "S4_disp", my_contract=mine) is not None,
   "C 上游 chain 未变 → 放行")
up_after = sc.write_contract(up, "S4_disp", params={"U": 3}, inputs={})
ck(up_after["chain_hash"] != up_before["chain_hash"], "上游 params 变 → 上游 chain 变")
ck_raises(lambda: sc.require_upstream("S5_fc", D, up, "S4_disp", my_contract=mine),
          sc.ContractError, "C 上游产物换过、本步旧产物失效 → 拒绝")
ck(sc.require_upstream("S5_fc", D, up, "S4_disp", my_contract=None) is not None,
   "C 本步没有旧契约（首次运行）→ 放行")

no_conf = D / "up_noconf"
no_conf.mkdir()
sc.write_contract(no_conf, "S4_disp", params={}, inputs={})
ck(sc.require_upstream("S5_fc", D, no_conf, "S4_disp") is not None,
   "上游契约缺 step_conf_sha256（如无 step.conf）时不做 B 比对，向后兼容")

print("")
if FAIL:
    print("suite_step_contract: FAILED %d 项" % len(FAIL))
    sys.exit(1)
print("suite_step_contract: ALL PASS")
