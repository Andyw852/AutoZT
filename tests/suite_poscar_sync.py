# -*- coding: utf-8 -*-
"""自测：v3 本地模式下"材料根 POSCAR 推送"与 contcar_to_poscar 的两个坑（2026-09-16）。

全程本地、不碰集群：直接用远端命令原样在临时目录里跑 bash，断言文件效应。

背景（MoS2 实测）：
  1. remote_gen 里 POSCAR 原来只在"远端缺这个文件"时才推 → 本地按结构复用流程改了
     材料根 POSCAR 后（TASKFLOW 结构复用），远端旧结构一直在用，第一步 VASP 拿旧结构跑。
     现在改成"内容（sha256）不一致就推，推完立刻校验"。
  2. contcar_to_poscar 原来用 [ -f CONTCAR ] 判断 → 崩掉的作业留下的 0 字节 CONTCAR
     会把刚 gen 好的 POSCAR 清成空文件，之后每次 retry/submit 都必然再崩。
     现在要求 CONTCAR 非空且 >= 8 行、首行非空。
"""
import base64
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from autozt.workflow import _contcar_to_poscar_line  # noqa: E402

FAILS = []
BASES = []


def ok(cond, msg):
    if cond:
        print("  ✓ %s" % msg)
    else:
        FAILS.append(msg)
        print("  ✗ %s" % msg)


def tmpdir(prefix):
    d = tempfile.mkdtemp(prefix="tf_%s_" % prefix)
    BASES.append(d)
    return d


def run(cmd, cwd):
    return subprocess.run(["bash", "-c", cmd], cwd=cwd,
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace")


POS = ("Mo S\n1.0\n  3.0 0.0 0.0\n -1.5 2.6 0.0\n 0.0 0.0 20.0\n"
       "Mo S\n1 2\nDirect\n0 0 0\n")


def push_cmd(data):
    """与 remote_gen 里 POSCAR 推送那段等价的 shell（同一套写法）。"""
    b64 = base64.b64encode(data).decode()
    psha = hashlib.sha256(data).hexdigest()
    return ("_tfp=$(sha256sum POSCAR 2>/dev/null | cut -d' ' -f1); "
            "[ \"$_tfp\" = %s ] || { _t=.tf_poscar.$$; "
            "echo %s | base64 -d > \"$_t\" && mv -f \"$_t\" POSCAR && "
            "echo \"[tf] POSCAR 已更新为本地版 %s（远端原为 $_tfp）\"; }; "
            "[ \"$(sha256sum POSCAR 2>/dev/null | cut -d' ' -f1)\" = %s ] || "
            "{ echo 'ERROR: POSCAR 推送后 sha256 校验失败' >&2; exit 1; }; "
            % (psha, b64, psha[:12], psha))


def test_contcar_guard():
    print("[1] contcar_to_poscar：空/半截 CONTCAR 不得覆盖 POSCAR")
    d = tmpdir("b2")
    with open(os.path.join(d, "POSCAR"), "w") as fh:
        fh.write(POS)
    with open(os.path.join(d, "CONTCAR"), "w"):     # 崩掉的作业留下的 0 字节占位
        pass
    run(_contcar_to_poscar_line(d), d)
    ok(open(os.path.join(d, "POSCAR")).read() == POS,
       "空 CONTCAR：POSCAR 保持原样（旧写法在这里会把它清成 0 字节）")

    with open(os.path.join(d, "CONTCAR"), "w") as fh:
        fh.write(POS.replace("3.0", "3.1"))
    run(_contcar_to_poscar_line(d), d)
    ok("3.1" in open(os.path.join(d, "POSCAR")).read(),
       "有效 CONTCAR：仍然覆盖 POSCAR（续跑语义没被破坏）")

    d2 = tmpdir("b2b")
    with open(os.path.join(d2, "POSCAR"), "w") as fh:
        fh.write(POS)
    run(_contcar_to_poscar_line(d2), d2)
    ok(open(os.path.join(d2, "POSCAR")).read() == POS, "没有 CONTCAR：POSCAR 不动")

    d3 = tmpdir("b2c")
    with open(os.path.join(d3, "POSCAR"), "w") as fh:
        fh.write(POS)
    with open(os.path.join(d3, "CONTCAR"), "w") as fh:
        fh.write("junk\n1.0\n3.0 0 0\n")            # 半截文件
    run(_contcar_to_poscar_line(d3), d3)
    ok(open(os.path.join(d3, "POSCAR")).read() == POS, "半截 CONTCAR（3 行）：POSCAR 不动")


def test_poscar_push():
    print("[2] 材料根 POSCAR：按内容差异推送 + 推后校验")
    src = open(os.path.join(ROOT, "autozt", "workflow.py")).read()
    ok("sha256sum POSCAR" in src and "[ -f POSCAR ] ||" not in src,
       "remote_gen 已改成按 sha256 差异推送（旧的「缺文件才推」写法已消失）")

    new, old = b"Mo S\nNEW-STRUCTURE\n", b"Mo S\nOLD-STRUCTURE\n"
    d = tmpdir("b1")
    r = run(push_cmd(new), d)
    ok(r.returncode == 0 and open(os.path.join(d, "POSCAR"), "rb").read() == new,
       "远端缺文件：推送成功（旧写法唯一能过的场景）")

    with open(os.path.join(d, "POSCAR"), "wb") as fh:
        fh.write(old)
    r = run(push_cmd(new), d)
    ok(r.returncode == 0 and open(os.path.join(d, "POSCAR"), "rb").read() == new
       and "已更新" in r.stdout,
       "远端是旧内容：覆盖为新内容（MoS2 那次漏掉的正是这条）")

    r = run(push_cmd(new), d)
    ok(r.returncode == 0 and "已更新" not in r.stdout, "远端已一致：跳过，不白写")

    d5 = tmpdir("b1b")
    with open(os.path.join(d5, "POSCAR"), "wb") as fh:
        fh.write(new)
    _ro = False
    try:
        os.chmod(d5, 0o500)                          # 目录只读 → 推送必然失败
        _ro = not os.access(os.path.join(d5, "POSCAR"), os.W_OK)
    finally:
        os.chmod(d5, 0o700)
    if not _ro:
        print("  ~ 跳过：目录只读未生效（root 或 WSL/NTFS 上 chmod 无效），无法测推送失败分支")
    else:
        os.chmod(d5, 0o500)
        r = run(push_cmd(old), d5)
        os.chmod(d5, 0o700)
        ok(r.returncode != 0 and "校验失败" in (r.stdout + r.stderr),
           "推送失败：gen 直接报错退出，不会带着旧结构继续跑")


def main():
    print("== suite_poscar_sync：POSCAR 推送 / contcar_to_poscar 自测 ==")
    test_contcar_guard()
    test_poscar_push()
    for d in BASES:
        shutil.rmtree(d, ignore_errors=True)
    print()
    if FAILS:
        print("FAIL %d 项：" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
