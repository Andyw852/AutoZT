#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""[0040] 本征载流子点（doping=0）的独立 AMSET 运行 —— 与主掺杂扫描分成两条线。

主 AMSET（掺杂扫描）已经跑完并 cp 成 transport.json 后，本脚本改 settings 重跑一次：
  doping = [0]（本征点，E_F 由 n = p 决定，amset_fermi_fix 求解）
  scattering_type 去掉 IMP（本征点无电离杂质，二维 IMP 屏蔽在 doping=0 会 nan）
  write_mesh = false（复用主运行的网格，不重写 mesh）
然后跑 Runner，把新产出的 transport_<mesh>.json 复制成 transport_undoped.json。

失败不影响主结果（本脚本是作业链尾部、尽力而为；失败只打印告警）。
"""
import glob
import os
import shutil
import sys


def main():
    try:
        import yaml
    except Exception as e:  # noqa: BLE001
        print("[WARN] run_undoped 缺 yaml：%s" % e, file=sys.stderr)
        return 1
    try:
        with open("settings.yaml", encoding="utf-8") as f:
            s = yaml.safe_load(f) or {}
        s["doping"] = [0.0]
        s["scattering_type"] = [m for m in s.get("scattering_type", [])
                                if str(m).upper() != "IMP"]
        s["write_mesh"] = False
        with open("settings.yaml", "w", encoding="utf-8") as f:
            yaml.safe_dump(s, f, allow_unicode=True)
    except Exception as e:  # noqa: BLE001
        print("[WARN] run_undoped 改 settings 失败：%s" % e, file=sys.stderr)
        return 1
    try:
        import amset2d_plugin  # noqa: F401
        import amset_fermi_fix  # noqa: F401
        import amset_ir_fix  # noqa: F401
        import amset_desym_fix  # noqa: F401
        from amset.log import initialize_amset_logger as L
        L()
        from amset.core.run import Runner
        Runner.from_directory(".").run()
    except Exception as e:  # noqa: BLE001
        print("[WARN] run_undoped AMSET 失败：%s" % e, file=sys.stderr)
        return 1
    tf = sorted(glob.glob("transport_*.json"), key=os.path.getmtime, reverse=True)
    if not tf:
        print("[WARN] run_undoped：没找到 transport_*.json", file=sys.stderr)
        return 1
    shutil.copy(tf[0], "transport_undoped.json")
    print("[OK] 本征点 AMSET 完成 -> transport_undoped.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
