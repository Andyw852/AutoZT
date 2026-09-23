#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""COMSOL 后端 —— 接口已留，未启用（需要 COMSOL 安装 + 授权）。

为什么是"留接口"而不是直接实现
------------------------------
COMSOL 是商业软件：安装包与 license 都在 COMSOL Access 账号后面，license 与
host ID 绑定，无法由本仓库自行获取或分发。因此这里实现同样的后端契约，并在
缺少 COMSOL 时给出明确报错与启用步骤，而不是静默失败。

后端契约（与 backend_skfem 完全一致）
------------------------------------
    solve(model, workdir=None, opts=None) -> dict
    dict 至少包含：
        backend, dT_peak_K, T_top_max_K,
        P_source_W_per_m, P_nominal_W_per_m, sink,
        mesh, y_nodes_m, T_mid_x_K
    其中 y_nodes_m/T_mid_x_K 是中轴剖面，必须与 device-thermal 的
    T_profile_mid_K 同坐标，S3 才能逐点比对。

启用步骤（Linux 无头为例）
--------------------------
1. 装 COMSOL：安装包放到目标机，装到 ~/software/AutoZT/comsol（见铁律 11），
   装完确认可执行文件存在：
       ls $COMSOL_ROOT/bin/comsol $COMSOL_ROOT/bin/comsolbatch
2. 配 license（三选一）：
       export LMCOMSOL_LICENSE_FILE=port@license_server     # 许可证服务器
       export LMCOMSOL_LICENSE_FILE=/path/to/license.dat    # 单机文件
       或把 license.dat 放到 $COMSOL_ROOT/license/
   用离线批处理自检：
       $COMSOL_ROOT/bin/comsol batch -help      # 能打印帮助即 license OK
3. 本技能 step.conf：
       FEM_BACKEND         = comsol
       COMSOL_BIN          = ~/software/AutoZT/comsol/bin/comsol
       COMSOL_LICENSE_FILE = ~/software/AutoZT/comsol/license/license.dat
       COMSOL_CORES        = 4
4. Python 侧（可选，二选一）：
   a) mph（推荐）：pip install mph
      - 用 MPh 按 model 建几何/材料/边界/源，直接拿 T 场；
   b) 纯 batch：本后端写一个 .java 或 .mph + 参数扫描脚本，交给
      comsol batch -input ... -output ... 执行，再解析导出的 CSV。

实现要点（接上时照此填 _solve_via_mph / _solve_via_batch）
----------------------------------------------------------
* 几何：矩形 [0,L] x [0,H]；两层材料（各向异性 kx/ky）+ 界面薄层（厚 delta，
  导热 delta*G），与 skfem 后端的物理严格一致；
* 边界：y=H 恒温 T_amb；x=0 与 x=L 依 sink 取恒温或绝热；
* 源：沟道层体热源 Q/t_ch，保证总功率 = Q*L（不要按单元重复计数！）；
* 网格：按 model["mesh"] 指定，输出中轴剖面供 S3 比对。
"""

REQUIRED = ("COMSOL_BIN", "FEM_BACKEND=comsol")


def available():
    """检查 COMSOL 是否可用；返回 (ok, 说明)。"""
    import os
    import shutil
    binp = os.environ.get("COMSOL_BIN") or ""
    if binp and os.path.exists(os.path.expanduser(binp)):
        return True, "COMSOL_BIN=%s" % binp
    found = shutil.which("comsol") or shutil.which("comsolbatch")
    if found:
        return True, "PATH 中找到 %s" % found
    return False, ("未找到 COMSOL：请设 COMSOL_BIN，或把 comsol/comsolbatch 放进 PATH；"
                   "并确认 license 可用")


def _not_ready():
    ok, why = available()
    return ("[ERROR] FEM_BACKEND=comsol 尚未就绪：%s\n"
            "        COMSOL 是商业软件，本仓库无法自带安装包与 license。\n"
            "        启用步骤见 backend_comsol.py 顶部文档（装软件 -> 配 license -> "
            "设 COMSOL_BIN/COMSOL_LICENSE_FILE -> 选 mph 或 batch 路径）。\n"
            "        在此之前请用 FEM_BACKEND=skfem（开源、pip 可装），二者物理一致。"
            % why)


def solve(model, workdir=None, opts=None):
    raise SystemExit(_not_ready())


def _solve_via_mph(model, workdir, opts):   # pragma: no cover - 待授权后实现
    raise SystemExit(_not_ready())


def _solve_via_batch(model, workdir, opts):  # pragma: no cover - 待授权后实现
    raise SystemExit(_not_ready())
