# =====================================================================
# incar_uniform_full_2d.tpl —— AMSET **全网格**波函数用的密网格自洽（2d）
#
# 【为什么单独一份模板】ISYM = -1。
#   0.4.19 / 0.5.1 的 WavefunctionOverlapCalculator.from_coefficients 里：
#       if np.prod(mesh_dim) == len(kpoints):  return cls.from_data(...)
#   即 **h5 里的 k 点若本身就构成完整网格，就完全不调用 desymmetrize_coefficients**。
#   ISYM = 2 时 VASP 按对称性把 47x47x3 约化成 416 个不可约点 -> AMSET 必须去对称化，
#   而这条路径在实测中可疑（见 step8.4_amset2d/VERIFICATION.md V22）。
#   ISYM = -1 让 VASP 把全部 47x47x3 = 6627 个 k 点都算出来 -> AMSET 直接走 from_data。
#   **必须 -1，不能 0**：ISYM = 0 仍然使用时间反演，只给半个网格，
#   AMSET 照样去对称化，而且走的正是实测失败率最高的时间反演分支。
#
# 【成本】k 点 x16；WAVECAR 约 x16；amsat wave 产出的 h5 约 85 MB -> 1.3 GB 量级。
# 【用途】诊断/裁决专用步骤，默认关闭（optional_steps.wavefunction_full.default = false），
#   由项目 project_setting 里写 wavefunction_full: true 打开；不影响原 S3/S4。
# =====================================================================

SYSTEM = {{SYSTEM}}

ISTART = 0
ICHARG = 2
GGA    = {{GGA}}
{{VDW_LINE}}

PREC   = Accurate
ENCUT  = {{ENCUT}}
LREAL  = .FALSE.
LASPH  = .TRUE.

ALGO   = Normal
EDIFF  = 1E-7            # AMSET 要求高精度波函数
NELM   = 200
NELMIN = 6
AMIN   = 0.01          # 2D 长真空层电子步稳定
ISMEAR = 0              # AMSET 一律高斯小展宽，别用 -5
SIGMA  = 0.01

IBRION = -1
NSW    = 0
ISYM   = -1          # ★ 全网格：禁止 VASP 按对称性约化 k 点

# ---- AMSET 关键：输出致密 DOS 与波函数 ----
LWAVE  = .TRUE.         # amset wave 读 WAVECAR
LCHARG = .TRUE.          # ★ step3c_uniform_offgrid 要读这份 CHGCAR 做非自洽
LORBIT = 11
NEDOS  = 5000
LOPTICS = .FALSE.       # amset 自己算跃迁，不需要 VASP 的 LOPTICS

NCORE  = 6
KPAR   = 2
