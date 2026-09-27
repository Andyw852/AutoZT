# phonopy 位移超胞单点取力（2D）。声子硬约束：ISYM=0（位移破缺对称，勿对称化力）。
# 精度（对齐 kl-dft-cpu S4，phonopy 官方超胞取力示例口径）：PREC=Accurate、EDIFF=1E-8、
# LREAL=.FALSE.。位移帧的力直接变成 fc2：LREAL=Auto 的实空间投影带格点噪声，是 fc 里
# 假软模的常见来源。三个值走 step.conf 的 FORCE_PREC / FORCE_EDIFF / FORCE_LREAL。
# 超胞很大（几百原子）时可用 FORCE_LREAL=Auto 换速度，但必须先在 3 帧上与 .FALSE.
# 比力（最大偏差 < 1 meV/Å）并显式设 ROPT。ADDGRID 降低力的格点噪声。
SYSTEM  = {{SYSTEM}}
ISTART  = 0
ICHARG  = 2
GGA     = {{GGA}}
{{VDW_LINE}}

PREC    = {{FORCE_PREC}}
ENCUT   = {{ENCUT}}
EDIFF   = {{FORCE_EDIFF}}
LREAL   = {{FORCE_LREAL}}
ADDGRID = .TRUE.
LASPH   = .TRUE.
ALGO    = Normal
NELM    = 200
NELMIN  = 6
AMIN    = 0.01
ISMEAR  = 0
SIGMA   = 0.05

IBRION  = -1
NSW     = 0
ISYM    = 0

LWAVE   = .FALSE.
LCHARG  = .FALSE.
LORBIT  = 0
NCORE   = 4
KPAR    = 2
