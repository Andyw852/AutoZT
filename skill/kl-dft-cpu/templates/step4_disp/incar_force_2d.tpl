# phono3py 位移超胞单点取力（2D）。声子硬约束：ISYM=0（位移破缺对称，勿对称化力）。
# 精度（phonopy 官方超胞取力示例的口径）：PREC=Accurate、EDIFF=1E-8、LREAL=.FALSE.。
#   位移帧的力直接变成 fc2/fc3：LREAL=Auto 的实空间投影带格点噪声，是 fc 里假软模的
#   常见来源，也会污染 ShengBTE/fourphonon 的散射矩阵。三个值走 step.conf 的
#   FORCE_PREC / FORCE_EDIFF / FORCE_LREAL（旧版注释写的是这套、正文却是
#   Normal/Auto/1E-7，已对齐）。
#   超胞很大（几百原子）时可用 FORCE_LREAL=Auto 换速度，但必须先在 3 帧上与 .FALSE.
#   比力的最大偏差（< 1 meV/Å 才采用），并显式设 ROPT（如 ROPT = 1E-3）逐步收紧。
#   ADDGRID 降低力的格点噪声，代价是略微变慢。
SYSTEM  = {{SYSTEM}}
ISTART  = 0
ICHARG  = 2
GGA     = {{GGA}}
{{VDW_LINE}}

PREC    = {{FORCE_PREC}}
ENCUT   = {{ENCUT}}
EDIFF   = {{FORCE_EDIFF}}
LREAL   = {{FORCE_LREAL}}
ADDGRID = .TRUE.          # 降低力的格点噪声（可选，phonopy 取力推荐）
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
# P2-3：偶极修正由 gen 从 S1/S2 的 INCAR 继承（三步必须同一静电边界条件）
{{DIPOLE_LINE}}

LWAVE   = .FALSE.
LCHARG  = .FALSE.
LORBIT  = 0
NCORE   = 4
KPAR    = 2
