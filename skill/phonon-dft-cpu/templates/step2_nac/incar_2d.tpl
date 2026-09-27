# step2_nac 2D：LEPSILON DFPT 介电 + Born 有效电荷（原胞）。
# ★ 2D 长真空层必加 AMIN=0.01，否则电子步收敛极慢。
# ★ 注意：phonopy 的 NAC 是 3D（Wang/Gonze）方案，对真 2D 是近似（LO-TO 在 q->0 应趋零）；
#   要严格 2D 处理请把本项目 nac 置 false，或对照 band_noNAC.yaml。
SYSTEM = {{SYSTEM}} (2D)
ISTART = 0
ICHARG = 2
GGA    = {{GGA}}
{{VDW_LINE}}

PREC   = Accurate
ENCUT  = {{ENCUT}}
EDIFF  = 1E-7
LREAL  = .FALSE.
LASPH  = .TRUE.
ALGO   = Normal
NELM   = 200
ISMEAR = 0
SIGMA  = 0.05

IBRION = -1
NSW    = 0
ISIF   = 2

LWAVE  = .TRUE.
LCHARG = .TRUE.
LORBIT = 11
KPAR   = 1

AMIN   = 0.01
