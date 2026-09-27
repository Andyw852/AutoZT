# step2_nac 3D：LEPSILON DFPT 介电 + Born 有效电荷（原胞；给 S3 声子谱的非解析项修正）。
# 关键：LEPSILON=.TRUE. 出 MACROSCOPIC STATIC DIELECTRIC TENSOR + Born；LPEAD 提升数值稳定；
#      NPAR/NCORE 与 LEPSILON 不兼容，gen 会移除它们；KPAR 保留。
SYSTEM = {{SYSTEM}}
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
KPAR   = 2
