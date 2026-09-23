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
NCORE  = 4
KPAR   = 1

AMIN    = 0.01          # 2D 长真空层必加，否则电子步收敛极慢（与 S1/S4 的 2D 模板一致；review 六）
