#!/bin/bash
# =====================================================================
# apply_new_dielectric.sh v2 —— 用新 DFPT 结果更新 AMSET settings.yaml
# 用法: bash apply_new_dielectric.sh <材料> <新DFPT目录> [--submit]
#
# v2 加固（2026-09-15）：
#   - POP 提取加 || true，避免 set -e 在 amset 失败时中止整个脚本
#   - 更新前**强制**检查 eps 数值有限（NaN/Inf 一律拒绝写入）
#   - 更新后用 yaml.safe_load 校验两份 settings.yaml
#   - 提交前校验迪特
# =====================================================================
set -e
M="$1"; NEWD="$2"; DO_SUBMIT="$3"
if [ -z "$M" ] || [ -z "$NEWD" ]; then echo "用法: $0 <材料> <DFPT目录> [--submit]"; exit 1; fi
[ -f "$NEWD/OUTCAR" ] || { echo "[ERROR] $NEWD/OUTCAR 不存在"; exit 1; }

source /public/home/.../miniconda3/etc/profile.d/conda.sh
conda activate amset051 >/dev/null 2>&1

echo "===== $M ====="
# ---- 0) 先用数值校验器把关新 DFPT 产物（NaN/单位矩阵/求和规则等）----
set +e
python3 /tmp/vd2.py --step-dir "$NEWD" --json /tmp/dielchk_$M.json > /tmp/dielchk_$M.txt 2>&1
RC=$?
set -e
echo "  [校验] validate_dielectric 退出码=$RC"
sed "s/^/    /" /tmp/dielchk_$M.txt | head -6
if [ $RC -ne 0 ]; then echo "  [中止] 新 DFPT 产物未通过数值校验，拒绝更新"; exit 2; fi

# ---- 1) 提取介电 ----
python3 /tmp/extract_diel.py "$NEWD" > /tmp/diel_$M.json
python3 - "$M" <<'PY'
import json, sys, math
M = sys.argv[1]
d = json.load(open("/tmp/diel_%s.json" % M))
bad = []
for name in ("eps_inf", "eps_static"):
    for r in d[name]:
        for v in r:
            if not math.isfinite(v): bad.append((name, v))
print("  NKPTS      :", d["nkpts"])
print("  eps_inf    :", [round(v, 4) for v in d["eps_inf_diag"]])
print("  eps_static :", [round(v, 4) for v in d["eps_static_diag"]])
for name in ("eps_inf_diag", "eps_static_diag"):
    if any(v < 1.0 for v in d[name]):
        bad.append((name, d[name]))
if bad:
    print("  [FATAL] 介电数值异常:", bad)
    sys.exit(7)
PY

# ---- 2) 重算 pop_frequency（失败不致命）----
POP=""
if [ -f "$NEWD/vasprun.xml" ]; then
  POP=$( (cd "$NEWD" && amset phonon-frequency -o OUTCAR -v vasprun.xml 2>/dev/null) \
         | grep -oE "pop_frequency[:= ]+[0-9.]+" | grep -oE "[0-9.]+$" | tail -1 ) || true
fi
if [ -z "$POP" ]; then
  echo "  [WARN] pop_frequency 未算出——保留原值"
  POP="-"
fi
echo "  pop_frequency: $POP"

# ---- 3) 备份 + 更新两套 settings.yaml ----
UPDATED=""
for S in /public/home/.../ke_work/$M/ke-dft-cpu/step8_amset/settings.yaml \
         /public/home/.../ke_soc_20260910/$M/ke-dft-cpu/step8_amset/settings.yaml; do
  if [ ! -f "$S" ]; then echo "  [跳过] $S 不存在"; continue; fi
  [ -f "$S.bak-9x9x1" ] || cp "$S" "$S.bak-9x9x1"
  python3 /tmp/update_amset_settings.py "$S" /tmp/diel_$M.json "$POP"
  UPDATED="$UPDATED $S"
done

# ---- 4) 校验更新后的 YAML ----
python3 /tmp/val_yaml.py $UPDATED

# ---- 5) 提交 ----
if [ "$DO_SUBMIT" = "--submit" ]; then
  for D in /public/home/.../ke_work/$M/ke-dft-cpu/step8_amset \
           /public/home/.../ke_soc_20260910/$M/ke-dft-cpu/step8_amset; do
    [ -f "$D/submit.sh" ] || continue
    ( cd "$D"
      if [ -f transport.json ]; then mv transport.json "transport.json.9x9x1-$(date +%m%d)"; fi
      jid=$(sbatch submit.sh 2>&1 | grep -oE "[0-9]+" | head -1)
      echo "  已提交 $D -> jobid $jid" )
  done
fi

