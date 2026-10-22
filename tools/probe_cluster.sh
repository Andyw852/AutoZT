#!/usr/bin/env bash
# =============================================================================
# probe_cluster.sh -- read-only site survey for AutoZT site-profile authoring.
# Collects facts about a compute site without submitting work (an optional
# --test-only mode performs scheduler acceptance checks that reserve no resources).
#
# Usage:
#   tools/probe_cluster.sh                          # current machine
#   tools/probe_cluster.sh <ssh-alias>              # remote site
#   tools/probe_cluster.sh <ssh-alias> --test-only  # also test per-partition acceptance
#   tools/probe_cluster.sh <ssh-alias> --test-only --max-partitions 8
#
# The script assumes nothing about the presence of SLURM, a module system or
# accelerators; unavailable tools are reported and skipped.
# =============================================================================
set -u
set -o pipefail

TARGET="${1:-local}"
[ $# -gt 0 ] && shift
TEST_ONLY=0
MAX_PART=24
while [ $# -gt 0 ]; do
  case "$1" in
    --test-only) TEST_ONLY=1 ;;
    --max-partitions) MAX_PART="${2:-24}"; shift ;;
    -h|--help) sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

if [ "$TARGET" = "local" ]; then
  SSH=""
else
  SSH="ssh -o BatchMode=yes -o ConnectTimeout=15"
  if ! $SSH "$TARGET" true 2>/dev/null; then
    echo "[ERROR] cannot connect to '$TARGET' in BatchMode; configure ~/.ssh/config and keys." >&2
    exit 1
  fi
fi

RUN() {
  # Execute through 'bash -s' on stdin: independent of the login shell and executed once.
  if [ -z "$SSH" ]; then printf '%s\n' "$1" | LC_ALL=C bash -s; else printf '%s\n' "$1" | $SSH "$TARGET" 'LC_ALL=C bash -s'; fi
}

echo "======================================================================"
echo " AutoZT site survey -- target=$TARGET   $(date '+%F %T')"
echo "======================================================================"

echo
echo "## 0. Host information"
RUN '
echo "hostname     : $(hostname 2>/dev/null)"
echo "user         : $(whoami 2>/dev/null)"
echo "login shell  : $(getent passwd $(whoami) 2>/dev/null | cut -d: -f7)"
echo "os           : $(grep PRETTY_NAME /etc/os-release 2>/dev/null | cut -d= -f2- | tr -d \"\\\"\")"
echo "kernel       : $(uname -r 2>/dev/null)  arch: $(uname -m 2>/dev/null)"
echo "cpus(login)  : $(nproc 2>/dev/null)"
echo "mem(login)   : $(grep MemTotal /proc/meminfo 2>/dev/null)"
echo "home         : $HOME"
'

echo
echo "## 1. Batch scheduler"
RUN '
for c in sbatch squeue sinfo sacct scontrol scancel; do command -v $c >/dev/null 2>&1 && echo "slurm  : $c -> $(command -v $c)"; done
command -v qsub >/dev/null 2>&1 && echo "PBS/SGE: qsub -> $(command -v qsub)"
command -v bsub >/dev/null 2>&1 && echo "LSF    : bsub -> $(command -v bsub)"
command -v sinfo >/dev/null 2>&1 && echo "slurm version: $(sinfo -V 2>/dev/null | head -1)"
if command -v sbatch >/dev/null 2>&1; then
  echo "SLURM available (AutoZT submits directly)"
elif command -v qsub >/dev/null 2>&1; then
  echo "non-SLURM scheduler (submit templates or an adapter layer required; see guide, Section 11.5)"
else
  echo "no batch scheduler: AutoZT local mode or tools/fakeslurm applies"
fi
'

echo
echo "## 2. SLURM partitions (available|partition|limit|nodes|cores/node|memory MB|GRES)"
RUN '
if command -v sinfo >/dev/null 2>&1; then
  sinfo -a -h -o "%P|%a|%l|%D|%c|%m|%G" 2>/dev/null | awk -F"|" "
    { key=\$1; cnt[key]+=\$4; cpu[key]=\$5; mem[key]=\$6; gres[key]=\$7; avail[key]=\$2; tl[key]=\$3 }
    END { for (k in cnt) printf \"%-22s %-4s %-12s nodes=%-5d %s-core %sMB %s\\n\", k, avail[k], tl[k], cnt[k], cpu[k], mem[k], gres[k] }" | sort
  echo
  echo "partition state counts:"
  sinfo -a -h -o "%t" 2>/dev/null | sed "s/[^a-z].*//" | sort | uniq -c | head -8
  echo
  echo "node samples (first 5):"
  sinfo -a -h -o "%N %P %t %c %G" 2>/dev/null | head -5
else
  echo "(no sinfo: not a SLURM site)"
fi
'

echo
echo "## 3. Account associations and limits (as far as readable)"
RUN '
if command -v sacctmgr >/dev/null 2>&1; then
  echo "--- assoc (account/user/partition/QOS/GrpTRES/MaxTRES/MaxJobs) ---"
  sacctmgr -n -P show assoc where user=$USER format=Account,User,Partition,QOS,GrpTRES,MaxTRES,MaxJobs 2>/dev/null | head -5
  echo "--- QOS ---"
  sacctmgr -n -P show qos where name=\"user_$USER\" 2>/dev/null | head -3
fi
if command -v sshare >/dev/null 2>&1; then
  echo "--- sshare ---"; sshare -P -u $USER 2>/dev/null | head -3
fi
echo "note: GrpTRES cpu= is the account concurrency limit; the per-job limit may be probed with --test-only (Section 10)."
'

echo
echo "## 4. Storage"
RUN '
echo "--- df -h ---"
df -h $HOME 2>/dev/null | tail -2
for d in /work /public /scratch /gpfs /lustre; do [ -d "$d" ] && df -h "$d" 2>/dev/null | tail -1; done
echo "--- quota (if available) ---"
( quota -s 2>/dev/null || lfs quota -h $HOME 2>/dev/null ) | head -6
echo "--- home writability ---"
[ -w "$HOME" ] && echo "$HOME is writable" || echo "$HOME is NOT writable"
'

echo
echo "## 5. Module system, compilers and MPI"
RUN '
if command -v module >/dev/null 2>&1; then
  echo "module system: $(module --version 2>&1 | head -2 | tr \"\\n\" \" \")"
  echo "available modules: $(module avail 2>&1 | wc -l) lines"
  echo "--- compiler/MPI/mathematics modules (excerpt) ---"
  module avail 2>&1 | grep -oE \"(compiler|mathlib|mpi)/[^ ]+\" | sort -u | head -25
else
  echo "(no module system: compilers and libraries come from PATH or /opt)"
fi
for c in gcc g++ gfortran icc icpc ifort ifx mpiicc mpicc mpif90 mpirun mpiexec cmake make; do
  p=$(command -v $c 2>/dev/null) && echo "$c -> $p ($($c --version 2>&1 | head -1 | cut -c1-60))"
done
'

echo
echo "## 6. Mathematics and FFT libraries"
RUN '
echo "MKLROOT=${MKLROOT:-(not set)}"
for d in /opt/intel /opt/hpc/software/compiler/intel /public/software/mathlib; do
  [ -d "$d" ] && echo "directory $d: $(ls "$d" 2>/dev/null | head -6 | tr \"\\n\" \" \")"
done
echo "--- dynamic-library visibility ---"
( ldconfig -p 2>/dev/null | grep -E \"libmkl_core|libopenblas|liblapack|libblas|libfftw3|libscalapack|libsymspg\" | head -12 ) || true
'

echo
echo "## 7. Accelerators"
RUN '
if command -v nvidia-smi >/dev/null 2>&1; then echo "NVIDIA: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -4 | tr \"\\n\" \"; \")"; else echo "NVIDIA: nvidia-smi not found"; fi
if [ -x /opt/hyhal/bin/hy-smi ]; then echo "Hygon DCU: /opt/hyhal/bin/hy-smi available"; /opt/hyhal/bin/hy-smi -L 2>/dev/null | head -4; fi
command -v rocm-smi >/dev/null 2>&1 && echo "AMD ROCm: $(command -v rocm-smi)"
[ -e /dev/kfd ] && echo "/dev/kfd present (ROCm/DCU class)"
'

echo
echo "## 8. Python environment (non-interactive PATH is required by AutoZT)"
RUN '
echo "python3 -> $(command -v python3 2>/dev/null || echo not found)"
python3 -V 2>/dev/null
python3 -c "import numpy;print(\"numpy\",numpy.__version__)" 2>/dev/null || echo "numpy: absent"
python3 -c "import scipy;print(\"scipy\",scipy.__version__)" 2>/dev/null || true
python3 -c "import pymatgen;print(\"pymatgen ok\")" 2>/dev/null || true
python3 -c "import spglib;print(\"spglib ok\")" 2>/dev/null || true
python3 -c "import h5py;print(\"h5py ok\")" 2>/dev/null || true
for c in conda mamba micromamba; do command -v $c >/dev/null 2>&1 && echo "$c -> $(command -v $c)"; done
echo "--- first writable directory on PATH (candidate for a python3 symlink) ---"
IFS=\":\"; for d in $PATH; do [ -d "$d" ] && [ -w "$d" ] && { echo "$d"; break; }; done; unset IFS
'

echo
echo "## 9. Existing AutoZT installation"
RUN '
if [ -d "$HOME/software/AutoZT" ]; then
  echo "~ /software/AutoZT exists:"
  du -sh "$HOME"/software/AutoZT/* 2>/dev/null | tail -12
else
  echo "~/software/AutoZT not found (create it according to the project layout convention)"
fi
echo "--- ~/.vaspkit ---"; [ -f "$HOME/.vaspkit" ] && head -4 "$HOME/.vaspkit" || echo "absent"
echo "--- candidate pseudopotential directories (first 3) ---"
find "$HOME/software" -maxdepth 5 -type d -name "potpaw*" 2>/dev/null | head -3
'

if [ "$TEST_ONLY" = "1" ]; then
  echo
  echo "## 10. Submission acceptance check (sbatch --test-only; reserves no resources)"
  RUN "
if ! command -v sbatch >/dev/null 2>&1; then echo '(no SLURM; skipped)'; exit 0; fi
tmp=\$(mktemp -d 2>/dev/null) || tmp=\$HOME
i=0
for p in \$(sinfo -a -h -o '%P' 2>/dev/null | tr -d '*' | sort -u); do
  i=\$((i+1)); [ \$i -gt $MAX_PART ] && break
  f=\$tmp/probe.slurm
  printf '#!/bin/bash\\n#SBATCH -p %s\\n#SBATCH -N 1\\n#SBATCH -n 1\\n#SBATCH -t 00:10:00\\nhostname\\n' \"\$p\" > \"\$f\"
  out=\$(sbatch --test-only \"\$f\" 2>&1 | head -1)
  echo \"P=\$p => \$out\"
  sleep 1
done
rm -f \"\$tmp/probe.slurm\" 2>/dev/null; rmdir \"\$tmp\" 2>/dev/null
"
fi

echo
echo "======================================================================"
echo " Next: fill setting/<site>.yaml (field reference: docs/DEPLOY-NEW-CLUSTER.md, Section 10)"
echo " Required: name / ssh_host / work_dir / partition / max_cpus / vasp.* / potcar_dir / conda_env"
echo " Acceptance: autozt check --cluster <site>, then the smoke calculation."
echo "======================================================================"
