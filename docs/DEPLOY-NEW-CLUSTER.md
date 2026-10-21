# Deploying AutoZT at a New HPC Site: Porting and Acceptance Guide

**Document version:** 1.0  
**Applies to:** AutoZT v2.0 and later  
**Audience:** site administrators and users who deploy AutoZT on a compute cluster or standalone node  
**Companion artifacts:** `tools/probe_cluster.sh` (read-only site survey), `docs/DEPLOY-RECORD-TEMPLATE.md` (per-site record), `docs/CONFIGURING.md` (configuration reference)  
**License:** MIT (see `LICENSE`)

---

## 1. Scope

AutoZT is distributed with templates rather than with a working site configuration. Everything
that identifies a machine -- scheduler, partition names, module system, compilers, mathematics
and FFT libraries, pseudopotential locations, scratch directories -- is site specific and must
be supplied by the operator. This guide describes a reproducible procedure for bringing up a new
site, from a first read-only survey to a documented acceptance test.

The guide deliberately makes no assumptions about the target machine. It covers:

- batch schedulers (SLURM, PBS/SGE/LSF, or none);
- compiler stacks (GNU, Intel, or mixed);
- mathematics libraries (MKL, OpenBLAS, vendor BLAS/LAPACK);
- accelerators (NVIDIA GPUs, Hygon DCUs, or none);
- module systems (Lmod/Tcl modules, Spack, or plain environment variables).

Each section states the acceptance criterion for that step and lists common site-to-site
variations. A worked example with measured numbers is given in Appendix C.

## 2. Terminology

| Term | Meaning in this document |
|---|---|
| Site / cluster | A compute resource accessed over SSH (possibly several accounts on one machine). |
| Site profile | The YAML file `setting/<cluster>.yaml` that describes a site to AutoZT. |
| Skill | An AutoZT workflow module (`skill/<name>/skill.yaml`) defining steps, criteria and artifacts. |
| Runtime environment | The Python environment (virtualenv or conda prefix) used by input-generation scripts. |
| Work directory | A writable file-system root on the site where AutoZT creates material directories. |
| Pseudopotential library | The POTCAR (or UPF) data set used by the electronic-structure engine. |
| Template | A scheduler submission script in `setting/<cluster>/templates/`, rendered at input generation. |
| Acceptance test | The evidence that the site is configured correctly (Sections 14 and 15). |

## 3. Relation to existing practice

AutoZT follows an established pattern in scientific software: the code is portable, while
machine-specific parameters are isolated in a declarative profile and validated by a dedicated
command. Representative examples:

| Project | Site-specific configuration | Validation / diagnostics |
|---|---|---|
| ASE | `siteconfig.py` (shipped template, generated once per site) | `ase info`, import test |
| AiiDA | `verdi computer setup` (hostname, scheduler, work directory, `mpirun` command, prepend/append text) | `verdi computer test` |
| Quantum ESPRESSO | `make.inc` plus `arch/` machine files | build of the test suite |
| CP2K | `arch/` files, one per platform | `make -j` and regression tests |
| LAMMPS | `src/MAKE/MACHINES/Makefile.*` | unit tests in `unittest/` |
| Snakemake | executor plugins and per-site profiles (`profiles/<site>/config.yaml`) | `--dry-run` |
| Nextflow | `nextflow.config` with named `profiles` selected by `-profile` | `-preview` |
| jobflow-remote | `jobflow.yaml` (runner, workers, scheduler, work directory, pre/post-run commands) | `jf runner check` |

The conventions shared by these projects, and adopted here, are:

1. **Ship a template, not a working configuration.** No site paths, host names or credentials
   are committed; `setting/*` is excluded from version control.
2. **One declarative file per site**, enumerating access details, scheduler, work directory,
   runtime environment and launch commands.
3. **Environment discovery over hard coding.** Compilers, MPI and libraries are selected through
   the site's own mechanism (modules, Spack, or environment variables), never assumed.
4. **A preflight command.** Configuration errors are reported before any job is submitted
   (`autozt check`, Section 14).
5. **Reference and tutorial are separate documents.** Field semantics live in
   `docs/CONFIGURING.md` and Section 11 here; the step-by-step procedure is this guide.
6. **Recorded acceptance evidence.** Version pins, build flags and smoke-test output are kept in
   a per-site record (`docs/<SITE>-DEPLOY-<date>.md`), preserving reproducibility.

## 4. Site survey (read-only)

Perform the survey before changing anything on the site. `tools/probe_cluster.sh` collects the
following items over SSH without submitting work (unless `--test-only` is given):

```bash
tools/probe_cluster.sh                       # current machine
tools/probe_cluster.sh <ssh-alias>           # remote site
tools/probe_cluster.sh <ssh-alias> --test-only   # also probe scheduler acceptance
```

The report covers host and OS information; the available batch scheduler; SLURM partitions with
node counts, cores, memory, GRES and time limits; account associations and QOS limits; file
systems and quotas; module system and compiler/MPI versions; mathematics and FFT libraries;
accelerators; the non-interactive Python environment; existing AutoZT installation; and
candidate pseudopotential directories.

### 4.1 Six quantities to determine manually

| Question | Read-only command | Recorded in |
|---|---|---|
| Which partitions accept jobs from this account? | per-partition `sbatch --test-only` | `partition` |
| Cores and memory per node? | `sinfo -o "%c %m"`, `scontrol show node` | `max_cpus` |
| Account-level concurrency limit? | `sacctmgr show assoc` (GrpTRES), or an oversized `--test-only` request | comment in the site profile |
| Maximum job wall time? | `sinfo -o "%l"`, `scontrol show partition` | `--time` in templates |
| Accelerators per node? | `nvidia-smi -L`, `hy-smi`, `sinfo` GRES | skills using GPUs/DCUs |
| Work and home file systems, sizes and quotas? | `df -h`, `quota` | `work_dir` |

When probing many partitions, pause for one to two seconds between requests; some schedulers
throttle rapid submissions (`sbatch: request too frequently`), which would otherwise make
available partitions appear unavailable.

## 5. Access and non-interactive execution

### 5.1 SSH configuration

```sshconfig
Host <alias>
    HostName <host>
    Port <port>
    User <account>
    IdentityFile ~/.ssh/id_rsa_<name>
    IdentitiesOnly yes
    ServerAliveInterval 30
    ServerAliveCountMax 6
    ControlMaster auto
    ControlPath ~/.ssh/cm-%r@%h:%p
    ControlPersist 4h
```

Acceptance criterion: `ssh -o BatchMode=yes <alias> true` exits with status 0. Private keys
must be readable only by the owner (`chmod 600`) and the parent directory must be `chmod 700`.

If the cluster is reachable only through a bastion host, use
`ProxyCommand ssh -W %h:%p <bastion-alias>` and verify each hop independently; a failure in
either hop is otherwise reported only as a timeout.

### 5.2 Requirements of non-interactive sessions

AutoZT invokes the remote host non-interactively (`ssh host 'python3 ...'`) for job collection
and submission. Two conditions therefore apply:

1. **The runtime Python must be on the non-interactive `PATH`.** Variables exported only from
   interactive start-up files are not visible. Two standard remedies are to write the `PATH`
   to `~/.ssh/environment` (requires `PermitUserEnvironment yes` in `sshd_config`), or to
   create a symbolic link to the interpreter in a directory that is already on the
   non-interactive `PATH` and is writable by the account; the survey script reports the first
   such directory.
2. **Module commands are not inherited.** Every job script must load its modules explicitly; a
   non-interactive SSH command does not run the login-shell module initialisation.

Acceptance criterion: `ssh <alias> 'python3 -V && command -v python3'` succeeds and reports
the intended runtime environment.

## 6. Directory layout and storage

The following layout is the project convention:

```
~/software/AutoZT/
├── vasp/                    # source trees and binaries, one directory per version
│   ├── vasp.6.4.3-pristine/
│   ├── vasp.6.4.3-optcell/  # only if cell-constrained relaxation is required
│   └── makefile.include     # site-specific build configuration, kept for provenance
├── vasp_pseudopotentials/   # e.g. potpaw_PBE_64/
├── venvs/                   # one Python environment per component
├── conda_envs/              # conda prefixes, if used
├── mlff_models/             # machine-learning interatomic potentials
├── qe-*/ pseudo-*/          # optional components
└── work/                    # work_dir for AutoZT
```

Sizing guidance: a VASP build tree occupies about 2 GB (sources, objects and six binaries for
the pristine and cell-constrained variants); a pseudopotential library 0.3 GB; Python
environments 0.6-2 GB each; machine-learning potentials 0.1-1 GB; optional components up to
3 GB. Check quotas (`df -h`, `quota`) before installing.

Acceptance criterion: `du -sh ~/software/AutoZT/*` remains within the available quota.

## 7. Toolchain selection

### 7.1 Compiler stack

| Available on the site | Recommended choice |
|---|---|
| System GCC only (often 4.8-7.x) | Load a newer module (`gcc/9` or later) or a developer toolset; older compilers may fail with an internal compiler error |
| Intel oneAPI / Parallel Studio | `mpiicc`/`mpiifort` with MKL; least configuration effort |
| Mixed (GNU compilers with Intel MKL) | Supported; link the MKL variants explicitly (`gnu_thread`) |

ABI rule: the MPI library determines the compiler and Fortran runtime to use. Mixing
`libgfortran` versions commonly produces binaries that link but fail at run time.

### 7.2 Mathematics and FFT libraries

| Source | Advantage | Caution |
|---|---|---|
| MKL | Fast and self-contained | With GNU compilers choose `intel_lp64` + `intel_thread` or `gnu_thread`; `libiomp5` resides in the compiler library directory |
| System BLAS/LAPACK | No installation | Link as `-l:liblapack.so.3 -l:libblas.so.3`; LAPACK before 3.8 provides `dgegv_`, later versions removed it (Section 8) |
| OpenBLAS | Good multithreaded performance | Account for thread oversubscription against MPI ranks |
| FFTW3 | Portable | Use the double-precision installation |

Acceptance criterion: `ldd <binary> | grep "not found"` produces no output, and a smoke run
completes without undefined-symbol errors.

## 8. Building the electronic-structure engine

The procedure below is written for VASP and applies unchanged to other source-distributed
engines.

```bash
module load <compiler> <mpi>            # also required inside job scripts
cd ~/software/AutoZT/vasp
# 1) unpack the source tree
# 2) write makefile.include for this site and keep it with the tree
# 3) start with a single target, then increase parallelism
make -k -j 16 std > build_std.log 2>&1 || true
```

### 8.1 Diagnostic table for build failures

| Symptom in the build log | Probable cause | Remedy |
|---|---|---|
| `internal compiler error: Segmentation fault` | Compiler too old for the source | Load a newer compiler |
| `cannot find -lxxx` | Wrong library path or name | Locate the library (`ldconfig -p`, `module show`); add `-L` and preferably `-Wl,-rpath` |
| `undefined reference to 'name_'` | Missing library, wrong link order, or ABI mismatch | Place libraries last; ensure MPI and Fortran runtime match the compiler |
| Many missing `.mod` files under parallel make | Module-dependency race in parallel builds | Repeat `make -k -j N` until convergence; each pass completes another group of objects |
| `invalid command 'bdist_wheel'` (Python packages) | `wheel` absent from the environment | `pip install wheel`, then retry |

After installation, verify dynamically: `ldd` for unresolved libraries, then a minimal
calculation (Section 8.3). A successful link is not sufficient evidence of a working binary.

### 8.2 Cell-constrained relaxation (optional)

Materials requiring in-plane relaxation with a fixed vacuum axis need a variant built from a
patched source tree. Keep both variants side by side and reference them through two VASP
profiles in the site profile (`standard`, `relax_2d`).

### 8.3 Smoke test and numerical baseline

Use a two-atom Si diamond cell with a minimal input set (`ENCUT=250`, `ISMEAR=0`, `SIGMA=0.05`,
`EDIFF=1e-5`, `NELM=20`) submitted through the site template. Acceptance criteria: exit status
0; `reached required accuracy` present in `OUTCAR`; a timing block written.

Record the total energy and compare it with an established value obtained on another machine
with identical inputs. A discrepancy larger than 1e-4 eV warrants checking the pseudopotential
file, the plane-wave cutoff and the engine version.

### 8.4 Parallel configuration

- `NCORE` (VASP 5 and later) is usually chosen between four and eight, or approximately the
  square root of the cores per node; use `KPAR` only when enough k-points are available.
- Whole-node requests are more efficient but queue longer than partial-node requests on busy
  systems; measure both policies (Section 12.2).
- The account concurrency limit bounds the product of job size and the number of simultaneous
  jobs; reserve one slot for debugging.

## 9. Pseudopotentials, post-processing and the Python runtime

### 9.1 Pseudopotential library layouts

Three layouts occur in practice; `potcar_dir` must point to the correct level:

```
A) potcar_dir/potpaw_PBE_64/<element>/POTCAR     (recommended; expected by defect-type skills)
B) potcar_dir/potpaw_PBE/<element>/POTCAR
C) potcar_dir/<element>/POTCAR                   (result of extracting the official archive)
```

The official archives extract as layout C (element directories directly beneath the current
directory). Create the parent directory, move the content into it, and make the read-only
directories writable first (`chmod -R u+w`), because the archives preserve restrictive modes.

Acceptance criterion: `autozt check` reports no missing pseudopotential data, and the
post-processing tool generates a `POTCAR` whose first line names the expected potential.

### 9.2 VASPKIT configuration

VASPKIT 1.x reads `~/.vaspkit` in a whitespace-separated `KEY VALUE` format. Writing
`KEY = VALUE` causes the value to be parsed as `=` and the pseudopotential is silently not
generated. Verify with a single structure:

```bash
printf '1\n103\n' | vaspkit      # should create POTCAR starting with 'PAW_PBE Si ...'
```

A related diagnostic issue exists in early releases: `autozt check` parsed only `KEY = VALUE`
lines and therefore reported `vaspkit_potcar` as missing although VASPKIT itself worked. Both
formats are now accepted; the native whitespace-separated format remains required by VASPKIT.

### 9.3 Python runtime environments

| Requirement | Recommended environment |
|---|---|
| Core AutoZT input generation | A virtualenv with numpy, scipy, matplotlib, spglib, pymatgen, h5py and seekpath |
| Component-specific packages (amset, phonon codes, MACE, ...) | One isolated environment per component; version requirements frequently conflict |

Version policy:

- `numpy`: legacy scientific packages (amset 0.5.x, pheasy) generally require 1.x, whereas
  `cogito-dft` 0.3.13 requires 2.x; such packages must not share an environment.
- `pymatgen`: many older tools are incompatible with recent major versions; verify through the
  component's own entry point, not only with `import pymatgen`.
- After installation, always test through a non-interactive invocation (Section 5.2).

Acceptance criterion:
`ssh <alias> 'python3 -c "import numpy,scipy,pymatgen,spglib,h5py; print(1)"'`.

#### 9.3.1 Environments created with `conda create -p`

AutoZT treats a runtime-environment path containing `/` as a virtualenv and activates it with
`source <env>/bin/activate`. Conda environments created by recent conda versions do not contain
a virtualenv-style `bin/activate`; activation therefore fails silently and the system
interpreter is used instead, which manifests as missing packages in `autozt check` and import
failures in generation scripts.

Remedy: supply a minimal activation shim:

```bash
cat > <env>/bin/activate <<'EOF'
_AZT_ENV="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
export VIRTUAL_ENV="$_AZT_ENV" CONDA_PREFIX="$_AZT_ENV"
export PATH="$_AZT_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$_AZT_ENV/lib:${LD_LIBRARY_PATH:-}"
EOF
chmod +x <env>/bin/activate
```

Verification: `bash -c 'source <env>/bin/activate; python -V; python -c "import <package>"'`
must use the interpreter of that environment.
## 10. Cluster profile reference (`setting/<cluster>.yaml`)

The site profile is the single source of truth for how AutoZT operates at a site. All keys are
optional unless marked required; missing keys fall back to package defaults or to the project
configuration, in this order: `project_setting/` > user configuration > package defaults
(see `docs/CONFIGURING.md`).

### 10.1 Field reference

| Key | Required | Meaning |
|---|---|---|
| `name` | yes | Cluster identifier, used by `setting/<name>.yaml`, `-hpc` and `check --cluster`. |
| `ssh_host` | yes | SSH alias from `~/.ssh/config`; an empty string selects local execution. |
| `template_map` | no | Logical-template-to-file mapping; normally `{}` because templates are resolved by name under `setting/<name>/templates/`. |
| `vasp.<profile>` | yes (DFT skills) | VASP profiles: `standard` (`cell_constraint: none`) and optionally `relax_2d` (`cell_constraint: ioptcell_tag`), each with absolute paths to `std`, `gam`, `ncl`. |
| `work_dir` | yes | Writable root on the site; must differ from other sites. |
| `partition` | yes | Default submission partition. |
| `max_cpus` | recommended | Cores per node; used by the preflight check to detect templates requesting more cores than a node provides. |
| `conda_sh` | no | Path to `profile.d/conda.sh`; leave empty for virtualenv-only sites. |
| `conda_env` | yes (generation) | Python environment for input-generation scripts; a path containing `/` is activated as a virtualenv. |
| `pheasy_env`, `kl_env`, `amset_env` | per skill | Component-specific environments; each is validated separately by `autozt check`. |
| `potcar_dir` | yes (DFT skills) | Pseudopotential root; see Section 9.1 for the three layouts. |
| `mlff_model_dir` | MLFF skills | Directory holding machine-learning potential weights. |
| `remote_path_prefix` | recommended | `PATH` prefix for remote commands (for example a `pybin` directory with tool links, followed by the environment's `bin`). |
| `qe_bin`, `qe_pseudo_dir`, `qe_setup` | eph-qe skill | Quantum ESPRESSO executables, UPF directory and pre-run environment statements. |
| Additional keys | as reported | `autozt check` names anything still missing; add keys as skills require them. |

### 10.2 Example profiles

```yaml
# (a) CPU cluster with SLURM and modules
name: examplecpu
ssh_host: examplecpu
work_dir: /work/<account>/software/AutoZT/work
partition: cpu64
max_cpus: 64
vasp:
  standard:
    cell_constraint: none
    std: /work/<account>/software/AutoZT/vasp/vasp.6.4.3-pristine/bin/vasp_std
    gam: .../vasp_gam
    ncl: .../vasp_ncl
conda_sh: ""
conda_env: /work/<account>/software/AutoZT/venvs/autozt
potcar_dir: /work/<account>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: /work/<account>/software/AutoZT/pybin:/work/<account>/software/AutoZT/venvs/autozt/bin
```

```yaml
# (b) GPU cluster: templates carry --gres=gpu:N
name: examplegpu
ssh_host: examplegpu
work_dir: /home/<account>/AutoZT/work
partition: gpu
max_cpus: 128
conda_env: /home/<account>/AutoZT/venvs/mace_gpu
mlff_model_dir: /home/<account>/AutoZT/mlff_models
```

```yaml
# (c) Scheduler-less node: submission is emulated (tools/fakeslurm)
name: exampleserver
ssh_host: exampleserver
work_dir: /home/<account>/AutoZT/work
conda_env: /home/<account>/AutoZT/venvs/autozt
```

```yaml
# (d) Local execution
name: local
ssh_host: ""
work_dir: /tmp/autozt_work
conda_env: ""
```

## 11. Submission templates and resource policy

### 11.1 Templates and rendering

Templates live in `setting/<cluster>/templates/` (for example `submit_std_3d.tpl`,
`submit_ncl_2d.tpl`, `step1_opt/submit_std_2d.tpl`). At input generation, AutoZT substitutes
the VASP executable path according to the active VASP profile, but it does **not** interpret
`module load` lines or scheduler directives. When copying templates from another site, verify
at least: partition name, module names, MPI launcher and flags, and wall-time limit.

Acceptance criterion: after `autozt -t <skill> -p <material> -j <step> init` (generation
without submission), inspect the generated `submit.sh` on the site and confirm partition,
cores and executable paths before executing `start`.

### 11.2 Resource sizing

1. Determine the partition limits (`sinfo`) and the account concurrency limit (GrpTRES).
2. Keep `job size x concurrent jobs` within the account limit; reserve one slot for debugging.
3. For multi-node jobs specify `-N` and `-n` consistently with `--ntasks-per-node`.
4. Set the wall-time limit below the partition maximum, otherwise submission is rejected.

Operating policy for partial-node versus whole-node requests: on busy systems, requesting whole
nodes forces the scheduler to wait for fully idle nodes and increases queueing time. A
configuration with more, smaller jobs is often more efficient in aggregate; the optimal policy
depends on the site and should be measured (Section 12.2).

### 11.3 MPI launch

| Observation | Probable cause | Remedy |
|---|---|---|
| Job fails within seconds; stderr contains `binding more processes than cpus` | OpenMPI default binding incompatible with the node or cgroup layout | `mpirun --bind-to none -np $SLURM_NTASKS ...`, after confirming on a single node |
| Intel MPI program launched with OpenMPI `mpirun` | ABI mismatch | Use the `mpiexec`/`mpirun` of Intel MPI, or load the matching `mpi/intelmpi` module |
| Several full copies of the program start | Template launches one full calculation per process | Ensure a single `mpirun -np <tasks> <binary>` invocation |
| Single-node communication issues | InfiniBand/UCX initialisation | For single-node jobs set `OMPI_MCA_btl=self,vader,tcp` and `OMPI_MCA_pml=ob1` |

### 11.4 Accelerators

- NVIDIA GPUs: request with `#SBATCH --gres=gpu:N`, install a CUDA-enabled runtime, and verify
  visibility with `nvidia-smi` inside the job.
- Hygon DCUs: request with `#SBATCH --gres=dcu:N`, use the DTK toolchain (`compiler/dtk/...`),
  and install the DTK-compatible PyTorch build rather than the upstream CUDA wheel.
- Acceptance criterion: the accelerator is listed inside the job and the runtime reports it as
  available (`torch.cuda.is_available()` or the vendor equivalent).

### 11.5 Sites without SLURM

AutoZT manages jobs through `sbatch`, `squeue` and `scancel` semantics. For a scheduler-less
server, the repository provides `tools/fakeslurm` (used automatically in local mode and
injected through `remote_path_prefix` on remote nodes). For PBS/SGE/LSF, a translation layer
from `sbatch` semantics to the local scheduler is required; until one exists, automatic job
progression is not supported and this limitation must be recorded in the site record.

## 12. Multi-cluster operation and project-level overrides

### 12.1 Step-level routing

A material may run different steps on different sites. This is expressed in the project
configuration (`tf_<project>.yaml`):

```yaml
step_hpc:
  S7_deform: examplegpu
  S7.1_read: examplegpu
```

Only the listed steps are routed; all others remain on the material default. Artifacts required
by downstream steps are transferred by AutoZT according to step dependencies. Switching the
whole material (`autozt -p <material> hpc <cluster>`) affects every step and may therefore
require re-running earlier steps; this decision should be recorded in the site record.

### 12.2 Project-level step configuration

Parameters can be overridden at project level, in increasing priority:
`project_setting/templates/step.conf` (shared) and
`project_setting/templates/<step>/step.conf` (step-specific). A single parameter can be changed
with

```bash
autozt -t <skill> -p <material> -j <step> conf --set params.<KEY>=<VALUE>   # writes project config
autozt -t <skill> -p <material> -j <step> conf                              # read-only, shows merged values
```

Writing project configuration is a controlled operation in AutoZT and must be requested from
the operator of the material.

### 12.3 Operational discipline

Job submission commands are not idempotent with respect to one another. `retry`, `init` and
`rerun` cancel queued jobs of the same material (`kill_if_queued`). If two operators or two
sessions act on the same material concurrently, the later action will cancel the earlier
jobs, so results are no longer attributable. The following rules avoid this:

1. One operator at a time per material; inspect the task graph (`autozt graph`) before acting.
2. Isolate site profiles between parallel worktrees. `setting/*` is excluded from version
   control, so each git worktree has its own copy of the site profiles.
3. Record the owner of a run in the site or project record.

## 13. Acceptance criteria

The site is considered operational when every item below passes. The evidence should be pasted
into the per-site record.

| # | Check | Command | Expected result |
|---|---|---|---|
| 1 | Connectivity | `ssh -o BatchMode=yes <alias> true` | exit status 0 |
| 2 | Non-interactive Python | `ssh <alias> 'python3 -V'` | intended version printed |
| 3 | Configuration preflight | `autozt doctor` | no errors; the new site is listed |
| 4 | Skill readiness | `autozt check --cluster <site>` | DFT-related skills pass |
| 5 | Template rendering | `autozt -p <material> -j <step> init` then inspect `submit.sh` | partition, cores and paths correct |
| 6 | Smoke calculation | minimal Si job (Section 8.3) | `reached required accuracy`; energy matches the baseline |
| 7 | Concurrency | submit two probe jobs | both accepted; no QOS or partition rejection |
| 8 | Clean-up | `squeue -u <account>` | no test jobs remain |
| 9 | Record | `docs/<SITE>-DEPLOY-<date>.md` | quotas, issues and acceptance numbers documented |

Interpreting `autozt check` output:

- **Hard failure**: the skill cannot run (missing engine binary, pseudopotential library, core
  Python package, or the environment specific to that skill).
- **Optional failure**: a secondary capability is unavailable (`optional exe:`/`optional py:`,
  GPU absent). The main workflow is unaffected.
- A missing package is reported **for the environment the skill uses** (suffixes such as
  `@amset_env`, `@kl_env`). Confirm the package is installed in that environment, not only in
  `conda_env`.

## 14. Optional components

Components are installed on demand. Each subsection gives the installation route, the
acceptance test and known issues.

### 14.1 Quantum ESPRESSO (skill `eph-qe-cpu`)

1. Prefer a site module (`module avail | grep -i quantum-espresso`). Point `qe_bin` to its
   `bin` directory and express the runtime environment in `qe_setup` (module loads, or explicit
   `LD_LIBRARY_PATH` for MKL and Intel MPI).
2. Otherwise use a conda-forge package (self-contained; requires network access and sufficient
   memory for the conda solver).
3. Acceptance: `pw.x -h` and `ph.x -h` run; a small Si SCF calculation reaches
   `convergence has been achieved`; norm-conserving pseudopotentials are present in
   `qe_pseudo_dir`.
4. Known issues: preloaded MPI modules may conflict with the QE module (unload the former
   first); MKL/Intel MPI runtime libraries sometimes need to be exported explicitly.

### 14.2 AMSET (skill `ke-dft-cpu`)

Use a dedicated environment: AMSET 0.5.x typically requires numpy 1.x and an older pymatgen,
which conflicts with machine-learning environments. A pip-based installation is:

```bash
python -m venv <env>
<env>/bin/pip install 'amset==0.5.1' 'pymatgen<2025' 'numpy<2' scipy matplotlib h5py numba interpolation sumo pyFFTW
```

Known issue: `BoltzTraP2` frequently fails to build from its PyPI source distribution (the
bundled spglib submodule lacks a CMakeLists file, and netCDF4 headers are required). A
conda-forge binary package can be unpacked and its module copied into the environment; the
interpolation feature is otherwise not needed by the transport calculation.
Acceptance: `import amset`, `amset deform create --help`, and the AMSET environment listed in
the site profile as `amset_env`.

### 14.3 MACE (MLFF skills)

Install a CPU or CUDA PyTorch build, then `mace-torch`, `e3nn`, `ase` and `matscipy`. Note the
PyTorch/e3nn compatibility: PyTorch 2.6 changed the default of `torch.load` to
`weights_only=True`, which breaks older e3nn releases with `UnpicklingError`; either use
PyTorch 2.4.x or a newer e3nn. The default model name expected by the skills is listed in
`skill/*/templates/step.conf` (`MACE_MODEL`); download the corresponding weights to
`mlff_model_dir` and verify with a single-point calculation.

### 14.4 COGITO (skill `cohp-cogito`)

Install `cogito-dft` (for example version 0.3.13) into a dedicated virtualenv, because it
requires numpy 2.x. Expose `COGITO`, `COGITOanalyze` and `COGITOpost` to the environment used
by the skill (symbolic links are sufficient). Prerequisites: a static VASP run with
`LWAVE=.TRUE.`, a symmetry-reduced k-point grid (`ISYM` not `-1`) and a sufficient number of
bands.

### 14.5 UniHamGNN (skill `unihamgnn`)

The component comprises the HamGNN source tree, the universal model weights (about 1 GB), the
OpenMX basis-set library `DFT_DATA19`, the helper executables `openmx_postprocess` and
`read_openmx`, and a Python environment with PyTorch (CPU), `torch_geometric`, `lightning` and
`e3nn`.

Building `openmx_postprocess` requires GSL, MPI, LAPACK and FFTW. Missing ScaLAPACK/BLACS
symbols can be satisfied by the MKL ScaLAPACK and BLACS libraries; the `libiomp5` runtime is
located in the Intel compiler library directory; and statically linking GSL avoids a runtime
dependency on `libcblas`. ELPA Fortran modules require several passes of `make -k -j N` because
of parallel module-dependency races. The skill parameters (`HAMGNN_DIR`, `UNI_MODEL`,
`DFT_DATA`, `OPENMX_POSTPROCESS`, `READ_OPENMX`) are project-level step parameters and belong in
`project_setting/templates/step.conf`, not in the site profile.

### 14.6 Phonon and force-constant fitting (phonopy, phono3py, ALM, pheasy)

- `phonopy` and `phono3py`: install binary wheels (`pip install --only-binary=:all: phonopy
  phono3py`); verify that the `phono3py` executable exists and both modules import.
- `ALM`: the PyPI package of this name is an empty 2.4 KB placeholder and must not be installed.
  conda-forge provides builds only for Python 3.6-3.8. For Python 3.10/3.11, extract the
  conda-forge Python 3.11 package, copy its `site-packages/alm/` into the environment, and
  provide the required `libsymspg.so.1` by symbolic link to the `libsymspg.so.2` shipped with
  the environment's `spglib` (the extension resolves libraries relative to `$ORIGIN/../../..`).
- `pheasy`: install from source (`pip install wheel` first, then
  `pip install --no-deps --no-build-isolation .`); verify that the `pheasy` executable runs.

## 15. Troubleshooting matrix

| Symptom | Probable cause | Remedy |
|---|---|---|
| `Invalid account or account/partition combination` | Account not authorised for the partition | Identify an authorised partition with `--test-only`; request access |
| `CPU count per node can not be satisfied` | Requested cores exceed those of one node (e.g. a template with a hard-coded 192) | Adjust the per-node core count in the template or step configuration; keep `max_cpus` consistent |
| `group max tres limit ...` | Account concurrency limit exceeded | Reduce the number of simultaneous jobs or their size; document the change |
| Job fails within seconds with `binding more processes than cpus` | OpenMPI default binding | `mpirun --bind-to none` |
| `error while loading shared libraries` | Library not available at run time | Use `ldd` to identify the library; rebuild with `-Wl,-rpath`, or load the corresponding module |
| `undefined symbol: dgegv_` | LAPACK 3.8 and later removed the routine | Add `dgegv.f` and `dggsvd.f` from netlib to the sources |
| `I REFUSE TO CONTINUE` at start-up | Inconsistent structure or corrupted POSCAR | Verify that the transferred POSCAR matches the local file (AutoZT performs a checksum) |
| VASPKIT does not create POTCAR | `~/.vaspkit` uses `KEY = VALUE` | Use `KEY VALUE`; verify the `PBE_PATH` root |
| `autozt check` reports missing packages although they are installed | Non-interactive `PATH`; a skill-specific environment; a conda prefix without `bin/activate` | Sections 5.2, 13 and 9.3.1 |
| `conda create` killed (exit 137) or extremely slow | Solver memory or index access | Use micromamba or `--solver=libmamba`; unpack conda packages directly; use a local mirror |
| `sbatch: request too frequently` | Submission throttling | Pause one to two seconds between submissions; submit in batches |
| Remote collection reports `Connection closed` or a timeout | Intermediate host, proxy, or ControlMaster failure | Re-establish the connection (`ssh -O exit` and reconnect); check key expiry on each hop |
| Jobs submitted by two sessions cancel each other | Both sessions operate on the same material (`kill_if_queued`) | One operator per material; isolate `setting/` in worktrees; inspect the task graph before acting |

## 16. Site record and reproducibility

For every site, create `docs/<SITE>-DEPLOY-<YYYYMMDD>.md` from
`docs/DEPLOY-RECORD-TEMPLATE.md`. The record should contain the measured resource limits, the
software versions and build flags, the acceptance output, and any remaining limitations. This
is the document that makes the site reproducible by another group.

Reproducibility requirements:

1. Pin versions of the engine, libraries, Python packages and model weights.
2. Keep site configuration out of version control (`setting/*` is ignored); templates and
   placeholders are distributed instead.
3. Cite the software version: see `README.md` and `CITATION.cff` for the version-specific DOI.
4. When reporting results, name the site and the site-record document alongside the input
   files.

---

## Appendix A: Minimal site profile template

```yaml
name: <SITE>
ssh_host: <ssh-alias>
template_map: {}
vasp:
  standard:
    cell_constraint: none
    std: <root>/software/AutoZT/vasp/<version>/bin/vasp_std
    gam: <root>/software/AutoZT/vasp/<version>/bin/vasp_gam
    ncl: <root>/software/AutoZT/vasp/<version>/bin/vasp_ncl
work_dir: <root>/software/AutoZT/work
partition: <partition>
max_cpus: <cores-per-node>
conda_sh: ""
conda_env: <root>/software/AutoZT/venvs/autozt
potcar_dir: <root>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: <root>/software/AutoZT/pybin:<root>/software/AutoZT/venvs/autozt/bin
# add for ke-dft-cpu:      amset_env: ...
# add for eph-qe-cpu:      qe_bin / qe_pseudo_dir / qe_setup
# add for MLFF skills:     mlff_model_dir: ...
```

## Appendix B: Procedure on one page

1. Survey the site: `tools/probe_cluster.sh <alias>`; record the six quantities of Section 4.1.
2. Configure SSH and the non-interactive runtime (Section 5).
3. Create `~/software/AutoZT/{vasp,vasp_pseudopotentials,venvs,work}` (Section 6).
4. Build the engine (pristine and, if required, the cell-constrained variant); verify with
   `ldd` (Sections 7-8).
5. Install the pseudopotential library, configure VASPKIT and create the Python environment
   (Section 9).
6. Write the site profile and templates (Sections 10-11).
7. Run `autozt doctor`, `autozt check --cluster <site>`, then the smoke calculation
   (Section 13).
8. Create the site record (Section 16).
9. Run one real material through a short chain to measure throughput and queueing.
10. Install optional components as required, re-running `autozt check` after each
    (Section 14).

## Appendix C: Worked example (measured values)

The following values were measured on SCNet Xi'an `clus_xa` (2026-10-19); the detailed record is
in `docs/CLUS_XA-DEPLOY-20261019.md`. They are given as an order-of-magnitude reference, not as
a configuration to copy.

| Item | Measured value | Consequence |
|---|---|---|
| Node specification | 64 cores, 255 GB | determines `max_cpus` and template size |
| Partition wall time | 3 days | `--time=3-00:00:00` |
| Account concurrency | 320 cores (also per-job) | e.g. 32-core jobs x 10 concurrent |
| Submission template | 32 cores per job | whole-node jobs queue for fully idle nodes |
| MPI issue | OpenMPI default binding failed on some nodes | `mpirun --bind-to none` |
| Compiler | GCC 9.3.0 (system 7.3.1 fails with an internal compiler error) | module load |
| MPI library | HPC-X 2.11 | matched to the compiler |
| BLAS/LAPACK | system 3.4.2 with `dgegv`/`dggsvd` sources | later LAPACK lacks the routines |
| Engine | VASP 6.4.3, pristine and IOPTCELL-patched | `standard` and `relax_2d` profiles |
| Pseudopotentials | potpaw_PBE_64, 346 elements | extraction requires mode adjustment |
| VASPKIT | 1.2.4, whitespace-separated `~/.vaspkit` | POTCAR generation |
| Quantum ESPRESSO | site module 6.8 (Intel MPI), MKL runtime, `Si_ONCV_PBE-1.2.upf` | SCF converged in 6 iterations |
| AMSET | 0.5.1 in an isolated environment (pymatgen 2024.11, numpy 1.26) | `ke-dft-cpu` ready |
| COGITO | cogito-dft 0.3.13 in an isolated environment (numpy 2.x) | `cohp-cogito` ready |
| MACE | mace-torch 0.3.16, PyTorch 2.4.1 CPU, e3nn 0.4.4; `MACE-matpes-pbe-omat-ft.model` | `mlff`/`opt-mlff-cpu` ready |
| UniHamGNN | HamGNN 2.3.0, 885 MB model, `DFT_DATA19`, self-built OpenMX tools | `unihamgnn` ready |

## Appendix D: References

- ASE, customising the installation (`siteconfig.py`): https://ase-lib.org/ase/installation.html
- AiiDA, setting up and testing a computer: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html
- Quantum ESPRESSO user guide (installation, `make.inc`, `arch/`): https://www.quantum-espresso.org/Doc/user_guide/
- CP2K, compiling with `arch/` files: https://manual.cp2k.org/
- LAMMPS, machine makefiles: https://docs.lammps.org/Install.html
- Snakemake executor plugins and profiles: https://snakemake.github.io/snakemake-plugin-catalog/
- Nextflow configuration and profiles: https://www.nextflow.io/docs/latest/config.html
- jobflow-remote configuration: https://matgenix.github.io/jobflow-remote/
- AutoZT repository and version-specific DOI: see `README.md` and `CITATION.cff`.
