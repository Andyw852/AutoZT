# <SITE> -- AutoZT Deployment and Acceptance Record (template)

Copy this file to `docs/<SITE>-DEPLOY-<YYYYMMDD>.md` and fill in the measured values.
Remove every `<placeholder>`. Site profiles (`setting/*`) are excluded from version control;
replace account names and absolute paths by `<account>` and `<root>` when publishing the record.

**Document version:** 1.0  
**Site profile:** `setting/<site>.yaml`  
**Related documents:** `docs/DEPLOY-NEW-CLUSTER.md` (general procedure), `docs/CONFIGURING.md`

---

## 1. Summary

| Item | Value |
|---|---|
| Site / organisation | <name> |
| Cluster identifier / scheduler | <id / SLURM version / PBS / none> |
| Account | <account> (key expiry <YYYY-MM-DD>) |
| SSH alias | <alias> (port <port>; direct / via <bastion>) |
| Login node | <hostname>, <n> cores / <memory> GB |
| Date / operator | <YYYY-MM-DD> / <name> |

## 2. Measured resources

### 2.1 Partition availability

| Partition | Accepts jobs | Nodes | Cores/node | Memory/node | GRES | Max wall time |
|---|---|---|---|---|---|---|
| <partition> | yes/no | <n> | <c> | <m> GB | <gpu:N / none> | <time> |

Method: per-partition `sbatch --test-only`. Record the rejection message for unavailable
partitions (for example `Invalid account or account/partition combination`).

### 2.2 Account limits (measured, not quoted)

| Item | Measured value | Method |
|---|---|---|
| Concurrency limit (CPU cores) | <N> | oversized request reports `group max tres limit ...` |
| Per-job limit | <N cores / M nodes> | `-N`/`-n` stepped with `--test-only` |
| Job-count limit | <none / N> | submit N single-core probe jobs |
| Maximum wall time | <X days> | partition limit, confirmed by `--time` |
| Storage quota | <work X TB / home Y GB> | `df -h`, `quota` |

### 2.3 Scheduling behaviour

- Whole-node jobs: <whether fully idle nodes are required; observed queueing>
- n simultaneous jobs: <how many started; reported pending reasons>
- Operating policy: <e.g. 32-core jobs, 10 concurrent>

## 3. Software and paths

| Component | Version | Path on site | Notes |
|---|---|---|---|
| Compiler | <GCC 9.3.0 / Intel 2021> | module <name> | |
| MPI | <HPC-X 2.11 / Intel MPI 2017> | <path> | ABI matched to compiler |
| BLAS/LAPACK | <MKL / system 3.4.2> | <path> | `dgegv`/`dggsvd` added if required |
| FFT | <FFTW3 double> | <path> | |
| Engine, pristine | <VASP 6.4.3> | <root>/.../bin | `vasp_std`, `vasp_gam`, `vasp_ncl` |
| Engine, cell-constrained | <version + IOPTCELL> | <root>/.../bin | `relax_2d` profile |
| Pseudopotentials | <potpaw_PBE_64, N elements> | <root>/vasp_pseudopotentials | layout A/B/C |
| VASPKIT | <1.2.4> | <path> | `~/.vaspkit` uses whitespace-separated keys |
| Python environment | <virtualenv/conda, version> | <path> | key package versions |
| Optional: Quantum ESPRESSO | <6.8 module> | <path> | |
| Optional: AMSET | <0.5.1, isolated env> | <path> | |
| Optional: MACE | <PyTorch x.y, mace-torch z> | <path> | model files |
| Optional: other | <COGITO / UniHamGNN> | <path> | |

## 4. Deployment issues encountered (site-specific only)

| # | Observation (key log line) | Cause | Remedy | Verification |
|---|---|---|---|---|
| 1 | <error line> | <cause> | <change> | <how it was confirmed> |
| 2 | | | | |

Common issues are listed in `DEPLOY-NEW-CLUSTER.md`, Section 15; record here only conditions
specific to this site (partition names, node models, module conflicts).

## 5. Site profile (redacted)

```yaml
name: <SITE>
ssh_host: <alias>
work_dir: <root>/software/AutoZT/work
partition: <partition>
max_cpus: <cores-per-node>
conda_env: <root>/software/AutoZT/venvs/autozt
potcar_dir: <root>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: <root>/software/AutoZT/pybin:<root>/software/AutoZT/venvs/autozt/bin
```

- Templates: `setting/<site>/templates/` (cores <N>, wall time <...>, MPI launcher <...>).
- Differences from other sites: <e.g. requires `--bind-to none`, specific modules, bastion host>.

## 6. Acceptance evidence

| # | Check | Result |
|---|---|---|
| 1 | `ssh -o BatchMode=yes <alias> true` | pass |
| 2 | non-interactive Python: `ssh <alias> 'python3 -V'` | pass, <version> |
| 3 | `autozt doctor` | pass |
| 4 | `autozt check --cluster <SITE>` | pass, <k>/<total> skills |
| 5 | smoke calculation (Si fcc primitive cell, 2 atoms, ENCUT 250) | electronic loop converged (`aborting loop because EDIFF is reached`), E = <value> eV (baseline delta <value>) |
| 6 | concurrency probe | <n> jobs x <m> cores accepted |
| 7 | clean-up | `squeue` empty |

## 7. Known limitations and outstanding work

- <e.g. skill default of 192 cores per node overridden through step configuration>
- <e.g. optional component not installed; partition unavailable>
- <e.g. bastion-host ControlMaster must be re-established daily>

## 8. Change log

| Date | Change | Effect |
|---|---|---|
| <YYYY-MM-DD> | initial deployment | -- |