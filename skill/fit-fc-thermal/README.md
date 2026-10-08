# fit-fc-thermal — force constants (fc2/fc3) + lattice thermal conductivity (κ)

Fit second- and third-order interatomic force constants from an existing
displacement + force dataset, then solve the three-phonon BTE for the lattice
thermal conductivity. Three interchangeable engines:

| `FIT_ENGINE` | Method | Dataset it accepts |
|---|---|---|
| `phono3py` | phono3py + symfc or alm (least squares) | finite-displacement **and** random-displacement |
| `pheasy` | pheasy compressive sensing (OLS / LASSO / ALASSO / RFE-OLS / RFE-OLS-TSQR / RIDGE / ARDR / RVM) | random-displacement |
| `hiphive` | hiphive cluster space + linear regression (ols / ridge / lasso / ard / bayes), optional Huang + Born-Huang projection | random-displacement |

The skill **never runs VASP** and **never generates displacements**. It is the
force-constant-fitting branch lifted out of `kl-dft-cpu` S5_fc
(`gen_step5_fc.py` + `kl_fc_backends.py` + `submit_fit_*.tpl` +
`templates/step5_fc/step.conf`), decoupled from the thermal-conductivity
pipeline so that any dataset can be fitted, and extended with a third engine
(hiphive).  It now **closes the loop with κ**: `S2_kappa` runs the phono3py BTE on
the fitted fc2/fc3 and writes `kappa_summary.json`, so a displacement/force dataset
goes in and the lattice thermal conductivity comes out.

## Pipeline

~~~text
S1_fit   (submitted job, compute node)
         prep  normalise the dataset  -> POSCAR / SPOSCAR / dataset_*.npy
                                         disp_matrix.pkl (pheasy input)
                                         phono3py_*.yaml passthrough, BORN
         fit   fc2 (+ fc3) with the selected engine   -> fc2.hdf5 (+ fc3.hdf5)
         post  ShengBTE export (optional) + imaginary-frequency gate + kl_bundle
                                         -> phonon_summary.json, fc_fit_summary.json
S2_kappa (submitted job, needs S1_fit)
         phono3py three-phonon BTE (RTA/LBTE) on the fitted fc2/fc3
                                         -> kappa_summary.json, kappa-m*.hdf5
S2_plot  (login node, optional)  phonon band figures from fc2.hdf5
                                 + ZA bending-branch exponent p
S3_sweep (optional, method_sweep, fanout m-*)  one job per
         engine x method x neighbour shell N: fit + BTE -> sweep_result.json
S4_compare (login node, needs S3_sweep)  -> fit_compare.json / fit_compare.csv
~~~

## Quick start

~~~bash
cd <material root>
autozt -tt fit-fc-thermal -p <material> -j step1_fit conf          # effective parameters
autozt -tt fit-fc-thermal -p <material> -j step1_fit \
   conf --set params.FIT_ENGINE=hiphive               # switch engine
autozt -tt fit-fc-thermal -p <material> start                     # generate inputs + submit
autozt -tt fit-fc-thermal -p <material> status                    # collect and inspect

# κ 参数（BTE 网格 / 温度 / 同位素 / 解法）走 step2_kappa 的 step.conf。
# 默认 MESH=auto + MESH_CONV=auto：按倒格矢自动定网格并自动做 q 网格收敛；
# 要固定网格（不做收敛）：
autozt -tt fit-fc-thermal -p <material> -j step2_kappa \
   conf --set params.MESH="16 16 16"
~~~

The gen step locates the dataset automatically. It searches `step4_disp`,
`step2_disp_force`, `step2_disp`, `step3_disp` in the skill directory, then the
sibling skills' step directories (`../kl-dft-cpu/step4_disp`,
`../kl-mlff-cpu/step2_disp_force`, `../kl-mlff-gpu/step2_disp_force`,
`../phonon-dft-cpu/step2_disp`, `../phonon-mlff-cpu/step2_disp_force`,
`../phonon-mlff-gpu/step2_disp_force`), then `step*disp*` / `step*force*`
directly under the material root (a hand-assembled dataset such as
`<material>/step4_disp`), and then every sibling skill directory whose step
starts with `step*disp*` or `step*force*`, so a producer that does not exist
yet is found without editing anything. The skill's own output directories
(`step1_fit/`, `step2_kappa/`) are never taken as input — a previous run's
normalised copy used to shadow the real dataset. Point it somewhere else explicitly:

~~~bash
autozt -tt fit-fc-thermal -p <material> -j step1_fit \
   conf --set params.FIT_INPUT_DIR=../<other-skill>/<step>
~~~

## Dataset contract

Six input layouts are recognised, in this priority order:

| # | Files in the dataset directory | Needs |
|---|---|---|
| 1 | `phono3py_params.yaml` (forces embedded) | phono3py |
| 2 | `phono3py_disp.yaml` + `FORCES_FC3` | phono3py |
| 3 | `phono3py_disp.yaml` (forces embedded) | phono3py |
| 4 | `dataset_disps.npy` + `dataset_forces.npy` | phonopy only |
| 5 | `disp_matrix.pkl` + `force_matrix.pkl` | phonopy only |
| 6 | `disp-*/vasprun.xml` (+ `phono3py_disp.yaml`) | phonopy + phono3py |

Layouts 1-3 and 6 are read with phono3py's own YAML/dataset machinery; layouts
4-5 are plain NumPy/pickle arrays. In every case prep normalises the data to

~~~text
POSCAR                                  unit cell
SPOSCAR                                 supercell
dataset_disps.npy, dataset_forces.npy   (n+1, natom_super, 3); Cartesian
                                        displacements, trailing frame = the
                                        zero equilibrium reference
disp_matrix.pkl, force_matrix.pkl       (n, natom_super, 3), equilibrium
                                        subtracted - this is what pheasy reads
fc_dataset.json                         frames, atoms, RMS displacement,
                                        coords_mode / reference_frame_index /
                                        equilibrium_source, provenance, matrices
~~~

`COORDS = auto` (default) inspects `dataset_disps.npy` and decides by itself:

* **absolute fractional coordinates** - every value in [0,1) and the frame
  closest to the ideal `SPOSCAR` is within ~1 A of it.  Each frame's
  minimum-image displacement from the ideal supercell is converted to Cartesian
  with the supercell lattice;
* **Cartesian displacements** - centred on zero (about half the values
  negative), so subtracting the ideal POSITIONS leaves a several-Angstrom
  residual.

The **equilibrium frame is located wherever it sits** (first, middle, last or
absent), not assumed to be the first/last frame: the frame with the smallest
RMS displacement from the ideal supercell becomes the reference and is dropped
from the training set.  When no frame is exactly ideal, every frame is kept and
the closest frame's forces are used as the equilibrium approximation when it is
within 0.25 A (otherwise set `EQUILIBRIUM_FORCES_NPY`).  Force the old behaviour
with `COORDS = fractional` / `COORDS = cartesian`; the auto choice and reference
index are recorded in `fc_dataset.json` (`coords_mode`,
`reference_frame_index`, `equilibrium_source`).

Layouts 4-6 do not always carry a supercell matrix (a bare `.npy`/`.pkl`
dataset has no YAML to record one), so the gen falls back to the `POSCAR` /
`SPOSCAR` edge ratio for `SUPERCELL`, and `prep` records the result in
`fc_dataset.json` together with the atom order (see the pitfalls).

Layouts 4-5 carry no phono3py YAML either.  For `FIT_ENGINE = phono3py`, `prep`
now **synthesises `phono3py_params.yaml`** from the normalised displacement/force
arrays, so a bare `.npy`/`.pkl` dataset fits with the phono3py engine too (it
used to require the source to ship a YAML).

### Equilibrium residual forces

The reference (perfect) supercell is rarely force-free. Subtract its residual
forces from every frame or the fit is biased:

* `SUBTRACT_EQUILIBRIUM = true` (default) does it automatically when the
  dataset carries a reference: `disp-00000/vasprun.xml` for the VASP layout, or
  the frame closest to the ideal supercell - **wherever it sits** (first,
  middle or last; in fractional mode any frame matching `SPOSCAR`, in Cartesian
  mode a near-zero displacement frame).  That frame is dropped from training.
* `EQUILIBRIUM_FORCES_NPY = <file.npy>` supplies the reference by hand
  (`(natom_super, 3)`), for datasets without one.
* With no reference available the driver prints a note and fits the raw forces.
* `EQ_FORCE_MAX = 0.2` (eV/A, default): prep stops when the reference frame's
  largest force exceeds it - subtraction only removes the linear term, a
  reference far from the minimum biases fc2 itself.  `0` disables the gate.

## Engines

### Choosing the fitting method(s): `FIT_METHODS`

One parameter in `step1_fit/step.conf` decides what is fitted:

| `FIT_METHODS` | effect |
|---|---|
| `auto` (default) | `FIT_ENGINE` + that engine's method key; the defaults are **pheasy ALASSO** |
| `pheasy:RIDGE` (or a unique bare name: `ALASSO`, `symfc`, `ard` ...) | switch to that one method (`ridge` is ambiguous: pheasy or hiphive) |
| `pheasy:OLS,ALASSO hiphive:ridge` | all listed methods are computed |
| `pheasy:all` / `all` | every pheasy method / all 13 methods |

~~~bash
autozt -tt fit-fc-thermal -p <material> -j step1_fit conf --set params.FIT_METHODS="pheasy:RIDGE"
autozt -tt fit-fc-thermal -p <material> -j step1_fit conf --set params.FIT_METHODS=all
~~~

With several methods the **first** is the primary: it is `step1_fit/` itself
(the S1 gate, the S2 headline kappa), so S1/S2 behave exactly as with one
method.  The others are complete S1 recipes in `step1_fit/methods/<tag>/`
that run **serially in the same job** (non-fatal: a failing method is logged
in `methods/<tag>/run.log` and skipped); `step1_fit/methods.json` lists them.
S2_kappa then computes kappa for every method in `step2_kappa/methods/<tag>/`,
each with its own cutoff scan and its own 2x2 figure (below);
`kappa_driver.py compare` collects the figures as
`step2_kappa/figures/kappa_vs_cutoff_<tag>.{png,pdf}` and the numbers
(kappa at each method's cutoff, force RMSE, gate) in `methods_compare.json`.  Wall time adds up - raise the
`[submit]` time, or compare many methods in parallel with S3_sweep.

ALASSO/LASSO note: when the cross-validated alpha lands on the edge of the
grid the quality gate fails on purpose (`.fit_gate_fail`); widen
`PHEASY_MU_MIN`/`PHEASY_MU_MAX` or add frames, or switch to `pheasy:RIDGE`.

### phono3py (FIT_ENGINE = phono3py)

`phono3py.produce_fc2` / `produce_fc3` with `FC_CALC = symfc | alm`, written as
`fc2.hdf5` / `fc3.hdf5`. Works for both finite-displacement and
random-displacement datasets, because phono3py rebuilds the finite-difference
set itself from the YAML + `FORCES_FC3`.

* `FC3_CUTOFF` - third-order cutoff in Angstrom.  It never exceeds the
  supercell safe cutoff (half the shortest supercell height - 0.1 A): empty or
  larger values are clamped to it by the gen (the same clamp applies to
  `PHEASY_C3_CUTOFF` and `HIPHIVE_CUTOFF2/3`).
* `fc2.hdf5` is full `(N,N,3,3)`; `fc3.hdf5` is **compact**
  `(n_prim,N,N,3,3,3)` with its `p2s_map` (phono3py's own format).  On the
  243-atom test cell the full fc3 needed >15 GB and 3.1 GB on disk per cutoff;
  compact is 43 s / 1.2 MB.  phono3py (S2_kappa, kl-dft-cpu import) reads it
  natively; the ShengBTE export and the fit residual expand it on demand
  (`compact_fc3_to_full_fc3`, or a translation-based fallback on older
  phono3py), still subject to `FC3_LOAD_GB_LIMIT`.  Plots use fc2 only.
* symfc is the fast path; if it is not installed, or it fails on a very large
  supercell, install it (`pip install symfc`) or set `FC_CALC = alm`.

### pheasy (FIT_ENGINE = pheasy)

Runs the pheasy CLI in four steps - cluster space (`-s`), symmetry constraints
(`-c`), sensing matrix (`-d --disp_file`), fit (`-f --full_ifc`) - with
`--hdf5` so the output is the same `fc2.hdf5` / `fc3.hdf5` pair. The environment
tuning that used to live in the submit template (`PHEASY_ASR_*`, `PHEASY_SM_*`,
celer, two-level solvers) now lives in the driver (`_pheasy_env`), so the same
job runs on any cluster.

* `PHEASY_FIT_METHOD` - `OLS` (most memory hungry), `LASSO`, `ALASSO`,
  `RFE-OLS`, `RFE-OLS-TSQR`, `RIDGE`, `ARDR` (automatic relevance
  determination), `RVM` (relevance vector machine) (default: `ALASSO`).  The
  pre-rename spelling `RFE` is not accepted; write `RFE-OLS`.  ARDR/RVM are
  sparse-Bayesian fitters: the driver always passes `--std` for them; the
  pheasy knobs are `PHEASY_ARD_STD` / `PHEASY_ARDR_MAX_ITER` / `PHEASY_RVM_*`.
* `PHEASY_TSQR_CRITERION` - `RFE-OLS-TSQR` only: `aic` (pheasy's default since
  the rename; n = the force-component rows), `bic`, or `cv` (reproduces
  `RFE-OLS` exactly).  Empty = pheasy's default.
* `PHEASY_C2_CUTOFF` / `PHEASY_C3_CUTOFF` - cutoffs in Angstrom (empty = all
  interactions). Keep both comfortably below half the smallest supercell edge,
  otherwise periodic images double-count interactions.
* `PHEASY_RASR = auto` (default) resolves to `BHH` for a 2D slab and `none`
  for a bulk crystal. RASR imposes Born-Huang rotational invariance and the
  Huang equilibrium conditions, which is what makes a truncated fit physical.
  It is passed to pheasy's **null-space construction step (`-c`)** — that is
  the only step that reads it (`-f` loads `ns_*.npz` and re-imposes nothing),
  and the driver aborts if `-c` does not log the constraint. For 2D slabs BHH
  is mandatory (without it ZA goes linear near Gamma). On a bulk crystal it is
  not neutral: the MnIn2Se4 189-atom 3x3x3 pheasy OLS fit went from **0.62 %**
  force error / min frequency -0.03 THz with `none` to **8.9 %** / -1.57 THz
  (spurious imaginary modes) with `BHH`, at the same cutoff — so `auto` never
  puts BHH on a bulk cell. An explicit `BH` / `H` / `BHH` / `none` still wins.
* `PHEASY_BIN = pheasy-gpu` (the CPU `pheasy` build was removed). It uses the
  GPU submit template; on a cluster with no GPU nodes the fit cannot run there.
* **The GPU only accelerates the two-level sparse matvec.** pheasy builds a
  dense/CSR sensing matrix for LASSO unless `PHEASY_LASSO_TWOLEVEL=1` (OLS
  defaults to the two-level path via `PHEASY_OLS_TWOLEVEL`), and
  `PHEASY_GPU_LASSO_RESIDENT=1` / `PHEASY_GPU_TWOLEVEL_LASSO=1` *require* a
  `TwoLevelSM` — with the dense matrix pheasy raises
  `NotImplementedError: Resident GPU LASSO/ALASSO requires TwoLevelSM input`
  after the matrix has already been written. `PHEASY_LASSO_TWOLEVEL` (empty =
  auto: on when the resident GPU is requested) is the knob that makes them
  consistent; the driver also drops a resident request it cannot honour, with a
  `[WARN]`. `PHEASY_NGPU = N` is the one number the job's allocation and the
  solver share: the gen renders `--gres=gpu:N` from it (unless `[submit]` sets
  `gres`) and the driver exports `PHEASY_GPU_SM_NGPU=N` plus
  `PHEASY_GPU_SM_DEVICES=0..N-1`, so "run on three cards" is
  `PHEASY_NGPU = 3` and nothing else. `PHEASY_GPU_LASSO_RESIDENT = true`
  requests the resident backend from step.conf rather than the submit template.

  **Accuracy.** None of this changes the answer: the GPU matvec/adjoint match
  the CPU ones to ~1e-16 (pheasy-gpu's acceptance suite), the two-level matvec
  computes the same `SM_prime @ (NS @ x)` product as the materialised `SM`
  (only the rounding order differs, and it skips the float32 dense product
  altogether), and the resident LASSO keeps the same CV grouping, alpha grid and
  debias flow, in float64. Treat a systematic change in
  `pheasy_relative_error` / `pheasy_worst_force_correlation` / `fc2.hdf5` as
  a bug report, not as cost of doing business.

  Whatever happens, the fit records what actually ran:
  `pheasy_gpu_used` plus `pheasy_gpu_evidence` (`gpu_sm_matvec`,
  `gpu_resident_lasso`, `gpu_cv_folds`, `twolevel_sm`) land in
  `fit_metrics.json` / `fc_fit_summary.json`, and a GPU request with no
  marker prints a `[WARN]` instead of passing as a GPU run.
* **Quality gate.** For `LASSO`/`ALASSO` the driver parses the selected
  `alpha`; if it sits on the edge of the cross-validation grid the selection is
  meaningless (more frames are needed), the run is flagged and the S1 gate fails
  on purpose. Widen `PHEASY_MU_MIN`/`PHEASY_MU_MAX` or add frames.
* **Force correlation.** pheasy prints a per-configuration correlation of the
  fitted forces and warns when it drops. The driver reads it and records it as
  `pheasy_worst_force_correlation`. Below ~0.5 the dataset and pheasy's
  supercell almost certainly disagree on the **atom order** (see
  "Physics notes and pitfalls"), and the run is failed on purpose; between 0.5
  and 0.98 it is ordinary model error and only a note is printed.

#### Environment tuning: `PHEASY_TUNING`

`kl-dft-cpu` hard-coded the author's production solver settings in
`submit_fit_pheasy.tpl`. Those settings are not neutral on every dataset:

| `PHEASY_TUNING` | what it sets | BaS 250-atom / 10-frame / fc2+fc3 / OLS |
|---|---|---|
| `safe` (default) | memory layout, threads, CV grouping only | 2.2 % relative error, correlation 0.9997, gate **stable** |
| `kl` | the above **plus** the production solver block (`ILP64`, `ATOL`, `BTOL`) | same as `safe` here — those knobs are inert on this dataset |

The production solver block used to carry one more variable,
**`PHEASY_OLS_RIDGE=1e-4`**, and that one is not inert. It is not a normalised
sklearn ridge: pheasy's OLS turns it into Tikhonov damping
`damp = sqrt(ridge * n_samples)` appended to the LSMR system
(`core/optimizer.py::_ols_lsmr`), so its strength grows with the number of
equations — 0.87 at ndata = 7500, which swamps the design matrix. One
dataset, one sensing matrix, one seed:

| `PHEASY_OLS_RIDGE` | LSMR iterations | relative error | correlation |
|---|---|---|---|
| 0 | 509 | **2.2 %** | 0.9997 |
| 1e-10 | 508 | 2.2 % | 0.9997 |
| 1e-8 | 342 | 3.1 % | 0.9995 |
| 1e-6 | 72 | 6.8 % | 0.9977 |
| 1e-4 | 18 | **58 %** | 0.8991 |

`ridge > 0` also silently disables pheasy's resident **GPU** OLS, which falls
back to the CPU path. The failure mode is what makes this worth knowing:
pheasy still exits 0 and still writes `fc2.hdf5`; the only hints are the
correlation line and the LSMR iteration count in the log. The four production
templates that carried it (`jzzn`, `a800`, `3090`, `hanhai25`,
`submit_fit_pheasy*.tpl`) were corrected to `PHEASY_OLS_RIDGE=0` on
2026-09-11, so neither tuning profile sets it any more; an explicit
`PHEASY_OLS_RIDGE` in `step.conf` still wins, and the driver never inherits a
stray value from the environment.

### hiphive (FIT_ENGINE = hiphive)

Builds a `ClusterSpace` on the primitive cell from `HIPHIVE_CUTOFF2`
(+ `HIPHIVE_CUTOFF3` when `ENABLE_FC = 3`), adds every frame as a structure
carrying the documented `displacements` and `forces` arrays, stacks the design
matrix and solves it.

* `HIPHIVE_FIT_METHOD = ols | ridge | lasso | ard | bayes`, `HIPHIVE_ALPHA` for
  ridge / lasso. Both hiphive generations are supported: the legacy
  `hiphive.fitting.Optimizer` when present, otherwise a scikit-learn regressor
  (hiPhive 1.5 dropped its fitting module).
* `HIPHIVE_ENFORCE_ASR = true` projects the Huang and Born-Huang rotational sum
  rules onto the fitted parameters (`hiphive.enforce_rotational_sum_rules`).
* `HIPHIVE_N_CONFIGS` subsamples the frames when the design matrix becomes too
  large to hold in memory.
* hiphive derives its own cutoff-bounded atom list, but it aligns an arbitrary
  supercell internally, so the dataset supercell can be larger than the
  cluster-space cell. It must however be **at least twice the cutoff** in every
  direction, otherwise the periodic images alias.
* hiphive builds its `ClusterSpace` on the **primitive cell** and needs
  `2 * cutoff < min periodic width` there.  An acute/layered primitive (e.g. a
  7-atom rhombohedral cell with a 3.5 A interplanar spacing) is too thin for any
  physical cutoff and hiphive fails deep in orbit construction with a cryptic
  `(0, N) is not in list`.  The driver now catches this up front and prints an
  actionable error (use phono3py/pheasy, lower the cutoff, or set
  `HIPHIVE_CELL = supercell` to use the dataset supercell as the reference -
  the parameter count grows a lot).
* `fc2` is written densely (small). `fc3` is written with hiphive's streaming
  writers (`write_to_phono3py`, `write_to_shengBTE`) because the dense
  `(N, N, N, 3, 3, 3)` array is N^3 * 27 * 8 bytes - already ~3.4 GB for a
  250-atom supercell.
* `EXPORT_SHENGBTE` writes `shengbte/FORCE_CONSTANTS_2ND _3RD POSCAR`.  The
  **pheasy engine reuses the files pheasy itself wrote** (directly from its
  compact cluster IFCs, so the dense fc3 is never materialised); the phono3py
  engine re-exports through hiphive, and `FC3_LOAD_GB_LIMIT` (default 8 GB)
  skips the fc3 text when the dense `fc3.hdf5` would need more than that.
* **ShengBTE layout (checked, not a bug):** `FORCE_CONSTANTS_2ND` is the phonopy
  `FORCE_CONSTANTS` of the **full supercell** — its header is
  `natoms x prod(scell)` (9 atoms x 3x3x1 → `81 81`; 7 x 3x3x3 → `189 189`), which is
  exactly what ShengBTE's `read2fc` requires; `FORCE_CONSTANTS_3RD`, `POSCAR` and
  `CONTROL` refer to the **unit cell** (3RD atom indices 1..natoms + R vectors).
  Do not fold the 2ND file to the unit cell — ShengBTE then stops with "wrong number of
  force constants for the specified scell".
* The export also writes a runnable `shengbte/CONTROL` (lattvec/positions of the unit
  cell, `scell` = the 2ND supercell, `ngrid` from `SHENGBTE_NGRID`, temperatures from
  `SHENGBTE_T`, Born charges + `nonanalytic=.TRUE.` when `BORN` exists) and
  `shengbte/shengbte_manifest.json` with the checks: 2ND header = natoms x prod(scell),
  3RD max index <= natoms, supercell ordering (the 2ND file is re-ordered into
  ShengBTE's x-fastest/atom-slowest order if the dataset's SPOSCAR is not), diagonal
  supercell (ShengBTE cannot express a non-diagonal one → `usable: false`).
  2D: ShengBTE divides by the full cell volume — rescale the in-plane κ by c/thickness.
* Re-export without refitting (e.g. in the fetched `result/step1_fit/`):
  `python fc_fit_driver.py shengbte fit_config.json` (needs POSCAR, SPOSCAR, fc2/fc3
  or the existing shengbte/ files; exits non-zero when the set is not runnable).
  It re-exports from the source (pheasy's native files or fc2/fc3.hdf5) instead of
  re-using shengbte/, so files written or mis-ordered by older versions are replaced.
* **2ND atom order is checked everywhere ShengBTE files are written** (S1 export, S2
  `SOLVER=shengbte`, kl-mlff S4 — `fc_common.shengbte_fc2_ensure_order`): the file's
  actual layout is decided from the force constants themselves (a correct labelling is
  invariant under a lattice translation), so the check is idempotent — a file already
  in ShengBTE order (pheasy native output, a previous fix) is kept, a file in the
  dataset's SPOSCAR order (ASE-repeat interleave) is re-ordered, and a file matching
  neither is refused (`fc2_order` in shengbte_manifest.json).

## Lattice thermal conductivity (S2_kappa)

`S2_kappa` takes the fc2/fc3 that `S1_fit` wrote and runs the phono3py
three-phonon BTE (RTA by default, `BTE_METHOD = lbte` for the full solver),
writing `kappa_summary.json`:

~~~json
{"KAPPA_DONE": true, "bte_method": "rta", "mesh": "18 18 18", "nac": false,
 "mesh_converged": true,
 "temperatures": [100, 200, 300, ...], "kappa_xx_yy_zz": [[...], ...],
 "kappa_300K_xx_yy_zz": [xx, yy, zz], "kappa_inplane_300K": 0.5*(xx+yy)}
~~~

* **q-mesh: automatic + converged by default.**  `MESH = auto` sizes the mesh
  from the reciprocal lattice of phono3py's primitive cell with the phonopy
  length convention `n_i = max(1, round(MESH_LENGTH * |b_i|))` (`|b_i|` in 1/A,
  no 2pi; phono3py makes it symmetry-consistent), so a small cell gets a dense
  mesh and a big cell a sparse one without editing anything.  Note it is the
  reciprocal length that matters: a rhombohedral primitive cell with 16 A edges
  and a 14 deg angle has `|b| = 0.29 1/A`, i.e. the density of a 3.5 A cube.
  `MESH_CONV = auto` then reruns the BTE with `L *= MESH_CONV_FACTOR` (1.25)
  until kappa at `MESH_CONV_T` (300 K) changes by less than
  `MESH_CONV_TOL_PCT` (3 %) against the latest earlier mesh that is **coarser
  along every axis** (not simply the previous mesh).  Each mesh records its
  per-axis counts (`axis_counts`; for a generalized grid the gcd of each
  `grid_matrix` row, i.e. the conventional-cell divisions) and the mesh it was
  compared with (`compared_with`); axes that cannot change below
  `MESH_CONV_MAX_LENGTH` (a 2D vacuum axis) are listed in `frozen_axes` and
  exempt.  Why: with `L *= 1.25` the long axis of an anisotropic cell only
  steps every two or three meshes.  On Mn2In2Se5 (R-3m) `[1,14,42]` and
  `[1,18,54]` share one k_z division, the old consecutive-pair test passed at
  2.4 % without k_z ever being densified, and kappa then moved by 12-18 % when
  k_z was refined.  **`MESH_CONV_MODE
  = per_component` (default)** judges each diagonal component against *itself*,
  `max_i |d kappa_ii| / max(|kappa_ii|, MESH_CONV_FLOOR * max_j |kappa_jj|)`
  (off-diagonals normalised by the largest diagonal), so an anisotropic cell
  whose out-of-plane `zz` converges more slowly than the in-plane `xx/yy` is not
  declared converged prematurely.  The floor (`MESH_CONV_FLOOR = 0.1`) keeps a
  negligible component — a 2D vacuum axis, or a zz below 10 % of xx — from
  blocking convergence on its own relative noise; `MESH_CONV_FLOOR = 0` is the
  strict per-component test.  `MESH_CONV_MODE = max_ii` is the legacy
  whole-tensor norm `max|d kappa_ij| / max|kappa_ii|`, which does mask a slow zz
  (measured on a rhombohedral R-3m primitive).  Any component whose own change
  is still above the tolerance when the mesh is declared converged is printed
  as a `[WARN]` and listed in `mesh_convergence.components_over_tol`.
  Capped by
  `MESH_CONV_MAX_LENGTH` (150 A) and `MESH_CONV_MAX_POINTS` (64000 q-points).
  The denser mesh of the converged pair is reported; every point goes to
  `mesh_convergence.json` (each record also carries
  `rel_change_per_component` = xx/yy/zz), and `kappa_summary.json` carries
  `mesh`, `mesh_converged` and `mesh_convergence`.  Hitting a cap leaves
  `mesh_converged: false` with a warning (the densest mesh is still reported).
  With a cutoff scan the order is q-mesh -> cutoff -> temperature: the mesh is
  converged on a reference cutoff (the nominal one, else the largest; report in
  `mesh_convergence_reference.json`), every candidate cutoff is then run at that
  converged mesh and the cutoff is picked there, and every BTE run carries all
  of T_MIN..T_MAX, so the chosen cutoff's result is the full kappa(T).  If the
  chosen cutoff is not the reference, its mesh is re-checked from the converged
  length (at least one denser mesh).  `MESH = "n n n"` keeps the old fixed-mesh
  behaviour (no convergence); `MESH = 60` is a fixed starting length.
  No single default length is safe for every material (high-kappa, light-atom
  or low-T work needs denser meshes than a low-kappa complex cell), which is
  why the convergence loop is on by default rather than a bigger fixed mesh.
* `kappa_summary.json` also carries the full tensor
  (`kappa_voigt_xx_yy_zz_yz_xz_xy`), its eigenvalues `kappa_principal_300K`
  and `kappa_avg_300K` (trace/3): the Cartesian xx/yy/zz follow the POSCAR
  orientation and are not the crystal axes for a non-standard cell.
* Keys live in `step2_kappa/step.conf`: `MESH`, `MESH_LENGTH`, `MESH_CONV*`,
  `T_MIN/T_MAX/T_STEP`,
  `ISOTOPE`, `NAC` (auto = on iff `S1_fit/BORN` exists and the cell is not 2D), `BTE_METHOD`,
  `P3PY_OMP_THREADS` (phono3py is OpenMP-only: 1 process x N threads, never
  `mpirun`), `SBATCH_QOS`, `MESH_2D_VACUUM`.
* The driver (`kappa_driver.py`) loads `phono3py_params.yaml` (or
  `phono3py_disp.yaml` - generated on the fly from `POSCAR` +
  `fc_dataset.json` when the fit engine produced none), sets fc2/fc3 from the
  hdf5 and runs the BTE.  The third-order interaction is the fitted fc3; this
  step never re-fits.
* 2D: `MESH_2D_VACUUM = on` (default) pins the vacuum-axis mesh to 1 (the
  primitive axis with the shortest reciprocal vector, resolved on the compute
  node in phono3py's primitive basis).
* **Cutoff chosen by kappa** (`CUT3_KAPPA_SCAN = auto`, default): when `S1_fit`
  ran with `CUT3_SCAN = auto` (default) it refits fc3 at every candidate cutoff
  into `cutoff_scan/cut3_<c>/`.  `S2_kappa` then runs the BTE once per candidate
  and applies the **same three criteria as `kl-dft-cpu` S6** - residual 1-SE,
  per-shell stability (`cut <= stable_upper_cut`) and the kappa plateau
  (`CUT3_KAPPA_TOL_PCT`) - via `cut3_select.select_cutoff`, writing
  `cutoff_selection.json` and promoting the chosen cutoff's kappa to the
  top-level `kappa_summary.json` (`chosen_cutoff_A`).  `CUT3_PICK = smallest`
  (default) takes the smallest cutoff satisfying all three; `not_converged`
  reports the largest data-determined cutoff.  A phono3py fit writes per-cut
  dirs too (no bootstrap, so only the kappa plateau discriminates).  With
  `CUT3_SCAN = off` (or a hiphive fit, which does not write per-cut dirs yet) it
  falls back to the nominal cutoff.
* **Without a bootstrap** (`CUT3_BOOTSTRAP = 0`) the per-shell stability
  (criterion 2) cannot be measured: those records carry
  `stability_measured: false`, criterion 2 is skipped and the cutoff is chosen
  by the kappa plateau alone - the same rule as a phono3py scan.
  `cutoff_selection.json` says so (`stability_gate_applied: false`).  Older
  `cutoff_scan.json` files are recognised by their `bootstrap: 0`.  With a
  bootstrap, a candidate whose `stable_upper_cut` is null is still unusable.
* **Plateau per component**: the kappa-plateau criterion compares xx, yy and
  zz separately (each against `max(|kappa_ii|, 0.1 * max_j |kappa_jj|)`), so a
  slowly converging out-of-plane zz is no longer hidden behind a flat in-plane
  kappa; a 2D vacuum axis (zz ~ 0) never blocks it.  Records without
  per-component kappa fall back to the in-plane value.
* **2D thickness normalisation** (`KAPPA_2D_THICKNESS = vdw | cell | <A>`,
  default `vdw`): phono3py divides by the whole cell volume, vacuum included, so
  the raw in-plane kappa of a slab is diluted by `h_perp/d`.  For a 2D cell
  `kappa_summary.json` keeps the raw fields and adds
  `kappa_2d_normalized_xx_yy_zz`, `kappa_2d_normalized_inplane_xx_yy`,
  `kappa_2d_normalized_inplane_300K` (= raw * `h_perp/d`) plus the geometry in
  `kappa_2d_norm`, and the thickness-free sheet conductance
  `sheet_conductance_W_per_K_xx_yy` / `sheet_conductance_inplane_300K_W_per_K`
  (= kappa * d = raw kappa * h_perp, W/K; Wu et al., arXiv:1607.06542) -
  the kl-dft-cpu S6 convention, computed by the same
  `_common/thickness_2d.py` (`d` = atomic span + top/bottom vdW radii; a number
  fixes `d`, e.g. 6.15 A, the MoS2 bulk interlayer spacing; `cell` = no
  normalisation).  Compare with the literature only after checking which `d`
  the paper used.  The zz component of a 2D layer is not physical (no
  dispersion across the vacuum); the out-of-plane transport of a stack is a
  bulk / multilayer calculation.
* **kappa vs cutoff**: every scan writes `kappa_vs_cutoff.json` (per cutoff:
  neighbour shells inside it, kappa xx/yy/zz and in-plane at 300 K, the
  per-component % change from the previous cutoff, the criteria flags) and
  `kappa_vs_cutoff.png` (300 dpi) + `kappa_vs_cutoff.pdf` (vector), journal
  style, always 2 x 2: (a) kappa_xx top left, (b) kappa_zz top right,
  (c) kappa_yy bottom left, (d) force RMSE bottom right; the method name sits
  above (a).  Every point is labelled with its value and every kappa segment
  with its % change (bold inside the plateau tolerance).  A 2D layer keeps
  panel (b) as an empty frame (kappa_zz is not physical) and plots xx/yy
  thickness-normalised; a 3D cell fills it.  Without a cutoff scan (e.g. the
  hiphive engine) the figure shows the single nominal cutoff, with the fit's
  own RMSE.  Panel (d): the force RMSE (meV/A) of each cutoff's
  fit - S1 evaluates every `cutoff_scan/cut3_<c>/` on `FIT_RMSE_FRAMES` frames
  and stores `force_rmse_eV_per_A` in `cutoff_scan.json` (pheasy's own
  relative error is shown when no RMSE was recorded).  The chosen cutoff is
  dashed in every panel; no axis zooms below 10 % of its values, so a tiny
  wiggle never looks like a trend.  Written also when no cutoff is usable, which is when the curve is
  most needed.  The points are at the converged q-mesh, i.e. the same mesh as the reported kappa.
* **Candidate spacing**: the candidates are the midpoints between *adjacent*
  neighbour shells of the structure (all atom pairs, shells closer than 0.05 A
  merged), so each step adds exactly one shell; a wide gap between two
  candidates means a wide gap between two shells of that structure (MoS2
  monolayer: 5.44 A = sqrt(3) a, 6.28 A = 2a -> candidates 5.86 / 6.50 A), not a
  skipped shell.
* **Cutoff cap**: fc3 can only be trusted up to half the supercell's smallest
  periodic width (5.27 A for the 3x3x3 / 189-atom cell).  Candidates beyond it
  are dropped.  To sweep larger cutoffs, regenerate the dataset with a larger
  supercell (4x4x4 -> ~7 A).
* **Per-cut phonon spectra**: `S2_plot` draws the dispersion for the nominal fc2
  and, when the fit wrote `cutoff_scan/cut3_<c>/`, one dispersion per candidate
  too (`phonon_band_plot/cut3_<c>/`); the per-cut min/max frequencies land in
  `phonon_band_summary.json` under `cut3_bands`.
* The compute host needs `phono3py` + `h5py` (the same env as the fit job).
## Method sweep (S3_sweep / S4_compare)

Comparing fitting methods on one dataset is a first-class step, not a script:
the optional group `method_sweep` (off by default, independent of S1/S2 —
`needs: []`) fans out one job per **engine x method x neighbour shell**, each
fitting the dataset and running the BTE on its own fc2/fc3.

~~~bash
autozt -tt fit-fc-thermal -p <material> -j S3_sweep conf \
   --set params.SWEEP_METHODS="phono3py:symfc pheasy:OLS,LASSO,RIDGE hiphive:ridge"
autozt -tt fit-fc-thermal -p <material> -j S3_sweep conf --set params.SWEEP_C3_SHELLS="4 5 6"
autozt -tt fit-fc-thermal -p <material> -j S3_sweep start      # enables the group for this run
# ... when every m-* is done, S4_compare writes step4_compare/fit_compare.{json,csv}
~~~

* `SWEEP_METHODS` - `all` (15: phono3py symfc/alm, pheasy x8, hiphive x5) or
  `engine[:m1,m2|all]` tokens; an engine alone means its default method.
* `SWEEP_C3_SHELLS` - the third-order cutoff as **"up to the N-th nearest
  neighbours"**, with the ShengBTE `thirdorder.py -n` definition
  (`thirdorder_common.calc_frange`): for every atom the distinct neighbour
  distances `u_1 < u_2 < ...` are listed and `cutoff(N) = max over atoms of
  (u_N + u_{N+1}) / 2`.  N therefore means the same as in a paper that fitted
  "fc3 up to the 5th neighbours".  `auto` = the largest N inside the supercell
  safe cutoff; `5`, `4 5 6` and `3-6` are accepted; an N beyond the safe cutoff
  is skipped with a warning (enlarge the supercell instead).  The N -> cutoff
  table and every atom's shells land in `step3_sweep/sweep_plan.json`.
  `SWEEP_SHELL_TOL` (1e-4 A, thirdorder's tolerance) merges nearly degenerate
  distances; raising it changes the N count.
* `SWEEP_C2` - `all` or a neighbour number for the pheasy/hiphive fc2
  (phono3py/symfc fits the full fc2).
* Every S1_fit/S2_kappa key (dataset, `SUPERCELL`, `HIPHIVE_CELL`, `MESH`,
  temperatures, conda env, `[submit]`) is accepted in `step3_sweep/step.conf`
  and applies to all variants; per variant the sweep only overrides engine,
  method and cutoffs, and turns the per-variant cutoff scan off.  A fixed `MESH
  = "n n 1"` is the cheap like-for-like choice; `auto` converges every variant.
* Each variant is the unmodified S1 recipe (`gen_step1_fit.main`) plus the S2
  recipe finished on the compute node (`sweep_driver.py kappa-prep`).  A fit
  with imaginary modes records `stable=false` and skips kappa
  (`SWEEP_KAPPA_IF_IMAG = true` runs it anyway); a tool error leaves no marker,
  so `retry` re-runs only the failed variants.  The gen is idempotent: a
  variant with `sweep_result.json` is never regenerated.
* `fit_compare.json` / `.csv`: per variant the gate (`stable`,
  `min_frequency_THz`), the fit residual (`fit_rmse_relative`,
  `pheasy_relative_error`) and kappa(300 K) (in-plane, trace/3, xx/yy/zz), plus
  the shell table.

## Artifacts

~~~text
step1_fit/
  fit_config.json          everything the compute-node driver reads
  submit.sh                rendered from the engine template
  POSCAR SPOSCAR           cells used for the fit and the gate
  dataset_*.npy            normalised dataset
  disp_matrix.pkl force_matrix.pkl   pheasy input
  phono3py_disp.yaml phono3py_params.yaml FORCES_FC3 BORN   passthrough
  fc2.hdf5 [fc3.hdf5]      fitted force constants (phono3py layout)
  FORCE_CONSTANTS          fc2 in phonopy text format (hiphive path)
  shengbte/FORCE_CONSTANTS_2ND _3RD POSCAR CONTROL shengbte_manifest.json
                           ready-to-run ShengBTE / fourphonon set (2ND = supercell)
  band-dft-cpu.yaml        phonopy band structure of the fit
  cutoff_scan.json         shell / cutoff determination (see below)
  cutoff_scan/cut3_<c>/    per-candidate fc2/fc3 + pheasy logs (scan mode)
  fc_dataset.json          dataset audit (frames, atoms, RMS displacement,
                           supercell, frozen atom order, ...)
  fit_metrics.json         engine-reported fit quality (RMSE, correlation, ...)
  phonon_summary.json      the gate verdict (this is the step's marker)
  fc_fit_summary.json      identical content, skill-facing name
  queue.out queue.err      job log
~~~

`phonon_summary.json` carries `tool_ok`, `stable` and `min_frequency_THz`, so autozt's
built-in `phonon` judge gives three outcomes:

* **stable** - the step is done, downstream steps may run;
* **imaginary frequency** - the fit finished but the mesh minimum is below
  `-IMAG_THR`; the step is *not* done and downstream steps stay held back (this
  is a physics result, not an error);
* **fc2 breaks the crystal symmetry** (`status: fc2_symmetry_broken`, autozt
  shows `fc2 breaks crystal symmetry (spread ... THz) -- check supercell`) -
  symmetry-equivalent q points differ by more than `FC2_SYM_TOL` (1e-3 THz;
  a symmetric fc2 gives ~1e-8).  Usually the supercell does not keep every
  point-group operation of the crystal.  Not a physics result: imaginary
  frequencies and kappa from this fc2 are not usable;
* **tool error** - the gate itself could not be evaluated; the job exits
  non-zero so the step shows as error and the log is worth reading.

How the gate samples (2026-10-07): the 3D mesh is the length
`IMAG_MESH_LENGTH` (100 A) evaluated on the **full** mesh with no symmetry
reduction, plus the band path.  The old `run_mesh(60.0)` with symmetry
reduction gave 17x17x4 on Mn2In2Se5 and evaluated 315 q points; the spurious
modes of its symmetry-broken fc2 sat at low-symmetry points inside the zone and
none were sampled, while the high-symmetry paths were all positive.  The gate
also checks the fc2 itself (`fc2_symmetry` in the summary: the largest
frequency spread inside 12 random symmetry stars) and, before evaluating
anything, reorders the fc2 from the dataset's SPOSCAR atom order into the order
phonopy rebuilds from POSCAR: `_phonopy_unitcell` wraps fractional coordinates
into [0, 1), and an atom written at exactly 1.0 (or a negative coordinate) then
has its periodic images permuted relative to the dataset.

**Supercell symmetry.**  `SUPERCELL_SYMMETRY = strict` (default) stops the fit
when the dataset supercell does not keep every point-group operation of the
crystal (checked at gen when the dataset ships a POSCAR, and again on the
compute node before the fit); `warn` only prints, `off` skips.  The same key
guards the displacement generators (kl-dft-cpu S4, phonon-dft-cpu S2,
phonon-mlff S2, kl-mlff S2) and S2_kappa (`kappa_driver` stops before the BTE).
A rhombohedral primitive cell with `a1`, `a2` in the plane and a tilted `a3`
repeated `n n 1` loses the 3-fold axis (Mn2In2Se5 3x3x1 keeps 4 of 12
operations); the message suggests the smallest symmetric diagonal supercell
(there: `3 3 3`).

`fc_fit_summary.json` adds the engine, frame count, atom count, the NAC/no-NAC
minima, the ShengBTE export status and - when `FIT_RMSE_FRAMES > 0` (default 10;
needs hiphive in the job environment) - the
training residual (`fit_rmse_eV_per_A`, `fit_rmse_relative`) evaluated by
predicting the forces of that many frames from the *fitted* force constants.
That number is engine independent and is the quickest way to tell a good fit
from a bad one; for a decent dataset it is well under 1 per cent relative.

The optional S2_plot step writes `phonon_band_plot/phonon_band_summary.json`
together with the two PNGs.  Besides the min/max frequencies it carries the ZA
bending-branch exponent fitted along the two inequivalent in-plane directions of
the cell - the same `za_power_law` fit and direction rule `kl-dft-cpu` uses:

* `za_exponent_q1`, `za_exponent_q2` - log-log least-squares exponent p of the
  lowest branch over q in (0, `za_qmax`]; `za_minfreq_q1_THz` /
  `za_minfreq_q2_THz` are the corresponding lowest frequencies;
* `za_qdir_q1` / `za_qdir_q2` - the reduced-coordinate directions, e.g.
  (1,0,0) and (1,1,0) for a hexagonal cell, plus `za_qmax`, `za_n_qpoints` and
  `za_p_range` recording the fit settings;
* `za_ok` - the **coarse machine verdict**: true when both directions satisfy
  1.7 < p < 2.3, the un-stressed 2D Born-Huang/Huang quadratic ZA branch.  It is
  a one-bit decision and must not be read on its own when a direction is close
  to a bound or the fit is noisy - use the fields below;
* `za_r2_q1` / `za_r2_q2` and `za_log_resid_rms_q1` / `za_log_resid_rms_q2` -
  the goodness of each log-log fit: R^2 and the RMS residual in log space.  A low
  R^2 (or a large residual) means the fitted p is uncertain whatever its value;
* `za_margin_q1` / `za_margin_q2` - `min(p - 1.7, 2.3 - p)`, how far each p sits
  from the nearer criterion bound.  A small margin means the verdict can flip on
  a tiny change of dataset or q-window;
* `za_needs_review` - true when some direction has `za_margin_qX < za_margin_min`
  (default 0.20) **or** `za_r2_qX < za_r2_min` (default 0.98); the two thresholds
  are recorded in the summary as `za_margin_min` / `za_r2_min`.  When true,
  `za_note` explicitly says 建议人工复核，不要仅凭 za_ok 下结论 (manual review, do
  not decide from `za_ok` alone) and the script prints a `[WARN]` line;
* `za_is_2d` / `za_vacuum_axis` - the dimension and vacuum axis, taken from the
  `dim` S1_fit already resolved and from the same `dim_common.detect_dimension`
  the fit gen uses; `za_note` explains the verdict.

The criterion is meaningful only for a 2D layer: for a 3D cell `za_is_2d` is
false, the exponents are informational, `za_ok` should not be read as a
stability verdict and no review is flagged.  A direction with no positive
frequency (a linearised or imaginary ZA branch) gets `null` for its exponent plus
a `za_note` saying so - it never aborts the plot step.

The engines also report their own quality numbers, which the driver scrapes out
of the log into `fit_metrics.json` and copies into both summaries:
`hiphive_backend`, `hiphive_parameters`, `hiphive_design_matrix`,
`pheasy_relative_error`, `pheasy_worst_force_correlation`,
`pheasy_free_ifcs`. pheasy's correlation is the fastest warning sign
available and the driver acts on it (see the pheasy section).

## Shell / third-order cutoff determination

Ported from `kl-dft-cpu` S4/S5 (2026-09-24 user 流程). Setting the fc3 cutoff by
hand is guesswork; instead the skill enumerates the neighbour shells of the unit
cell and takes the **midpoints between adjacent shells** as candidate cutoffs
(never on a shell distance), then lets the fit decide how far the data can
actually determine the fc3.

* The gen always enumerates the shells and stores the candidates in
  `fit_config.json` (it prints them). Candidates beyond the supercell safe
  cutoff (`0.5 x` the shortest cell height, minus 0.1 A) are dropped with a
  `[WARN]` - a larger cutoff makes periodic images double-count.
* `CUT3_CANDIDATES = auto | off | "3.6 4.2 4.8 ..."` chooses them;
  `CUT3_MAX`, `CUT3_MIN_SHELLS`, `CUT3_GAP_TOL`, `CUT3_MIN_GAP` tune the
  enumeration (see `fc_common.neighbor_shells` / `cut3_candidates`).
* `CUT3_SCAN = off` by default. Set it to `auto`/`on` to run the determination:

  * **pheasy** (the full port): every candidate is refitted (`-s -c -d -f`) and
    then refit `CUT3_BOOTSTRAP` times on bootstrap-resampled frames. For each
    shell of the supercell the driver records the mean |Phi^3| and its scatter
    across the bootstrap fits; walking outwards while `sigma/|mean| <
    CUT3_STABILITY_THR` gives **`stable_upper_cut`**, the largest cutoff the
    data can determine. A candidate is usable when `cut <= stable_upper_cut`;
    `recommended_cut` is the largest usable candidate. Each candidate's
    `fc2.hdf5` / `fc3.hdf5` and pheasy logs land in `cutoff_scan/cut3_<c>/`.
  * **phono3py / hiphive**: a single-fit report - the per-shell mean |Phi^3| of
    the nominal fit. No bootstrap, so `stable_upper_cut` stays `null`; use it
    to see which shell the nominal cutoff falls between, not as a verdict.

  The result is `cutoff_scan.json` (`mode`, `shells`, `candidates`, `records`
  with `shell_stats` / `stable_upper_cut`, `recommended_cut`, `note`).
* The final cutoff choice is deliberately left to the caller: `fit-fc-thermal` never
  computes kappa, so it cannot apply the kl S6 plateau criterion. The
  `stable_upper_cut` is the data-side bound; a downstream kappa step decides
  between the usable candidates (kl's `select_cutoff`).

`tests/suite_fcfit_shell.py` covers the algorithms without touching a cluster.

## Physics notes and pitfalls

* **Supercell vs cutoff.** A cutoff larger than half the smallest supercell
  edge makes periodic images contribute twice. Either enlarge the dataset
  supercell or lower the cutoff.
* **Atom order (pheasy only).** pheasy does not read `SPOSCAR`: it rebuilds the
  supercell from `POSCAR` + `SUPERCELL` as *per-primitive-atom blocks* (all
  images of primitive atom 0, then atom 1, ...). phonopy and phono3py agree with
  that layout, so `kl-*` datasets are fine. A dataset assembled with ASE's
  `Atoms.repeat()` is *not*: it interleaves the images, so every atom except
  the first is scrambled. The fit still converges, still exits 0 and still
  writes `fc2.hdf5` -- it just lands near 50 % residual with a force correlation
  around 0.7. `prep` therefore compares the dataset supercell against pheasy's
  convention and freezes the answer in `fc_dataset.json`; when they differ the
  pheasy engine permutes `disp_matrix.pkl` / `force_matrix.pkl` into pheasy's
  order (the npy files and every other engine keep the dataset's order).
  Measured on the 16-atom test dataset: relative error 0.70 -> 0.0046,
  correlation 0.67 -> 1.0000.
* **pheasy rewrites `SPOSCAR`.** Its CLI writes the supercell in its own atom
  order into the step directory, so `SPOSCAR` describes pheasy's order after a
  pheasy run and the dataset's order after `prep`. `prep` always refreshes it
  from the dataset, and everything order-sensitive reads the frozen record
  instead of the file.
* **Training residual is not a validation score.** A random-displacement fit
  with enough parameters can reproduce its own training forces while
  extrapolating badly. Raise `FIT_RMSE_FRAMES` and, ideally, hold frames out.
* **2D materials.** For `DIM = 2d` the gate verdict uses the no-NAC minimum: the
  3D Coulomb kernel produces spurious imaginary frequencies near Gamma in a
  strictly 2D material. `DIM = auto` inherits the dimension from the dataset's
  `kl_params.txt` / `workflow_method.txt` when a sibling kl/phonon skill wrote
  one, and otherwise detects the vacuum axis.
* **NAC.** A `BORN` file in the dataset (from `kl-dft-cpu` step3_nac) is picked
  up automatically; both the NAC and no-NAC minima are reported, and only the
  3D verdict uses the NAC one.
* **fc3 memory**: an empty `FC3_CUTOFF` used to mean a full-supercell fc3 and
  was the classic way to exhaust a compute node; it is now clamped to the safe
  cutoff and fc3 is written compact.  For the ShengBTE export,
  `FC3_LOAD_GB_LIMIT` still skips the fc3 text when the expanded array would
  exceed it (the pheasy engine is unaffected — it writes ShengBTE from compact
  IFCs).

## Relation to the kl skills

| | `kl-dft-cpu` S5_fc | `fit-fc-thermal` S1_fit |
|---|---|---|
| Input | hard-wired to its own `step4_disp` | any dataset, `FIT_INPUT_DIR` or auto-search |
| Engines | phono3py, pheasy | phono3py, pheasy, **hiphive** |
| pheasy logic | inside the submit template | inside `fc_fit_driver.py` (testable, cluster independent) |
| Shell/cutoff determination | S4 candidates + S5 scan + S6 kappa selection | candidates + S5-equivalent scan ported; **no** S6 (no kappa step) |
| Output | `step5_fc/phono3py/fc2.hdf5` (+ shengbte/) | `step1_fit/fc2.hdf5` (+ shengbte/, cutoff_scan.json, **kl_bundle/**) |
| Marker | `phonon_summary.json` | `phonon_summary.json` (+ `fc_fit_summary.json`) |
| Next step | kappa from the BTE solver | kappa: `kl-dft-cpu` S5_fc auto-imports `kl_bundle/` |

`fit-fc-thermal` S1_fit also writes a small `kl_bundle/` (`fc2/fc3.hdf5` + `shengbte/` +
`phonon_summary.json` + `fc_dataset.json` + `POSCAR`; `EXPORT_KL_BUNDLE=true`). It is
fetched locally with the rest of the step, and `kl-dft-cpu`'s S5_fc declares a
`push_paths` entry so tf **pushes it to whichever cluster κ runs on** and the kl gen
auto-imports it (no manual scp, no step.conf edit). "Fit on 3090, κ on jzzn/hanhai"
is therefore the default flow.

Nothing here replaces the kl skills: `kl-dft-cpu` / `kl-mlff-*` still own the
end-to-end thermal-conductivity workflow. `fit-fc-thermal` is for fitting force
constants on their own - from a dataset another skill produced, from a
hand-assembled dataset, or as a sandbox for comparing engines on the same data.

## AutoZT integration notes

Two things are specific to how this skill plugs into `autozt`:

* **`submit_required`.** `autozt/workflow.py::_remote_submit_preflight` checks,
  before every `sbatch`, that the step's inputs exist locally, and its
  historical default is `INCAR, POSCAR, KPOINTS, submit.sh` for anything whose
  skill name does not contain *mace*. A fitting skill runs no VASP and its gen
  writes neither INCAR nor KPOINTS (the cells are copied by `prep` on the
  compute node), so the default would block every submission. The skill
  therefore declares what it really produces:

  ~~~yaml
  submit_required: [fit_config.json, submit.sh]
  ~~~

  `autozt/bootstrap.py` carries the key through `_MANIFEST_TYPE_KEYS`, and the
  preflight falls back to the old default when a skill does not declare it,
  so every existing skill behaves exactly as before.
* **The dataset is not part of the skill.** `autozt` pushes only the gen, its
  `gen_need` files and the templates; the displacement+force dataset stays
  where it is. Point `FIT_INPUT_DIR` at it (an absolute path on the cluster
  works) or let the auto-search find a sibling `kl-*`/`phonon-*` step.

## Verification status

The rows below were run with the `atomate2_p_a` environment, pinned to what
the jzzn login node has (phonopy 2.47.1, phono3py 3.24.0, symfc 1.6.0,
hiphive 1.5, numpy 2.3.5), by invoking the compute-node driver directly — the
same `prep` / `fit` / `post` sequence the submit script runs -, except the
last row, which is a real `autozt start` submission. Datasets: the 250-atom BaS
supercell produced by `kl-dft-cpu` (10 frames, `phono3py_params.yaml`), a
16-atom / 100-frame model dataset in the three array layouts, and a 128-atom
Si 4x4x4 supercell (90 frames) on the cluster.

| Case | Engine | Result |
|---|---|---|
| BaS 250 atoms, fc2 | hiphive (cutoff 6 A, 34 parameters) | stable, mesh min -2.5e-07 THz, training RMSE 6.4 % |
| BaS 250 atoms, fc2+fc3 | phono3py + symfc, c3 = 4 A | stable, min -8.8e-07 THz, RMSE 1.3 %, `fc2.hdf5` + `fc3.hdf5` + ShengBTE |
| BaS 250 atoms, fc2+fc3 | pheasy OLS, c3 = 4 A (189 free IFCs) | stable, mesh min -0.000 THz, relative error 2.2 %, correlation 0.9997 |
| 16 atoms, fc2+fc3 | pheasy LASSO (40-alpha grid) | relative error 0.46 %, correlation 1.0000 |
| 16 atoms, fc2+fc3 | hiphive (streaming fc3 writers, `shengbte/FORCE_CONSTANTS_3RD` written) | stable, RMSE 0.14 % |
| In2MnSe4 189 atoms (arrays only, 45 frames, 3x3x3) | phono3py + symfc, c3 = 5 A | stable, min -0.030 THz; `fc2.hdf5` 0.4 MB + `fc3.hdf5` 8.9 MB + ShengBTE + `kl_bundle/` (11 files) |
| In2MnSe4 189 atoms (arrays only, 45 frames) | pheasy RIDGE, c3 = 5 A (1272 free IFCs) | stable, min -0.031 THz, RMSE 9.5e-4 eV/A, relative error 0.63 %, `kl_bundle/` (10 files) |
| In2MnSe4 189 atoms | pheasy OLS, c3 = 5 A | stable, RMSE 9.5e-4 eV/A, relative error 0.63 % |
| In2MnSe4 189 atoms | pheasy LASSO, c3 = 5 A | stable, RMSE 1.66e-3 eV/A, relative error 1.10 % |
| In2MnSe4 189 atoms | pheasy ALASSO, c3 = 5 A | fit relative error 0.66 % but **stable=False** (significant imaginary frequency) - use RIDGE/OLS/RFE-OLS here |
| In2MnSe4 189 atoms | pheasy RFE-OLS / RFE-OLS-TSQR, c3 = 3.5 A (CPU smoke) | both stable, relative error 0.94 % each |
| In2MnSe4 189 atoms, phono3py fc2/fc3 | **S2_kappa** (phono3py BTE RTA, mesh 8x8x8, 200-400 K) | kappa300K in-plane 3.68, zz 1.63 W/mK -> `kappa_summary.json` KAPPA_DONE |
| In2MnSe4 189 atoms | hiphive | not applicable: primitive min periodic width 3.51 A < 2xcutoff; driver now errors clearly (see the hiphive section) |
| 16 atoms, fc2 | hiphive, pkl-only dataset (auto-detected, supercell deduced) | stable, mesh min -9.9e-08 THz |
| BaS 250 atoms | plot step (optional S2_plot, login node) | both PNGs + `phonon_band_summary.json`, bands 0-21 THz |
| Si 128 atoms, 90 frames, fc2 | pheasy OLS, c2 = 5 A, submitted through `autozt start` to jzzn | ran on cu18, gate **stable**, but 80 % fit error -- the cutoff is far too tight, see below |
| Si 128 atoms, 90 frames, fc2+fc3 | pheasy OLS, c2 = inf, c3 = 5 A, same dataset | 147 free IFCs and 69.5 % relative error -- **identical to the user's own bench run**; gate correctly reports **imaginary** |

The three-way judge was exercised in all three states: **stable** (cases
above), **imaginary** (the deliberately under-fitted pheasy run, mesh min
-6.6 THz, job exits 0, marker unsatisfied) and **tool error** (the phonopy
route was checked against an ASE-Atoms supercell, which made the gate exit
non-zero instead of reporting a false verdict).

### Cluster run (jzzn, through `autozt`)

A real submission, `autozt -tt fit-fc-thermal -p Si_fcfit -j step1_fit start`, on a
2-atom Si cell / 128-atom 4x4x4 supercell / 90-frame pkl dataset on the
shared filesystem: job `3828909` ran on `cu18` and produced
`fc2.hdf5`, `shengbte/FORCE_CONSTANTS_2ND`, `phonon_summary.json`
(`stable=true`, min -4.2e-08 THz) and `fc_fit_summary.json`, with
`pheasy` and `celer` both present in the environment. The whole path is
exercised: login-node gen -> push -> `sbatch` -> compute node
(prep/fit/post) -> judge -> fetch.

Two submissions on the same dataset bracket the gate's behaviour:

* `PHEASY_C2_CUTOFF = 5.0` gives Si a **3-cluster** cluster space (27 IFCs, 6
  free) — a severe truncation — and the fit reports 80 % relative error with
  correlation 0.60. The verdict is still `stable`, and the gate prints the
  right diagnosis (*ordinary model error -- widen the cutoff or fit fc3 as
  well*) instead of passing silently.
* With the production cutoffs (`c2` = inf, `c3` = 5 A, fc2+fc3) the skill
  builds exactly the cluster space the user's own bench builds (HARM 297 IFCs
  / 138 free, ANHARM3 81 / 9) and lands on **69.5 % relative error, 147 free
  IFCs, correlation 0.698** — the bench's own run reports 69.45 % and the
  same 147. The independent RMSE check over 4 frames agrees (70.8 %).
  With that much residual the resulting spectrum is genuinely unstable
  (min -1.51 THz) and the judge says so: `stable=false`, status `imaginary`.

So the pipeline, the order check (`is_pheasy_order: true`, identity
permutation on this dataset), the quality gate and the phonon judge all
behaved correctly; the ~70 % residual is a property of that dataset. Both my
run and the user's own log report `Space group: R-3m (166), 12 symmetry
operations` for a 2-atom Si cell (diamond Si is Fd-3m, 192 operations), which
is worth a look before that dataset is used for anything beyond timing.

Still not exercised:

* **`PHEASY_BIN = pheasy-gpu`** is exercised as of 2026-09-12 (Mg4C60,
  512-atom supercell / 128 frames, 3090, `--gres=gpu:4`): cluster space, null
  space and the sensing matrix all ran, and the fit engine was the GPU build —
  but the *cards were idle*, because LASSO builds a dense/CSR matrix and only
  the two-level path has a GPU implementation. A GPU request on that path
  (`PHEASY_GPU_LASSO_RESIDENT=1` without a two-level matrix) fails hard with
  `NotImplementedError` ~2 minutes into the fit, after a 10.8 GB
  `sm_dense.npy`; the driver now reconciles the two flags and records
  `pheasy_gpu_used`. A non-zero `PHEASY_OLS_RIDGE` also disables the resident
  GPU OLS path.
* **Least-squares `phono3py` with `FC_CALC = alm`** — only `symfc` was run
  (the local phono3py 3.24.0 / phonopy 2.47.1 pair matches the cluster).
* **A dask/multi-node `phono3py` fit** on a supercell much larger than 250
  atoms.
