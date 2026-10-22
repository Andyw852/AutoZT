# Changelog

All notable changes to AutoZT are documented here.
This project adheres to Semantic Versioning; versions before 1.0.0 are development snapshots.

## [1.2.0] - 2026-09-16

### Changed
- Renamed the software to AutoZT: package autozt, commands autozt (short form az), environment
  variables AUTOZT_*, configuration directory ~/.config/autozt. The previous command names
  (pa, phonoagent) remain as compatibility wrappers pointing at the same program, so existing
  scripts keep working.
- Version banner and metadata report AutoZT 1.2.0.

## [Unreleased]

### Changed (2026-10-21: [0065] ke-dft-cpu S6_elastic cheaper defaults; S8.1 no longer waits for S6)
- `S8.1_bt2` drops `step6_elastic` from its needs: it never reads the elastic constants (C only enters the DPT τ of
  `S8.2_dpt`, which S8.1 still needs, so semiconductors wait for S6 through S8.2 as before). Metals (dpt group off)
  no longer sit behind a multi-day IBRION=6 run.
- `S6_elastic` factory defaults: NFREE 4 -> 2, EDIFF 1E-7 -> 1E-6, KSPACING 0.02 -> 0.03. Mg2C60 with the old values
  was 144 SCF steps x ~40 min (~4 days) and IBRION=6 cannot restart. Only newly generated S6 inputs change; finished
  S6 runs are not recomputed. Per material, `[incar.final] NFREE = 4 / EDIFF = 1E-7` in step.conf restores the old
  accuracy.

### Changed (2026-10-21: ke-dft-cpu two-stage HSE band structure)
- New step `S2.28_hsescf` (`step2_bandgap/step2.28_hse_scf`): one HSE SCF without KPOINTS_OPT, with NKRED on the
  axes whose k mesh is even (`HSE_SCF_NKRED = auto|off|"a b c"` in step.conf), saving WAVECAR.
- `S2.3_hse` now depends on it: every slice restarts from that HSE WAVECAR (ISTART=1, no NKRED, KPOINTS and
  NBANDS taken from stage 1) and only computes its own path points. Previously each of the four slices redid the
  full HSE SCF from the PBE WAVECAR; for large cells (Mg4C60, 128 atoms) one slice took 40-50 h, beyond the 48 h
  wall clock.
- Materials whose S2.3 already finished the old way: the new `hse_scf` check marks S2.28 complete, so no SCF job is
  submitted. Regenerating S2.3 for them without a stage-1 directory falls back to the old per-slice SCF with a
  warning; `rerun S2.28_hsescf` switches them to the two-stage route (S2.3 is archived and regenerated).

### Added (2026-10-21: [0061-m] ADP_DEGEN strain operator from DFT wavefunctions, plan A1)
- `skill/ke-dft-cpu/deform_wave.py` (numpy + h5py): streams WAVECAR (standard / ISPIN=2 / non-collinear), extracts the
  band-edge k-points, an energy window of bands and coefficients up to 300 eV into `deform_wave.h5`, then builds
  `deformation_matrix.h5` in S7.1. For each k-point it projects every strained configuration's states onto the
  unstrained band basis (Löwdin) and takes the central difference: X^(P)_mn = ⟨m|∂H/∂ε_P|n⟩, i.e. the whole first-order
  strain Hamiltonian at fixed fractional coordinates, including interband elements and their signs.
  `StrainGrid` interpolates the operator itself trilinearly on the S7 mesh (gauge-invariant; time reversal fills the
  ISYM=0 half mesh, spinors included).
- The band window is fixed by the unstrained configuration only, and is widened over near-degenerate groups within
  the bands every strained configuration still covers (`MARGIN` 6, `PAD` 2). Per-configuration energy windows used to
  cascade-shrink the window down to the valence bands.
- S7 (`DEFORM_WAVE`, default on; `DEFORM_WAVE_KEEP`, `DEFORM_WAVE_MARGIN`): WAVECAR is written only where the
  eigenvalues come from — the ionrelax static stage, or the top-level single point when there is no ionrelax. It is
  extracted in the job and then deleted. The INCAR templates are unchanged. LWAVE is ignored when deciding whether
  an ionrelax result can be reused. A finished ionrelax without `deform_wave.h5` reruns only its static stage
  (archive suffix `nowave`). Extraction failures only warn.
- S7.1 runs `deform_wave.py build`; missing files or a failed build only warn. The signed h5 now stores `ibands_<spin>`.
- S8 / S8.4 (`ADP_DEGEN_MATRIX = auto|off`): link `deformation_matrix.h5` and ship `deform_wave.py`.
- Runtime behaviour:
  - Reference-level alignment is taken from non-degenerate signed potentials; the matrix is rejected when the
    spread is too large.
  - At the valley centre, Z = Bᴴ D̂(k0) B. This step checks the reciprocal-lattice representation by coverage,
    symmetrizes, and cross-checks against the signed potentials; on any failure it falls back to the eigenvalue fit.
  - Pairwise replacement uses the precomputed Zk = Bᴴ D̂(k) B. This already includes the kinematic term, so T is no
    longer subtracted.
  - With cached coefficients, a grid formula replaces the ADP matrix element for every final state of the carrier type
    that has a degenerate valley. It uses AMSET's coefficients × the DFT operator, interpolates the vertex bilinears
    with AMSET's own `get_cross_section_values`, self-checks |⟨f|i⟩|² against AMSET's cross-section overlap, and keeps
    k·p states within 3 DFT mesh spacings of the valley centre. Large arrays go through AMSET's shared memory.
- Fixed after the first Si run (2026-10-22): `build()` crashed on patterns with only one strain sign (IndexError, wrong
  index into the common-k eigenvalues). Such patterns (cubic shear after `amset deform create`'s symmetry reduction)
  now get the missing side from the space-group operation that flips the strain (wavefunctions transformed, central
  difference, < 1% on the EPM bench vs ~10% for a forward difference). `load_dft_nodes` set `efermi` from the file
  when pymatgen leaves it None with `parse_dos=False` (AMSET's `get_band_structure` raised TypeError, so ADP_DEGEN had
  silently fallen back to plain AMSET in real runs).
- `build()` now reconstructs every strain component AMSET's symmetry reduction left out, not only the missing sign:
  for cubic Si `amset deform create` computes only xx± and yz+, so `deformation_matrix.h5` had 2 of 6 patterns and the
  runtime rejected it ("缺图样"). yy/zz now come from xx and xy/xz from yz via space-group operations
  (strain a·E_P -> a·RᵀE_P R), each with both signs, still a central difference (< 1% on the EPM bench).
  k-points where neither side's symmetry image was extracted are stored as NaN and treated as missing grid points.
- The matrix path's check at the valley centre now compares sums over the degenerate subspace. Per-band values of
  S7.1's energy-sorted ±δ central differences carry an O(δ) second-order error at a degeneracy (EPM Si Γ, ε_xx:
  ±0.7 eV per band, 0.1% on the sum), so the old per-band check could reject a correct DFT matrix; the per-band
  residual is kept as a diagnostic. `tools/adp_k0_probe.py` replays both runtime checks from the three h5 files in seconds.
- No published precedent for feeding a DFT strain-Hamiltonian matrix into AMSET's tetrahedron integration; this is
  the project's own construction.
- Tests: `tests/test_0061m_deform_wave.py`, `test_0061m_matrix_mode.py`, `test_0061m_grid_formula.py` and
  `test_0061m_wiring.py`. They use an empirical-pseudopotential Si test bench (`tests/_epm_si.py`), where ∂H/∂ε is
  known exactly. Results:
  - Projection is within 1% of the exact value.
  - Off-grid interpolation error falls at second order with mesh refinement.
  - Pairwise and grid paths are bitwise equal to the vertex sums.
  - The S7/S7.1/S8/S8.4 wiring is covered.

### Changed (2026-10-21: [0061-p] ADP_DEGEN replaces matrix elements pairwise inside AMSET's integration)
- Default (`ADP_DEGEN_PAIR = auto`): instead of scaling each initial state's ADP rate by a model-computed ratio, the
  plugins hook three intermediates of AMSET's `calculate_rate` (identical in 0.4.19 and 0.5.1): `tbs.get_masks` (final
  band per tetrahedron), the cross-section overlap from `get_cross_section_values`, and `get_fine_mesh_qpoints` (q-point
  → tetrahedron mapping). For (initial, final) pairs inside the degenerate subspace, ADP's `factor()` returns
  F_full / AMSET's own overlap at that point. AMSET multiplies by exactly that overlap, so the full matrix element
  remains. Tetrahedron weights, momentum-relaxation factors, DOS and the final-state set are AMSET's own.
- Falls back to the per-initial-state ratio when the source structure does not match, or with `ADP_DEGEN_PAIR = off`.
  `adp_degen_fit.json` records the `method`.
- Fixed: the ratio table's m_e v⊗v term assumed Å; AMSET's runtime reciprocal lattice is in bohr (`kunit`).
- Tests: pairwise and ratio paths agree to 1e-8 on the same final states (3D Si, 2D doublet), and a sign error in the
  kinetic term breaks it. An AMSET-shaped `calculate_rate` checks the hooks, the exact overlap cancellation and the
  untouched out-of-subspace states.

### Changed (2026-10-21: [0061-g] ADP_DEGEN generalized; 2D support)
- The ADP interband correction no longer assumes a cubic Γ triplet. `adp_degen.py` now finds band-edge valleys
  (local extrema of the edge band on the full DFT mesh, grouped into stars by point group + time reversal; holes and
  electrons), takes the bands within 50 meV at the valley centre as the subspace, and reads everything from the
  material's own data:
  - k·p matrices from projecting AMSET's wavefunction coefficients onto the k0 subspace (Löwdin; gauge-invariant).
  - The little-group representation matrices are solved from those matrices (unitary and antiunitary ops, closed into a
    group), and the subspace is split into common invariant blocks when it decouples.
  - The k·p model and strain operator are fitted in the symmetry-allowed subspace. The k·p model is refined against
    AMSET's dense bands. The strain operator reproduces S7.1's ±δ sorted central differences at k0 exactly and uses
    the non-degenerate signed potentials near k0.
  - The ratio table is energy-dependent, with 3 acoustic modes in 3D and 2 in-plane modes in 2D.
- Works for any crystal system, non-Γ valley stars, SOC spinor manifolds (e.g. Γ8 + Γ7) and near-degenerate pairs.
  Kramers pairs and non-degenerate valleys are left unchanged. When the S7 mesh cannot pin a sign, solutions whose
  factors agree within 10% are averaged; otherwise the valley is skipped with a reason.
- S8.4 (2D): `ADP_DEGEN` wired into `amset2d_plugin.AcousticDeformation2D` and `gen_step14_amset2d.py`. The signed file
  matches the chosen h5 (`deformation_vac_signed.h5` for the vacuum reference). It needs real overlaps
  (`UNITY_OVERLAP = false`). The S8.3 summary covers both runs.
- S8 no longer refuses SOC band structures.
- The rate ratio now mirrors AMSET's elastic-rate formula on both sides: the momentum-relaxation factor
  1 − v·v′/|v|², and the AMSET form |D(fixed fractional k)| + m_e v⊗v. The full element includes the same
  intraband kinetic terms. Only the interband strain coupling differs. Si factors change by a few % (0.1 eV:
  1.81/2.78/2.04 vs 1.93/3.21/2.26 before). `amset_terms=False` keeps the first-version convention for tests.
- Validated on synthetic DFT data (known k·p + strain operator, random wavefunction gauges): cubic, SOC, 2D hexagonal
  Γ/K, orthorhombic near-degenerate, 2D p_z + p_xy blocks, remote-band coupling (≤ 4%). Matches the first cubic
  version to machine precision. Tests: `tests/test_0061_adp_degen.py` (20 tests, pure numpy).

### Added (2026-10-21: [0061] ADP interband correction for degenerate valence-band tops)
- AMSET's ADP matrix element (diagonal deformation potential of the initial band × wavefunction overlap) misses the shear
  coupling between heavy and light holes of a cubic Γ-point p-like triplet (Si/Ge/diamond/zincblende holes, non-SOC).
  `skill/ke-dft-cpu/adp_degen.py` fits a 3-band k·p (to AMSET's own bands) and Bir–Pikus l, m, n (from S7.1's signed/abs
  deformation potentials at Γ; sign of n from nearby non-degenerate k), and scales each initial state's ADP rate by
  full-matrix-element / AMSET-form (Si, relaxed ions: ADP-only mobility 1.6–2.1× too high without it).
- `amset_dpo3d.py` applies the factor in the ADP class (after the ir_fix star average); S8 `ADP_DEGEN = auto|on|off`,
  report `adp_degen_fit.json`, S8.3 summary line. Warns when S7.1 `deform_geometry` is not ionrelax (Si clamped-ion shear
  potential is ~1/4 of relaxed; the 2026-10-19 Si hole ADP-only 7610 matches the clamped-ion AMSET form).
- Conventions checked against local QE 7.4 strained Si (k·p L/M/N, b, |d|, sign of d). Tests: `tests/test_0061_adp_degen.py`.

### Fixed (2026-10-20: [0059-r3] DFPT intravalley ODP)
- S8 / S8.4 put the intravalley ODP into the rates twice when both DFPT (`ivs_dfpt_params.json`) and frozen-phonon/manual
  parameters existed (e.g. Si holes 7.77 + 6.0). DFPT intravalley now wins; frozen/manual intravalley channels are reference
  only (`odp_reference`), frozen is recorded as `superseded_dfpt` and no longer triggers a re-run.
- eph-qe-cpu S13: the 3 acoustic modes at q=Γ are excluded explicitly (QE gives 0.6–0.8 meV, above the 0.5 meV threshold).
- eph-qe-cpu S10: single-valley carriers (Si holes at Γ) are kept for ODP when `IVS_ODP` is on; electron/hole sharing a k
  share one nscf k-point.

### Added (2026-10-20: [0059-r3])
- eph-qe-cpu S12: per-q AHC window widened from S11's valley bands when the band-edge degeneracy exceeds `IVS_AHC_NB`
  (Si Γ holes → 3); other q's are untouched. A stale S11 (no valley bands for the new initial k) is an error.
- `IVS_ODP_POLAR = exclude|include` (+ `IVS_ODP_POLAR_TOL`): polar modes at Γ are judged from Born charges computed by
  `epsil` in the q-odp DFPT run (compounds only), consistent with frozen-phonon `ODP_POLAR`.
- S8.3 summary prints an `[ODP 参数来源]` line. Local QE 7.4 check on Si: hole Γ ODP D = 7.89 eV/Å @ 63.9 meV vs frozen 7.77.
- Tests: `test_r3_*` in `tests/test_0059_ivs_dfpt.py`.

### Fixed (2026-10-19: [0059-r2] 打印/找谷口径)
- `postprocess_intrinsic.py` printed "intrinsic mobility" from `mobility[0][t300]`, i.e. the **first doping row** (−1e21),
  not the intrinsic point → now prints the **lowest non-zero doping** row ("本征近似"), Si = 1459 cm²/Vs at n = 1e10.
  (The exact 0-doping row stays unusable in AMSET and is documented as such.)
- S8.3 `[AMSET 机制]` summary printed AMSET's "inactive mechanism" mobilities (1e6–1e9, meaning scattering rate ≈ 0)
  as real numbers → values above `INACTIVE_MU_CM2_VS = 1e5` now print "—" and never become the limiting mechanism.
- `odp_frozen.find_valleys_3d` merged time-reversal partners into one plateau: on an 11³ grid Si's ±Δ valleys
  (t = 0.909 / 1.091) are adjacent and degenerate, so 6 valleys became 3 (V1×3). The plateau flood fill now checks
  each edge and refuses to absorb q ≡ −c (mod 1).
- Tests: `tests/test_0060_ivs_fixes.py` (3 tests); `_SKILL_REV` bumped to 2026-10-19-0060, postprocess script v5.

### Added (2026-10-19: IVS-DFPT Si cluster run and acceptance)
- `docs/RESULT-ivs-dfpt-Si-20261019.md`: first cluster run of the `ivs_dfpt` group on Si (hfeshell, QE 7.4, SG15 NC) and the
  ke-dft-cpu S8 integration on jzzn (AMSET 0.5.1). Per-pair D / ħω, star channels, comparison with the Si intervalley
  deformation potentials compiled in Pop's thesis (Table 3.2; per-pair D compared directly with ∆), and the acceptance
  number: electron μ(300 K) at n = 1e10 cm⁻³ = 1459 cm²/Vs (ADP-only 3883), consistent with the experimental
  lattice-limited ~1400–1450; IVS (2747) is the limiting phonon mechanism. Selection rules (TA/LA = 0) reproduce the
  local smoke run.
- Known issues recorded for follow-up: (1) `postprocess_intrinsic.py` prints `mobility[0][t300]`, whose index 0 is the
  -1e21 doping row rather than the intrinsic point; (2) the exact 0-doping row is unusable in AMSET (all mechanism
  subsets return the same 7610); (3) the IVS single-mechanism mobility is doping dependent (2747 @1e10 → 483 @1e21);
  (4) summaries print "inactive" mechanism values (1e6–1e9) as if they were physical; (5) S8.21's 3D valley finder
  merges Si's ±Δ valleys into one plateau on an 11³ grid (V1×3 instead of 6).

### Fixed (2026-10-19: [0059-r1])
- `ivs_dfpt.jsonable` rounded floats to 10 decimals, so `tr2_ph = 1e-16` was recorded as 0.0 in ivs_qlist.json /
  qrun.json / ivs_dfpt_params.json (metadata only; ph.x inputs were correct). It now keeps 12 significant digits;
  S12 regeneration treats an old 0.0 record as "not recorded", so finished q directories are not archived.
- RESULT-ivs-dfpt-Si: the literature comparison used D = ∆·√(Z_f/2), which doubles the rate; per-pair D equals ∆.
  The "2381 → ~1400" baseline came from the −1e21 row printed by postprocess_intrinsic and is replaced by the
  experimental lattice-limited mobility. Conclusions unchanged.

### Added (2026-10-18: [0059] DFPT intervalley scattering)
- eph-qe-cpu optional group `ivs_dfpt` (S10_ivsq, S11_ivsscf, S12_ivsph fan-out `q-*`, S13_ivs): zero-order intervalley
  deformation potentials and phonon energies from QE single-q DFPT plus `ph.x electron_phonon='ahc'` matrix elements at
  the valley k and k+q; valley pairs reduced by point group, time reversal and the reverse process; output
  `ivs_dfpt_params.json` (schema `ivs_dfpt/1`) folded into dp_optical star channels. Shared math in
  `skill/_common/ivs/ivs_dfpt.py`. Generic submit template `submit_ivs_qe.tpl`; QE paths from the cluster yaml
  (`qe_bin`, `qe_pseudo_dir`, `qe_setup`, injected as `QE_BIN` / `QE_PSEUDO_DIR` / `QE_SETUP`).
- ke-dft-cpu S8 (3D) and S8.4 (2D): `IVS_DFPT = auto|on|off` (+ `IVS_DFPT_DMIN`, `IVS_DFPT_KTOL`). DFPT channels take
  priority over frozen-phonon K↔K′ and hand-written intervalley channels (those are kept as references); stars are
  matched by k centres; the record `dp_optical_ivs_dfpt` goes into `2d_correction.json` / `dp_optical_3d.json`.
  S8.3 writes the IVS parameter source into comparison_summary.txt.

### Changed (2026-10-18: [0059])
- eph-qe-cpu: the Perturbo chain S0–S6 is now the default-on group `perturbo` (`steps: []`); materials that only need the
  intervalley chain set `perturbo: false`. Wannier90/Perturbo moved to `optional_exe`.
- Skill loader accepts an empty main step list when an optional group provides steps.

### Added (2026-10-17: ke-dft-cpu metal mode)
- S7.1: a metallic undeformed band structure no longer fails the step. It writes `metal.json`, skips
  `band_edges.json` / `deformation_vac.h5`, and still moves `deformation.h5` into the step directory.
- S8.1: metal mode. Skips the band-type gate, the scissor and the DPT τ. Reports S, σ/τ, κ_e/τ and the Lorenz
  number at the charge-neutral E_F(T); absolute σ, ρ, κ_e only when step.conf gives `METAL_LAMBDA_TR` (+ source;
  high-temperature limit ħ/τ = 2πλ_tr k_B T, Allen 1978) or an assumed `METAL_TAU_FS`. Writes `metal_transport.json`.
- S8.2 / S8.3 stop with "not applicable to metals" and list the groups to switch off in tf_<project>.yaml.
- METHODOLOGY: section on how metals are treated without experimental data.

### Fixed (2026-10-17)
- `cores: N` normalisation also keeps NCORE=1 for finite-difference INCARs (IBRION=5/6): Mg4C60 S6 elastic
  failed with "only possible if NPAR=number of nodes".

### Fixed (2026-10-17: DFPT NCORE, S7/S7.1 AMSET interpreter, relax watchdog)
- `cores: N` normalisation no longer rewrites NCORE in linear-response INCARs (IBRION=7/8 or
  LEPSILON=.TRUE.): VASP only accepts NCORE=1 there (Mg4C60 S5: template NCORE=1 became 4 and the job failed).
- ke S7 / S7.1: after `conda activate` the scripts call `$CONDA_PREFIX/bin/python` / `$CONDA_PREFIX/bin/amset`
  instead of whatever `python3` / `amset` comes first on PATH, and log the environment they used
  (Mg2C60 / Mg4C60 S7.1 failed with "No module named 'numpy'" although the amset env has numpy).
- run_relax.sh watchdog: `STALL_MINUTES = 0` means "off" again (since FREEZE_MIN was added the stall branch
  used 0 as its threshold and killed the segment after one minute); `FREEZE_MINUTES` can be set in step.conf.
  The FREEZE_MIN default goes from 120 to 240 minutes so it no longer fires before STALL_MIN (180), which was
  raised on 2026-09-17 because file-system pauses of 60-180 minutes killed healthy runs.

### Fixed (2026-10-15: fail counter counted state flapping as retries; conf/dir polluted state)
- autozt/progress.py: _history_stats now counts a failure only when the FAIL came from an actual
  attempt (previous state R / PD / TODO / PREP). OK->FAIL and WAIT->FAIL are flapping (the remote step
  directory is briefly invisible, an upstream file was refreshed, or a read-only command did not collect
  the material) and used to inflate the count until the step was promoted to human_review and stopped
  being advanced. Measured on Si_hanhai/ke-dft-cpu/S2.15_discr: 71 -> 2.
- autozt/collect.py: the "not collected" placeholder from conf / dir (CONFIG_ONLY) is marked
  config_only; autozt/history.py and autozt/progress.py skip it. Previously every conf/dir call
  recorded a WAIT transition for every step of every material and the next real collect flipped it
  back to FAIL - corrupting both the failure count and the progress view.
- Regression tests: tests/test_ai_integration.py::test_fail_count_counts_attempts_not_flapping and
  ::test_config_only_round_is_not_recorded_as_state.
### Fixed / Added (2026-10-15: ke S7.2 frozen-phonon ODP follow-up to the kz=1 change)
- S7.2 `ODP_KZ_VACUUM` default is now `auto`: if the step already has single points that ran (or are
  running) on a kz>1 grid it behaves as `keep`, so regenerating an in-flight material no longer changes
  KPOINTS and archives finished points; fresh materials get `reduce`. Explicit `reduce` / `keep` still
  win. The manifest records requested vs effective mode and the reason.
- S7.2 `ODP_KZ_VACUUM=reduce` now records what actually happened: only automatic-mesh KPOINTS are
  reduced; an explicit k-point list is left alone and the manifest no longer claims `kz_vacuum_reduced`.
  The vacuum axis is taken from the structure instead of assuming the third component; the manifest has a
  `kz_vacuum` record (mode, axis, status, S7 / ODP grids). Shared parser: `odp_frozen.kpoints_mesh` /
  `reduce_vacuum_kz`.
- S8.21 accepts the "in-plane equal, vacuum = 1" KPOINTS difference only when the manifest says the grid
  was really reduced, and compares along the recorded vacuum axis.
- S8.21 checks that every displaced configuration stays on the same magnetic solution as odp-ref (total
  moment from OSZICAR, per-site moments from OUTCAR, 0.1 uB tolerance) and stops otherwise; a flipped or
  high/low-spin configuration would otherwise turn into a spurious large D. Skipped when not spin-polarised.

### Fixed (2026-10-15: ke S7.1 vacuum alignment with a missing ionrelax LOCPOT)
- S7.1: when a configuration's eigenvalues come from `ionrelax/vasprun.xml`, its LOCPOT (and POSCAR) are
  read only from `ionrelax/`; a missing `ionrelax/LOCPOT` raises `RuntimeError` ("retry S7_deform to add
  the static segment") instead of silently using the rigid configuration's LOCPOT. The same rule applies
  to the embedded vacuum-h5 script.
- S7.1: if `deformation_vac.h5` cannot be built for a 2D material (`DEFORM_REF=auto`) or with explicit
  `DEFORM_REF=vacuum`, the step fails instead of falling back to the core reference. `DEFORM_REF` can now be
  set in step.conf (reserved key; write `DEFORM_REF = core` in the material-wide step.conf to opt out).
- S7: `retry S7_deform` re-runs only the ionrelax static segment for 2D vacuum-reference materials whose
  ionrelax finished but has no LOCPOT and no `deformation_vac.h5` yet: the static OUTCAR/vasprun.xml/OSZICAR
  are renamed `*.nolocpot-<time>`, submit.sh is rewritten to run the static segment only (from
  `CONTCAR.relax`, `ISTART=0`, `ICHARG=2`); the rigid single point and the relaxation are not repeated.
- `autozt graph` text view: a WAIT step now lists its already-finished non-main dependencies
  (the S8.22 -> S8.4 edge was missing from the text graph although the JSON had it).

### Fixed (2026-10-14: review follow-up to the band/functional settings)
- SOC impact gate: reports `inconclusive` (warning) instead of `ok` when either run's real VBM/CBM
  is not on the k points shared by the S2.2 SCF mesh and the S3 grid (PbTe without L, Bi2Te3).
- Hybrid vs semilocal m* check: the semilocal bands are fitted on the k points the downsampled
  hybrid path kept, with the same window; no more false "valley order changed" alarms for edges
  between hybrid path points, and no ratio bias for non-parabolic bands.
- FUNC=auto: autozt now pushes the local `result/step1*/workflow_method.txt` of both ke and kl into
  `.autozt_hints/` before S1 gen (`push_paths` entries may be `{from, to}`; mirrored, deleted when
  absent), so `rerun S1`, a cluster switch, or ke/kl on different clusters keep / align the
  functional instead of silently re-deciding. FUNC_REASON keeps the original reason across retries,
  names the script default correctly, a blank `FUNC_ALIGN` means true, and backup `step1*.xxx`
  directories are ignored.
- zt S20: falls back to kl S1's workflow_method.txt (when kappa_L comes from kl-dft-cpu) if
  `kl_method.json` is missing; a non-object `kl_method.json` no longer crashes.
- kl S5.1: the band cache fingerprint includes the segment end points and `PLOT_BAND_POINTS`, so the
  60-degree K fix and a new point count take effect on `retry`; a bad step.conf value is reported
  instead of falling back to defaults; a seekpath failure is no longer reported as "missing seekpath".
- ke S8/S8.4: the scissor-gap source label comes from band_summary.json `hybrid.label`
  (HSE06 / HSEsol / HSE03 / PBE0-type / custom) instead of "HSE06" for every hybrid run.
  Docs: changing `HYBRID` on a finished material needs `rerun S2.3`.

### Fixed (2026-10-14: ke/kl band and functional settings)
- kl-dft-cpu / phonon-dft-cpu / kl-mlff / phonon-mlff 2D phonon path: on a hexagonal primitive
  cell with gamma = 60 deg the "K" point (1/3, 1/3) sat on the Gamma-M' line at 0.577 |K|; it is
  now chosen from the in-plane angle ((2/3, 1/3) for 60 deg, (1/3, 1/3) for 120 deg).
- Shared fitting engine (`fc_fit_driver._stability_gate`, used by kl-dft-cpu S5 and fit-fc-thermal)
  and fit-fc-thermal S2_plot: 2D band yaml/plots now use the kz=0 2D path (`za_2d.band_path_2d`,
  new pool copy kept identical to the skill copies by tests) instead of seekpath, which added
  Gamma-A segments along the vacuum direction. S2_plot records `band_path` in its summary.
- kl-dft-cpu lattice_kappa (reference engine) and phonon-mlff-cpu S5 cutoff_conv: 2D phonon bands use
  the kz=0 2D path; cutoff_conv also uses the same explicit 2D q-mesh as S3 (vacuum 1, in-plane 60)
  so its screening matches S3. lattice_kappa's pymatgen fallback now refuses non-standard primitive
  cells (gamma = 60 deg put K at 0.58 |K|) and imports json (S5 phonon_summary.json was never read).
- 2D band paths (all four band_path_2d copies) are returned as one connected path; as separate
  two-point segments phonopy silently dropped the Gamma/M/K labels.
- phonon-mlff-cpu S3 referenced an undefined `dim` (NameError before fitting); supercell size is now
  |det(scm)| (np.prod of the matrix was 0, so the underdetermined-fit check always fired).
  phonon-mlff cpu/gpu S3: `vac_axis or 2` turned a vacuum along a (index 0) into c.
- relax_common 2D layer thickness (`_layer_factor`, `write_material_thickness`): same `vac_axis or 2`
  bug (vacuum along a computed as c; thickness_2d.json recorded axis 0 but used c). Latent in the
  S1 flow, which rejects or rotates a non-c vacuum before these run.
- kl-dft-cpu S5.1: no longer falls back silently to a hexagonal Gamma-M-K-Gamma path when seekpath
  is missing (3D) or the 2D path fails; the step now errors. seekpath is listed under
  `optional_python` (only 3D S5.1 needs it; as a hard requirement it blocked every kl step on the local backend).
  New S5.1 key `PLOT_BAND_POINTS` (default 31) replaces the environment-only setting.
- ke-dft-cpu global `FUNC` changed from a hard-coded `pbesol` to `auto` (same as kl-dft-cpu).
  `FUNC=auto` in relax_common now resolves in order: this skill's existing step1 FUNC (retry keeps
  the functional) -> the same material's paired skill (ke <-> kl, `FUNC_ALIGN`, default true) ->
  method_select bond-network rule. The reason is written to workflow_method.txt (`FUNC_REASON=`).
  Explicit hybrid FUNC values are rejected with a pointer to S2.3 `HYBRID`.
- zt-dft-cpu S20: checks that kappa_e (ke) and kappa_L (kl) come from the same functional
  (`functional_check` in zt_summary.json); kl-dft-cpu S6 writes `kl_method.json` for it.

### Added (2026-10-14)
- ke-dft-cpu S2.3 `HYBRID = auto | hse06 | hsesol | hse03 | pbe0 | custom` (custom:
  `HYBRID_GGA`, `HYBRID_AEXX`, `HYBRID_HFSCREEN`); unscreened choices get HFRCUT=-1 automatically.
  The keys are reserved, so they can sit in the material-wide step.conf.
- ke-dft-cpu S8/S8.4 warning-only band checks (`band_consistency.json`): SOC impact gate (S2.2 SOC
  SCF levels vs S3 non-SOC levels on shared k points: edge splitting > 0.05 eV, VBM/CBM k change)
  and hybrid-vs-semilocal band-edge effective-mass comparison (ratio outside 1 +/- 0.25).
  The transport chain still runs without SOC by default; nothing is switched automatically.

### Added
- Supercell symmetry guard (`SUPERCELL_SYMMETRY = strict | warn | off`, default strict):
  the displacement generators (kl-dft-cpu S4, phonon-dft-cpu S2, phonon-mlff S2, kl-mlff S2,
  zt-dft-cpu SK4), the fit-fc-thermal S1 fit and S2_kappa stop when the supercell does not
  keep every point-group operation of the crystal, and suggest the smallest symmetric
  diagonal supercell.  Mn2In2Se5 (R-3m) on a 3x3x1 supercell kept 4 of 12 operations; every
  fitting engine then produced a symmetry-broken fc2 (spurious zone-interior imaginary
  modes, kappa jumping by 20-30 % between q meshes).  `symmetry_audit.supercell_symmetry_gate`.
- fit-fc-thermal S1 gate: fc2 symmetry check (`FC2_SYM_TOL`, frequency spread between
  symmetry-equivalent q points; new status `fc2_symmetry_broken`, shown by autozt as
  "fc2 breaks crystal symmetry") and a denser full mesh (`IMAG_MESH_LENGTH = 100`, no
  symmetry reduction) instead of `run_mesh(60.0)` with reduction.
- docs/mcp.md: the MCP interface (tool table, risk tiers, measured safety numbers).
- AI integration (from a 295-material production run where every problem was in discovery,
  name resolution or status retrieval, not in scheduling):
  - `.tf_progress.json` + `autozt progress` / `autozt agent progress` / MCP `get_progress`:
    per-material state, structured FAIL codes (`class`/`code`/`action`/`fail_count`),
    `since`/`eta_s`, per-project counts and monitor liveness, read from a local file with
    no collection and no ssh. Written by the monitor every round and by any collecting
    command; partial scopes merge instead of truncating the full view.
  - `autozt doctor` / `agent doctor` / MCP `doctor`: ssh-free configuration preflight
    (effective `max_jobs` per project, tuning-knob sources, blocked projects, cross-project
    duplicate names and shadowing by blocked projects, monitor/cron keep-alive, agent gate).
  - `--project NAME` scope on every command, agent op and MCP state tool; segments are
    pruned before collection.
  - Stable material id `<project>/<full name>` (`qualified_name`, `id` in agent JSON),
    accepted by `-p` and the `material` argument.
  - `auto on --dry-run`; `auto on` skips files already at the target value.
  - Structured diagnosis codes for segmented relaxations (`segment_zbrent`,
    `segment_watchdog`, `segment_interrupted`, `segments_incomplete`, `segments_exhausted`),
    pressure, walltime, OOM, disk-full, missing outputs, stale inputs, plus a coarse
    `diag_class` and an `action_hint` for parameter tuning.
  - Tuning knobs `collect_chunk`, `collect_workers`, `op_workers`, `init_workers`,
    `cache_ttl` can live in tf.yaml (environment variables still override).

### Changed
- `_converge_mesh` compares each mesh with the latest earlier mesh that is coarser along every
  axis (axes frozen below `MESH_CONV_MAX_LENGTH` exempt), not with the previous mesh: with
  `L *= 1.25` the long axis of an anisotropic cell only steps every few meshes, and two meshes
  with the same k_z division were declared converged without k_z ever being tested.
  `mesh_convergence.json` records `axis_counts`, `compared_with` and `frozen_axes`.
- Project-level `max_jobs` in `tf_*.yaml` is now enforced as a per-project cap in addition
  to the skill-level cap. It was previously ignored silently. **Check `autozt doctor`
  before restarting a monitor: a project that sets `max_jobs: 16` now really runs at 16.**
- The agent approval gate is on by default (`agent_gate: auto`): besides `AUTOZT_ACTOR`,
  known AI-agent environment markers and non-interactive stdin make destructive commands
  require an approval token. `agent_gate: off` / `AUTOZT_AGENT_GATE=off` restores the old
  behaviour; `AUTOZT_AGENT_STRICT=0` still disables tokens for your own scripts.
- Monitor output: one `[round] {json}` line per round plus at most 20 step changes; the full
  table is printed only with `monitor_table: true`.
- A basename that also exists in a blocked (stale) project no longer aborts the command when
  an unblocked material with that name exists; errors now name the blocking project and
  its reason, and list all rejected targets at once.
- Legacy skill-name warnings are aggregated by rename pair once there are more than ten.

### Fixed
- fit-fc-thermal S1 gate: the fc2 (dataset SPOSCAR order) is reordered into the supercell order
  phonopy rebuilds from the wrapped POSCAR; an atom at fractional coordinate 1.0 had its
  periodic images permuted, so the gate evaluated a mislabelled fc2.
- E2BIG during collection: the collector payload is sent on stdin instead of as a
  `--config64` argument (Linux caps one argument at 128 KiB). The sbatch guard uses the
  same transport.
- `autozt agent …` and `autozt mcp` did not work through the pip-installed `autozt`
  command (`agent` was rejected; `mcp` collected all state first). All entry points now
  share one router. Agent/MCP subprocess calls fall back to `python -m autozt` when
  `bin/autozt` is not installed.
- The monitor never wrote `history.jsonl` (`history_record` was not imported and the
  NameError was swallowed).
- A single failing monitor round (ssh timeout, `sys.exit` inside collection) no longer kills
  the daemon; the error is logged and the next round retries.
- `monitor -d` checks the pid file before scanning project configurations, so a cron
  keep-alive no longer spends minutes rescanning a slow 9p mount when the daemon is alive.
- Blocked-project conflict checks cache material discovery per project (previously
  re-scanned for every target × blocked config) and see three-level layouts.
- stdout is line-buffered when redirected, so `autozt … > log` shows progress.
- `skills`, `schema`, `skill` and untargeted `history` no longer scan every project
  configuration before running (1–2 minutes per call on a WSL 9p mount). This also speeds
  up `autozt agent skills/contract/research_plan` and MCP `list_skills`/`describe_skill`/
  `research_plan`, which all call `schema --json`. With an explicit `-p` target the
  blocked-project check still runs first.
- Loading project configurations with many blocked projects was slow again (>10 minutes on
  a 9p mount with ~30 blocked projects): the three-level blocked-material discovery used a
  `*/*/*/POSCAR` glob that listed the contents of every material directory. It now walks
  level by level and never descends into a material directory (≈18× fewer filesystem
  operations on a 30×295 synthetic layout, same materials found).
- `history_record` merges its state file across scopes instead of overwriting it, so two
  differently scoped monitors (or a monitor plus a scoped CLI call) no longer erase each
  other's baselines and lose transitions and FAIL counts.
- Monitor `[round]` lines carry `timing_s` (config check / collect / cache+history / fetch /
  hang check / advance) and the global `auto_advance` flag, so "nothing submitted" is not
  mistaken for "stuck".
- Monitor rounds do less local filesystem work (first real round on a 9p mount: 48.7 min,
  of which `cfg_check` 14.4 min and `collect` 32.3 min):
  - the configuration-change walk no longer descends into `step*` directories of materials
    and their skill subfolders (only computation outputs and fanout `disp-*` folders live
    there); ~5.5× fewer filesystem operations on a synthetic tree, same config files tracked;
  - segments sharing a `local_root` are discovered once per round;
  - `fill_local_dim` caches the POSCAR-based dimension guess by file stat across rounds;
  - `[round]` gains `collect_detail` (local resolve vs parallel ssh vs post-processing
    seconds, segment/material/ssh-call counts) and a separate `fill_dim` timing, so the
    local-9p and ssh parts of `collect` can be told apart.
- Agent calls no longer collect the whole skill for every read or action (on the 328-material
  production set each `inspect`/`cycle` or `act start -p X` took tens of minutes on 9p):
  - `agent snapshot/inspect/plan/propose/cycle/run` and the MCP `get_snapshot`/`inspect`/
    `cycle`/`propose_actions`/`apply_actions` tools take `source` (`auto` default,
    `progress`, `live`). `auto` reads `.tf_progress.json` when it is fresh (live monitor, or
    not older than `narrow_max_age`) and covers the requested scope, and falls back to live
    collection otherwise. Every result says which source it used and how old it is.
  - Execution groups actions by (action, skill, step) into one `act -p A,B,C` call: one
    collection and one shared max_jobs gate per group instead of one full collection per
    action. `start`/`retry` carry `--expect-state`, so each target is re-checked against the
    freshly collected state and skipped if it changed (a plan from the progress file cannot
    resubmit a finished step or `retry` away a running job). Each action reports an
    `outcome`: `submitted`, `deferred` (max_jobs reached), `skipped_stale`, `ok`, `failed`.
  - Commands with explicit `-p` targets (`start`, `retry`, `fetch`, `advance`, `list`, `json`,
    `diagnose`, `conf`, `dir`, `prove`) resolve and probe only the matching materials of a
    shared-layout segment when the progress file recorded a whole-skill collection within
    `narrow_max_age` (default 7200 s; 0 disables). The matcher is a superset of `-p`
    resolution, so ambiguity errors are unchanged. The max_jobs gate adds the running/queued
    jobs of the materials that were not collected, from the progress file plus a submission
    ledger (`.tf_submit_ledger.jsonl`) that covers jobs submitted after the last collection.
    `autozt doctor` shows per skill whether this fast path is available.
- `autozt agent setup` prints a ready-to-paste MCP server entry (absolute interpreter and
  entry paths, `AUTOZT_CONFIG`, `workflow` profile) and a short rules card for the model;
  `--rules` prints only the card.
- Progress views carry `fail_summary` (FAILs bucketed by suggested action with codes and
  ids), and jobs submitted after the last collection show as `PD`
  (`submitted_after_collect`) instead of being proposed again.
- max_jobs was not enforced for explicit `start -p X -j S` (the path every agent
  `start_step` takes), and `start -p A,B,C` built a fresh gate per material from the same
  pre-submission counts, so a batch could overshoot the cap. Both now share one gate; `-f`
  still lets a human exceed the cap on purpose. The gate also counts jobs of materials
  outside the collected scope (`--project`, skill-subdir `-p` pruning) when the progress
  file is fresh.
- `act --project P …` was classified as an unknown (destructive) command because
  `--project`'s value was taken for the command word.

### Fixed
- Cluster-switch self-checks, each backed by a real failure observed in testing:
  work_dir provenance is printed by status; resource requests above the cluster's
  declared max_cpus are flagged; a partition the target cluster does not use is
  reported with the partitions its own templates use; and conda environments or
  conda.sh paths that do not exist on the target cluster are flagged before
  submission (this was the cause of a silent stage death: activation failed, the
  preparation stage still passed and the fit stage exited without a message).
- start now prints the exact command needed to resubmit a FAIL step
  (autozt -tt SKILL -p MATERIAL -j STEP start -f) instead of claiming that all
  steps are running, finished or dependency-blocked.
- docs/CONFIGURING.md section 5 documents the four cluster-switch traps.

## [1.1.0] - 2026-09-15

### Added
- MCP interface: autozt mcp exposes the tool as a JSON-RPC 2.0 stdio server
  (initialize / tools/list / tools/call). The tool table holds 14 generic verbs
  (list_skills, list_materials, get_summary, get_status, get_json, conf_get, start_step,
  retry_step, fetch_results, conf_set, stop_step, rerun_step, clean_material,
  session_export) and is capped at 20: skills are discovered from skill.yaml, so adding a
  skill never adds a tool.
- MCP policy boundary: read-only tools call the CLI directly; mutating and destructive
  tools are routed through the action gateway (autozt act), so an agent session gets
  needs-approval instead of execution, and approvals only work from a real TTY.
- work_dir provenance: every status output names the configuration layer that supplied
  work_dir, and a failing generation prints that layer plus the hint that
  project_setting/setting.yaml outranks hpc.yaml after a cluster switch.
- tests/test_mcp.py: tool-table shape and cap, risk tiers, protocol handshake, and the
  gateway refusal for destructive calls.

### Fixed
- Status did not reveal which configuration layer supplied work_dir (a cluster switch
  could silently keep writing to the old cluster's path).

## [1.0.0] - 2026-09-15

First public release, derived from the taskflow code base (v1.0 refactor).

### Added
- Package layout autozt/ with the autozt CLI (zero third-party dependency core).
- Criterion-driven step engine: every step declares a physical convergence criterion and
  the engine advances, retries or reports accordingly.
- Provenance per step plus autozt prove --verify, and autozt session export for
  reproducible archival (manifest with per-file sha256, history timeline, replay hints).
- Action gateway (autozt act / approve) with risk tiers, audit log and TTL-bounded
  human approvals for agent-driven sessions.
- Unified core count control: one "cores" setting normalises submit scripts and INCAR
  (NCORE/KPAR) for every skill, whatever template names it uses.
- Unified pre-submission input check: required inputs are discovered from the remote step
  directory instead of per-skill hard-coded lists.
- Supercell specification accepts both diagonal notation ("4 4 4") and a general 3x3
  integer matrix ("2 1 0 -1 2 0 0 0 1", phonopy/phono3py --dim semantics).
- 18 skills covering VASP (band, elastic, defect, ke, kl, opt, phonon), MACE
  (kl/opt/phonon CPU and GPU, mlff), plus te-screen and unihamgnn.
- Drop-in skill interface with autozt schema --strict validation.

### Fixed
- vasp_ncl submit template loaded an Intel module environment while executing the
  locally built (AOCC/AOCL) binary, failing with a missing libscalapack; the template now
  uses the matching environment, plus a single-node MPI guard for nodes with a broken
  mlx5/UCX stack.
- defect-dft-cpu shipped ENCUT = auto which VASP rejects (it now uses the documented 370 eV).
- unihamgnn could never submit its graph step (non-VASP inputs were checked against the
  VASP default list).
