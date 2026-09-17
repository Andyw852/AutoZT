# Configuring AutoZT for your own clusters

AutoZT ships with templates, not with a working site configuration. Cluster names,
ssh aliases, module systems, pseudopotential libraries and scratch directories are site
specific and must be filled in by you.

## 1. Where configuration lives

- Package defaults: setting/ in the repository
- User configuration: ~/.config/autozt/
- Per project: <material>/<skill>/project_setting/{setting.yaml,hpc.yaml,step.conf}

Precedence, highest first: project_setting/ > user configuration > package defaults.

## 2. Adding a cluster

Copy setting/template-cluster.yaml to setting/<your-cluster>.yaml and fill in: name,
ssh_host, work_dir, remote_path_prefix (the python used by generation scripts),
conda_sh / conda_env / amset_env, the vasp.<profile> binary paths, and templates/ for
your queue system.

## 3. Keeping site-specific paths out of the code

    grep -rn "/home/\|/public/home/" skill/ setting/ | grep -v template

## 4. Running without a cluster

    pip install -e ".[yaml,mace]"
    autozt -c my.yaml -tt te-screen -p Si_demo start

te-screen is a pure-python surrogate screener and needs no cluster. Everything else
needs a cluster entry, a pseudopotential library and (for the MACE skills) a MACE model
that you obtain yourself.

## 5. Switching a material between clusters: the four traps

Switching a material with autozt ... hpc <cluster> only rewrites
project_setting/hpc.yaml. Three other things stay behind, and every one of them has
bitten a real run:

| # | Symptom | Cause | Fix |
|---|---|---|---|
| 1 | generation dies with mkdir: cannot create directory '/public': Permission denied | project_setting/setting.yaml still pins the old cluster's work_dir, and it outranks hpc.yaml | rewrite work_dir there as well; autozt status now prints "work_dir source: ..." |
| 2 | job sits in PD(PartitionConfig) forever, log says nothing | the submit template requests more resources than the machine provides (measured: cpus-per-task=48 on a 24-core server) | lower cpus-per-task in a project-level copy of the template; autozt warns when a cluster declares max_cpus |
| 3 | job sits in PD(PartitionConfig) with a partition the cluster does not know | template hardcodes the source cluster's partition name | copy the template to project_setting/templates and use a partition the target cluster accepts; autozt prints the partitions the cluster's own templates use |
| 4 | the fit/relax stage dies silently: no traceback, no products | the template activates a conda env or conda.sh path that only exists on the source cluster (e.g. atomate2_p_a and /public/home/... on a machine that has mace-gpu) | point CONDA_SH / CONDA_ENV in that step's step.conf at the target cluster's environment (values belong in the cluster's setting file) |

Trap 4 is the nastiest: the shell activation block fails, the script keeps going with
whatever python is on PATH, the preparation stage may still succeed, and the compute
stage exits without a message.
## 6. Two-dimensional cells: the constrained-relaxation interface

A 2D slab relaxed with a plain `ISIF=3` will relax the vacuum axis too, and a slab
relaxed with `ISIF=2` keeps whatever in-plane lattice the input file had. Both give a
structure whose in-plane stress is non-zero, and that stress is not cosmetic: a
tensile in-plane strain adds a linear term to the ZA flexural branch, which breaks the
Huang zero-stress condition the 2D phonon/thermal-conductivity steps rely on. Measured
on a screening batch (2026-09-16, 10 structures, PBE-derived lattices relaxed with
PBEsol): every step-1 directory was `ISIF=2` with no constraint, and the in-plane
stress in the layer convention (sigma x h_perp / d) ranged from -2 to -62 kbar, all
tensile.

Declare how your cluster constrains the cell in `vasp.relax_2d.cell_constraint`:

| value | what it does | requirement |
|---|---|---|
| `ioptcell_tag` | keeps an `IOPTCELL` line in the INCAR | a VASP built with the IOPTCELL patch (the binary is `vasp.*-optcell`) |
| `optcell_file` | converts the `IOPTCELL` line into an `OPTCELL` file | same patched build |
| `lattice_constraints` | writes `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` | official VASP >= 6.5.0, and c parallel to cartesian z |
| `none` | neither | see below |

For `lattice_constraints` the generator reads the VASP version from the binary path
(e.g. `.../vasp.6.5.0/bin/vasp_std`). If the path carries no version, put
`export AUTOZT_VASP_VERSION=6.5.0` in that step's submit template. An unreadable or
too-old version is a hard error, never a silent fallback: older VASP ignores the tag
completely, so `ISIF=3` would relax the vacuum.

`none` is not a neutral default. On a 2D structure it means the in-plane lattice is
never relaxed, and the energy-versus-area scan that would replace it does not exist in
the pipeline. Since 2026-09-16 the 2D generator therefore **stops with an error** in
both cases:

- `cell_constraint: none`, and
- `ioptcell_tag` with no valid `IOPTCELL` line in the INCAR template (the usual cause
  is a project-level copy of `incar_2d.tpl` shadowing the skill template).

If you really do want a fixed cell, say so explicitly with `ALLOW_2D_FIXED_CELL = true`
in that step's `step.conf`. The generator then warns instead of stopping. Do not use it
for anything whose phonons or thermal conductivity you intend to report.

## 7. Overridden templates go stale silently

`project_setting/templates/` and `setting/<cluster>/templates/` outrank the skill's own
templates. A copy made months ago still renders, still submits, and still runs - it
simply lacks every placeholder the skill template has gained since.

Two measures exist:

- at generation time autozt compares the template it resolved against the skill's own
  copy and prints `a warning naming the missing {{PLACEHOLDERS}}`;
- for the placeholders whose absence changes the physics, the 2D kappa step refuses to
  render a submit script that lacks `kappa_2d_normalized`, and the relaxation steps
  force the force-precision tags and dipole-correction lines back into the INCAR.

When the warning appears, copy the skill template over the override (keep a backup) and
re-run `autozt ... retry` so the new file is pushed.
