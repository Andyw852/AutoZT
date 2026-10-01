# Quick start

## Without a cluster (local backend)

Everything below is typed in an ordinary terminal; no AI agent is involved.
`setting/local.yaml` runs jobs on this machine. If `sbatch` is not installed,
AutoZT puts the bundled `tools/fakeslurm` on PATH, which queues jobs on a CPU
budget (`FAKESLURM_MAX_CPUS`, default: all cores).

    pip install -e ".[yaml]"            # CLI; add the extras your skills need (see installation.md)
    autozt check                        # which skills can run here, and what is missing

`autozt check` probes every skill on this machine: Python packages, programs
(vasp_std, vaspkit, phono3py, ShengBTE...), the POTCAR library, MACE models and GPUs.
For each missing item it prints how to fix it. Use `--cluster jzzn` for a cluster
(over ssh), `--cluster all` for every configured machine, and `--json` for scripts.

Point `setting/local.yaml` at what you have:

    conda_sh: ""            # both empty = use the python on PATH (pip/venv users)
    conda_env: ""
    vasp: {standard: {cell_constraint: none, std: vasp_std, gam: vasp_gam, ncl: vasp_ncl}}
    potcar_dir: /path/to/potpaw_PBE      # VASP skills that read POTCAR_DIR
    mlff_model_dir: /path/to/mace_models # MLFF skills

Skills that build POTCAR with vaspkit read `PBE_PATH` from `~/.vaspkit`.

Then register a material and run it:

    autozt -tt fit-fc-thermal -p MoS2 register --dataset /data/MoS2/step4_disp --cluster local
    autozt -tt band-dft-cpu  -p Si   register --poscar /data/Si/POSCAR --cluster local
    autozt -tt band-dft-cpu -p Si start          # generates inputs, submits locally
    autozt -tt band-dft-cpu -p Si status         # or: autozt summary / autozt monitor -d

`start` refuses to submit on this machine when a hard dependency is missing,
and it prints the same hints as `check`; `-f` overrides this.
Generic submit templates for the local backend live in `setting/local/templates/`.
To adapt one for a single project, copy it to the project's
`project_setting/templates/`.

The thermoelectric screener needs only numpy:

    autozt -tt te-screen -p Si_demo register --poscar POSCAR --cluster local
    autozt -tt te-screen -p Si_demo start

## Using it from an AI client (MCP)

    autozt agent setup --write          # writes .mcp.json and .cursor/mcp.json in the AutoZT folder

Start the client (Claude Code, Cursor, ...) in the AutoZT folder and the `autozt`
MCP server is picked up automatically. Other clients: paste the `mcp_server` block
printed by `autozt agent setup`. The CLI behaves the same with or without an agent:
errors go to stderr with a nonzero exit code. Agents additionally get them on stdout.

## With a cluster

1. Copy setting/template-cluster.yaml to setting/<your-cluster>.yaml and fill it in.
2. Create a project configuration with autozt -tt <skill> -p <material> init.
3. Generate inputs without submitting: autozt -tt <skill> -p <material> -j 2 init.
4. Inspect them, then submit: autozt -tt <skill> -p <material> -j 2 start.
5. Watch the pipeline: autozt -tt <skill> -p <material> status.

Step advancement, retries and stuck-job handling are described in README.md.
