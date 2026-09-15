# Quick start

## Without a cluster

The thermoelectric screener runs locally and needs only numpy:

    python -m pip install -e ".[yaml,mace]"
    phonoagent -c demo.yaml -tt te-screen -p Si_demo start

with demo.yaml containing at least:

    project_roots:
      - /path/to/materials
    task_types:
      te-screen:
        max_jobs: 1

The two generation steps write te_features.json (14 DFT-free features) and
te_screen_summary.json (predicted ZT_e and log10(PF)).

## With a cluster

1. Copy setting/template-cluster.yaml to setting/<your-cluster>.yaml and fill it in.
2. Create a project configuration with phonoagent -tt <skill> -p <material> init.
3. Generate inputs without submitting: phonoagent -tt <skill> -p <material> -j 2 init.
4. Inspect them, then submit: phonoagent -tt <skill> -p <material> -j 2 start.
5. Watch the pipeline: phonoagent -tt <skill> -p <material> status.

Step advancement, retries and stuck-job handling are described in README.md.
