# Configuring PhonoAgent for your own clusters

PhonoAgent ships with templates, not with a working site configuration. Cluster names,
ssh aliases, module systems, pseudopotential libraries and scratch directories are site
specific and must be filled in by you.

## 1. Where configuration lives

- Package defaults: setting/ in the repository
- User configuration: ~/.config/phonoagent/
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
    phonoagent -c my.yaml -tt te-screen -p Si_demo start

te-screen is a pure-python surrogate screener and needs no cluster. Everything else
needs a cluster entry, a pseudopotential library and (for the MACE skills) a MACE model
that you obtain yourself.

## 5. Switching a material between clusters: the four traps

Switching a material with phonoagent ... hpc <cluster> only rewrites
project_setting/hpc.yaml. Three other things stay behind, and every one of them has
bitten a real run:

| # | Symptom | Cause | Fix |
|---|---|---|---|
| 1 | generation dies with mkdir: cannot create directory '/public': Permission denied | project_setting/setting.yaml still pins the old cluster's work_dir, and it outranks hpc.yaml | rewrite work_dir there as well; phonoagent status now prints "work_dir source: ..." |
| 2 | job sits in PD(PartitionConfig) forever, log says nothing | the submit template requests more resources than the machine provides (measured: cpus-per-task=48 on a 24-core server) | lower cpus-per-task in a project-level copy of the template; phonoagent warns when a cluster declares max_cpus |
| 3 | job sits in PD(PartitionConfig) with a partition the cluster does not know | template hardcodes the source cluster's partition name | copy the template to project_setting/templates and use a partition the target cluster accepts; phonoagent prints the partitions the cluster's own templates use |
| 4 | the fit/relax stage dies silently: no traceback, no products | the template activates a conda env or conda.sh path that only exists on the source cluster (e.g. atomate2_p_a and /public/home/... on a machine that has mace-gpu) | point CONDA_SH / CONDA_ENV in that step's step.conf at the target cluster's environment (values belong in the cluster's setting file) |

Trap 4 is the nastiest: the shell activation block fails, the script keeps going with
whatever python is on PATH, the preparation stage may still succeed, and the compute
stage exits without a message.
