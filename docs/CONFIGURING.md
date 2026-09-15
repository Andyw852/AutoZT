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
