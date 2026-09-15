# Installation

PhonoAgent's command-line core has no third-party dependency.

    git clone https://github.com/Andyw852/PhonoAgent
    cd PhonoAgent
    pip install -e .            # CLI only
    pip install -e ".[yaml]"    # + PyYAML (read yaml configuration)

Optional extras pull the scientific stack used by the generation scripts on the
computing side:

    pip install -e ".[vasp]"    # numpy, ase, pymatgen, phonopy, phono3py, spglib
    pip install -e ".[mace]"    # + mace-torch
    pip install -e ".[hiphive]" # hiphive + ase
    pip install -e ".[dev]"     # pytest, ruff, build

PhonoAgent does not redistribute pseudopotentials (POTCAR), MACE model weights or
third-party training data. Obtain them from their licensed sources and point the
configuration at them.
