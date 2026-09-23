# Table 1 | AutoZT software and workflow overview

| Domain | Skill | Main function | Main packages / engines | Execution |
|:--|:--|:--|:--|:--:|
| Structure | `opt-dft-cpu` | DFT structure optimisation and total energy / formation energy | VASP | CPU |
| Structure | `opt-mlff-cpu` | MACE structure optimisation and formation energy | MACE, ASE | CPU |
| Structure | `opt-mlff-gpu` | MACE structure optimisation and formation energy | MACE, ASE | GPU |
| Electrons | `band-dft-cpu` | PBE/HSE band structures and k-path generation | VASP, VASPKIT, SeeK-path | CPU |
| Electrons | `ke-dft-cpu` | Electronic transport: conductivity, Seebeck coefficient and electronic thermal conductivity | AMSET, BoltzTraP2, CRTA check | CPU |
| Electrons | `unihamgnn` | ML Hamiltonian generation and band calculation | Uni-HamGNN | CPU / GPU |
| Electrons | `cohp-cogito` | Chemical-bonding and COHP/ICOHP analysis | COGITO, VASP | CPU |
| Electrons | `eph-qe-cpu` | Electron–phonon coupling and transport workflow | Quantum ESPRESSO, Perturbo | CPU |
| Defects & mechanics | `defect-dft-cpu` | Point defects, charge states and formation energies | VASP | CPU |
| Defects & mechanics | `elastic-dft-cpu` | Elastic constants and mechanical response | VASP | CPU |
| Phonons | `phonon-dft-cpu` | Harmonic phonons and dynamical stability | VASP, Phonopy | CPU |
| Phonons | `phonon-mlff-cpu` | MACE phonon spectrum and stability screening | MACE, Phonopy | CPU |
| Phonons | `phonon-mlff-gpu` | GPU-accelerated MACE phonon spectrum and stability screening | MACE, Phonopy | GPU |
| Force constants | `fc-fit` | Harmonic / anharmonic force-constant fitting | Phono3py, hiPhive, symfc, ALM | CPU / local |
| Thermal transport | `kl-dft-cpu` | Lattice thermal conductivity from DFT forces | Phono3py, ShengBTE, VASP | CPU |
| Thermal transport | `kl-mlff-cpu` | Lattice thermal conductivity from MACE forces | MACE, Phono3py | CPU |
| Thermal transport | `kl-mlff-gpu` | GPU-accelerated MACE force evaluation for lattice thermal conductivity | MACE, Phono3py | GPU |
| ML potentials | `mlff` | Random-displacement MLFF training and MACE potential generation | MACE, VASP | CPU / GPU |
| Thermoelectrics | `te-screen` | Fast surrogate screening of thermoelectric performance | Python | Local |
| Thermoelectrics | `zt-dft-cpu` | End-to-end DFT thermoelectric figure of merit (ZT) workflow | VASP, AMSET, BoltzTraP2, Phono3py | CPU |

*Notes:* The table lists the 20 current production skills in `skill/`; `_common` and `_template` are internal support directories and are excluded. CPU/GPU indicates the execution target supported by the skill, while Local means a login-node or workstation task. The shared structure stack is pymatgen, ASE, spglib and SeeK-path. The older labels `opt-mlip`, `phonon-mlip`, `kl-mlip`, `mlff-mlip` and `zt` map to the corresponding MACE or `zt-dft-cpu` entries above.
