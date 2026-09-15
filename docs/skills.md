# Skills

A skill is a directory under skill/ with a skill.yaml declaring its steps, generation
scripts, convergence criteria and default cluster. AutoZT ships 18 skills:

| Family | Skills |
|---|---|
| VASP (CPU) | band-dft-cpu, defect-dft-cpu, elastic-dft-cpu, ke-dft-cpu, kl-dft-cpu, opt-dft-cpu, phonon-dft-cpu |
| MACE | kl-mace-cpu, kl-mace-gpu, opt-mace-cpu, opt-mace-gpu, phonon-mace-cpu, phonon-mace-gpu, mlff-mace |
| Auxiliary | te-screen (thermoelectric surrogate screening), unihamgnn (graph data generation) |
| Fitting | fc-fit (phono3py / pheasy / hiphive force-constant fitting) |

Adding a skill means adding a directory; autozt schema --strict validates it against
the skill schema.
