#!/usr/bin/env python3
import json
from pathlib import Path

root = Path.cwd()
checks = {
    'preflight': root/'step0_preflight'/'qe_env.json',
    'scf': root/'step1_scf'/'scf.ok',
    'wannier': root/'step2_wannier'/'wannier.ok',
    'phonon': root/'step3_phonon'/'phonon.ok',
    'qe2pert': root/'step4_qe2pert'/'si_epr.h5',
    'perturbo': root/'step5_perturbo'/'perturbo.ok',
}
missing = [k for k, p in checks.items() if not p.is_file() or p.stat().st_size == 0]
if missing:
    raise SystemExit('[错误] 上游工件缺失：' + ', '.join(missing))
env = json.loads(checks['preflight'].read_text(encoding='utf-8'))
params_path = root / 'step5_perturbo' / 'eph_parameters.json'
params = json.loads(params_path.read_text(encoding='utf-8')) if params_path.is_file() else {}
result = {
    'workflow': 'QE-Wannier90-Perturbo electron-phonon matrix-element smoke test',
    'prefix': 'si', 'status': 'complete', 'environment': env,
    'artifacts': {key: str(path) for key, path in checks.items()},
    'parameters': params,
    'convergence_scope': {
        'official_reference': 'Perturbo 3.0.0 tests/epr_computation/epr9',
        'requires_wannier_disentanglement_converged': True,
        'requires_grid_cutoff_smearing_comparison': True,
    },
    'scientific_scope': 'Si is a semiconductor: this validates e-ph data plumbing, not superconducting lambda or Tc.'
}
(root/'step6_summary').mkdir(exist_ok=True)
(root/'step6_summary'/'eph_summary.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print('[DONE] eph_summary.json written')
