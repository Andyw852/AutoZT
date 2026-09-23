#!/usr/bin/env python3
import json
from pathlib import Path

root = Path.cwd()
manifest_candidates = [root/'step5_perturbo'/'eph_parameters.json',
                       root/'step4_qe2pert'/'eph_parameters.json',
                       root/'step0_preflight'/'eph_parameters.json']
manifest_path = next((path for path in manifest_candidates if path.is_file()), None)
if manifest_path is None:
    raise SystemExit('[错误] 上游工件缺失：eph_parameters.json')
manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
if not isinstance(manifest, dict):
    raise SystemExit('[错误] eph_parameters.json 根对象不是 JSON 映射。')
prefix = str(manifest.get('prefix', '')).strip()
if not prefix:
    raise SystemExit('[错误] eph_parameters.json 缺少 prefix。')
checks = {
    'preflight': root/'step0_preflight'/'qe_env.json',
    'scf': root/'step1_scf'/'scf.ok',
    'wannier': root/'step2_wannier'/'wannier.ok',
    'phonon': root/'step3_phonon'/'phonon.ok',
    'qe2pert': root/'step4_qe2pert'/'eph.h5',
    'perturbo': root/'step5_perturbo'/'perturbo.ok',
    'ephmat': root/'step5_perturbo'/'ephmat',
    'ephmat_yml': root/'step5_perturbo'/'ephmat.yml',
}
missing = [k for k, p in checks.items() if not p.is_file() or p.stat().st_size == 0]
if missing:
    raise SystemExit('[错误] 上游工件缺失：' + ', '.join(missing))
env = json.loads(checks['preflight'].read_text(encoding='utf-8'))
params = manifest
result = {
    'workflow': 'QE-Wannier90-Perturbo electron-phonon matrix-element workflow',
    'prefix': prefix, 'status': 'complete', 'environment': env,
    'artifacts': {key: str(path) for key, path in checks.items()},
    'parameters': params,
    'convergence_scope': {
        'official_reference': 'Perturbo 3.0.0 tests/epr_computation/epr9',
        'requires_wannier_disentanglement_converged': True,
        'requires_grid_cutoff_smearing_comparison': True,
    },
    'scientific_scope': (
        'This validates the configured e-ph data plumbing and interpolation only; '
        'lambda, alpha2F, omega_log and Tc require a separately converged metallic '
        'workflow and must not be inferred from this summary.'
    )
}
(root/'step6_summary').mkdir(exist_ok=True)
(root/'step6_summary'/'eph_summary.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print('[DONE] eph_summary.json written')
