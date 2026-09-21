"""Offline regression tests: actual discovery and preflight, no cluster access."""
from pathlib import Path
import shlex
import sys
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt
from autozt import workflow as wf

class FanoutPreflightTests(unittest.TestCase):
    def test_stale_archive_filter(self):
        for name in ('IBZKPT.stale-grid-15x15x1', 'REPORT.stale-old',
                     'WAVECAR.stale-grid-46x46x3', 'ionrelax/INCAR.stale-old'):
            self.assertFalse(wf._step_input_name_ok(name), name)
        self.assertTrue(wf._step_input_name_ok('ionrelax/INCAR.relax'))
        self.assertTrue(wf._step_input_name_ok('custom-input.dat'))

    def run_case(self, names, declared=None, missing=None, fetch=False):
        base = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'submit.sh']
        ir = ['ionrelax/' + x for x in ('INCAR.relax', 'INCAR.static', 'KPOINTS', 'POSCAR', 'POTCAR')]
        files = {'inplane': base + ir + ['custom-input.dat'], 'clamped': base}
        checked, commands, transferred = [], [], []
        available = {'/mirror/step/' + n + '/' + f for n in names for f in (declared or files[n])}
        if missing: available.discard('/mirror/step/' + missing)
        if fetch: available.clear()
        def ssh(cfg, host, args):
            return ['SSH_MOCK', args[0]]
        def run(cmd, **kwargs):
            if cmd[0] == 'tar':
                available.update('/mirror/step/' + f for f in transferred)
                return types.SimpleNamespace(returncode=0, stderr='')
            self.assertEqual(cmd[0], 'SSH_MOCK')
            command = cmd[1]; commands.append(command)
            self.assertIn("! -name '*.stale-*'", command)
            name = shlex.split(command)[1].split('/')[-1]
            # Include archives even though real remote find excludes them: also test Python defense.
            return types.SimpleNamespace(returncode=0, stdout='\n'.join(files[name] + ['REPORT.stale-old', 'IBZKPT.stale-grid-15x15x1']), stderr='')
        def popen(cmd, **kwargs):
            import io
            transferred.extend(shlex.split(cmd[1].split('tar --ignore-failed-read -cf - ', 1)[1]))
            return types.SimpleNamespace(stdout=io.BytesIO(), stderr=io.BytesIO(), wait=lambda: 0)
        def isfile(p):
            checked.append(str(p)); return str(p) in available
        cfg = {'task_types': {'test': {'submit_required': declared}}}
        m = {'tt': 'test', 'host_eff': 'mock', 'hpc_name': 'mock', 'result_dir': '/mirror'}
        s = {'name': 'step', 'dir': '/remote/step', 'fanout': '*', 'fan_todo': names}
        with patch.object(autozt, '_ssh_cmd', ssh), patch.object(autozt, '_load_yaml_file', return_value={}), \
             patch.object(autozt, 'run_remote', return_value=(0, '')), \
             patch('subprocess.run', side_effect=run), patch('subprocess.Popen', side_effect=popen), \
             patch('os.path.isfile', side_effect=isfile), patch('os.makedirs'):
            ok, reason = wf._remote_submit_preflight(cfg, m, s, {})
        self.assertTrue(ok, reason)
        self.assertFalse(any('.stale-' in p for p in checked + transferred))
        self.assertFalse(any('/clamped/ionrelax/' in p for p in checked + transferred))
        for n in names:
            for f in declared or files[n]: self.assertIn('/mirror/step/' + n + '/' + f, checked)
        if declared: self.assertEqual(commands, [])
        else: self.assertEqual(len(commands), len(names))

    def test_ionrelax_first(self): self.run_case(['inplane', 'clamped'])
    def test_clamped_first(self): self.run_case(['clamped', 'inplane'])
    def test_fetch_per_child_without_archives(self): self.run_case(['inplane', 'clamped'], fetch=True)
    def test_declared_priority(self): self.run_case(['inplane', 'clamped'], declared=['explicit.dat'])

if __name__ == '__main__': unittest.main(verbosity=2)
