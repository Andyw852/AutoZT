import contextlib
import copy
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import autozt
from autozt import bootstrap as b

class Conflicts(unittest.TestCase):
    def setUp(self):
        b.reset_config_conflicts()
        scratch = Path(__file__).resolve().parents[1] / "tmp"
        scratch.mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(dir=scratch)
        self.root = Path(self.tmp.name).resolve()
    def tearDown(self):
        self.tmp.cleanup()
        b.reset_config_conflicts()
    def config(self, group, name='dup', skill='band', text=''):
        d = self.root / group
        (d / 'project_setting').mkdir(parents=True, exist_ok=True)
        p = d / 'project_setting' / ('tf_' + name + '.yaml')
        p.write_text('task_types:' + chr(10) + '  '+skill+':' + chr(10) + '    local_root: ".."' + chr(10) +text)
        return p
    def mat(self, name):
        d = self.root / name
        d.mkdir(parents=True, exist_ok=True)
        (d/'POSCAR').write_text('fixture')
        return d
    def fixture(self):
        ps = [self.config('a/Clash'), self.config('b/Clash', text='    desc: different'+chr(10))]
        for group in ['a/Clash', 'b/Clash', 'a/Good']:
            self.mat(group+'/M' if group.endswith('Clash') else group)
        return ps
    def test_conflict_all_paths_once_and_reverse(self):
        ps = self.fixture()
        self.config('a/Good', 'safe')
        err=io.StringIO()
        with contextlib.redirect_stderr(err):
            first=b.scan_project_configs([self.root/'b', self.root/'a'])
            again=b.scan_project_configs([self.root/'a', self.root/'b'])
        self.assertEqual(first, again)
        self.assertEqual([n for n,p,d in first], ['safe'])
        self.assertEqual(err.getvalue().count('警告：'),1)
        for p in ps: self.assertIn(str(p),err.getvalue())
        self.assertIn('已屏蔽',err.getvalue())
    def test_overlap_and_symlink_silent(self):
        p=self.config('a/M','only')
        alias=self.root/'alias'; alias.symlink_to(self.root/'a',target_is_directory=True)
        err=io.StringIO()
        with contextlib.redirect_stderr(err):
            found=b.scan_project_configs([self.root,self.root/'a',alias])
        self.assertEqual(len(found),1); self.assertEqual(err.getvalue(),'')
    def test_identical_files_are_conflicts(self):
        self.config('a','same');self.config('b','same')
        self.assertEqual(b.scan_project_configs([self.root]),[])
    def test_shared_root_main_and_segment_keep_good(self):
        self.fixture(); b.scan_project_configs([self.root])
        for root in [self.root, self.root/'a', self.root/'b']:
            mats=b.discover_local(str(root),tt='band')[1]
            self.assertTrue(all('Clash' not in m['lpath'] for m in mats))
        self.assertEqual([m['name'] for m in b.discover_local(str(self.root/'a'),tt='band')[1]],['Good'])
        self.assertTrue(b.discover_local(str(self.root/'a'),tt='other')[1])
    def test_targets_and_context(self):
        self.fixture(); b.scan_project_configs([self.root])
        for target in ['M','Clash/M','a/Clash/M',str(self.root/'a/Clash/M'),'Good,M','dup/Clash/M']:
            with self.assertRaisesRegex(SystemExit,'已屏蔽'):
                b.reject_config_conflict_targets(target,'band')
        b.reject_config_conflict_targets('Good','band')
        b.reject_config_conflict_targets('M','other')
    def test_cached_and_legacy_materials_removed(self):
        self.fixture(); b.scan_project_configs([self.root])
        data={'types':[{'key':'band','materials':[{'name':'Clash/M','lpath':str(self.root/'a/Clash/M')},{'name':'Good','lpath':str(self.root/'a/Good')},{'name':'M'}]}]}
        self.assertEqual([m['name'] for m in b.filter_config_conflicts(data)['types'][0]['materials']],['Good'])
    def test_ops_fallback_rejects(self):
        self.fixture(); b.scan_project_configs([self.root])
        with self.assertRaisesRegex(SystemExit,'已屏蔽'):
            autozt.resolve_mat_dir({'project_roots':[str(self.root)]},[], 'band','M')
    def test_cli_rejects_before_watch_collect_filters_each_invocation(self):
        self.fixture()
        cfg={'project_roots':[str(self.root)],'task_types':{'band':{'local_root':str(self.root/'a'),'steps':[]}}}
        for cmd in ['list','start','history','init','auto','monitor']:
            for run in range(2):
                err=io.StringIO()
                with patch.object(autozt,'load_config',return_value=(copy.deepcopy(cfg), None)), patch.object(autozt,'agent_direct_gate',return_value=None), patch.object(autozt,'apply_skills',side_effect=lambda c,**kw:c), patch.object(autozt,'_watch_ensure',side_effect=AssertionError('watch called')), patch.object(autozt,'collect_data',side_effect=AssertionError('collect called')), patch('sys.argv',['autozt','-tt','band','-p','Good,M','-x','M','-status','done',cmd]), contextlib.redirect_stderr(err):
                    with self.assertRaisesRegex(SystemExit,'已屏蔽'):
                        autozt.cli.main()
                self.assertEqual(err.getvalue().count('警告：'),1)

    def test_other_skill_and_ops_cache_context(self):
        self.fixture(); b.scan_project_configs([self.root])
        ts=[{'key':'other','local_root':str(self.root/'a')}]
        self.assertIn('Clash/M',autozt._skill_local_mats({},ts,'other'))
        self.assertEqual(autozt.resolve_mat_dir({},ts,'other','M'),str(self.root/'a/Clash/M'))
        with self.assertRaisesRegex(SystemExit,'已屏蔽'):
            autozt.resolve_mat_dir({},ts,'band','M')
    def test_reverse_scandir(self):
        self.fixture(); self.config('a/Good','safe')
        expected=b.scan_project_configs([self.root])
        real=os.scandir
        class Reversed:
            def __init__(self,p):
                with real(p) as it: self.entries=list(it)[::-1]
            def __iter__(self): return iter(self.entries)
            def __enter__(self): return self
            def __exit__(self,*args): pass
        b.reset_config_conflicts()
        with patch.object(os,'scandir',side_effect=Reversed):
            self.assertEqual(b.scan_project_configs([self.root]),expected)

    def test_init_cannot_reintroduce_conflict(self):
        self.fixture(); b.scan_project_configs([self.root])
        for target,name in [(self.root/'a/Clash/M',None),(self.root/'a/Good','dup')]:
            with self.assertRaisesRegex(SystemExit,'屏蔽'):
                autozt._init_one_skill({'task_types':{'band':{'steps':[]}}},[],str(target),name=name,tt='band',known_names={})
        self.assertFalse((self.root/'a/Good/project_setting').exists())
    def test_monitor_daemon_reject_before_control(self):
        self.fixture()
        cfg={'project_roots':[str(self.root)]}
        for flag in ['--daemon','--restart']:
            with patch.object(autozt,'load_config',return_value=(copy.deepcopy(cfg),None)), patch.object(autozt,'agent_direct_gate',return_value=None), patch.object(autozt,'_watch_daemon',side_effect=AssertionError('daemon')), patch.object(autozt,'_watch_stop',side_effect=AssertionError('stop')), patch('sys.argv',['autozt','-tt','band','-p','M','monitor',flag]):
                with self.assertRaisesRegex(SystemExit,'已屏蔽'): autozt.cli.main()

    def test_qualified_name_with_default_local_root(self):
        for group in ['treeA','treeB']:
            p=self.config(group,'batch')
            p.write_text('task_types:' + chr(10) + '  band: {}' + chr(10))
            self.mat(group+'/Si')
        b.scan_project_configs([self.root])
        with self.assertRaisesRegex(SystemExit,'已屏蔽'):
            b.reject_config_conflict_targets('batch/Si','band')
    def test_invalid_conflict_yaml_is_conservatively_blocked(self):
        self.fixture()
        with patch.object(b,'_load_yaml_file',return_value=['invalid']):
            self.assertEqual(b.scan_project_configs([self.root]),[])
        self.assertTrue(b.config_material_blocked(str(self.root/'a/Clash/M'),'anything'))

if __name__=='__main__': unittest.main(verbosity=2)
