"""Offline alias/safety contracts; fixtures never leave this worktree's tmp/."""
import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

import autozt
from autozt import bootstrap as b
from autozt import data as data_api


ALIASES = {
    "opt-mace-cpu": "opt-mlff-cpu",
    "opt-mace-gpu": "opt-mlff-gpu",
    "kl-mace-cpu": "kl-mlff-cpu",
    "kl-mace-gpu": "kl-mlff-gpu",
    "phonon-mace-cpu": "phonon-mlff-cpu",
    "phonon-mace-gpu": "phonon-mlff-gpu",
    "mlff-mace": "mlff",
}
KNOWN = "band-dft-cpu"
UNKNOWN = "typo-no-such-skill"
BLOCKED = r"已屏蔽|未知|不支持|迁移|unknown|blocked|unsupported"


class SkillAliases(unittest.TestCase):
    def setUp(self):
        b.reset_config_conflicts()
        scratch = Path(__file__).resolve().parents[1] / "tmp"
        scratch.mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(prefix="skill-alias-tests-", dir=scratch)
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(b.reset_config_conflicts)
        self.root = Path(self.tmp.name).resolve()
        # Deliberately independent of HEAD skill directories. These manifests have
        # no renamed steps: test name aliases separately from migration safety.
        self.registry = {
            key: {"skill_subdir": True,
                  "steps": [{"name": "step1_relax", "label": "S1_relax", "seq": 1}]}
            for key in [*ALIASES.values(), KNOWN]
        }
        for module in (b, autozt):
            p = patch.object(module, "discover_skills",
                             side_effect=lambda *a, **kw: copy.deepcopy(self.registry))
            p.start()
            self.addCleanup(p.stop)
        # Any accidental subprocess/network launch is a test failure, not an
        # opportunity to touch a user's remote project or daemon.
        p = patch("subprocess.Popen", side_effect=AssertionError("offline test launched process"))
        p.start()
        self.addCleanup(p.stop)

    def material(self, name):
        d = self.root / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "POSCAR").write_text("offline fixture\n", encoding="utf-8")
        return d

    def project(self, owner, name, keys, extra=""):
        d = self.root / owner / "project_setting"
        d.mkdir(parents=True, exist_ok=True)
        p = d / ("tf_" + name + ".yaml")
        p.write_text("task_types:\n" + "".join(
            "  " + key + ":\n    local_root: '..'\n" + extra for key in keys),
            encoding="utf-8")
        return p

    def cfg(self):
        return {"project_roots": [str(self.root)], "_config_dir": str(self.root),
                "task_types": {KNOWN: {"local_root": str(self.root)}}}

    def assemble(self, cfg=None):
        return b.merge_project_configs(b.apply_skills(cfg or self.cfg()))

    def snapshot(self):
        return {str(p.relative_to(self.root)): p.read_bytes()
                for p in self.root.rglob("*") if p.is_file()}

    def test_actual_new_manifests_enforce_layout_boundary(self):
        skill_root = Path(__file__).resolve().parents[1] / "skill"
        for old, new in ALIASES.items():
            manifest = b._load_yaml_file(str(skill_root / new / "skill.yaml"))
            self.registry[new]["steps"] = manifest["steps"]
            self.material(new)
            self.project(new, new, [old])
        with contextlib.redirect_stderr(io.StringIO()):
            cfg = self.assemble()
        for old, new in ALIASES.items():
            if old == "mlff-mace":
                b.reject_config_conflict_targets(new, new)
                types = b.get_types(cfg, tt=new)
                self.assertEqual(types[0]["dir_name"], old)
                self.assertEqual(autozt.ops._physical_skill_key(types, new, str(self.root / new)), old)
            else:
                with self.assertRaisesRegex(SystemExit, "迁移"):
                    b.reject_config_conflict_targets(new, new)

    def test_canonical_key_with_explicit_old_directory_is_blocked(self):
        self.registry["opt-mlff-gpu"]["steps"][0]["name"] = "step1_mlff_relax"
        self.material("Legacy")
        self.project("Legacy", "legacy", ["opt-mlff-gpu"], "    dir_name: opt-mace-gpu\n")
        with contextlib.redirect_stderr(io.StringIO()):
            self.assemble()
        with self.assertRaisesRegex(SystemExit, "迁移"):
            b.reject_config_conflict_targets("Legacy", "opt-mlff-gpu")

    def test_global_alias_does_not_abort_unrelated_type_selector(self):
        self.registry["opt-mlff-gpu"]["steps"][0]["name"] = "step1_mlff_relax"
        cfg = self.cfg()
        cfg["task_types"]["opt-mace-gpu"] = {"local_root": str(self.root)}
        with contextlib.redirect_stderr(io.StringIO()):
            cfg = self.assemble(cfg)
        self.assertTrue(b.get_types(cfg, tt=KNOWN))
        with self.assertRaisesRegex(SystemExit, "迁移"):
            b.get_types(cfg, tt="opt-mlff-gpu")

    def test_global_legacy_inheritance_blocks_canonical_project(self):
        self.registry["opt-mlff-gpu"]["steps"][0]["name"] = "step1_mlff_relax"
        cfg = self.cfg()
        cfg["task_types"]["opt-mace-gpu"] = {}
        self.material("Legacy")
        self.project("Legacy", "legacy", ["opt-mlff-gpu"])
        with contextlib.redirect_stderr(io.StringIO()):
            self.assemble(cfg)
        with self.assertRaisesRegex(SystemExit, "迁移"):
            b.reject_config_conflict_targets("Legacy", "opt-mlff-gpu")

    def test_all_seven_aliases_map_in_memory_and_warn_once_without_writes(self):
        paths = []
        for i, (old, new) in enumerate(ALIASES.items()):
            owner = "M" + str(i)
            self.material(owner)
            paths.append(self.project(owner, owner, [old]))
        before = self.snapshot()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            cfg = self.assemble()
            types = b.get_types(cfg)
            # Repeated scan/assembly in one invocation must not spam warnings.
            self.assemble()
        for old, new in ALIASES.items():
            self.assertNotIn(old, cfg["task_types"])
            matching = [t for t in types if t["key"] == new and t.get("_from")]
            self.assertEqual(len(matching), 1, (old, matching))
            self.assertEqual(matching[0]["steps"][0]["name"], "step1_relax")
            self.assertIn(old, err.getvalue())
            self.assertIn(new, err.getvalue())
        for p in paths:
            self.assertIn(str(p), err.getvalue())
        notices = [line for line in err.getvalue().splitlines()
                   if line.startswith(("警告：", "提示："))]
        self.assertEqual(len(notices), 1, err.getvalue())
        self.assertEqual(before, self.snapshot())

    def test_old_type_selector_matches_canonical_type(self):
        self.material("Training")
        self.project("Training", "training", ["mlff-mace"])
        cfg = self.assemble()
        old_types = b.get_types(cfg, tt="mlff-mace")
        new_types = b.get_types(cfg, tt="mlff")
        self.assertEqual(old_types, new_types)
        self.assertTrue(old_types)
        self.assertEqual({t["key"] for t in old_types}, {"mlff"})

    def test_multiskill_project_keeps_each_canonical_segment(self):
        self.material("Multi")
        p = self.project("Multi", "multi", ["opt-mace-cpu", "mlff-mace", KNOWN])
        cfg = self.assemble()
        types = [t for t in b.get_types(cfg) if t.get("_from") == str(p)]
        self.assertEqual({t["key"] for t in types}, {"opt-mlff-cpu", "mlff", KNOWN})
        self.assertTrue(all(t["steps"] for t in types))

    def unknown_fixture(self, mixed=False):
        self.material("Bad")
        self.material("Good")
        p = self.project("Bad", "bad", [UNKNOWN, KNOWN] if mixed else [UNKNOWN])
        self.project("Good", "good", [KNOWN])
        return p

    def test_unknown_single_project_blocks_only_owner_not_shared_root(self):
        p = self.unknown_fixture()
        before = self.snapshot()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            cfg = self.assemble()
            types = b.get_types(cfg)
            self.assemble()
        self.assertNotIn(UNKNOWN, {t["key"] for t in types})
        self.assertEqual([m["name"] for m in b.discover_local(str(self.root), tt=KNOWN)[1]], ["Good"])
        self.assertEqual([m["name"] for m in b.discover_local(str(self.root), tt="mlff")[1]], ["Good"])
        self.assertIn(str(p), err.getvalue())
        self.assertIn(UNKNOWN, err.getvalue())
        self.assertEqual(err.getvalue().count("警告："), 1)
        self.assertEqual(before, self.snapshot())

    def test_unknown_in_multiskill_project_fails_closed_for_whole_material(self):
        self.unknown_fixture(mixed=True)
        cfg = self.assemble()
        self.assertTrue(b.get_types(cfg, tt=KNOWN))
        self.assertEqual([m["name"] for m in b.discover_local(str(self.root), tt=KNOWN)[1]], ["Good"])
        with self.assertRaisesRegex(SystemExit, BLOCKED):
            autozt.resolve_mat_dir(cfg, b.get_types(cfg), KNOWN, "Bad")

    def test_unknown_cached_materials_and_ops_fallback_cannot_reappear(self):
        self.unknown_fixture()
        cfg = self.assemble()
        types = b.get_types(cfg)
        data = {"types": [{"key": KNOWN, "materials": [
            {"name": "Bad", "lpath": str(self.root / "Bad")},
            {"name": "Bad"}, {"name": "Good", "lpath": str(self.root / "Good")},
        ]}]}
        filtered = b.filter_config_conflicts(copy.deepcopy(data))
        self.assertEqual([m["name"] for m in filtered["types"][0]["materials"]], ["Good"])
        # Patch only the unrelated signature encoding (JSON list vs tuple);
        # run the real cache loader and its safety filter against stale data.
        sig = ["offline-fixture"]
        cache = self.root / ".tf_state_cache.json"
        cache.write_text(json.dumps({"ts": time.time(), "sig": sig, "data": data}))
        with patch.object(data_api, "_state_cache_sig", return_value=sig):
            cached = data_api._state_cache_load(cfg, types, None, None, 60)
        self.assertIsNotNone(cached)
        self.assertEqual([m["name"] for m in cached["types"][0]["materials"]], ["Good"])
        for target in ("Bad", str(self.root / "Bad")):
            with self.subTest(target=target), self.assertRaisesRegex(SystemExit, BLOCKED):
                autozt.resolve_mat_dir(cfg, [], KNOWN, target)
        self.assertEqual(autozt.resolve_mat_dir(cfg, types, KNOWN, "Good"), str(self.root / "Good"))
        self.assertNotIn("Bad", autozt._skill_local_mats(cfg, types, KNOWN))

    def test_cli_unknown_explicit_target_rejected_before_side_effects_and_early_returns(self):
        self.unknown_fixture()
        commands = ["list", "start", "history", "init", "auto", "monitor", "monitor -d",
                    "watch", "watch -d", "probe",
                    "level", "hpc", "conf", "skills", "schema", "skill"]
        traps = ["_watch_ensure", "_watch_daemon", "_watch_cron", "_watch_stop",
                 "collect_data", "cmd_start", "cmd_history", "cmd_init", "cmd_auto",
                 "cmd_auto_project", "cmd_auto_skill", "cmd_watch", "cmd_level",
                 "cmd_hpc", "cmd_conf", "cmd_skills", "cmd_schema", "cmd_skill_show"]
        before = self.snapshot()
        for command in commands:
            for run in range(2):
                with self.subTest(command=command, run=run), contextlib.ExitStack() as stack:
                    err = stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                    stack.enter_context(patch.object(autozt, "load_config", return_value=(self.cfg(), None)))
                    stack.enter_context(patch.object(autozt, "agent_direct_gate", return_value=None))
                    stack.enter_context(patch.object(autozt, "set_active_cfg"))
                    for name in traps:
                        stack.enter_context(patch.object(autozt, name, side_effect=AssertionError(name + " reached")))
                    stack.enter_context(patch("sys.argv", ["autozt", "-p", "Good,Bad", "-tt", KNOWN,
                                                          "-x", "Bad", "-status", "done", *command.split()]))
                    with self.assertRaisesRegex(SystemExit, BLOCKED):
                        autozt.cli.main()
                    self.assertIn(UNKNOWN, err.getvalue())
        self.assertEqual(before, self.snapshot())

    def test_renamed_steps_fail_closed_even_without_local_results(self):
        for old, new in ALIASES.items():
            if old == "mlff-mace":
                continue
            self.registry[new]["steps"][0].update(name="step1_mlff_relax", gen="gen_step1_mlff_relax.py")
            self.material(old)
            self.project(old, old, [old])
        self.material("Good")
        self.project("Good", "good", [KNOWN])
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            cfg = self.assemble()
            types = b.get_types(cfg)
        self.assertEqual([m["name"] for m in b.discover_local(str(self.root))[1]], ["Good"])
        self.assertRegex(err.getvalue(), r"迁移|migration")
        for old, new in ALIASES.items():
            if old != "mlff-mace":
                with self.subTest(old=old), self.assertRaisesRegex(SystemExit, BLOCKED):
                    autozt.resolve_mat_dir(cfg, types, new, old)

    def test_existing_mace_step_is_not_silently_repointed_to_mlff(self):
        old, new = "opt-mace-cpu", "opt-mlff-cpu"
        self.registry[new]["steps"][0].update(name="step1_mlff_relax", gen="gen_step1_mlff_relax.py")
        mat = self.material("Legacy")
        self.project("Legacy/" + old, "legacy", [old])
        result = mat / old / "result" / "step1_mace_relax"
        result.mkdir(parents=True)
        (result / "CONTCAR").write_text("valuable existing result\n")
        before = self.snapshot()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            cfg = self.assemble()
            types = b.get_types(cfg)
        with self.assertRaisesRegex(SystemExit, BLOCKED):
            autozt.resolve_mat_dir(cfg, types, new, "Legacy")
        self.assertRegex(err.getvalue(), r"迁移|migration")
        self.assertEqual(before, self.snapshot())
        self.assertFalse((mat / new).exists())

    def test_mlff_unchanged_steps_preserve_old_physical_skill_directory(self):
        mat = self.material("Training")
        p = self.project("Training/mlff-mace", "training", ["mlff-mace"])
        result = mat / "mlff-mace" / "result" / "step1_relax"
        result.mkdir(parents=True)
        (result / "CONTCAR").write_text("existing training result\n")
        before = self.snapshot()
        cfg = self.assemble()
        types = [t for t in b.get_types(cfg, tt="mlff") if t.get("_from") == str(p)]
        self.assertEqual(len(types), 1)
        t = types[0]
        self.assertEqual(t.get("dir_name"), "mlff-mace")
        m = {"name": "Training", "lpath": str(mat)}
        with patch.object(b, "pkg_setting_path", return_value=None):
            b.resolve_material_local(t, str(self.root), m)
        self.assertEqual(m["_skill_dir_local"], str(mat / "mlff-mace"))
        self.assertEqual(t["steps"][0]["name"], "step1_relax")
        self.assertEqual(before, self.snapshot())
        self.assertFalse((mat / "mlff").exists())


if __name__ == "__main__":
    unittest.main()
