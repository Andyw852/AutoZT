#!/usr/bin/env python3
"""Smoke tests for the thermal-sampling → DFT-label → NEP89 fine-tune route."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TMP_ROOT = ROOT / "tmp"
GPUMD = ROOT / "skill" / "_common" / "gpumd"
KL_TEMPLATES = ROOT / "skill" / "kl-gpumd-3090" / "templates"
STEPCONF = ROOT / "skill" / "_common" / "opt" / "stepconf.py"
sys.path.insert(0, str(GPUMD))
import thermal_md_sampler


def fake_poscar():
    lines = ["synthetic Mn2In2Se5", "1.0", "10 0 0", "0 10 0", "0 0 20",
             "Mn In Se", "18 18 45", "Direct"]
    for i in range(81):
        lines.append("%.8f %.8f %.8f" % ((i % 9) / 10, (i % 7) / 9, (i % 5) / 7))
    return "\n".join(lines) + "\n"


def sample_manifest(root):
    frames = []
    split = thermal_md_sampler.stratified_validation_indices([13, 12, 12, 13], 10)
    (root / "frames").mkdir(parents=True)
    serial = 0
    for ti, (temp, count) in enumerate(zip((300, 500, 650, 800), (13, 12, 12, 13))):
        for i in range(count):
            serial += 1
            name = "POSCAR-%03d" % serial
            (root / "frames" / name).write_text(fake_poscar())
            frames.append({"id": name, "file": "frames/" + name,
                           "temperature_K": temp,
                           "split": "validation" if i in split[ti] else "train",
                           "source_frame_index": i})
    (root / "thermal_manifest.json").write_text(json.dumps({
        "total": 50, "train": 40, "validation": 10, "natoms": 81,
        "frames": frames,
    }))


class ThermalNEPPipelineTests(unittest.TestCase):
    def setUp(self):
        TMP_ROOT.mkdir(exist_ok=True)

    def run_script(self, script, cwd):
        env = dict(os.environ)
        env["PYTHONPATH"] = str(cwd)
        return subprocess.run([sys.executable, str(script)], cwd=str(cwd), env=env,
                              capture_output=True, text=True)

    def test_exact_temperature_stratification(self):
        indices = thermal_md_sampler.stratified_validation_indices([13, 12, 12, 13], 10)
        self.assertEqual([len(indices[i]) for i in range(4)], [3, 2, 2, 3])
        self.assertEqual(sum(len(v) for v in indices.values()), 10)
        self.assertEqual(indices[0], {3, 6, 9})
        self.assertEqual(indices[3], {3, 6, 9})

    def test_periodic_dedup_ignores_translation_and_rejects_near_duplicates(self):
        import numpy as np
        from ase import Atoms
        base = Atoms("Mn2", positions=[[0.2, 0.3, 0.4], [1.1, 1.2, 1.3]],
                     cell=[10, 10, 10], pbc=True)
        shifted = base.copy()
        shifted.positions += [9.7, 0.0, 0.0]
        near = base.copy()
        near.positions[1, 0] += 0.01
        novel1 = base.copy()
        novel1.positions[1, 0] += 0.4
        novel2 = base.copy()
        novel2.positions[1, 1] += 0.5
        self.assertLess(thermal_md_sampler.periodic_rms_displacement(base, shifted), 1e-8)
        self.assertLess(thermal_md_sampler.periodic_rms_displacement(base, near), 0.05)
        self.assertGreater(thermal_md_sampler.periodic_rms_displacement(base, novel1), 0.05)
        chosen, stats = thermal_md_sampler.select_novel_frames(
            [base, near, novel1, novel2], [base], 2, 0.05)
        self.assertEqual([i for i, _ in chosen], [3, 2])
        self.assertEqual(stats["candidate_count"], 4)
        self.assertGreaterEqual(stats["rejected_near_duplicate"], 1)

    def test_sample_generator_prepares_job_without_running_md(self):
        with tempfile.TemporaryDirectory(dir=TMP_ROOT) as temp:
            root = Path(temp)
            for name in ("gen_thermal_sample.py", "thermal_md_sampler.py",
                         "gpumd_common.py", "stepconf.py", "submit_gpumd.tpl"):
                src = (STEPCONF if name == "stepconf.py" else
                       (KL_TEMPLATES / name if name == "submit_gpumd.tpl" else GPUMD / name))
                shutil.copy2(src, root / name)
            (root / "POSCAR").write_text(fake_poscar())
            (root / "seed.txt").write_text("nep4 3 Mn In Se\n")
            (root / "step.conf").write_text(
                "[params]\nGPUMD_BIN=/fake/gpumd\nSEED_NEP=seed.txt\n"
                "TEMPERATURES=300,500,650,800\nCOUNTS=13,12,12,13\n")
            run = self.run_script(root / "gen_thermal_sample.py", root)
            self.assertEqual(run.returncode, 0, run.stderr + run.stdout)
            submit = (root / "step1_thermal_sample" / "submit.sh").read_text()
            self.assertIn("thermal_md_sampler.py", submit)
            self.assertIn("--production-steps 100000", submit)
            self.assertTrue((root / "step1_thermal_sample" / "seed.txt").is_file())

    def test_dft_label_generator_writes_50_inputs_and_refuses_stale_force_labels(self):
        with tempfile.TemporaryDirectory(dir=TMP_ROOT) as temp:
            root = Path(temp)
            (root / "step1_thermal_sample" / "sampled").mkdir(parents=True)
            sample_manifest(root / "step1_thermal_sample" / "sampled")
            (root / "INCAR_thermal_source").write_text(
                "ENCUT=410\nISPIN=2\nMAGMOM=18*5.0 18*0.6 45*0.6\n"
                "LDAU=.TRUE.\nLDAUTYPE=2\nLDAUL=2 -1 -1\n"
                "LDAUU=4 0 0\nLDAUJ=0 0 0\nNCORE=12\nKPAR=4\n")
            (root / "KPOINTS_thermal_source").write_text("mesh\n0\nGamma\n2 2 6\n")
            (root / "POTCAR_thermal_source").write_text(
                "TITEL = PAW_PBE Mn test\nTITEL = PAW_PBE In test\n"
                "TITEL = PAW_PBE Se test\n")
            (root / "submit_std_3d.tpl").write_text(
                "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=24\n"
                "echo {{JOBNAME}}\n")
            (root / "step.conf").write_text(
                "[params]\nINCAR_SOURCE=INCAR_thermal_source\n"
                "KPOINTS_SOURCE=KPOINTS_thermal_source\n"
                "POTCAR_SOURCE=POTCAR_thermal_source\nEXPECTED_FRAMES=50\n"
                "[submit]\nntasks_per_node=48\n")
            shutil.copy2(STEPCONF, root / "stepconf.py")
            script = root / "gen_thermal_label.py"
            shutil.copy2(GPUMD / script.name, script)
            run = self.run_script(script, root)
            self.assertEqual(run.returncode, 0, run.stderr + run.stdout)
            self.assertEqual(len(list((root / "step1_thermal_label").glob("cfg-*"))), 50)
            first = root / "step1_thermal_label" / "cfg-001"
            self.assertIn("2 2 6", (first / "KPOINTS").read_text())
            self.assertIn("ntasks-per-node=48", (first / "submit.sh").read_text())
            self.assertIn("LWAVE      = .FALSE.", (first / "INCAR").read_text())
            self.assertIn("ISPIN=2", (first / "INCAR").read_text())

            (first / "OUTCAR").write_text(
                "General timing and accounting informations\n")
            source = root / "step1_thermal_sample" / "sampled" / "frames" / "POSCAR-001"
            source.write_text(source.read_text().replace("0.00000000", "0.01000000", 1))
            stale = self.run_script(script, root)
            self.assertNotEqual(stale.returncode, 0)
            self.assertIn("旧力标签", stale.stderr + stale.stdout)

    def test_nep_generator_uses_relayed_dataset_and_matching_base_pair(self):
        with tempfile.TemporaryDirectory(dir=TMP_ROOT) as temp:
            root = Path(temp)
            for name in ("gen_step2_nep.py", "gpumd_common.py", "nep_train.py",
                         "vasp_to_xyz.py", "submit_nep.tpl", "stepconf.py"):
                src = (STEPCONF if name == "stepconf.py" else
                       (KL_TEMPLATES / name if name == "submit_nep.tpl" else GPUMD / name))
                shutil.copy2(src, root / name)
            (root / "step1_struct").mkdir()
            (root / "step1_struct" / "gpumd_params.json").write_text(json.dumps({
                "DATA_DIR": "/nonexistent", "DATA_GLOB": "POSCAR-*/OUTCAR",
                "TRAIN_FRAC": 0.9, "SPLIT_SEED": 1,
            }))
            dataset = root / "step1_thermal_dataset"
            dataset.mkdir()
            (dataset / "dataset_summary.json").write_text(
                json.dumps({"DATASET_DONE": True, "train_frames": 1, "test_frames": 1}))
            frame = "1\nLattice=\"10 0 0 0 10 0 0 0 10\" Properties=species:S:1:pos:R:3:force:R:3 energy=0\nMn 0 0 0 0 0 0\n"
            (dataset / "train.xyz").write_text(frame)
            (dataset / "test.xyz").write_text(frame)
            (dataset / "nep89_20250409.txt").write_text("nep4 1 Mn\n")
            (dataset / "nep89_20250409.restart").write_text("matched restart\n")
            (root / "step.conf").write_text(
                "[params]\nNEP_BIN=/fake/nep\nPRETRAINED_NEP=nep89_20250409.txt\n"
                "PRETRAINED_RESTART=nep89_20250409.restart\nELEMENTS=Mn In Se\n"
                "GENERATION=1000\nBATCH=1\n")
            run = self.run_script(root / "gen_step2_nep.py", root)
            self.assertEqual(run.returncode, 0, run.stderr + run.stdout)
            out = root / "step2_nep"
            self.assertEqual((out / "train.xyz").read_text(), frame)
            self.assertTrue((out / "nep89_20250409.txt").is_file())
            self.assertTrue((out / "nep89_20250409.restart").is_file())
            self.assertIn("python nep_train.py --step-dir .", (out / "submit.sh").read_text())
            self.assertNotIn("vasp_to_xyz.py", (out / "submit.sh").read_text())


if __name__ == "__main__":
    unittest.main()
