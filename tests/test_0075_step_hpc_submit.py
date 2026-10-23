# -*- coding: utf-8 -*-
"""[0075] 步骤级换集群（项目 step_hpc）的提交参数（纯本地、不需要 VASP / 集群）。

  · S2.15_discr 的 submit.sh 按本步集群推过来的 submit_std_*.tpl 渲染、套本步 step.conf 的 [submit]；
    以前复制 S2.1 的 submit.sh，S2.1 在 hfeshell、本步在 jzzn 时带着 hfacexclu08 -> sbatch: invalid partition。
    没有模板（旧 gen_need）才退回复制 S2.1 的。
  · resolve_cores：步骤换了集群时，hpc.yaml（材料默认集群）那一档的 cores 换成 setting/<步骤集群>.yaml 的。
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KE = REPO / "skill" / "ke-dft-cpu"
GEN = KE / "step2_bandgap" / "step2.15_discriminant" / "gen_step2.15_discriminant.py"

POSCAR = "MgC\n1.0\n 5.0 0 0\n 0 5.0 0\n 0 0 6.0\nMg C\n1 2\nDirect\n 0 0 0\n 0.5 0.5 0.3\n 0.5 0.5 0.7\n"
S21_SUBMIT = ("#!/bin/bash\n#SBATCH --job-name=old\n#SBATCH --partition=hfacexclu08\n#SBATCH --nodes=1\n"
              "#SBATCH --ntasks-per-node=48\nmpirun vasp_std\n")
JZZN_TPL = ("#!/bin/bash\n#SBATCH --job-name={{JOBNAME}}\n#SBATCH --partition=cpu192\n#SBATCH --nodes=1\n"
            "#SBATCH --ntasks-per-node=192\nmpirun vasp_std\n")


def _material(td, tpl=True, conf=""):
    m = Path(td) / "Mg2C60"
    s21 = m / "step2_bandgap" / "step2.1_static"
    s21.mkdir(parents=True)
    (s21 / "INCAR").write_text("ISTART = 0\nNCORE = 4\n")
    (s21 / "POSCAR").write_text(POSCAR)
    (s21 / "POTCAR").write_text("fake\n")
    (s21 / "KPOINTS").write_text("auto\n0\nGamma\n 4 4 3\n 0 0 0\n")
    (s21 / "submit.sh").write_text(S21_SUBMIT)
    if tpl:                                   # autozt 按本步集群（jzzn）推过来的模板
        (m / "submit_std_3d.tpl").write_text(JZZN_TPL)
    shutil.copy(REPO / "skill" / "_common" / "opt" / "stepconf.py", m)
    if conf:
        (m / "step.conf").write_text(conf)
    return m


def _gen(m):
    p = subprocess.run([sys.executable, str(GEN)], cwd=str(m), capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return (m / "step2_bandgap" / "step2.15_discriminant" / "submit.sh").read_text(), p.stderr


def test_s215_renders_cluster_template(tmp_path):
    m = _material(tmp_path, conf="[submit]\nqos = regular\ntime = 12:00:00\n")
    sub, log = _gen(m)
    assert "--partition=cpu192" in sub and "hfacexclu08" not in sub
    assert "#SBATCH --job-name=Mg2C60-ke-dft-cpu-S2.15_discr" in sub and "{{" not in sub
    assert "#SBATCH --qos=regular" in sub and "#SBATCH --time=12:00:00" in sub   # 本步 step.conf 的 [submit] 生效
    assert "本集群模板 submit_std_3d.tpl" in log


def test_s215_falls_back_to_s21_copy(tmp_path):
    m = _material(tmp_path, tpl=False)
    sub, log = _gen(m)
    assert "--partition=hfacexclu08" in sub and "--job-name=Mg2C60-ke-dft-cpu-S2.15_discr" in sub
    assert "由 step2.1_static 复制" in log


def test_skill_yaml_pushes_templates():
    txt = (KE / "skill.yaml").read_text(encoding="utf-8")
    seg = txt[txt.index("label: S2.15_discr"):txt.index("label: S2.155_discr")]
    assert "submit_std_2d.tpl" in seg and "submit_std_3d.tpl" in seg


def test_resolve_cores_uses_step_cluster(tmp_path, monkeypatch):
    import autozt
    jz = tmp_path / "jzzn.yaml"
    jz.write_text("name: jzzn\ncores: 192\n")
    plain = tmp_path / "plain.yaml"
    plain.write_text("name: plain\n")
    monkeypatch.setattr(autozt, "pkg_setting_path",
                        lambda n: {"jzzn.yaml": str(jz), "plain.yaml": str(plain)}.get(n))
    t = {"key": "ke-dft-cpu", "steps": [{"name": "step2_bandgap/step2.15_discriminant", "hpc": "jzzn"},
                                       {"name": "step2_bandgap/step2.2_pbe", "hpc": "plain"},
                                       {"name": "step7_deform"}]}
    m = {"ps": {"setting": {}, "hpc": {"cores": 48}}, "tt": "ke-dft-cpu", "hpc_name": "hfeshell"}
    assert autozt.resolve_cores({}, t, m, "step2_bandgap/step2.15_discriminant") == 192
    assert autozt.resolve_cores({}, t, m, "step2_bandgap/step2.2_pbe") is None        # 步骤集群没写 cores：不套 hfeshell 的 48
    assert autozt.resolve_cores({"cores": 64}, t, m, "step2_bandgap/step2.2_pbe") == 64  # 全局 cores 照旧兜底
    assert autozt.resolve_cores({}, t, m, "step7_deform") == 48                         # 没换集群的步骤：照旧
    assert autozt.resolve_cores({}, t, m) == 48
    m2 = dict(m, ps={"setting": {"cores": 16}, "hpc": {"cores": 48}})
    assert autozt.resolve_cores({}, t, m2, "step2_bandgap/step2.15_discriminant") == 16   # 项目 setting.yaml 仍最高


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
