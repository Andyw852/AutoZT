# -*- coding: utf-8 -*-
"""scan_project_configs：步骤目录剪枝 + 进程内缓存。

-p 命令在 9p 上慢，主因是扫描 project_roots 时把材料/技能目录下的步骤目录
（step4_disp/disp-*、step5_label/cfg-* …）整棵列一遍。剪枝只在父目录确认是
材料（有 POSCAR）或技能目录时发生，名字以 step 开头的项目目录照常发现。
"""
import os

from autozt import bootstrap as b


def _mk(root, rel, fname=None):
    d = os.path.join(root, rel)
    os.makedirs(d, exist_ok=True)
    if fname:
        open(os.path.join(d, fname), "w").close()


def _tree(root):
    _mk(root, "ProjA/project_setting", "tf_A.yaml")
    _mk(root, "Si", "POSCAR")
    _mk(root, "Si/project_setting", "tf_Si.yaml")
    _mk(root, "Si/kl-dft-cpu/project_setting", "tf_Si_kl.yaml")
    for i in range(50):
        _mk(root, "Si/kl-dft-cpu/step4_disp/disp-%03d" % i)
    _mk(root, "Si/phonon-dft-cpu/phonon_band_plot/x")
    _mk(root, "Legacy", "POSCAR")
    for i in range(50):
        _mk(root, "Legacy/step4_disp/disp-%03d" % i)
    _mk(root, "step2_proj/project_setting", "tf_step2.yaml")


def test_prune_keeps_results_and_skips_step_dirs(tmp_path, monkeypatch):
    root = str(tmp_path)
    _tree(root)
    b._note_skill_dirs({"kl-dft-cpu": {"steps": [{"name": "step4_disp"}]},
                        "phonon-dft-cpu": {"steps": [{"name": "phonon_band_plot"}]}})
    b.invalidate_project_scan()
    seen = []
    real = os.scandir

    def counting(p):
        seen.append(p)
        return real(p)

    monkeypatch.setattr(os, "scandir", counting)
    got = b.scan_project_configs([root])
    assert [n for n, _, _ in got] == ["A", "Si", "Si_kl", "step2"]
    assert not any("disp-" in p or "phonon_band_plot" in p for p in seen)
    assert any(p.endswith("step2_proj") for p in seen)      # step* 项目目录照常下探


def test_cache_and_invalidate(tmp_path, monkeypatch):
    root = str(tmp_path)
    _tree(root)
    b.invalidate_project_scan()
    first = b.scan_project_configs([root])
    calls = []
    real = os.scandir
    monkeypatch.setattr(os, "scandir", lambda p: calls.append(p) or real(p))
    assert b.scan_project_configs([root]) == first and not calls   # 同进程命中缓存
    _mk(root, "NewProj/project_setting", "tf_New.yaml")
    b.invalidate_project_scan()
    assert "New" in [n for n, _, _ in b.scan_project_configs([root])]
