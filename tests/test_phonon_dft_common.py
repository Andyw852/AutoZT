"""回归：phonon-dft-cpu 的 S2/S3 数据桥接（2026-09-27 集成测试发现的两个真 bug）。

  1) phonopy 的 phonopy_disp.yaml 里 displacements[k] 是【单个 dict】
     （不是 phono3py 的 list-of-dict），旧移植代码 rec[0] 直接 KeyError；
  2) yaml 里的 atom 是【1-based】（实测 phonopy 2.47：atom:1 → 0-based 下标 0），
     旧代码当 0-based → 真实数据逐帧错位，S3 抽帧门禁会误报失败。
本测试锁死这两点。不需要 phonopy（只造文本）。
"""
import importlib.util
import os
import sys

import numpy as np
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILL = os.path.join(ROOT, "skill", "phonon-dft-cpu")
COMMON = os.path.join(ROOT, "skill", "_common")
OPT = os.path.join(COMMON, "opt")


def _load_kc():
    for p in (SKILL, OPT, COMMON):
        if p not in sys.path:
            sys.path.insert(0, p)
    spec = importlib.util.spec_from_file_location(
        "kc_phonon_dft", os.path.join(SKILL, "kl_common.py"))
    mod = importlib.util.module_from_spec(spec)
    old = os.getcwd()
    os.chdir(SKILL)
    try:
        spec.loader.exec_module(mod)
    finally:
        os.chdir(old)
    return mod


def _write_case(root, move_atom0=True):
    """造一个 2 原子超胞 + 1 帧位移；move_atom0=True 时真实移动 0-based 原子 0。"""
    lat = np.eye(3) * 3.0
    pts = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    disp0 = np.array([0.01, 0.0, 0.0])          # yaml 说 atom:1（1-based）动了
    dd = {"supercell": {"lattice": lat.tolist(),
                        "points": [{"coordinates": c} for c in pts]},
          "displacements": [{"atom": 1, "displacement": disp0.tolist()}]}
    (root / "phonopy_disp.yaml").write_text(yaml.safe_dump(dd, sort_keys=False),
                                            encoding="utf-8")
    cart = (np.array(pts) @ lat).copy()
    if move_atom0:
        cart[0] += disp0                        # 正确：0-based 0 = 1-based 1
    else:
        cart[1] += disp0                        # 错误：移动 0-based 1
    frac = cart @ np.linalg.inv(lat)
    d = root / "disp-00001"
    d.mkdir(parents=True, exist_ok=True)
    rb = "".join("<v> %.10f %.10f %.10f </v>\n" % tuple(v) for v in lat)
    rp = "".join("<v> %.10f %.10f %.10f </v>\n" % tuple(v) for v in frac)
    (d / "vasprun.xml").write_text(
        '<?xml version="1.0"?>\n<modelization>\n<structure name="initialpos">\n'
        '<crystal><varray name="basis" >\n%s</varray></crystal>\n'
        '<varray name="positions" >\n%s</varray>\n</structure>\n</modelization>\n'
        % (rb, rp), encoding="utf-8")


def test_phonopy_disp_dict_format_and_1based_atom(tmp_path):
    kc = _load_kc()
    root = tmp_path
    _write_case(root, move_atom0=True)
    ok, note = kc.check_frames_match_displacements(str(root), disp_yaml="phonopy_disp.yaml")
    assert ok is True, note                       # dict 格式不再 KeyError，且 1-based 对齐

    # 反例：坐标对不上（移动了错误的原子）必须被抓住
    _write_case(root, move_atom0=False)
    ok2, note2 = kc.check_frames_match_displacements(str(root), disp_yaml="phonopy_disp.yaml")
    assert ok2 is False, note2
    assert "坐标最大偏差" in note2


if __name__ == "__main__":
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as t:
        test_phonopy_disp_dict_format_and_1based_atom(Path(t))
    print("ALL PASS")
