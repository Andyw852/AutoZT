# -*- coding: utf-8 -*-
"""V166：band-dft-cpu 的 HSE gen 写死 HFSCREEN=0.11（bohr⁻¹ 当成 Å⁻¹；标准 HSE06 是 0.2 Å⁻¹）。
Si 的带隙因此是 1.340 eV（按 0.2 是 1.091）。改回 0.2，并在写完 INCAR 后核对杂化参数。离线、只读源码与临时 INCAR。
"""
import importlib.util
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "skill" / "band-dft-cpu" / "gen_step4_HSE.py"


def _load_helper():
    """只取 hybrid_warning（gen 顶层 import 了 dim_common/stepconf，不整体导入）。"""
    src = GEN.read_text(encoding="utf-8")
    i = src.index("def hybrid_warning(")
    j = src.index("\ndef main():", i)
    mod = types.ModuleType("band_hybrid")
    exec(compile("from pathlib import Path\n" + src[i:j], str(GEN), "exec"), mod.__dict__)
    return mod.hybrid_warning


class BandHfscreenTests(unittest.TestCase):
    def test_default_is_hse06(self):
        src = GEN.read_text(encoding="utf-8")
        m = re.search(r'"HFSCREEN":\s*"([^"]+)"', src)
        self.assertEqual(m.group(1), "0.2")

    def test_warning(self):
        hw = _load_helper()

        def inc(t):
            p = Path(tempfile.mkdtemp()) / "INCAR"
            p.write_text(t)
            return p
        self.assertIn("bohr", hw(inc("LHFCALC = .TRUE.\nHFSCREEN = 0.11\nAEXX = 0.25\n")))
        self.assertIsNone(hw(inc("LHFCALC = .TRUE.\nHFSCREEN = 0.2\nAEXX = 0.25\n")))
        self.assertIsNone(hw(inc("LHFCALC = .TRUE.\nHFSCREEN = 0.3\n")))
        self.assertIsNone(hw(inc("ISMEAR = 0\n")))

    def test_wired_after_incar_write(self):
        src = GEN.read_text(encoding="utf-8")
        i = src.index('with open(os.path.join(out_dir, "INCAR"), "w") as f:')
        self.assertLess(i, src.index("hybrid_warning(os.path.join(out_dir", i))


if __name__ == "__main__":
    unittest.main(verbosity=2)
