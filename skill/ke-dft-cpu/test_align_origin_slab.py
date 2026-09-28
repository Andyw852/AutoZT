"""align_origin 的 2D slab 回归测试 + S8/S8.4 层厚读取的跨界测试。

用法（需要 pymatgen + spglib 的环境）：  python skill/ke-dft-cpu/test_align_origin_slab.py
背景：520ef11 的 align_origin 对不在 z=0.5 的层（例如 Mo 在 z=0.13）会把 σh 镜面
放到 z=0，层被劈到周期边界两侧，S8/S8.4 的 max(z)-min(z) 层厚读成 ≈c。
"""
import sys, importlib.util, tempfile
from pathlib import Path
import numpy as np
from pymatgen.core import Structure, Lattice
HERE = Path(__file__).resolve().parent                 # skill/ke-dft-cpu
for _p in (HERE, HERE.parent / "_common" / "opt", HERE.parent / "_common"):
    sys.path.insert(0, str(_p))
import ke_common as kc

def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(sp)
    sys.path.insert(0, str(Path(path).parent)); sp.loader.exec_module(m); return m

def zinfo(st):
    M = st.lattice.matrix
    h = abs(np.linalg.det(M)) / np.linalg.norm(np.cross(M[0], M[1]))
    zc, span, slab = kc._slab_center_span(st.frac_coords[:, 2], h)
    fz = np.mod(st.frac_coords[:, 2], 1.0)
    return zc, span * h, (fz.max() - fz.min()) > span + 1e-6

fails = 0
def check(cond, msg):
    global fails
    print(("  PASS " if cond else "  FAIL ") + msg); fails += (not cond)

a, c = 3.19, 20.0
hexL = Lattice.hexagonal(a, c)
dz = 1.56 / c
cases = {
  "MoS2 Mo@(0.1,0.2,0.5)":   ([.1, .2, .5], ["Mo", "S", "S"]),
  "MoS2 Mo@(1/3,2/3,0.3)":   ([1/3, 2/3, .3], ["Mo", "S", "S"]),
  "MoS2 Mo@(2/3,1/3,0.13) [旧版会劈层]": ([2/3, 1/3, .13], ["Mo", "S", "S"]),
  "MoS2 Mo@(0.05,0.4,0.97)": ([.05, .4, .97], ["Mo", "S", "S"]),
  "Janus MoSSe P3m1 @(0.2,0.1,0.12)": ([.2, .1, .12], ["Mo", "S", "Se"]),
}
with tempfile.TemporaryDirectory() as td:
    for name, (mo, sp) in cases.items():
        mo = np.array(mo); s_xy = np.array([mo[0] + 1/3, mo[1] + 2/3]) % 1
        fr = np.array([mo, [*s_xy, mo[2] + dz], [*s_xy, mo[2] - dz]]) % 1.0
        st = Structure(hexL, sp, fr)
        f = Path(td) / "POSCAR"; st.to(filename=str(f), fmt="poscar")
        before = kc._sym_ops_max_tau(st)[0]
        r = kc.align_origin(f)
        st2 = Structure.from_file(str(f))
        after = kc._sym_ops_max_tau(st2)[0]
        zc, span_A, wrapped = zinfo(st2)
        print(f"{name}: tau {before:.4f} -> {after:.4f}  ret={r}  zc={zc:.3f} span={span_A:.3f}Å wrapped={wrapped}")
        check(after < 1e-4, "tau 全为 0")
        check(not wrapped, "层没被周期边界劈开")
        check(abs(zc - 0.5) < 1e-3, "层心在 0.5")
        check(abs(span_A - 3.12) < 1e-3, "层跨度不变")
        check(st2.lattice == st.lattice and sorted(st2.species) == sorted(st.species), "晶格/元素不变")

    # 3D：GaAs（F-43m，原点偏移可消），Si（Fd-3m 非简单，应放弃且不改文件）
    for name, lat, spc, fr in [
        ("GaAs F-43m 原胞 原点偏 (0.1,0.2,0.3)", Lattice([[0,2.825,2.825],[2.825,0,2.825],[2.825,2.825,0]]), ["Ga","As"],
         np.array([[0,0,0],[.25,.25,.25]]) + [.1,.2,.3]),
        ("Si Fd-3m 原胞 (非简单)", Lattice([[0,2.715,2.715],[2.715,0,2.715],[2.715,2.715,0]]), ["Si","Si"],
         np.array([[0,0,0],[.25,.25,.25]]) + [.1,.2,.3]),
    ]:
        st = Structure(lat, spc, fr % 1.0); f = Path(td) / "POSCAR"; st.to(filename=str(f), fmt="poscar")
        txt0 = f.read_text(); r = kc.align_origin(f); st2 = Structure.from_file(str(f))
        after = kc._sym_ops_max_tau(st2)[0]
        print(f"{name}: ret={r} after={after:.4f}")
        if "Si" in name:
            check(r is not None and r[1] is None and f.read_text() == txt0, "非简单群：返回复核失败且文件不动")
        else:
            check(after < 1e-4, "3D 简单群 tau 全为 0")

    # 层厚读取：人为造一个跨界的 POSCAR，两个 gen 的读取函数都应给出 3.12 Å
    fr = np.array([[0, 0, 0.0], [1/3, 2/3, dz], [1/3, 2/3, 1 - dz]])
    st = Structure(hexL, ["Mo", "S", "S"], fr)
    wd = Path(td) / "mat"; (wd / "step3_uniform").mkdir(parents=True)
    st.to(filename=str(wd / "step3_uniform/POSCAR"), fmt="poscar")
    for g in ["step8_amset/gen_step10_amset.py", "step8.4_amset2d/gen_step14_amset2d.py"]:
        m = load(HERE / g, g.split("/")[-1][:-3])
        c_len, span = m.get_slab_geometry(wd)
        t, info = m.vdw_thickness(wd)
        print(f"{g}: c={c_len:.3f} span={span:.3f} t_vdw={t}")
        check(abs(span - 3.12) < 1e-3, f"{g} 跨界层 z 跨度 = 3.12 Å")
        check(abs(t - (3.12 + 2 * 1.80)) < 1e-3, f"{g} 跨界层 vdW 厚度 = 6.72 Å")
print("\nFAILS =", fails); sys.exit(1 if fails else 0)
