#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""defect-dft-cpu 自检（通用化正确性回归）。

不依赖外部材料文件：用最小合成结构 + 一段内嵌的 A2B2Te5 9 原子原胞做检查。
优先用 spglib；没有 spglib 时跳过积分测试（打印 SKIP）；生产不允许 legacy 回退。

用法：  python3 selftest.py            # 自动
        python3 selftest.py --strict   # 没有 spglib 直接失败
"""
import sys, os
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import defects_common as D

FAILS = []
PASSES = []


def check(name, got, want):
    ok = got == want
    (PASSES if ok else FAILS).append(name)
    print("%-4s %-42s got=%s want=%s" % ("OK" if ok else "FAIL", name, got, want))


def have_spglib():
    try:
        import spglib  # noqa
        return True
    except Exception:
        return False


# ---- 合成 1：NaCl 型 MgC（Fm-3m），4 Mg + 4 C -------------------------------
def nacl():
    lat = [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]]
    mg = [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]
    c = [[.5, .5, .5], [.5, 0, 0], [0, .5, 0], [0, 0, .5]]
    return {"lat": lat, "atoms": ["Mg"] * 4 + ["C"] * 4, "coords": mg + c}


# ---- 合成 2：低对称（把 NaCl 里 1 个 C 挪 0.05）-----------------------------
def nacl_distorted():
    st = nacl()
    st["coords"] = [c[:] for c in st["coords"]]
    st["coords"][4][0] += 0.05
    return st


# ---- 内嵌 A2B2Te5 9 原子原胞（P-3m1）---------------------------------------
A2B2TE5 = """Pb2Sb2Te5_B test (real primitive, P-3m1)
1.0
4.3949485444699379 0.0 0.0
-2.1974742722349689 3.8061370877757201 0.0
0.0 0.0 17.4367173625845666
Sb Te Pb
2 5 2
Direct
0.0000000000 0.0000000000 0.3341525125
0.0000000000 0.0000000000 0.6658474875
0.6666666870 0.3333333430 0.5722598601
0.3333333130 0.6666666570 0.4277401399
0.6666666870 0.3333333430 0.2248345367
0.3333333130 0.6666666570 0.7751654633
0.0000000000 0.0000000000 0.0
0.3333333430 0.6666666870 0.1128108416
0.6666666570 0.3333333130 0.8871891584
"""


def main():
    strict = "--strict" in sys.argv
    spg = have_spglib()
    print("spglib available:", spg)

    # 1) element_roles 自动判定
    cats, ans = D.element_roles({}, ["Mg", "C"])
    check("roles auto Mg/C -> cations", cats, ["Mg"])
    check("roles auto Mg/C -> anions", ans, ["C"])
    cats, ans = D.element_roles({}, ["Pb", "Sb", "Te"])
    check("roles auto Pb/Sb/Te -> cations", sorted(cats), ["Pb", "Sb"])
    check("roles auto Pb/Sb/Te -> anions", ans, ["Te"])

    # 2) spin_seed_magmom 分量数
    check("magmom ncomp=1 token count", len(D.spin_seed_magmom(3, ncomp=1).split()), 3)
    check("magmom ncomp=3 token count", len(D.spin_seed_magmom(3, ncomp=3).split()), 9)

    # 3) charge_states window / valence
    check("charge window -1..1", D.charge_states("v_Mg", {"QMIN": "-1", "QMAX": "1"}), [-1, 0, 1])
    vconf = {"CHARGE_MODE": "valence", "CATIONS": "Mg", "ANIONS": "C"}
    check("valence v_Mg (cation)", D.charge_states("v_Mg", vconf, ["Mg", "C"]), [0, -1, -2])
    check("valence v_C (anion)", D.charge_states("v_C", vconf, ["Mg", "C"]), [0, 1, 2])
    check("valence Mg_C", D.charge_states("Mg_C", vconf, ["Mg", "C"]), [0, -1, -2])
    check("valence C_Mg", D.charge_states("C_Mg", vconf, ["Mg", "C"]), [0, 1, 2])
    check("valence pair", D.charge_states("A_B__B_A_pair", vconf, ["Mg", "C"]), [0])

    # 4-7) Occupied-site symmetry and deduplicated candidate pools.
    # Explicit dependency failure is preferable to legacy z-coordinate grouping.
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory(prefix="defect-selftest-") as td:
        pp = Path(td) / "POSCAR"
        pp.write_text(A2B2TE5)
        st_a = D.parse_poscar(str(pp))
    if spg:
        st = nacl()
        check("NaCl occupied orbit count", len(D.symmetry_orbits(st, mode="spglib")), 2)
        ds = D.enumerate_defects(st, {"CATIONS":"Mg", "ANIONS":"C", "SYMMODE":"spglib"})
        cs = Counter(m["type"] for *_, m in ds)
        check("NaCl vacancy/substitution counts", (cs["vacancy"],cs["antisite"],cs["pair"]), (2,2,0))
        check("NaCl interstitial budget is upper bound", 0 < cs["interstitial"] <= 4, True)
        seeds = [[.125,.125,.125],[.125,.625,.625]]
        reduced = D.interstitial_orbits(st, seeds, mode="spglib")
        check("translation equivalent voids merged", len(reduced), 1)
        check("merged void tracks raw seeds", reduced[0][1]["seed_count"], 2)
        check("void orbit multiplicity positive", reduced[0][1]["site_count"] > 0, True)
        check("distorted NaCl symmetry runs", len(D.symmetry_orbits(nacl_distorted(), mode="spglib")) >= 2, True)
        check("A2B2Te5 occupied orbit count", len(D.symmetry_orbits(st_a, mode="spglib")), 5)
        big = D.supercell(st_a, (3,3,1))
        ds = D.enumerate_defects(big, {"CATIONS":"Pb,Sb", "ANIONS":"Te", "SYMMODE":"spglib"})
        cs = Counter(m["type"] for *_, m in ds)
        check("A2B2Te5 occupied defect counts", (cs["vacancy"],cs["antisite"],cs["pair"]), (5,10,1))
        check("A2B2Te5 interstitial budget upper bound", 0 < cs["interstitial"] <= 6, True)
        check("A2B2Te5 vacancy multiplicity sums to81", sum(m["site_count"] for *_,m in ds if m["type"]=="vacancy"),81)
    else:
        print("SKIP symmetry integration (spglib unavailable; --strict fails)")
    try:
        D.symmetry_orbits(nacl(), mode="legacy")
        check("unsafe legacy symmetry is rejected", False, True)
    except (SystemExit, ValueError):
        check("unsafe legacy symmetry is rejected", True, True)

    # 8) 凸包通用化：二元 Mg-C（Mg4C60 同型）+ 三元已知顶点
    try:
        import convex_hull as CH
        els = {"Mg": -1.5, "C": -9.0}
        tgt = {"formula": {"Mg": 8, "C": 120}, "E": 8*-1.5 + 120*-9.0 - 2.0}
        phs = [{"name": "Mg2C3", "formula": {"Mg": 2, "C": 3},
                "E": 2*-1.5 + 3*-9.0 - 0.2}]
        verts, dH = CH.convex_hull_window(tgt, els, phs)
        check("hull binary Mg-C vertex count", len(verts), 2)
        check("hull binary all dmu<=0",
              all(x <= 1e-9 for v in verts for x in v["dmu"].values()), True)
        check("hull binary equilibrium",
              all(abs(sum(tgt["formula"][e]*v["dmu"][e] for e in els) - dH) < 1e-4
                  for v in verts), True)
        els3 = {"Sn": -3.0, "Sb": -4.0, "Te": -5.0}
        tgt3 = {"formula": {"Sn": 2, "Sb": 2, "Te": 5},
                "E": 2*-3.0 + 2*-4.0 + 5*-5.0 - 1.10}
        phs3 = [{"name": "SnTe", "formula": {"Sn": 1, "Te": 1}, "E": -3.0-5.0-0.3},
                {"name": "Sb2Te3", "formula": {"Sb": 2, "Te": 3}, "E": 2*-4.0+3*-5.0-0.5},
                {"name": "SnSb", "formula": {"Sn": 1, "Sb": 1}, "E": -3.0-4.0-0.1}]
        v3, _ = CH.convex_hull_window(tgt3, els3, phs3)
        got3 = sorted(tuple(round(v["dmu"][e], 4) for e in ("Sn", "Sb", "Te")) for v in v3)
        want3 = sorted([(-0.1333, 0.0, -0.1667), (-0.3, -0.25, 0.0)])
        check("hull ternary known vertices", got3, want3)
    except ImportError:
        print("SKIP convex_hull (module missing)")

    # 9) 金属形成能：位点多重度优先取 manifest meta
    try:
        import formation_energy as FE
        stb = {"lat": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
               "atoms": ["Mg", "C"], "coords": [[0, 0, 0], [.5, .5, .5]]}
        check("site_mult uses meta.site_count",
              FE.site_multiplicity(stb, "v_Mg", "v_Mg(z=0.0)", {"site_count": 7}), 7)
        check("site_mult fallback without meta",
              FE.site_multiplicity(stb, "v_Mg", "v_Mg(z=0.0)"), 1)
    except ImportError:
        print("SKIP formation_energy (module missing)")

    # 10) submit 模板：build_job 按 SOC 自动选 vasp_std / vasp_ncl
    import tempfile, shutil
    td = tempfile.mkdtemp(prefix="dsub_")
    try:
        for nm, txt in (("submit_std_3d.tpl", "vasp_std {{JOBNAME}}\n"),
                        ("submit_ncl_3d.tpl", "vasp_ncl {{JOBNAME}}\n")):
            with open(os.path.join(td, nm), "w") as f:
                f.write(txt)
        os.makedirs(os.path.join(td, "potpaw_PBE", "Mg"))
        with open(os.path.join(td, "potpaw_PBE", "Mg", "POTCAR"), "w") as f:
            f.write("Mg\n")
        st1 = {"lat": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
               "atoms": ["Mg"], "coords": [[0, 0, 0]]}
        oldcwd = os.getcwd()
        os.chdir(td)
        try:
            conf = {"POTCAR_DIR": td, "SOC": "0", "KMESH": "1 1 1",
                    "NCORE": "1", "KPAR": "1", "ENCUT": "300"}
            D.build_job("j_std", st1, conf, {}, "std")
            conf["SOC"] = "1"
            D.build_job("j_ncl", st1, conf, {}, "ncl")
            std_txt = open("j_std/submit.sh").read()
            ncl_txt = open("j_ncl/submit.sh").read()
            check("build_job SOC=0 -> vasp_std", "vasp_std" in std_txt, True)
            check("build_job SOC=1 -> vasp_ncl", "vasp_ncl" in ncl_txt, True)
            lso0 = open("j_std/INCAR").read()
            check("SOC=0 INCAR has no LSORBIT", "LSORBIT" in lso0, False)
            check("SOC=1 INCAR has LSORBIT", "LSORBIT" in open("j_ncl/INCAR").read(), True)
        finally:
            os.chdir(oldcwd)
    finally:
        shutil.rmtree(td)

    # 11) 凸包：energy_above_hull（稳定~0 / 亚稳>0）与亚稳窗口诊断
    try:
        import convex_hull as CH2
        els_ab = {"A": 0.0, "B": 0.0}
        ph_ab = [{"name": "AB", "formula": {"A": 1, "B": 1}, "E": -1.0}]
        eah_s = CH2.energy_above_hull({"formula": {"A": 1, "B": 1}, "E": -1.0}, els_ab, ph_ab)
        check("energy_above_hull stable ~0", abs(eah_s) < 1e-6, True)
        eah_m = CH2.energy_above_hull({"formula": {"A": 2, "B": 1}, "E": -0.5}, els_ab, ph_ab)
        check("energy_above_hull metastable >0", eah_m > 0.1, True)
        raised = ""
        try:
            CH2.convex_hull_window({"formula": {"A": 2, "B": 1}, "E": -0.5}, els_ab, ph_ab)
        except SystemExit as ex:
            raised = str(ex)
        check("metastable window -> 亚稳 diagnostic", "亚稳" in raised, True)
    except ImportError:
        print("SKIP energy_above_hull (module missing)")

    # 12) 类型汇总 _aggregate_types：按浓度（含位点多重度）汇总，而非只看 E_f
    try:
        import formation_energy as FE2
        conc = [
            {"defect": "v_Mg", "conc_cm3": 7.0e5, "Ef_min": 0.95, "m": 8},
            {"defect": "v_C",  "conc_cm3": 1.5e6, "Ef_min": 1.00, "m": 120},
            {"defect": "v_C",  "conc_cm3": 1.0e5, "Ef_min": 1.01, "m": 4},
            {"defect": "Mg_i", "conc_cm3": 5.4e2, "Ef_min": 1.10, "m": 4},
        ]
        nt = {"v_Mg": "vacancy", "v_C": "vacancy", "Mg_i": "interstitial"}
        bn, bt = FE2._aggregate_types(conc, nt)
        check("agg by_name dominant = v_C", bn[0]["name"], "v_C")
        check("agg v_C m_total", [x for x in bn if x["name"] == "v_C"][0]["m_total"], 124)
        check("agg v_C conc summed", round([x for x in bn if x["name"] == "v_C"][0]["conc_cm3"]), 1600000)
        check("agg by_type dominant = vacancy", bt[0]["type"], "vacancy")
        check("agg vacancy nsites", bt[0]["nsites"], 3)
    except ImportError:
        print("SKIP _aggregate_types (module missing)")

    # Shared charge states include the pristine state once, never overfill a site.
    import formation_energy as FC
    probs = FC.site_probabilities([0.0]*5, [-2,-1,0,1,2], 0.0, 1.0)
    check("five equal charge states occupancy 5/6", abs(sum(probs)-5.0/6.0)<1e-12, True)
    probs = FC.site_probabilities([-10000.0,10000.0], [0,1],0.0,40.0)
    check("extreme formation energies stay finite bounded", all(0 <= p <= 1 for p in probs) and sum(probs)<=1.0+1e-12, True)

    # Regression: a one-point chemical-potential region is not empty.
    import convex_hull as HC
    for label, target, elements in [("unary", {"formula":{"Si":1},"E":-5.0}, {"Si":-5.0}),
                                     ("binary point", {"formula":{"A":1,"B":1},"E":0.0}, {"A":0.0,"B":0.0})]:
        vs, _ = HC.convex_hull_window(target, elements, [])
        check(label + " window valid", len(vs), 1)
    try:
        HC.convex_hull_window({"formula":{"Si":1},"E":-5.0}, {"Si":-5.0},
                             [{"formula":{"Si":1},"E":-6.0}])
        check("unstable unary rejected", False, True)
    except SystemExit:
        check("unstable unary rejected", True, True)

    # Regression: fractional wrapping is not MIC in a skew lattice.
    G = D._metric([[1.,0.,0.],[.9,.2,0.],[0.,0.,2.]])
    mic = D._frac_sep([.49,.49,0.], [0.,0.,0.], G)
    check("skew minimum image", abs(mic - (0.031**2 + 0.102**2)**0.5) < 1e-10, True)

    print()
    print("PASS %d / FAIL %d" % (len(PASSES), len(FAILS)))
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    if strict and not spg:
        print("strict: spglib missing")
        return 2
    return 0


# Persistent workflow regression tests, isolated fixtures retained under repository tmp/.
import ast, json, os, runpy, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
ROOT = Path(__file__).resolve().parents[2]
SKILL = ROOT / 'skill/defect-dft-cpu'
TMP = ROOT / 'tmp'
TMP.mkdir(exist_ok=True)

def functions(path, names):
    tree = ast.parse(path.read_text())
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    scope = {'os': os, 'json': json}
    exec(compile(tree, str(path), 'exec'), scope)
    return scope

class Repair(unittest.TestCase):
    def setUp(self):
        self.old = Path.cwd()
        self.root = Path(tempfile.mkdtemp(prefix='defect-s34-', dir=TMP))
        os.chdir(self.root)
        self.s3 = runpy.run_path(str(SKILL/'gen_step3_charged.py'))
        self.s4 = runpy.run_path(str(SKILL/'gen_step4_analysis.py'))
        Path('step2_defects/def-X').mkdir(parents=True)
        Path('step2_defects/defects_manifest.json').write_text(json.dumps([{'dir':'def-X','name':'V_X','disp':'X','counts':{'X':1}}]))
        Path('step2_defects/def-X/CONTCAR').write_text('X\n1\n1 0 0\n0 1 0\n0 0 1\nX\n1\nDirect\n0 0 0\n')
        Path('step2_defects/def-X/OUTCAR').write_text('fixture neutral reached required accuracy')
        Path('step.conf').write_text('[params]\nCHARGE_MODE=window\nQMIN=0\nQMAX=0\nIS_METAL=  \n')
    def tearDown(self):
        os.chdir(self.old)
    def test_neutral_no_jobs_and_stale(self):
        with patch.object(self.s3['D'], 'build_job', side_effect=AssertionError('VASP generated')):
            self.s3['main']()
        p = Path('step3_charged')
        self.assertEqual(list(p.glob('def-*')), [])
        ck = functions(ROOT/'autozt/_collector_remote.py', ['empty_fanout_ok'])['empty_fanout_ok']
        self.assertTrue(ck(str(p.resolve()), 'charged_manifest.json'))
        Path('step.conf').write_text('[params]\nQMIN=-1\nQMAX=1\n')
        self.assertFalse(ck(str(p.resolve()), 'charged_manifest.json'))
        import hashlib
        data=json.loads((p/'charged_manifest.json').read_text())
        data['input_files']['step.conf']=hashlib.sha256(Path('step.conf').read_bytes()).hexdigest()
        (p/'charged_manifest.json').write_text(json.dumps(data))
        self.assertFalse(ck(str(p.resolve()), 'charged_manifest.json'))
    def test_vertical_and_relax_no_q0(self):
        for mode in ['vertical','relax']:
            with patch.object(self.s3['D'], 'charge_states', return_value=[-1,0,1]), patch.object(self.s3['D'], 'read_nelect_from_potcar', return_value=[4]), patch.object(self.s3['D'], 'load_stepconf', return_value={'CHARGED_GEOMETRY':mode}), patch.object(self.s3['D'], 'build_job', side_effect=lambda directory, *args: Path(directory).mkdir(exist_ok=True)) as build:
                self.s3['main']()
                self.assertEqual(build.call_count, 2)
                self.assertEqual(build.call_args.args[3]['NSW'], '60' if mode=='relax' else '0')
                ck=functions(ROOT/'autozt/_collector_remote.py', ['empty_fanout_ok'])['empty_fanout_ok']
                self.assertTrue(ck(str(Path('step3_charged').resolve()), 'charged_manifest.json', require_empty=False))
    def test_failure_archives_old_marker(self):
        Path('step4_analysis').mkdir()
        marker=Path('step4_analysis/analysis_complete.json')
        marker.write_text('{"status":"done"}')
        with patch.object(self.s4['FE'], 'read_energy', return_value=None):
            with self.assertRaises(SystemExit): self.s4['main']()
        self.assertFalse(marker.exists())
        self.assertEqual(len(list(marker.parent.glob('analysis_complete.json.previous-*'))), 1)
    def test_analysis_failure_and_success(self):
        Path('step0_references').mkdir()
        Path('step0_references/references_energy.json').write_text('{}')
        Path('energies.json').write_text('{}')
        Path('step1_bulk').mkdir()
        Path('step1_bulk/OUTCAR').write_text('fixture bulk reached required accuracy')
        marker=Path('step4_analysis/analysis_complete.json')
        for stage in ['CH','EB','FE','ok']:
            def fail(): raise SystemExit('expected fixture failure')
            with patch.object(self.s4['FE'], 'read_energy', return_value=-1.0), patch.object(self.s4['CH'], 'run_from_references', side_effect=fail if stage=='CH' else None), patch.object(self.s4['EB'], 'main', side_effect=fail if stage=='EB' else None), patch.object(self.s4['FE'], 'analyze', side_effect=fail if stage=='FE' else None, return_value={'status':'done','defects':{'def-X':{}}}):
                if stage=='ok':
                    self.s4['main']()
                    self.assertTrue(marker.exists())
                    self.assertEqual(json.loads(Path('step4_analysis/formation_energy_results.json').read_text())['energy_interpretation'], 'neutral_only')
                    ck=functions(ROOT/'autozt/_collector_remote.py', ['ck_plot'])['ck_plot']
                    sc={'done_marker':'analysis_complete.json','strict_done_marker':True}
                    self.assertTrue(ck(str(marker.parent.resolve()), sc)[0])
                    Path('step.conf').write_text('changed')
                    self.assertFalse(ck(str(marker.parent.resolve()), sc)[0])
                else:
                    with self.assertRaises(SystemExit): self.s4['main']()
                    self.assertFalse(marker.exists())
    def test_unconverged_bulk_rejected(self):
        Path('step1_bulk').mkdir()
        Path('step1_bulk/OUTCAR').write_text('finite energy but no ionic convergence')
        with patch.object(self.s4['FE'], 'read_energy', return_value=-1.0):
            with self.assertRaises(SystemExit): self.s4['main']()
        self.assertFalse(Path('step4_analysis/analysis_complete.json').exists())
    def test_explicit_metal_skips_band(self):
        Path('step1_bulk').mkdir()
        Path('step1_bulk/OUTCAR').write_text('reached required accuracy')
        Path('step0_references').mkdir()
        Path('step0_references/references_energy.json').write_text('{}')
        Path('energies.json').write_text('{}')
        Path('step.conf').write_text(Path('step.conf').read_text().replace('IS_METAL=  ', 'IS_METAL=1'))
        with patch.object(self.s4['FE'], 'read_energy', return_value=-1.0), patch.object(self.s4['CH'], 'run_from_references'), patch.object(self.s4['EB'], 'main', side_effect=AssertionError('metal must not require band')), patch.object(self.s4['FE'], 'analyze', return_value={'status':'done','defects':{'def-X':{}}}):
            self.s4['main']()
        self.assertIs(json.loads(Path('energies.json').read_text())['is_metal'], True)
    def test_guard_neutral_never_sbatch(self):
        import base64, subprocess, sys
        self.s3['main']()
        tree = ast.parse((ROOT/'autozt/workflow.py').read_text())
        guard = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id=='_SBATCH_GUARD' for t in n.targets))
        config = {'dir':str(Path('step3_charged').resolve()), 'fanout':'def-*', 'empty_fanout_manifest':'charged_manifest.json'}
        for stale in (False, True):
            if stale: Path('step.conf').write_text('changed')
            result = subprocess.run([sys.executable, '-B', '-', '--config64', base64.b64encode(json.dumps(config).encode()).decode()], input=guard, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            payload=json.loads(result.stdout.split('__TF_RESULT__ ')[1])
            self.assertEqual(payload['ok'], not stale)
            self.assertEqual(payload['jobids'], [])
    def test_old_q0_rejected(self):
        Path('step3_charged/def-X_q+0').mkdir(parents=True)
        with self.assertRaises(SystemExit): self.s3['main']()
        self.assertFalse(Path('step3_charged/charged_manifest.json').exists())
    def test_strict_marker_no_png(self):
        import glob
        ck=functions(ROOT/'autozt/_collector_remote.py', ['ck_plot'])
        ck['glob']=glob
        Path('stale.png').write_text('old')
        self.assertFalse(ck['ck_plot']('.', {'done_marker':'analysis_complete.json','strict_done_marker':True})[0])
        self.assertTrue(ck['ck_plot']('.', {})[0])


if __name__ == "__main__":
    rc = main()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Repair)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(rc or (0 if result.wasSuccessful() else 1))
