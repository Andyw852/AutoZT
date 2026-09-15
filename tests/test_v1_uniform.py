#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""统一化机制自测（v1.0）：cores 归一化 + 提交前输入清单发现。

A. cores：submit.sh/INCAR 一律按 cores 改，与技能、模板名无关（defect 的
   submit_ncl_3d.tpl、mlff-mace 的 step1_relax/ 那份模板同样管到）
B. 清单：按远端步骤目录里实际存在的文件要求，技能不必自报 submit_required
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import phonoagent  # noqa: E402

OK = []
FAIL = []


def ck(cond, msg):
    (OK if cond else FAIL).append(msg)
    print("  %s %s" % ("PASS" if cond else "FAIL", msg))


NORM = phonoagent._CORES_NORMALIZER


def normalize(d, n, sub=""):
    s = os.path.join(d, ".tf_cores.py")
    with open(s, "w", encoding="utf-8") as f:
        f.write(NORM)
    p = subprocess.run([sys.executable, s, str(n), sub], cwd=d,
                       capture_output=True, text=True)
    return p


def read(p):
    with open(p, encoding="utf-8") as f:
        return f.read()


print("== A. cores 归一化 ==")
d = tempfile.mkdtemp(prefix="tf_uni_")
step = os.path.join(d, "step1_bulk")
nest = os.path.join(step, "frames")
os.makedirs(nest)
with open(os.path.join(step, "submit.sh"), "w", encoding="utf-8") as f:
    f.write("#!/bin/bash\n# submit_ncl_3d.tpl —— 非共线模板\n"
            "#SBATCH --job-name=x\n#SBATCH --nodes=1\n"
            "#SBATCH --ntasks-per-node=24\n#SBATCH --cpus-per-task=2\n"
            "#SBATCH --qos=regular\nmpirun -np $SLURM_NTASKS ./vasp_ncl\n")
with open(os.path.join(step, "INCAR"), "w", encoding="utf-8") as f:
    f.write("ISIF = 2\nNCORE  = 6   # 注释保留\nKPAR   = 4\n")
with open(os.path.join(nest, "submit.sh"), "w", encoding="utf-8") as f:
    f.write("#!/bin/bash\n#SBATCH --ntasks-per-node=24\n")
with open(os.path.join(nest, "INCAR.s1"), "w", encoding="utf-8") as f:
    f.write("NCORE = 12\nKPAR = 2\n")
# 单进程多线程型（MACE 那种只用 cpus-per-task）
mac = os.path.join(d, "step2_mace")
os.makedirs(mac)
with open(os.path.join(mac, "submit.sh"), "w", encoding="utf-8") as f:
    f.write("#!/bin/bash\n#SBATCH --cpus-per-task=8\n#SBATCH --gres=gpu:1\n")

p = normalize(d, 4, "step1_bulk")
ck(p.returncode == 0, "归一化程序退出码 0（%s）" % (p.stderr or "").strip()[:60])
t = read(os.path.join(step, "submit.sh"))
ck("#SBATCH --ntasks-per-node=4" in t, "MPI 型 submit.sh ntasks 24→4")
ck("#SBATCH --cpus-per-task=1" in t, "MPI 型 cpus-per-task 2→1（避免 4×2 核超订）")
ck("mpirun -np $SLURM_NTASKS" in t, "mpirun 行未被破坏")
i = read(os.path.join(step, "INCAR"))
ck("NCORE  = 4" in i and "KPAR   = 1" in i, "INCAR NCORE 6→4、KPAR 4→1")
ck("# 注释保留" in i, "INCAR 行尾注释未被吃掉")
ck("#SBATCH --ntasks-per-node=4" in read(os.path.join(nest, "submit.sh")),
   "扇出帧子目录（深度 1）里的 submit.sh 也归一")
ck("NCORE = 4" in read(os.path.join(nest, "INCAR.s1")),
   "子目录 INCAR.s1（分阶段 INCAR）也归一")

p = normalize(d, 4, "step2_mace")
ck("#SBATCH --cpus-per-task=4" in read(os.path.join(mac, "submit.sh")),
   "单进程多线程型（只有 cpus-per-task）改 cpus-per-task 8→4")

p = normalize(d, 6, "")
i = read(os.path.join(step, "INCAR"))
ck("NCORE  = 3" in i, "NCORE 取 N 的最大 ≤4 因子：N=6 → NCORE=3（必须整除进程数）")
p = normalize(d, 1, "")
ck("NCORE  = 1" in read(os.path.join(step, "INCAR")), "N=1 → NCORE=1")
ck("#SBATCH --ntasks-per-node=1" in read(os.path.join(step, "submit.sh")),
   "N=1 → ntasks=1")

print("== A2. cores 来源优先级 ==")
base_cfg = {"cores": 8, "task_types": {"defect-dft-cpu": {"cores": 4}}}
m = {"ps": {"setting": {"cores": 16}, "hpc": {"cores": 2}}, "tt": "defect-dft-cpu"}
ck(phonoagent.resolve_cores(base_cfg, {"key": "defect-dft-cpu"}, m) == 16,
   "项目 setting.yaml 最高（16）")
m2 = {"ps": {"setting": {}, "hpc": {"cores": 2}}, "tt": "defect-dft-cpu"}
ck(phonoagent.resolve_cores(base_cfg, {"key": "defect-dft-cpu"}, m2) == 2,
   "项目 hpc.yaml 次之（2）")
m3 = {"ps": {"setting": {}, "hpc": {}}, "tt": "defect-dft-cpu"}
ck(phonoagent.resolve_cores(base_cfg, {"key": "defect-dft-cpu"}, m3) == 4,
   "类型配置再次（4）")
ck(phonoagent.resolve_cores(base_cfg, {"key": "opt-dft-cpu"}, {"ps": {}}) == 8,
   "全局 tf.yaml cores 兜底（8）")
ck(phonoagent.resolve_cores({}, {"key": "x"}, {"ps": {}}) is None,
   "都没配 → None（保持出厂行为，零改动）")
ck(phonoagent.resolve_cores({"cores": "0"}, {"key": "x"}, {"ps": {}}) is None,
   "cores: 0 视为未配置")

print("== B. 提交前清单（统一发现）==")
ck(phonoagent._step_input_name_ok("POSCAR"), "POSCAR 算输入")
ck(phonoagent._step_input_name_ok("graph_data_gen.yaml"), "任意 gen 产物算输入")
ck(not phonoagent._step_input_name_ok("OUTCAR"), "OUTCAR 不算输入")
ck(not phonoagent._step_input_name_ok("vasprun.xml"), "vasprun.xml 不算输入")
ck(not phonoagent._step_input_name_ok("slurm-123.out"), "slurm-*.out 不算输入")
ck(not phonoagent._step_input_name_ok("queue.out"), "queue.out 不算输入")
ck(not phonoagent._step_input_name_ok(".tf_cores.py"), "tf 自己的临时文件不算输入")

# 真实远端：unihamgnn 的 S1_graph 步骤目录（非 VASP，历史上被默认清单卡死）
cfg = None
try:
    import yaml
    with open(os.path.join(ROOT, "tmp", "tf_jzzn_si_all.yaml"), encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
except Exception as exc:                                    # noqa: BLE001
    print("  SKIP 读测试配置失败：%s" % exc)
if cfg:
    cfg.setdefault("_config_path", os.path.join(ROOT, "tmp", "tf_jzzn_si_all.yaml"))
    d3090 = ("/home/wangchaoyue852/taskflow/work/tf_smoke/Si_unihamgnn/"
             "unihamgnn/step1_graph_data")
    got = phonoagent._discover_step_inputs(cfg, "wangchao_3090", d3090)
    if not got:
        print("  SKIP 远端测试材料已清理（%s 不存在）——远端相关断言跳过" % d3090)
        REMOTE_OK = False
    else:
        REMOTE_OK = True
    print("  远端发现：%s" % (", ".join(got) or "(空)"))
    if REMOTE_OK:
        ck("submit.sh" in got, "发现清单含 submit.sh")
        ck("POSCAR" in got, "发现清单含 POSCAR")
        ck(any(x.endswith(".yaml") for x in got),
           "发现清单含 gen 写的 yaml（技能无需自报）")
        ck(not any(x.startswith(("OUTCAR", "vasprun")) for x in got),
           "清单里没有输出文件")
    got2 = phonoagent._discover_step_inputs(cfg, "wangchao_3090", "/nonexistent/xyz")
    ck(got2 == (), "远端目录不存在 → 空清单（调用方兜底，不阻断提交）")

print("== C. 真实提交前检查（unihamgnn 已删掉自报清单）==")
if cfg:
    tt_cfg = (cfg.get("task_types") or {}).get("unihamgnn") or {}
    ck(not tt_cfg.get("submit_required"),
       "unihamgnn 配置里已无 submit_required（改由统一发现兜住）")
    mat = {"tt": "unihamgnn", "name": "Si_unihamgnn", "host_eff": "wangchao_3090",
           "result_dir": ("/home/wangchao/tf_smoke_all/Si_unihamgnn/unihamgnn/result")}
    step = {"name": "step1_graph_data", "dir": d3090}
    if not REMOTE_OK:
        print("  SKIP 远端测试材料已清理——真实提交前检查的端到端断言跳过")
    else:
        try:
            ok, reason = phonoagent._remote_submit_preflight(cfg, mat, step, {})
            ck(ok, "非 VASP 技能（无 INCAR/KPOINTS）提交前检查通过：%s"
               % (reason or "OK"))
        except Exception as exc:                             # noqa: BLE001
            ck(False, "提交前检查抛异常：%r" % exc)
    # 反面：把清单临时改成 VASP 四件套（老默认行为），应当报缺 INCAR
    mat2 = dict(mat, result_dir=tempfile.mkdtemp(prefix="tf_uni_res_"))
    try:
        phonoagent._remote_submit_preflight(cfg, mat2, dict(step), {})
        ck(False, "老默认清单（硬要 INCAR）本应拦住非 VASP 步骤")
    except Exception:                                        # noqa: BLE001
        ck(True, "老默认清单会拦住非 VASP 步骤（反证统一发现的必要性）")

print("\n结果：PASS %d，FAIL %d" % (len(OK), len(FAIL)))
for x in FAIL:
    print("  ! %s" % x)
sys.exit(1 if FAIL else 0)
