"""提交前自检的纯函数单测（CI 可跑，不需要集群）。

fixture 取自真实踩坑现场：
  · kl-mace-gpu S3（3090 上跑通）：cpu192 / 24 核
  · fc-fit S1（切到 3090 后静默死亡）：cpu192 / 48 核 + atomate2_p_a + /public/home/...
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from phonoagent import preflight as pf  # noqa: E402

OK_SUBMIT = """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR
source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh
conda activate mace-gpu
python fc_fit_driver.py 2>&1 | tee -a fc_build.log
"""

BAD_SUBMIT = """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
cd $SLURM_SUBMIT_DIR
if [ -n "atomate2_p_a" ] && [ -x "atomate2_p_a/bin/python" ]; then
    source "atomate2_p_a/bin/activate"
else
    source /public/home/wangchao/miniconda3/etc/profile.d/conda.sh
    conda activate atomate2_p_a
fi
python fc_fit_driver.py fit fit_config.json
"""


def test_parse_submit_reads_partition_and_resources():
    r = pf.parse_submit(OK_SUBMIT)
    assert r == {"partition": "cpu192", "cpus_per_task": 24, "ntasks_per_node": 1}


def test_resource_check_flags_overcommit_only():
    assert pf.check_resources(pf.parse_submit(OK_SUBMIT), 24) is None
    msg = pf.check_resources(pf.parse_submit(BAD_SUBMIT), 24)
    assert msg and "48" in msg and "24" in msg
    assert pf.check_resources(pf.parse_submit(BAD_SUBMIT), 0) is None  # 未声明上限不提示


def test_partition_check_messages():
    res = pf.parse_submit(OK_SUBMIT)
    assert pf.check_partition(res, "cpu192") is None
    assert "cpu192" in pf.check_partition(res, "")
    assert "gpu" in pf.check_partition(res, "", ["gpu"])
    assert "other" in pf.check_partition(res, "other")


def test_conda_activations_extracted_from_both_branches():
    acts = pf.parse_conda_activations(BAD_SUBMIT)
    assert ("env", "atomate2_p_a") in acts
    assert ("sh", "/public/home/wangchao/miniconda3/etc/profile.d/conda.sh") in acts
    assert acts == pf.parse_conda_activations(BAD_SUBMIT)  # 幂等、无重复


def test_probe_command_and_missing_detection():
    acts = pf.parse_conda_activations(BAD_SUBMIT)
    cmd = pf.conda_probe_command(acts)
    assert "MISS-0" in cmd and "MISS-1" in cmd
    assert pf.missing_activations(acts, "OK-0\nMISS-1") == [acts[1][1]]
    assert pf.missing_activations(acts, "OK-0\nOK-1") == []
    assert "atomate2_p_a" in pf.conda_missing_message(["atomate2_p_a"], "3090", "mace-gpu")


def test_good_template_produces_no_warnings():
    res = pf.parse_submit(OK_SUBMIT)
    assert pf.check_resources(res, 24) is None
    assert pf.check_partition(res, "", ["cpu192"]) is not None  # 未声明分区时给提示（不阻断）
    acts = pf.parse_conda_activations(OK_SUBMIT)
    assert acts == [("sh", "/home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh"),
                    ("env", "mace-gpu")]
    assert pf.missing_activations(acts, "OK-0\nOK-1") == []
