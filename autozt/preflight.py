# -*- coding: utf-8 -*-
"""preflight.py —— 提交前自检的纯函数部分（可单测、不碰集群）。

这些检查都是真机踩出来的坑（2026-09-15 一连串"切集群后静默失败"）：
  · 资源超配 → 作业永远 PD(PartitionConfig)，日志不提示
  · 分区名与目标集群不符 → 同上
  · 提交模板激活的 conda 环境/路径只在源集群存在 → 计算阶段无日志死亡
把"从 submit.sh 里解析什么、发现什么情况该说什么话"抽成纯函数，
远端探测（存在性）留在 workflow 里，便于 CI 覆盖。
"""
import re

from autozt import i18n as _i18n

_PART_RE = re.compile(r"^#SBATCH\s+--partition=(\S+)", re.M)
_CPT_RE = re.compile(r"^#SBATCH\s+--cpus-per-task=(\d+)", re.M)
_NPN_RE = re.compile(r"^#SBATCH\s+--ntasks-per-node=(\d+)", re.M)
_ACT_RE = re.compile(r"^\s*conda\s+activate\s+(\S+)", re.M)
_SH_RE = re.compile(r"^\s*source\s+(\S*profile\.d/conda\.sh)", re.M)


def parse_submit(text):
    """从 submit.sh 文本里取出分区与资源声明（缺省 0/空串）。"""
    cpt = _CPT_RE.search(text or "")
    npn = _NPN_RE.search(text or "")
    part = _PART_RE.search(text or "")
    return {"partition": part.group(1) if part else "",
            "cpus_per_task": int(cpt.group(1)) if cpt else 0,
            "ntasks_per_node": int(npn.group(1)) if npn else 0}


# 英文版提示（按"检查种类"取用；缺失时回退中文，绝不空输出）
EN_HINTS = {
    "resources": ("hint: the submit template asks for too many cores; the job then waits "
                  "in PD(PartitionConfig) forever and no log line says why. Lower "
                  "cpus-per-task in a project-level copy of the template."),
    "partition": ("hint: the submit template requests a partition the target cluster does "
                  "not use. Copy the template into project_setting/templates and set a "
                  "partition that cluster accepts."),
    "conda": ("hint: the template activates a conda environment or conda.sh path that does "
              "not exist on the target cluster. Activation fails silently, the preparation "
              "stage may still pass and the compute stage exits without a message. Point "
              "CONDA_SH/CONDA_ENV in that step's step.conf at the cluster's environment."),
    "fail_step": ("hint: a FAIL step is not in the DAG's active set, so -f alone does not "
                  "pick it up. Resubmit it explicitly with -j STEP start -f."),
}


def hint(kind, zh_text="", cluster=""):
    """按语言返回提示文本（英文缺失时回退中文，绝不空输出）。"""
    if not _i18n.is_en():
        return zh_text or EN_HINTS.get(kind, "")
    return EN_HINTS.get(kind, zh_text)


def check_resources(res, max_cpus):
    """资源申请超过集群上限时返回提示文本，否则 None。"""
    if not max_cpus:
        return None
    need = int(res.get("cpus_per_task") or 0) * max(int(res.get("ntasks_per_node") or 0), 1)
    if need <= int(max_cpus):
        return None
    return ("提示：提交模板申请 %d 核（cpus-per-task=%d × ntasks-per-node=%d）超过集群上限 %d。"
            % (need, int(res.get("cpus_per_task") or 0),
               int(res.get("ntasks_per_node") or 0), int(max_cpus)))


def check_partition(res, declared, cluster_partitions=()):
    """分区与集群声明不一致（或集群未声明）时返回提示文本，否则 None。"""
    want = (res.get("partition") or "").strip()
    declared = (declared or "").strip()
    if not want:
        return None
    if declared and want != declared:
        return "提示：提交模板要求 --partition=%s，集群声明的是 %s。" % (want, declared)
    if not declared:
        extra = ("该集群自带模板用的是：%s。" % ", ".join(sorted(cluster_partitions))
                 if cluster_partitions else "")
        return "提示：提交模板要求 --partition=%s，但集群未声明分区。%s" % (want, extra)
    return None


def parse_conda_activations(text):
    """取出模板里引用的 conda 环境名与 conda.sh 路径（保持出现顺序、去重）。"""
    hits = []
    for m in _ACT_RE.finditer(text or ""):
        hits.append((m.start(), ("env", m.group(1))))
    for m in _SH_RE.finditer(text or ""):
        hits.append((m.start(), ("sh", m.group(1))))
    hits.sort(key=lambda x: x[0])          # 保持模板里的出现顺序
    out, seen = [], set()
    for _pos, key in hits:
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def conda_probe_command(acts):
    """给远端的存在性探测命令（conda env 用 $HOME/miniconda3/envs/<env> 判定）。"""
    parts = []
    for i, (kind, target) in enumerate(acts):
        if kind == "sh":
            parts.append("test -f %s && echo OK-%d || echo MISS-%d" % (target, i, i))
        else:
            parts.append("ls -d \"$HOME\"/miniconda3/envs/%s >/dev/null 2>&1 && echo OK-%d"
                         " || echo MISS-%d" % (target, i, i))
    return " ; ".join(parts)


def missing_activations(acts, probe_output):
    """根据探测输出挑出缺失项。"""
    return [t for i, (_k, t) in enumerate(acts)
            if ("MISS-%d" % i) in (probe_output or "")]


def conda_missing_message(missing, cluster, suggested=""):
    return ("提示：提交模板引用的环境/路径在集群 %s 上不存在：%s。"
            "激活失败后脚本会继续跑，计算阶段可能静默死亡；"
            "把该步 step.conf 的 CONDA_SH/CONDA_ENV 改成该集群的值%s。"
            % (cluster, ", ".join(missing),
               ("（配置里写的是 %s）" % suggested) if suggested else ""))
