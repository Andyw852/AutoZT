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
_ACT_RE = re.compile(r"^[ \t]*conda[ \t]+activate[ \t]+(\S+)", re.M)
_SH_RE = re.compile(r"^[ \t]*source[ \t]+(\S*profile\.d/conda\.sh)", re.M)


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
    "qos": ("hint: the job is queued by a QOS per-user job/submit limit, not by cluster load. "
            "That QoS only allows a fixed number of jobs for one user; the slots are currently "
            "held by other materials. Check 'sacctmgr show qos format=Name,MaxJobsPU,MaxSubmitPU' "
            "and either wait for slots, move this step back to a QoS with a larger limit, or "
            "submit fewer configurations at a time. Do not lower the core count for this."),
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


# ---- QOS 作业数 / 提交数上限造成的排队（不是故障）--------------------------------
# 2026-10-23 实测（LS 的 S7）：jzzn 的 premium QoS MaxJobsPU=5，把 2D VASP 模板全局切
#   premium 后，LS 的 6 个 S7 帧全部 PD(QOSMaxJobsPerUserLimit)、--start 直接 N/A；
#   而 regular 的 MaxJobsPU=50。同一个用户的多个项目在互相挤占 premium 的 5 个槽位。
# 以前程序只把原因原样显示成 PD(...)，不说"为什么"和"怎么办"，只能靠人记。
QOS_LIMIT_REASONS = (
    "QOSMaxJobsPerUserLimit",
    "QOSMaxSubmitJobPerUserLimit",
    "QOSMaxJobsPerAccountLimit",
    "QOSMaxSubmitJobPerAccountLimit",
    "QOSMaxGRESPerUser",
)


def is_qos_limit_reason(reason):
    """排队原因是不是 QOS 的作业数/提交数上限（槽位占满；不是配置错误，也不是没资源）。"""
    r = (reason or "").strip().strip("()").strip()
    return any(r.startswith(p) for p in QOS_LIMIT_REASONS)


def qos_limit_message(reason="", qos="", in_use=0, limit=0):
    """QOS 槽位占满时的中文说明（英文走 hint('qos')）。

    reason：squeue 的 %R（如 QOSMaxJobsPerUserLimit）
    qos：撞限的 QoS 名（采集带 %q 时才有）
    in_use/limit：本账号在该 QoS 已占用的作业数 / 该 QoS 的 MaxJobsPU 或 MaxSubmitPU
    """
    r = (reason or "").strip().strip("()").strip()
    kind = "提交数" if "Submit" in r else "作业数"
    who = "账号" if "Account" in r else "用户"
    head = "排队原因 %s：该 QoS 的每%s%s上限已满" % (r, who, kind)
    if qos:
        head += "（QoS=%s" % qos
        if in_use or limit:
            head += "，本账号已占 %s/%s" % (in_use or "?", limit or "?")
        head += "）"
    body = ("\n      这不是集群没资源，也不是配置错误——同一个 QoS 里的其它材料占满了槽位。"
            "可以：① 等槽位（作业结束就释放）；② 把该步换到限额更大的 QoS"
            "（如 regular 的 MaxJobsPU 通常大得多）；③ 少提交几帧、分批跑。"
            "\n      别为这个降核数：核数不解决槽位上限，还会改变并行设置（铁律 9）。")
    return head + body


def qos_limit_groups(jobs):
    """把 QOS 限流导致的排队按 QoS 归并 -> {qos: [job, ...]}（qos 可能为空串）。"""
    out = {}
    for j in jobs or []:
        if (j.get("state") or "") != "PD":
            continue
        if not is_qos_limit_reason(j.get("info")):
            continue
        out.setdefault((j.get("qos") or "").strip(), []).append(j)
    return out


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


# conda env 存在的几个常见根。2026-09-27 修：原来只探 $HOME/miniconda3/envs/<env>，
#   把 $HOME/.conda/envs/<env>（conda create 的默认位置之一）与绝对 venv 路径全判成
#   "不存在" —— 例如 hanhai25 的 atomate2_p_a 就在 ~/.conda/envs，导致持续误报。
_CONDA_ENV_ROOTS = ("$HOME/miniconda3/envs", "$HOME/.conda/envs",
                    "$HOME/anaconda3/envs", "/opt/miniconda3/envs")


def conda_probe_command(acts):
    """给远端的存在性探测命令。

    conda env 名：依次探 _CONDA_ENV_ROOTS 下的 <env>；
    目标是绝对路径或 ~ 开头（venv）：直接按路径判定。
    """
    parts = []
    for i, (kind, target) in enumerate(acts):
        if kind == "sh":
            parts.append("if [ -f %s ]; then echo OK-%d; else echo MISS-%d; fi"
                         % (target, i, i))
        elif str(target).startswith(("/", "~")):
            parts.append("if [ -e %s ]; then echo OK-%d; else echo MISS-%d; fi"
                         % (target, i, i))
        else:
            _tests = " || ".join('[ -d "%s/%s" ]' % (r, target)
                                 for r in _CONDA_ENV_ROOTS)
            parts.append("if %s; then echo OK-%d; else echo MISS-%d; fi"
                         % (_tests, i, i))
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
