# B（AutoZT）提交批次方案（2026-09-18）

> 目的：把工作区未提交的改动按功能分批提交，便于审阅与回退，并为论文引用的
> 版本可追溯（每批一条独立提交信息）。本文件只是方案 + 可直接复制执行的命令。

## 0. 当前状态（生成时刻）

- 仓库：`/home/wangchao/software/AutoZT`，HEAD `7e2f7ff`（master，ahead 5）
- 未提交：**44 个已跟踪文件被修改 + 10 个未跟踪**（含另一个会话正在开发的
  `skill/eph-qe-cpu/`、`skill/zt-dft-cpu/`、`docs/figures/fig1_autozt_positioning_plugin.svg`）
- ⚠️ **另一个会话此刻仍在同一工作区写代码**（01:25 建 eph-qe-cpu、01:37 改 TASKFLOW.md、
  01:4x 建 zt-dft-cpu）。开工前先确认它停下了，否则会把半成品一起提交。

## 1. 建议批次

| 批次 | 主题 | 文件 |
|---|---|---|
| 0 | 清理误入库的历史残留（可选，建议单独一条提交；文件已归档到 `tmp/_archive_20260918/`） | `hfeshell_tmp`（根目录 POSCAR 残留）、`skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl.bak2-20260906` |
| 1 | 内核/agent 层：纠错、报告、MCP、工作流 | `autozt/agent_cli.py` `autozt/agent_service.py` `autozt/corrections.py` `autozt/mcp.py` `autozt/report.py` `autozt/workflow.py` `bin/autozt` |
| 2 | 公共弛豫池 + ke-dft-cpu 物理修复 | `skill/_common/opt/relax_common.py` `skill/_common/opt/checks_relax.py` `skill/_common/_corrections/README.md`；`skill/ke-dft-cpu/step3_uniform/incar_uniform_{2d,3d}.tpl` `step3b_uniform_full/*` `step3c_uniform_offgrid/`（新） `step8.3_output/gen_step13_output.py` `step8.3_output/overlap_ratio_guard.py`（新） `step8_amset/gen_step10_amset.py` `step8.4_amset2d/*`（含新 `overlap_preflight.py`）；`skill/kl-dft-cpu/templates/step1_std_opt/incar_2d.tpl`；**本次从 A 移植 P49**：`skill/kl-dft-cpu/gen_step6_kappa.py`（完成判据只认 CONV/RTA，去掉 small-grain limit 误判，并清掉重复代码块）+ `skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl`（OMP_STACKSIZE=1G 防 SIGSEGV；_sg 改记诊断字段） |
| 3 | schema 2 元数据与参数说明一致性（含本次修复与守卫测试） | 11 个 `skill/*/skill.yaml` + `tests/suite_io_schema.py`（新） `tests/suite_s6_marker.py`（新，P49 完成判据守护） `tests/test_suites.py` `scripts/export_params.py`（新） `docs/PARAMETERS.md`（新） |
| 4 | 文档 | `README.md` `TASKFLOW.md` `mkdocs.yml` `docs/{DELIVERY,agent-cli,index,mcp,skills}.md` `docs/HANDOVER-2D-kl-20260917.md`（新） |
| 5 | 并行会话的新技能（**建议等它写完再提交**） | `skill/eph-qe-cpu/` `skill/zt-dft-cpu/` `docs/figures/fig1_autozt_positioning_plugin.svg` |

## 2. 可直接执行的命令（逐批）

```bash
cd /home/wangchao/software/AutoZT

# 批次 0（清理；两个文件已移到 tmp/_archive_20260918/，这里只把删除记进版本库）
git add -u hfeshell_tmp skill/kl-dft-cpu/templates/step6_kappa
git commit -m 'chore: 清理误入库的历史残留（hfeshell_tmp POSCAR、ShengBTE 模板 .bak2）'

# 批次 1
git add autozt/agent_cli.py autozt/agent_service.py autozt/corrections.py autozt/mcp.py autozt/report.py autozt/workflow.py bin/autozt
git commit -m 'feat(agent): 纠错/报告/MCP/工作流层更新（corrections、agent_cli、report）'

# 批次 2（★ 提交前建议先跑：python3 tests/suite_ke_common.py && python3 tests/suite_uniform.py）
git add skill/_common/opt/relax_common.py skill/_common/opt/checks_relax.py skill/_common/_corrections/README.md \
        skill/ke-dft-cpu/step3_uniform/incar_uniform_2d.tpl skill/ke-dft-cpu/step3_uniform/incar_uniform_3d.tpl \
        skill/ke-dft-cpu/step3b_uniform_full skill/ke-dft-cpu/step3c_uniform_offgrid \
        skill/ke-dft-cpu/step8.3_output skill/ke-dft-cpu/step8_amset/gen_step10_amset.py \
        skill/ke-dft-cpu/step8.4_amset2d skill/kl-dft-cpu/templates/step1_std_opt/incar_2d.tpl \
        skill/kl-dft-cpu/gen_step6_kappa.py skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl
git commit -m 'fix(physics): 变胞段能量判据+2D 面内应力口径；ke 真空口径形变势/Voigt 重排/AMSET 二维核；kl S6 完成判据不再接受 small-grain limit + OMP 线程栈修复'

# 批次 3
git add skill/*/skill.yaml tests/suite_io_schema.py tests/suite_s6_marker.py tests/test_suites.py scripts/export_params.py docs/PARAMETERS.md
git commit -m 'docs(skillspec): io_schema 声明默认值与 step.conf 对齐 + 自动核对测试 + 参数表导出'

# 批次 4
git add README.md TASKFLOW.md mkdocs.yml docs/DELIVERY.md docs/agent-cli.md docs/index.md docs/mcp.md docs/skills.md docs/HANDOVER-2D-kl-20260917.md docs/dev/COMMIT-PLAN-20260918.md
git commit -m 'docs: 2D κ 交接说明、CLI/MCP/技能文档更新'

# 批次 5（确认并行会话写完后）
git add skill/eph-qe-cpu skill/zt-dft-cpu docs/figures/fig1_autozt_positioning_plugin.svg
git commit -m 'feat(skills): 新增 eph-qe-cpu（QE/Wannier90/Perturbo）与 zt-dft-cpu'
```

## 3. 提交前必做的三件小事

1. `for t in tests/suite_*.py; do python3 $t >/dev/null || echo FAIL $t; done`（本机无 pytest，直接跑脚本）
2. `python3 scripts/export_params.py --out docs/PARAMETERS.md` 重跑一次，确保参数表与代码同步
3. `git status --porcelain` 再确认一遍没混入别的会话的半成品

## 4. 论文引用

投稿时引用**打完标签的版本**（建议 `git tag -a v1.3.0-paper -m ...`），并 push 到远端 + Zenodo 存档取 DOI；
不要引用带 `+dirty` 的工作目录。参数表的「生成时的仓库版本」行会写明当时提交号，可直接进附录。
