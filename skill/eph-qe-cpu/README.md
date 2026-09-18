# eph-qe-cpu：QE–Perturbo 电子–声子耦合

此技能以 Quantum ESPRESSO、Wannier90 与 Perturbo 计算并检查电子–声子矩阵元数据链：

`SCF → NSCF/Wannier90` 与 `DFPT phonon` 并行完成，随后 `qe2pert.x → perturbo.x ephmat`。

## 适用范围和边界

- v0.1 是 **Si 两原子金刚石结构的小核端到端验证**；生成器会拒绝非 Si/非两原子的 POSCAR，避免把未经验证的元素、质量或赝势映射伪装成通用计算。
- Si 是半导体。本技能完成只表示 `si_epr.h5` 和点位电子–声子矩阵元链条有效，**不**产生或证明金属超导的 `lambda`、`alpha2F`、`omega_log`、`Tc`。
- 常规超导筛选必须另用金属体系并增加经收敛的 Fermi 面积分和 EPW/Eliashberg 后处理；普通 phonon DOS 不可替代 `alpha2F`。

## 必备环境

Perturbo 3.0.0 必须在兼容 QE 源码树内编译：其配置文件明确要求 QE 7.3。Si 冒烟测试使用官方 QE 7.3 源码在 jzzn 私有目录 `~/software/AutoZT/qe-7.3-perturbo/` 中构建。jzzn 的 HDF5 Fortran module 与 Intel 2021 ABI 不兼容，因此本技能固定使用已验证的 `gcc/10.2.0 + openmpi/4.0.5 + hdf5/1.14.6`，并加载 `aocl/4.1.0`、`fftw/3.3.8` 及 AOCC 编译器运行库；QE 与 qe2pert 使用 4 MPI 进程，Perturbo 使用官方 GNU/OpenMP 配置和 4 线程。Wannier90 版本须不低于 3，并在作业环境提供：

- `QE_BIN`：`pw.x`、`ph.x`、`pw2wannier90.x` 所在目录；
- `WANNIER_BIN`：`wannier90.x` 所在目录；
- `PERTURBO_BIN`：`qe2pert.x`、`perturbo.x` 所在目录；
- `PSEUDO_DIR`：含命名为 `Si.upf` 的兼容 Si UPF 文件的目录。

`S0_env` 是唯一环境探针：它在 4 核 Slurm 作业中检查以上可执行文件和赝势，并只在全部满足时写 `preflight.ok`。缺失时后续步骤会保持等待，`qe_env.json` 给出精确缺项；不要绕过它手动提交。

## 步骤与恢复

| 步骤 | 产物 | 完成条件 |
| --- | --- | --- |
| S0_env | `qe_env.json` | 全部软件与 `Si.upf` 可用 |
| S1_scf | QE `tmp/si.save` | `JOB DONE` |
| S2_wannier | Wannier 矩阵 | 三个程序均成功 |
| S3_ph | DFPT 2×2×2 数据 | `ph.x` 成功 |
| S4_qe2pert | `si_epr.h5` | HDF5 文件非空 |
| S5_ephmat | `pert.out` | Perturbo 成功 |
| S6_summary | `eph_summary.json` | 所有上游工件存在 |

环境/模块或路径错误先修安装和测试项目模板，再对失败步骤执行 `autozt -tt eph-qe-cpu -p <材料> -j <步骤> retry`，检查生成输入后再 `start`。`retry` 保留已有产物；`rerun` 会删除步骤目录，只应在用户明确要求时使用。

`qe2pert.x` 按官方接口将 `si_epr.h5` 写到其 QE `outdir`；本技能在 S4
创建步骤目录内的软链接，确保 AutoZT 的声明产物、fetch 和 S5 输入一致。

## Si 测试资源

`submit_qe_test.tpl` 是本技能专用的 4 核、`cpu192` 模板；它不触碰全局 VASP 提交模板。网格为 4×4×4 k 和 2×2×2 q，仅用于安装与接口冒烟测试，不能用于收敛后的材料结论。
