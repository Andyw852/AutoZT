# eph-qe-cpu：QE–Perturbo 电子–声子耦合

此技能以 Quantum ESPRESSO、Wannier90 与 Perturbo 计算并检查电子–声子矩阵元数据链：

`SCF → NSCF/Wannier90` 与 `DFPT phonon` 并行完成，随后 `qe2pert.x → perturbo.x ephmat`。

## 适用范围和边界

- 当前生成器支持任意元素数和原子数的 **3D 周期 VASP5 POSCAR**，包括 Direct/Cartesian 坐标和 Selective dynamics；QE、Wannier90、DFPT、qe2pert 和 Perturbo 的前缀、原子种类、质量、赝势和文件名均由材料上下文动态生成。
- `PSEUDO_FILES`、`ATOMIC_MASSES` 和 `WANNIER_PROJECTIONS` 可在项目级 `step.conf` 显式覆盖。Si 两原子金刚石仍自动使用 `si`、`Si.upf`、Si 标准质量和 `Si:sp3`，仅作为回归案例；其他体系若不能安全推断 Wannier 投影会明确报错，不会默默套用 Si 参数。
- 这是 **3D、非 SOC、非自旋极化** 的基础数据链。金属展宽、极性体系的介电/NAC、磁性、SOC、强关联和超导后处理需要额外参数与独立收敛验证，不能因为流程成功就视为物理结果已通用。
- 流程输出只证明电子–声子矩阵元数据链和插值工件存在，**不**自动产生或证明 `lambda`、`alpha2F`、`omega_log`、`Tc`。普通 phonon DOS 不能替代 `alpha2F`，金属超导筛选还需要经收敛的 Fermi 面积分和 EPW/Eliashberg 后处理。

多元素项目在项目级 `templates/step.conf` 中按元素配置，例如：

```ini
QE_PREFIX = c2si
PSEUDO_FILES = C:C.pbe-n-kjpaw_psl.1.0.0.UPF Si:Si.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_MASSES = C:12.011 Si:28.085
WANNIER_PROJECTIONS = C:sp2 Si:sp3
NUM_BANDS = 24
NUM_WANN = 12
DFT_BAND_MIN = 1
DFT_BAND_MAX = 24
```

`NUM_WANN` 必须与所选投影和能窗相容，且不能大于 `NUM_BANDS`；生成器只负责语法和流程契约，不替用户判断投影是否物理合理。

## 必备环境

Perturbo 3.0.0 必须在兼容 QE 源码树内编译：其配置文件明确要求 QE 7.3。Si 冒烟测试使用官方 QE 7.3 源码在 jzzn 私有目录 `~/software/AutoZT/qe-7.3-perturbo/` 中构建。jzzn 的 HDF5 Fortran module 与 Intel 2021 ABI 不兼容，因此本技能固定使用已验证的 `gcc/10.2.0 + openmpi/4.0.5 + hdf5/1.14.6`，并加载 `aocl/4.1.0`、`fftw/3.3.8` 及 AOCC 编译器运行库；QE 使用 4 MPI 进程，jzzn 小测试中的 `qe2pert.x`/Perturbo 固定串行运行（`OMP_NUM_THREADS=1`），并预加载 MKL 的 LP64 BLAS 符号绕开现有 AOCL BLIS 的 `zdotc` 崩溃。Wannier90 版本须不低于 3，并在作业环境提供：

- `QE_BIN`：`pw.x`、`ph.x`、`pw2wannier90.x` 所在目录；
- `WANNIER_BIN`：`wannier90.x` 所在目录；
- `PERTURBO_BIN`：`qe2pert.x`、`perturbo.x` 所在目录；
- `PSEUDO_DIR`：包含 manifest/`step.conf` 中列出的所有 UPF 文件的目录。

`S0_env` 是唯一环境探针：它在 4 核 Slurm 作业中检查以上可执行文件和赝势，并只在全部满足时写 `preflight.ok`。缺失时后续步骤会保持等待，`qe_env.json` 给出精确缺项；不要绕过它手动提交。

## 步骤与恢复

| 步骤 | 产物 | 完成条件 |
| --- | --- | --- |
| S0_env | `qe_env.json` | 全部软件与 manifest 列出的赝势可用 |
| S1_scf | QE `tmp/<prefix>.save` | `JOB DONE` |
| S2_wannier | Wannier 矩阵 | 三个程序均成功 |
| S3_ph | DFPT `PH_QGRID` 指定网格数据 | `ph.x` 成功且保存的 q 网格与配置一致 |
| S4_qe2pert | `eph.h5` | HDF5 文件非空 |
| S5_ephmat | `pert.out` | Perturbo 成功 |
| S6_summary | `eph_summary.json` | 所有上游工件存在 |

环境/模块或路径错误先修安装和测试项目模板，再对失败步骤执行 `autozt -tt eph-qe-cpu -p <材料> -j <步骤> retry`，检查生成输入后再 `start`。`retry` 保留已有产物；`rerun` 会删除步骤目录，只应在用户明确要求时使用。

`qe2pert.x` 按官方接口生成 `<prefix>_epr.h5`；本技能在 S4 将其复制为稳定别名
`eph.h5`，确保 AutoZT 的声明产物、fetch 和 S5 输入不依赖材料前缀。

## 参数与收敛测试

`submit_qe_test.tpl` 是本技能专用的 4 核、`cpu192` 模板；它不触碰全局 VASP 提交模板。所有电子结构、Wannier、DFPT 和 Perturbo 点位参数都从合并后的 `step.conf` 读取，并在每个步骤写入 `eph_parameters.json`。因此不同项目只需提供 POSCAR、覆盖项目级 `templates/step.conf`，不用改生成器。

本技能在 `skill.yaml` 中关闭 AutoZT 的通用 `hang_check`：该检测器使用 VASP 的 `OUTCAR/OSZICAR` 指纹，而 QE 的 `ph.x` 主要更新 `ph.out`、`dyn*.xml` 和 `dvscf`。对长时间 DFPT 作业沿用 VASP 指纹会把正常计算误判为无输出并取消重交；S3 的状态应以 AutoZT 作业状态和 QE 日志/产物判读。

默认值仍是安装/接口冒烟档：`SCF/NSCF=4^3`、`DFPT=2^3 q + 4^3 k`、Perturbo 单个 Gamma 点和 1 条带。这些默认值只能验证程序链路。新材料应在项目级配置中给出赝势、投影、能窗和适合的带数。

官方 Perturbo 3.0.0 `tests/epr_computation/epr9` 对照档为：

```text
ECUTWFC=40, ECUTRHO=320
SCF_KGRID=12 12 12, NSCF_KGRID=2 2 2, Q2PERT_KGRID=2 2 2
PH_QGRID=2 2 2, PH_KGRID=6 6 6
NUM_BANDS=16, NUM_WANN=8
DIS_NUM_ITER=500, DIS_WIN=[-100,17.2], DIS_FROZ=[-100,9], WANNIER_NUM_ITER=10000
PERT_LIST_MODE=official, PERT_BAND_MIN=1, PERT_BAND_MAX=1, PHFREQ_CUTOFF=1
```

可调参数包括：

- `NUM_BANDS` 与 `DIS_WIN_*`：扩大 Wannier 电子带范围；`NUM_WANN` 保持与投影一致。
- `SCF_KGRID`、`NSCF_KGRID`、`PH_KGRID`、`PH_QGRID`、`Q2PERT_KGRID`：电子/声子网格。
- `Q2PERT_KGRID` 是 `qe2pert.x` 读取的电子 DFT k 网格，必须与 `NSCF_KGRID` 一致；同时每个方向都必须满足 `NSCF_KGRID[i] % PH_QGRID[i] == 0`（例如 q=2 可配 k=4，q=8 可配 k=8 或 16）。不满足整除关系时 q/k 网格不相容，`qe2pert.x` 可能将多个 q 点映射到同一索引并报 `nc already filled`。
- `PH_SEARCH_SYM` 是 `ph.x` 的模式对称性分析开关，不控制 q 网格的不可约化；它与 QE-PW 的 `nosym/noinv` 独立。S3 正常产物是“不可约 q 文件 + 每个文件的 q-star”，`qe2pert.x` 会将它们合并为完整网格。不要为了让 `dyn0` 中的文件数等于 `nq1*nq2*nq3` 而关闭 SCF 对称性：在 QE 7.3/Perturbo 3.0.0 组合下，密网格会重复写入 q/-q 并触发 `nc already filled`。SCF 保持对称性；NSCF/Wannier 使用显式完整 k 网格并关闭其 k→−k 约简。
- `ECUTWFC/ECUTRHO`：截断能扫描，例如 40/50/60 Ry（电荷密度按 8 倍）。
- `SCF_CONV_THR`、`NSCF_CONV_THR`、`TR2_PH`：QE/DFPT 自洽阈值；默认是小核冒烟值，严格对照官方 epr9 时分别设为 `1e-15`、`3e-12`、`1e-17`。
- `OCCUPATIONS=smearing`、`SMEARING`、`DEGAUSS`：QE 电子展宽扫描；官方 Si 档保持 `fixed`。
- `PERT_DELTA_SMEAR`：Perturbo 的 delta 函数展宽（meV）；`PHFREQ_CUTOFF` 是声子频率 cutoff，不要混为同一参数。
- `PERT_LIST_MODE=grid` 配合 `PERT_KGRID/PERT_QGRID` 生成致密均匀点；`official` 复现官方 epr9 的 25 个随机 k/q 点。

建议在 jzzn 小核上至少做以下独立 case，并比较 `ephmat`/`ephmat.yml`、Wannier 解缠警告、声子频率和文件摘要：

| 轴 | case | 目的 |
| --- | --- | --- |
| 截断能 | 40/50/60 Ry，`ECUTRHO=8×ECUTWFC` | 总能、声子和矩阵元稳定性 |
| 展宽 | `DEGAUSS=0.01/0.02/0.04 Ry`（`OCCUPATIONS=smearing`），`PERT_DELTA_SMEAR=5/10/20 meV` | 电子占据与 Perturbo 展宽敏感性 |
| k/q 网格 | `NSCF/PH_KGRID=4^3/6^3/8^3`，`PH_QGRID=2^3/3^3/4^3` | Wannier、DFPT 与插值收敛 |
| 电子带窗 | `NUM_BANDS=16/24`，同时检查 `dis_win_max` | 解缠是否稳定、插值误差是否下降 |

收敛判据不是“作业退出码为 0” בלבד：Wannier `wout` 中不能再出现 `Disentanglement convergence criteria not satisfied`；Γ 点声学残差应接近数值零；矩阵元和网格/展宽变化应达到预先设定的相对容差。Perturbo `ephmat` 结果只证明电子–声子数据链和插值，不可解释为超导 `lambda/Tc`。
