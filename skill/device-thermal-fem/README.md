# device-thermal-fem —— 器件级传热 FEM 独立复核（技能介绍，v0.1）

> 用**第二套独立实现**复算 `device-thermal` 的器件传热，交叉验证其数值结果；
> 求解后端可插拔，**COMSOL 接口已留**（见第 5 节）。

---

## 1. 为什么需要它

- `device-thermal` 是手写的「固定两层 + 2D + 各向异性」有限体积法。**同一份代码无法自证**：
  SOR 求解器与独立 scipy 原型共用同一套源项/边界约定，两者「一致」只说明求解器没写错，
  说明不了物理建模没错。
- 一个独立的 FEM 实现能抓到**系统性**错误。本技能首次运行就抓到一处：
  `device-thermal` 的源项把每单位面积功率 `Q*dx` 加到了沟道的**每个**单元上，
  而沟道被离散成 2 层，导致**实际注入功率是名义值 Q*L 的 2.05 倍**，dT 因此偏高约 2.02 倍。
  该 bug 已经用户批准修复，修正后两套实现差 0.89%。详见第 10 节。
- 同时为将来换更强引擎（任意几何 / 3D / 电-热耦合 / COMSOL）留出可插拔位置。

---

## 2. 它做什么

**不重复定义器件**：直接复用 `device-thermal` 的产物

- `step1_device_spec/device_spec.json` —— 几何 / 焦耳热源 / 边界条件
- `step2_props/thermal_props.json` —— 各层 κx/κy/厚度/热容 + 界面热导 G

然后解**同一个** PDE：`div(k(x,y) grad T) + Q = 0`（各向异性 k，界面 G，底面恒温，
左右接触恒温或绝热），再做四类独立检查。

---

## 3. 四类检查（这是本技能的核心）

| 检查 | 判据 | 为什么能抓错 |
|---|---|---|
| `backend_vs_analytic_1d` | 侧壁绝热时与**解析严格解**相对差 ≤ `FEM_1D_TOL` | 绝热极限下 1D 有闭式解（含热源沿沟道厚度均匀分布），是外部真值，不依赖任何求解器 |
| `source_power_audit` | 实际注入功率 / 名义 `Q*L` 相对差 ≤ `FEM_POWER_TOL` | 直接审计功率口径，**独立于温度场**；单元重复计数在这里立刻现形 |
| `mesh_convergence` | 最细网格 vs 原网格 dT 相对差 ≤ `FEM_MESH_TOL` | 排除「结果依赖网格」 |
| `fem_vs_device_thermal_fvm` | dT 比值与 1 的偏差 ≤ `FEM_FVM_TOL`，并逐点比对中轴剖面 | 两个独立实现对同一物理的交叉校验 |

`verdict = pass`（全通过）/ `review`（有 FAIL），FAIL 时给出 `findings` 定位建议。

---

## 4. 三步流水线（全部 `run: gen`，CPU 秒级~十几秒）

| 步骤 | 名称 | 干什么 | 产物 | 判据 |
|---|---|---|---|---|
| S1 | step1_fem_model | 复用 device-thermal 产物，装配与后端无关的 FEM 模型 | `fem_model.json` | `FEM_MODEL_DONE` |
| S2 | step2_fem_solve | 按 `FEM_BACKEND` 求解（主算例 + 侧壁绝热 1D 自检算例） | `fem_summary.json`, `T_field_fem.json` | `FEM_SOLVE_DONE` |
| S3 | step3_fem_verify | 四类检查 + 结论 | `fem_verification.json` | `FEM_VERIFY_DONE` |

---

## 5. 后端（可插拔）

统一契约：`solve(model, workdir, opts) -> dict`，字段见 `backend_skfem.py` 顶部。

### 5.1 `FEM_BACKEND = skfem`（默认，已实测）

- 开源、**免 license**，纯 Python：`pip install scikit-fem`
- P1 三角元，各向异性 `(kx, ky)`，界面用**等效薄层**（厚度 δ、导热 δ·G，热阻恒为 1/G，δ 对结果不敏感）
- 网格：x 方向 `nx` 单元；y 方向分三层（沟道 / 界面薄层 / 氧化层）
- 实测：与解析 1D 差 **5.5e-6**（相对），网格收敛差 **0.011%**

### 5.2 `FEM_BACKEND = comsol`（**接口已留，未启用**）

`backend_comsol.py` 已按同一契约写好骨架 + 完整启用文档。**COMSOL 是商业软件**：
安装包与 license 都在 COMSOL Access 账号后面、license 与 host ID 绑定，本仓库无法自带。
启用需要：

1. 在目标机装 COMSOL（建议装到 `~/software/AutoZT/comsol`），确认 `bin/comsol` 存在；
2. 配 license（许可证服务器 `LMCOMSOL_LICENSE_FILE=port@host`，或单机 `license.dat`），
   用 `comsol batch -help` 自检；
3. step.conf 设 `FEM_BACKEND=comsol`、`COMSOL_BIN`、`COMSOL_LICENSE_FILE`、`COMSOL_CORES`；
4. Python 侧二选一：`pip install mph`（推荐）走 MPh 直接建模取场；或纯 batch 走
   `.java`/`.mph` + `comsol batch -input ... -output ...` + 导出 CSV 解析。

接上时只需实现 `backend_comsol._solve_via_mph` 或 `_solve_via_batch`，
**S1/S3 与全部检查逻辑无需改动**——它们只认契约里的 `dT_peak_K` / `y_nodes_m` / `T_mid_x_K` 等字段。

未就绪时该后端会给出明确报错并列出上述步骤，**不会静默失败**。

---

## 6. 输入

| 来源 | 文件 | 必需 | 说明 |
|---|---|---|---|
| `skill:device-thermal` | `device_spec.json` | ✅ | 器件几何 / 热源 / 边界 |
| `skill:device-thermal` | `thermal_props.json` | ✅ | 各层 κ/厚度/热容 + 界面 G |
| `skill:device-thermal` | `device_thermal_summary.json` | ⭕ | 有则做二维 FVM 交叉验证，无则只做解析/自洽检查 |

查找链：`params.FEM_FROM_DEVICE_THERMAL_DIR` > `<材料>/device-thermal/{stepdir, result/stepdir}/`。

---

## 7. 输出

| 文件 | 关键字段 |
|---|---|
| `fem_model.json` | `geometry`、`layers`、`interface_G_W_m2K`、`Q_W_m2`、`sink`、`mesh`、`analytic`（`dT_series_K` / `dT_distributed_K`） |
| `fem_summary.json` | `main`（`dT_peak_K`、`P_source_W_per_m`、`source_ratio`、`mesh`、`y_nodes_m`、`T_mid_x_K`）、`ref_1d`（含 `analytic_check.rel_err`） |
| `T_field_fem.json` | 中轴温度剖面 `y_nodes_m` / `T_mid_x_K` |
| `fem_verification.json` | `verdict`、`checks[]`、`findings[]`、`fvm.ratio_fvm_over_fem`、`fvm.profile_diff` |

---

## 8. 参数

| 参数 | 默认 | 含义 |
|---|---|---|
| `FEM_BACKEND` | `skfem` | `skfem` 或 `comsol` |
| `FEM_FROM_DEVICE_THERMAL_DIR` | 空 | 显式指定 device-thermal 产物目录 |
| `FEM_INTERFACE_LAYER_FRAC` | `1e-3` | 界面等效薄层厚度 / 氧化层厚度（热阻恒为 1/G） |
| `FEM_NX` / `FEM_N_CHANNEL` / `FEM_N_INTERFACE` / `FEM_N_OXIDE` | `60`/`8`/`4`/`200` | 网格分辨率 |
| `FEM_DO_1D_CHECK` | `true` | 是否额外解侧壁绝热算例做解析自检 |
| `FEM_MESH_CHECK` | `true` | 是否做网格收敛扫描（30/60/120/240） |
| `FEM_1D_TOL` / `FEM_POWER_TOL` / `FEM_FVM_TOL` / `FEM_MESH_TOL` | `0.01`/`0.02`/`0.05`/`0.01` | 四类检查容差 |
| `COMSOL_BIN` / `COMSOL_LICENSE_FILE` / `COMSOL_CORES` | 空/空/`4` | 仅 COMSOL 后端用 |

---

## 9. 用法与依赖

```bash
# 依赖（本机一次即可；scikit-fem 免 license）
pip install scikit-fem

# 前置：device-thermal 至少跑到 S2_props，做 FVM 交叉验证则要跑到 S3_solve
autozt -tt device-thermal   -p <材料> start

# 本技能
autozt -tt device-thermal-fem -p <材料> -j S1_fem_model init
autozt -tt device-thermal-fem -p <材料> start
```

> PEP 668 环境（externally-managed）需 `pip install --user --break-system-packages scikit-fem`，
> 或把 scikit-fem 装进该集群的 python 环境。

---

## 10. 实测结果（MoS2 参考算例，2026-09-23）

器件：MoS2 0.672 nm / SiO2 300 nm / Si 热沉；L=200 nm、W=1 µm、Q=1e8 W/m²、接触 sink；
κ_MoS2 = (21.425, 4.5) W/mK，κ_SiO2 = 1.4 W/mK，G = 14 MW/m²K。

| 检查 | 结果 | 判据 |
|---|---|---|
| 后端 vs 解析 1D（侧壁绝热） | FEM 28.5576 K vs 解析 28.5575 K，相对差 **5.5e-6** | PASS |
| 源功率审计 | 实际 20.0 W/m = 名义 `Q*L` 20.0 W/m，差 **0** | PASS |
| 网格收敛 | 9.7250 / 9.7301 / 9.7314 / 9.7317 K（30/60/120/240），相对差 **0.011%** | PASS |
| 与 device-thermal FVM 比对（**修正前**） | FVM 19.6355 K vs FEM 9.7307 K，比值 **2.018** | **FAIL** |
| 与 device-thermal FVM 比对（**修正后**） | FVM 9.8177 K vs FEM 9.7307 K，比值 **1.0089** | PASS |

**结论**：独立 FEM 定位到 `device-thermal` 的源项口径错误 ——
它在沟道层的**每个**单元都加了一份 `Q*dx`，而沟道被离散成 2 层，
使实际注入功率变成名义值的约 2.05 倍。

- 证据 1：源功率审计（与温度场无关）显示 FVM 口径下 `P_source/P_nominal = 2.05`；
- 证据 2：FEM 与解析 1D 差 5.5e-6，说明 FEM 侧无系统误差；
- 证据 3：把 FVM 源项乘 0.5 后，其 dT 降到 9.80 K，与 FEM 的 9.73 K 只差 0.7%。

**已修复**（2026-09-23，经用户批准）。`device_common.py` 的源项改为按沟道各单元厚度分摊：

```python
q[:, lay == 0] = (spec["Q_joule_W_m2"] * dx
                  * dy[lay == 0] / ch["thickness_m"])
```

修正后 FVM dT = **9.8177 K**（原 19.6355 K，正好减半），与 FEM 差 **0.89%**，本技能四类检查全 PASS；
`power_uW` = 20.0 µW 现在也与实际注入功率一致。

> **残留 0.89% 的来源**（不是错）：FVM 的有效求解域长 `L+dx`（细胞中心在 0..L，通量面在半格处），
> 以及界面通量面 `gy = G*dx` 未含两侧半单元材料热阻。两者都随网格细化消失。
>
> **仍待确认**：`device_thermal_summary.json` 的 `R_th_eff_m2K_per_W` / `G_eff_W_m2K` 报的是
> **解析 1D 串联网络**阻值（2.859e-7），并非由 dT/Q 得到的有效值（修正后应为 9.82e-8）。
> 字段语义如此，改名前需先确认没有下游消费者。

---

## 11. 局限

1. **本技能只复核**：源项修复是在 `device-thermal` 源码里做的（2026-09-23 经用户批准）。
   修复后两套实现的 dT 仍差约 0.89%（FVM 有效域长 `L+dx` + 界面半单元处理），不是零；
   该残差随网格细化收敛。
2. 界面用**等效薄层**表示，与零厚度界面热处理有 O(δ/t_ox) 差异；δ 默认为 t_ox/1000，
   已实测对结果不敏感（δ 变 30 倍时 dT 变化 <0.2%）。
3. 目前只做**稳态**、两层、2D 截面；瞬态、3D、电-热耦合需要扩展模型（后端契约已预留）。
4. `skfem` 后端用 scipy 稀疏直接解，网格约 2.5 万单元时耗秒级；再大需要换迭代解。
5. COMSOL 后端**未实现**，只留契约与启用文档。

---

## 12. 自检

```bash
python3 bin/autozt skills                                  # device-thermal-fem 应出现
python3 bin/autozt schema device-thermal-fem --strict       # 0 问题
```
