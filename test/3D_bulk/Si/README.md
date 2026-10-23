# Si（金刚石结构，3D bulk）—— ke / kl / zt 关键结果

- 材料：Si，金刚石结构（Fd-3m），3D bulk
- 技能链：ke-dft-cpu（电子输运）→ kl-dft-cpu（晶格热导）→ zt-dft-cpu（ZT 汇总）
- 导出日期：2026-10-23（关键文件；中间/大文件已省略）
- 溯源：
  - 集群 jzzn：`/public/home/wangchao/Fullerene_Network/work/Si/{ke-dft-cpu,kl-dft-cpu,zt-dft-cpu}`
  - 本地完整数据：`/mnt/d/tf_data/work_AutoZT/Si/`
  - Windows 归档：`D:\joint_research_project\AutoZT\基准测试\3D_bulk\Si\`

## 关键数值

| 项 | 值 | 出处 |
|---|---|---|
| κ_L(300 K) | 109.4816 W/m·K | kl S6_kappa（原始元胞口径） |
| κ_L 胞高 Lz | 3.8652 Å | kappa_summary.json（[0063] Lz_ang） |
| ke 胞高 c | 3.865 Å | 元胞口径闸门通过（差 0.0%） |
| 电子段散射机制 | ADP+IMP+ODP+IVS | ke S8_kappa（AMSET） |
| 带隙 | 1.1209 eV | ke S8_kappa |
| ZT 峰值 n 型 | 0.0934 @ 800 K, n=−1e20 cm⁻³ | κ_L=38.03 W/m·K |
| ZT 峰值 p 型 | 0.0564 @ 800 K, p=+1e20 cm⁻³ | κ_L=38.03 W/m·K |
| ZT(300 K, n=1e20) | 0.009 | zt_summary.txt |

Z = S²σT/(κ_e+κ_L)；κ_L 按温度线性插值到电子段网格（不外推）。

## 文件

```
Si/
├── POSCAR
├── ke/   transport.json            AMSET 3D 输运（S8_kappa，含各散射机制）
│         dpt_result.json           形变势拟合（S8.2_dpt）
│         comparison_300K.csv       300 K 多温度/掺杂对照表（S8.3）
│         comparison_vs_T.csv       温度依赖对照（S8.3）
├── kl/   kappa_summary.json        κ_L(T)、Lz_ang、截断/网格信息（S6_kappa）
└── zt/   zt_summary.txt/json       ZT 汇总（峰值 + 逐温度/掺杂表）
          zt_vs_T.png / zt_vs_doping.png / zt_heatmap.png
```

注意：旧版（2026-09-18）ZT 峰值 n 0.3163 / p 0.2386 已被覆盖（旧版仅 ADP+IMP、κ_L 不同），勿引用。

## 2026-10-23 更新：[0061-d] DPT 简并带边 E1 修正（Si 空穴）

- **代码**：`skill/ke-dft-cpu/step8.2_dpt/gen_step12_dpt.py` —— 用 S8 的 `deformation_matrix.h5`（S7.1 DFT 应变算符）
  在带边简并子空间取 `Y_αβ`，复用 `adp_degen.directions/acoustic_modes`（3D 三支声学模，弹性/均分），
  对带边 ±6 k_BT 内的 DFT 态做 Boltzmann 热平均：`r=⟨Σ_ab|M_ab|²⟩/⟨Σ_a|M_aa|²⟩`，
  折算 `E1_eff = E1_naive·√r`（DPT μ∝1/E1²）。d=1 时原路径逐位不变；手填 E1 优先（manual）。
- **结果（300 K）**：
  - 空穴：E1 **3.4700 → 3.7872 eV**（d=3，k=Γ，bands [1,2,3]，ratio=1.1912，热平均 42 个 DFT 态）；
    DPT μ **4878.586 → 4095.589 cm²/V·s（−16.1%）**
  - 电子：E1=8.72 eV、μ=2433.354 cm²/V·s（不变）
  - BT2 τ（来自 DPT）：τ_hole **393.4 fs**（原约 469 fs），τ_electron 128.7 fs（不变）；
    BT2 文献口径 ZT_xx 峰值（600 K，空穴支）0.973
  - S8.3 对照表 dptμ：p 侧 **4096**（n 侧 2433）
- **与 AMSET ADP-only 对比**（p 1e18/1e19/1e20 = 1441/1378/1019）：DPT 仍高 2.8–4.0×。
  原因：AMSET 的修正是费米海内逐对矩阵元替换（含动量弛豫因子），DPT 是带边单带闭式公式，只能吸收带边热平均的 1.19×。
  两者方法口径不同，不能互相校准。
- **溯源**：`dpt_result.json` 的 `E1_provenance` 前缀 `adp_degen(...)`，并有结构化 `E1_degen`；
  `_SKILL_REV = 2026-10-23-0061d-dptdegen`。
