# 技能参数表（自动生成）

> 由 `python3 scripts/export_params.py --out docs/PARAMETERS.md` 生成；来源：`skill/*/skill.yaml`、`skill/*/templates/**/step.conf`、gen 脚本的 SPEC。
> 生成时的仓库版本：`7e2f7ff+dirty`。请勿手工编辑本文件，改 step.conf 后重跑即可。

## 1. 技能与步骤

| 技能 | 说明 | 步骤数 | 步骤 |
|---|---|---|---|
| band-dft-cpu | 能带(DFT/CPU) | 3 | S1_opt, S2_static, S3_WAVECAR |
| cohp-cogito | COHP/ICOHP/ICOBI 成键分析 (COGITO) | 2 | S1_COHP, S2_summary |
| defect-dft-cpu | 缺陷形成能(DFT/CPU) | 5 | S0_refs, S1_bulk, S2_def, S3_chg, S4_anlys |
| elastic-dft-cpu | 弹性常数(DFT/CPU) | 3 | S1_opt, S2_elastic, S3_post |
| eph-qe-cpu | QE–Wannier90–Perturbo 电子–声子耦合（CPU） | 7 | S0_env, S1_scf, S2_wannier, S3_ph, S4_qe2pert, S5_ephmat, S6_summary |
| fc-fit | 力常数拟合 fc2/fc3 (pheasy/hiphive/phono3py) | 1 | S1_fit |
| ke-dft-cpu | 电子热导率(DFT/CPU) | 9 | S1_opt, S3_uniform, S4_wave, S5_dielect, S5.1_dievalid, S6_elastic, S7_deform, S7.1_read, S8_kappa |
| kl-dft-cpu | 晶格热导率(DFT/CPU) | 5 | S1_opt, S2_static, S4_disp, S5_fc, S6_kappa |
| kl-mace-cpu | 晶格热导率(MACE/CPU) | 4 | S1_relax, S2_force, S3_fc, S4_kappa |
| kl-mace-gpu | 晶格热导率(MACE/GPU) | 4 | S1_relax, S2_force, S3_fc, S4_kappa |
| mlff-mace | 随机位移法MLFF训练 | 9 | S1_relax, S2_cell, S3_calib, S4_gen, S5_label, S6_data, S7_ft, S8_bench, S9_pub |
| opt-dft-cpu | 结构优化(DFT/CPU) | 3 | S1_opt, S2_static, S3_energy |
| opt-mace-cpu | 结构优化(MACE/CPU) | 3 | S1_relax, S2_static, S3_energy |
| opt-mace-gpu | 结构优化(MACE/GPU) | 3 | S1_relax, S2_static, S3_energy |
| phonon-dft-cpu | 声子谱(DFT/CPU) | 3 | S1_opt, S2_disp, S3_phonon |
| phonon-mace-cpu | 声子谱(MACE/CPU) | 3 | S1_relax, S2_force, S3_phonon |
| phonon-mace-gpu | 声子谱(MACE/GPU) | 3 | S1_relax, S2_force, S3_phonon |
| te-screen | 热电快速筛选(替代模型) | 2 | S1_feat, S2_pred |
| unihamgnn | 能带(Uni-HamGNN) | 3 | S1_graph, S2_predict, S3_band |
| zt-dft-cpu | 热电优值 ZT 全流程(电子输运 + 晶格热导 + ZT 汇总) | 15 | S1_opt, S3_uniform, S4_wave, S5_dielect, S5.1_dievalid, S6_elastic, S7_deform, S7.1_read, S8_kappa_e, SK1_opt, SK2_static, SK4_disp, SK5_fc, SK6_kappa, S20_zt |

## 2. io_schema 声明的参数（schema 2，机器可读）

| 技能 | 参数 | 声明默认值 | 取值域 | 位置 | 说明 |
|---|---|---|---|---|---|
| band-dft-cpu | BANDGAP | hse | pbe, hse | step.conf | 带隙层级；pbe = 只到 PBE 带隙+画图（跳过整段 HSE），hse = 继续 HSE 段 |
| band-dft-cpu | plot_steps | True | True, False | tf.yaml/项目配置 | 画图步开关（S3.1_plot / S4.1_plot） |
| band-dft-cpu | vacuum_align | True | True, False | tf.yaml/项目配置 | PBEsol 真空对齐步开关（S3.2_vac / S4.2_vac） |
| band-dft-cpu | bandgap_hse | True | True, False | tf.yaml/项目配置 | HSE 段开关（S4_HSE / S4.1_plot / S4.2_vac）；也可由 BANDGAP=pbe 关掉 |
| cohp-cogito | NUM_OUTER | 4 | - | step.conf | COGITO 轨道优化外循环次数 |
| cohp-cogito | DENSIFY | - | - | step.conf | COHP 曲线加密级别；留空 = 用 COGITO 自身默认（1），2/3 更平滑但更慢 |
| eph-qe-cpu | QE_BIN | $PATH | - | job environment | 包含 pw.x、ph.x、pw2wannier90.x 的目录或 PATH。 |
| eph-qe-cpu | PERTURBO_BIN | $PATH | - | job environment | 包含 qe2pert.x、perturbo.x 的目录或 PATH。 |
| eph-qe-cpu | WANNIER_BIN | $PATH | - | job environment | 包含 wannier90.x 的目录或 PATH。 |
| eph-qe-cpu | PSEUDO_DIR | - | - | job environment | UPF 赝势目录；Si 测试需要其中的 Si.upf。 |
| fc-fit | FIT_ENGINE | phono3py | phono3py, pheasy, hiphive | step.conf | 力常数拟合引擎 |
| ke-dft-cpu | BANDGAP | hse | pbe, hse | step.conf | PBE 或 HSE 带隙分支 |
| kl-dft-cpu | METHOD | alm | alm, findiff | step.conf | 位移生成方式：alm = 随机位移+压缩感知（帧数少）；findiff = phono3py 对称有限位移 |
| kl-dft-cpu | FIT_ENGINE | auto | auto, phono3py, pheasy | step.conf | 力常数拟合引擎；auto = 按 DIM/帧数选（多为 phono3py），pheasy 需要 METHOD=alm（随机位移压缩感知） |
| kl-dft-cpu | FC_CALC | symfc | symfc, alm | step.conf | phono3py 的拟合后端：symfc（对称性约束）或 alm |
| kl-dft-cpu | SOLVER | phono3py | phono3py, shengbte | step.conf | BTE 求解器（shengbte 走 FORCE_CONSTANTS_2ND/3RD 双格式导出） |
| kl-dft-cpu | BTE_METHOD | rta | rta, lbte | step.conf | 单模弛豫时间近似（快）或直接解线性化 BTE（准，内存 O(N_mode^2)） |
| kl-dft-cpu | KAPPA_MESH | auto | - | step.conf | BTE q 网格：auto = 2D 按倒空间长度估（Q_LEN_2D/MESH_MIN，真空轴=1）、3D 用 15 15 15；也可写 auto3d 或 "N N N" |
| kl-dft-cpu | IMAG_THR | 0.1 | - | step.conf | 虚频闸（THz）：最小频率低于 -IMAG_THR 判结构不稳定，S6_kappa 不启动 |
| kl-dft-cpu | MAX_DISP | 500 | - | step.conf | 位移帧数硬闸：超过就报错停步（防超胞设太大白烧机时） |
| kl-dft-cpu | OVERSAMPLE | 3 | - | step.conf | 随机位移过采样系数（帧数 = ALM 自由力常数个数 × 它） |
| kl-mace-cpu | DEVICE | cpu | cpu | step.conf | MACE 推理设备 |
| kl-mace-cpu | SOLVER | phono3py | phono3py, shengbte | step.conf | BTE 求解器 |
| kl-mace-gpu | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | step.conf | MACE 势：内置名或本地 .model 路径（换势 = 换这一项）；技能默认用 Materials Project 微调势 |
| kl-mace-gpu | DEVICE | auto | auto, cuda, cpu | step.conf | MACE 推理设备：auto = 按集群/模板解析（GPU 机器上即 cuda；S2_disp 步骤默认写死 cuda）；CPU 版见 kl-mace-cpu |
| kl-mace-gpu | METHOD | random | random, findiff | step.conf | random = hiphive MC-rattle 随机位移（默认）；findiff = 对称有限位移 |
| kl-mace-gpu | N_RANDOM | auto | - | step.conf | 随机位移帧数：auto 按 ALM 自由力常数个数反推，也可写死 |
| kl-mace-gpu | OVERSAMPLE | 3 | - | step.conf | 随机位移过采样系数 |
| kl-mace-gpu | DISP_DISTANCE | 0.03 | - | step.conf | MC-rattle 目标位移 RMS（Å）；虚频可疑时用 0.01 / 0.05 各跑一次看敏感性 |
| kl-mace-gpu | RANDOM_SEED | 2025 | - | step.conf | MC-rattle 随机种子（可复现） |
| kl-mace-gpu | RESIDUAL_TOL | 0.002 | - | step.conf | 残余力闸（eV/Å）：S1 弛豫后残余力超它就判 FAIL——这才是声子能不能算的真闸门 |
| kl-mace-gpu | IMAG_THR | 0.1 | - | step.conf | 虚频闸（THz）：最小频率低于 -IMAG_THR 判不稳定，S4 不启动 |
| kl-mace-gpu | SOLVER | phono3py | phono3py, shengbte | step.conf | BTE 求解器 |
| kl-mace-gpu | BTE | rta | rta, lbte | step.conf | rta = --br（默认，快）；lbte = --lbte 直接解（准，内存 O(N_mode^2)） |
| kl-mace-gpu | KAPPA_MESH | 24 24 24 | - | step.conf | BTE q 网格（step2 写进 klmace_params 供 step4；2D 自动压真空方向） |
| mlff-mace | GENERATION | 0 | - | step.conf | 迭代代数 |
| mlff-mace | FORCE_CONTINUE | False | True, False | step.conf | 停机闸后的显式继续开关 |
| opt-mace-cpu | DEVICE | cpu | cpu | step.conf | MACE 推理设备 |
| opt-mace-gpu | DEVICE | cuda | cuda | step.conf | MACE GPU 推理设备 |
| phonon-dft-cpu | SUPERCELL | - | - | step.conf | phonopy 位移超胞设置 |
| phonon-mace-cpu | DEVICE | cpu | cpu | step.conf | MACE 推理设备 |
| phonon-mace-gpu | DEVICE | cuda | cuda | step.conf | MACE GPU 推理设备 |
| unihamgnn | SOC | True | True, False | step.conf | 是否生成 SOC 图数据（技能模板默认 true = 同时出 non-SOC 与 SOC 两份） |
| zt-dft-cpu | BANDGAP | hse | pbe, hse | step.conf | 电子段带隙层级：pbe 只算 PBE，hse 继续算 HSE（透传给上游 ke 段） |
| zt-dft-cpu | KTEMP_MODE | interp | interp, const300 | step.conf | ZT 汇总取 κ_L 的方式：interp = 按温度插值 κ_L(T)；const300 = 固定用 300 K 的 κ_L |

## 3. 出厂默认值总表（step.conf 与 SPEC 内建）

| 技能 | 来源 | 参数 | 默认值 | 类型 | 备注 |
|---|---|---|---|---|---|
| band-dft-cpu | templates/step.conf（技能级） | BANDGAP | hse | - | - |
| band-dft-cpu | templates/step.conf（技能级） | FUNC | pbesol | - | - |
| cohp-cogito | templates/step.conf（技能级） | CONDA_ENV | atomate2_p_a | - | - |
| cohp-cogito | templates/step.conf（技能级） | CONDA_SH | ~/miniconda3/etc/profile.d/conda.sh | - | - |
| cohp-cogito | templates/step.conf（技能级） | DENSIFY | (空) | - | - |
| cohp-cogito | templates/step.conf（技能级） | NUM_OUTER | 4 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | BAND_SUMMARY | (空) | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | ENCUT | 370 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | EPSILON | 100.0 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | FUNC | pbe-d3 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | KMESH | 2 2 1 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | LVHAR_CHARGED | 1 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | MSTAR_E | 0.2 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | MSTAR_H | 0.2 | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | POTCAR_VARIANT | Pb:Pb_d, Sn:Sn_d, Sb:Sb, Bi:Bi_d, Te:Te | - | - |
| defect-dft-cpu | templates/step.conf（技能级） | SUPERCELL | 3 3 1 | - | - |
| elastic-dft-cpu | templates/step.conf（技能级） | FUNC | pbesol | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | CONDA_ENV | atomate2_p_a | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | CONDA_SH | /public/home/wangchao/miniconda3/etc/profile.d/conda.sh | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | COORDS | cartesian | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | DIM | auto | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | ENABLE_FC | 3 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | EQUILIBRIUM_FORCES_NPY | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | EXPORT_SHENGBTE | true | - | also write shengbte/FORCE_CONSTANTS_2ND / _3RD |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FC3_CUTOFF | (空) | - | fc3 cutoff in A; empty = no cutoff (slow, big) |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FC3_LOAD_GB_LIMIT | 8.0 | - | skip materialising fc3 above this (ShengBTE/RMSE) |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FC_CALC | symfc | - | symfc / alm |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FIT_ENGINE | phono3py | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FIT_INPUT_DIR | auto | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | FIT_RMSE_FRAMES | 0 | - | 0 = off; else N training frames for the residual |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_ALPHA | 1e-10 | - | regularisation for ridge / lasso |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_CUTOFF2 | 6.0 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_CUTOFF3 | 6.0 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_ENFORCE_ASR | true | - | project Huang / Born-Huang sum rules |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_FIT_METHOD | ridge | - | ols / ridge / lasso / ard / bayes |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_N_CONFIGS | 0 | - | 0 = use every frame, else a subsample |
| fc-fit | templates/step1_fit/step.conf（步骤级） | HIPHIVE_SYMPREC | 1e-5 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | IMAG_THR | 0.10 | - | min frequency below -IMAG_THR THz counts as unstable |
| fc-fit | templates/step1_fit/step.conf（步骤级） | NULL_SPACE_EPS | 0.001 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_BIN | pheasy | - | pheasy / pheasy-gpu (GPU build, needs a GPU template) |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_C2_CUTOFF | (空) | - | fc2 cutoff in A; empty = all interactions |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_C3_CUTOFF | 6 | - | fc3 cutoff in A; empty = all interactions |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_CV | 5 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_FIT_METHOD | RFE | - | OLS/LASSO/ALASSO/RFE/RFE-OLS-TSQR/RIDGE |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_GPU_LASSO_RESIDENT | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_LASSO_TWOLEVEL | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_MAX_ITER | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_MU_MAX | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_MU_MIN | -8 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_NGPU | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_NMU | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_OLS_MAXITER | (空) | - | OLS LSMR iteration cap; empty = pheasy default 5000 |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_OLS_RIDGE | (空) | - | OLS ridge override; empty = pheasy default |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_RASR | BHH | - | BH / H / BHH / none (rotational sum rules) |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_SEED | 666666 | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_STD | false | - | --std: standardise the training data.  pheasy_fit.sh applies |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_TOL | (空) | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | PHEASY_TUNING | safe | - | safe / kl |
| fc-fit | templates/step1_fit/step.conf（步骤级） | SUBTRACT_EQUILIBRIUM | true | - | - |
| fc-fit | templates/step1_fit/step.conf（步骤级） | SUPERCELL | (空) | - | - |
| fc-fit | gen_step1_fit.py | BAND_POINTS | 51 | int | SPEC 内建默认 |
| fc-fit | gen_step1_fit.py | PHEASY_CV_MAX_ITER | - | str | SPEC 内建默认 |
| fc-fit | gen_step1_fit.py | PHEASY_CV_TOL | - | str | SPEC 内建默认 |
| ke-dft-cpu | step.conf（技能级） | BANDGAP | hse | - | - |
| ke-dft-cpu | step.conf（技能级） | FUNC | pbesol | - | - |
| ke-dft-cpu | step8_amset/gen_step10_amset.py | LAYER_THICKNESS | LAYER_THICKNESS | str | SPEC 内建默认 |
| ke-dft-cpu | step8_amset/gen_step10_amset.py | NWORKERS | NWORKERS | int | SPEC 内建默认 |
| ke-dft-cpu | step3b_uniform_full/gen_step5b_uniform_full.py | VACUUM_KZ_MIN | VACUUM_KZ_MIN | int | SPEC 内建默认 |
| kl-dft-cpu | templates/step.conf（技能级） | FUNC | auto | - | - |
| kl-dft-cpu | templates/step2_static/step.conf（步骤级） | FUNC | pbesol | - | 与 step1 一致；auto=继承 step1 的 workflow_method.txt |
| kl-dft-cpu | templates/step2_static/step.conf（步骤级） | KSPACING | 0.03 | - | 静态密网格 |
| kl-dft-cpu | templates/step3_nac/step.conf（步骤级） | FUNC | pbesol | - | - |
| kl-dft-cpu | templates/step3_nac/step.conf（步骤级） | KSPACING | 0.03 | - | - |
| kl-dft-cpu | templates/step3_nac/step.conf（步骤级） | SKIP_IF_METAL | true | - | 仅提示；实际跳过用 nac: false |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | ALM_CUT2 | (空) | - | 二阶截断(Å)，空=不截断 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | ALM_CUT3 | 6.0 | - | 三阶截断(Å)。帧数超闸时先调小这个 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | DISP_RMS | 0.03 | - | 目标位移模长 RMS(Å)（每分量 RMS≈0.0173） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FC3_CUTOFF_PAIR | (空) | - | 三阶对距离上限(Å)；空=全对称集(最准最贵)。大超胞务必设 4~6 控制帧数 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FD_DISTANCE | 0.03 | - | 有限位移幅度(Å) |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FORCE_EDIFF | 1E-8 | - | 原为 1E-7（与模板注释不符，已对齐） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FORCE_LREAL | .FALSE. | - | 倒空间投影，力无格点噪声；原为 Auto。超胞几百原子时可换 Auto， |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FORCE_PREC | Accurate | - | phonopy 官方超胞取力示例口径 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | FUNC | pbesol | - | - |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | KAPPA_MESH | auto | - | - |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | KSPACING | 0.03 | - | - |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MAX_DISP | 500 | - | ★ 位移帧数硬闸：超过就报错停步（正常应 ≤ 几百帧） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MAX_MULTIPLE | 6 | - | - |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MC_DMIN_SCALE | 0.75 | - | MC 最小间距 d_min = 平衡最近邻 × 此系数 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MC_NITER | 10 | - | MC 迭代数（实际位移 ∝ √n_iter，脚本内自动标定 rattle_std） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MC_SEED | 2025 | - | - |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MESH_MIN | 20 | - | auto 时每方向下限（BZ 采样下限） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | METHOD | alm | - | alm / findiff |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | MIN_SC_LEN | 12.0 | - | 超胞每个非真空方向最小胞长(Å) |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | OVERSAMPLE | 3 | - | 过采样系数（帧数公式里的 ×3） |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | Q_LEN_2D | 275.0 | - | auto 时的目标倒空间长度(Å)，文献 2D 起步值 250~300 |
| kl-dft-cpu | templates/step4_disp/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "3 3 3"；空=自动。2D 会把真空方向压 1 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | BAND_POINTS | 51 | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | EXPORT_SHENGBTE | true | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | FC3_CUTOFF | 6 | - | fc3 截断 A；留空或 None = 不截断 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | FC_CALC | symfc | - | symfc / alm |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | FIT_ENGINE | auto | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | IMAG_THR | 0.10 | - | 最小声子频率 < -IMAG_THR(THz) 判不稳定，S6 不启动 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | MIN_SUCCESS_FRAMES | 0 | - | 成功帧绝对下限，0 = 不限 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | MIN_SUCCESS_RATIO | 0.9 | - | 成功帧占比下限，低于它报错（0 = 不限） |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | NULL_SPACE_EPS | 0.001 | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | PHEASY_C3_CUTOFF | 6 | - | pheasy fc3 截断 A；None = 不截断 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | PHEASY_ENABLE_FC | 3 | - | 拟合最高阶：热导率需 >=3 |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | PHEASY_FIT_METHOD | auto | - | auto / LASSO / RFE(内存友好) / OLS(最吃内存) |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | PHEASY_RASR | auto | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | ZA_CHECK | auto | - | - |
| kl-dft-cpu | templates/step5_fc/step.conf（步骤级） | ZA_QMAX | 0.05 | - | 拟合 q 上限（倒格子约化单位）；太小主要反映插值误差 |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | BTE_METHOD | rta | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | COMPARE_LBTE | auto | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | FOURPHONON_CPUS_PER_GPU | 8 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | FOURPHONON_NGPU | 4 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | ISOTOPE | True | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | KAPPA_2D_THICKNESS | vdw | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | KAPPA_NAC | auto | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | MESH_MIN | 20 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | MESH_SCAN | auto | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | SCALEBROAD | 1.0 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | SHENGBTE_CORES_PER_NUMA | 8 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | SHENGBTE_TOTAL_CORES | 96 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | SOLVER | phono3py | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | T_MAX | 800 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | T_MIN | 100 | - | - |
| kl-dft-cpu | templates/step6_kappa/step.conf（步骤级） | T_STEP | 100 | - | - |
| kl-dft-cpu | gen_step5.1_plot_phonon.py | CONDA_ENV | atomate2_p_a | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step5.1_plot_phonon.py | CONDA_SH | - | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step2_static.py | ENCUT | None | int | SPEC 内建默认 |
| kl-dft-cpu | gen_step2_static.py | ENCUT_FACTOR | 1.5 | float | SPEC 内建默认 |
| kl-dft-cpu | gen_step6_kappa.py | FOURPHONON_EXE | - | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step6_kappa.py | KAPPA_CONVERGENCE | True | bool | SPEC 内建默认 |
| kl-dft-cpu | gen_step2_static.py | KSCHEME | 2 | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step6_kappa.py | MESH_OVERRIDE | None | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step5_fc.py | PHEASY_BIN | pheasy | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step5.1_plot_phonon.py | PLOT_CORES | 16 | int | SPEC 内建默认 |
| kl-dft-cpu | gen_step5.1_plot_phonon.py | SBATCH_QOS | regular | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step6_kappa.py | SHENGBTE_EXE | ShengBTE | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step6_kappa.py | SHENGBTE_NTASKS | auto | str | SPEC 内建默认 |
| kl-dft-cpu | gen_step4_disp.py | STRESS_2D_THR | 0.5 | float | SPEC 内建默认 |
| kl-dft-cpu | gen_step2_static.py | VASPKIT_EXE | vaspkit | str | SPEC 内建默认 |
| kl-mace-cpu | templates/step.conf（技能级） | DEVICE | cpu | - | - |
| kl-mace-cpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| kl-mace-cpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | none | - | primitive = 先 spglib 取标准原胞；none = 原样用（2D 恒为 none） |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | auto / 2d / 3d。2D 会锁死真空方向的胞长和相关剪切 |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | 挂 ASE FixSymmetry。关掉的话数值噪声会把空间群降到 P1， |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-4 | - | eV/A。机器学习势很便宜，别学 DFT 只收到 1e-2 |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | ase.optimize 里的名字：FIRE / LBFGS / BFGS |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | false = 不优化只测残余力（你确信结构已在本势极小点时才用） |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | 晶格也放开。做定容/加压对比时置 false |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | eV/A。残余力超过它就判本步 FAIL —— 这才是声子能不能算的真闸门 |
| kl-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | ALM_CUT2 | (空) | - | 二阶截断(Å)；空=不截断 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | ALM_CUT3 | 6.0 | - | 三阶截断(Å) |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | CKPT | 10 | - | 每多少帧落一次断点。CPU 单帧慢，落密一点，被墙钟砍了少赔 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | DEVICE | cpu | - | - |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | DISP_DISTANCE | 0.03 | - | MC-rattle 目标位移 RMS(Å)；三阶用 0.03 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | FC2_SUPERCELL | (空) | - | 二阶专用大超胞（--dim-fc2）。CPU 上这也是性价比最高的 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | KAPPA_MESH | 15 15 15 | - | 起步网格；2D 自动将真空轴压成 1 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MAX_MULTIPLE | 6 | - | - |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MC_DMIN_SCALE | 0.85 | - | MC-rattle d_min = 最近邻 × 此系数 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MC_NITER | 10 | - | MC-rattle 迭代数 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | METHOD | random | - | random = hiphive MC-rattle 随机位移（帧数按 ALM nfree 反推） |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MIN_SC_LEN | 12.0 | - | 超胞每个非真空方向的最小胞长(A)。与 kl-dft-cpu/GPU 版一致 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | N_RANDOM | auto | - | random 帧数：auto=按 ALM nfree 反推（N=max(10, ceil(Σnfree/DOF)*OVERSAMPLE)） |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | N_RANDOM_FC2 | auto | - | fc2 专用超胞随机帧数（auto=按 ALM 反推，只算二阶） |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | OVERSAMPLE | 3 | - | 随机位移过采样系数（越大越准越慢） |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | RANDOM_SEED | 2025 | - | MC-rattle 随机种子 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | SUBTRACT_RESIDUAL | true | - | 扣掉未位移超胞的残余力。不扣会破坏声学求和规则 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "3 3 3"；空 = 按 MIN_SC_LEN 自动。2D 真空方向恒 1 |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | nodes | 1 | - | - |
| kl-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | ntasks_per_node | 1 | - | - |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | BAND_POINTS | 101 | - | - |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | FIT | auto | - | phono3py 拟合器：auto（findiff->sym-fc，random->symfc）/ symfc / alm |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | FIT_SOFTWARE | pheasy | - | - |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | IMAG_THR | 0.10 | - | 虚频阈值(THz)：最小频率 < -IMAG_THR 判不稳定，step4 不启动 |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | NAC_BORN | (空) | - | 外部 BORN 文件路径。MACE 给不出 Born 电荷和 eps_inf，极性 |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | PHEASY_C3_CUTOFF | 6.0 | - | pheasy 三阶截断(Å)；None/空=不截断 |
| kl-mace-cpu | templates/step3_fc/step.conf（步骤级） | PHEASY_METHOD | LASSO | - | pheasy 方法：OLS / LASSO / RFE / RFE_TSQR |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | BTE | rta | - | rta = --br（默认）；lbte = --lbte 直接解，内存 O(N_mode^2) |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | EXTRA_ARGS | (空) | - | 原样附加给 phono3py-load，如 "--boundary-mfp 1e6" |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | ISOTOPE | true | - | 同位素散射（自然丰度） |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | KAPPA_2D_THICKNESS | vdw | - | 2D 层厚：vdw=原子z跨度+两侧vdW半径（照文献口径，勿写死；超晶格在项目 step.conf 锁 6.73） |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | MESH_OVERRIDE | (空) | - | 空 = 用 step2 写进 klmace_params 的 MESH |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | MESH_SCAN | (空) | - | 默认单网格；收敛测试时显式填写多套网格 |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | SCALEBROAD | 1.0 | - | shengbte 高斯展宽系数（=ShengBTE 默认；0.1 更快但会漏过程） |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_CORES_PER_NUMA | 8 | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_EXE | ShengBTE | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_TOTAL_CORES | 96 | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | SOLVER | phono3py | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | T_MAX | 800 | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | T_MIN | 100 | - | - |
| kl-mace-cpu | templates/step4_kappa/step.conf（步骤级） | T_STEP | 100 | - | - |
| kl-mace-gpu | templates/step.conf（技能级） | DEVICE | auto | - | - |
| kl-mace-gpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| kl-mace-gpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | primitive | - | primitive = 先 spglib 取标准原胞；none = 原样用（2D 恒为 none） |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | auto / 2d / 3d。2D 会锁死真空方向的胞长和相关剪切 |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | 挂 ASE FixSymmetry。关掉的话数值噪声会把空间群降到 P1， |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-4 | - | eV/A。机器学习势很便宜，别学 DFT 只收到 1e-2 |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | ase.optimize 里的名字：FIRE / LBFGS / BFGS |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | false = 不优化只测残余力（你确信结构已在本势极小点时才用） |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | 晶格也放开。做定容/加压对比时置 false |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | eV/A。残余力超过它就判本步 FAIL —— 这才是声子能不能算的真闸门 |
| kl-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | ALM_CUT2 | (空) | - | 二阶截断(Å)；空=不截断 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | ALM_CUT3 | 6.0 | - | 三阶截断(Å) |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | CKPT | 50 | - | 每多少帧落一次断点（GPU 快，落盘可以稀一点） |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | DEVICE | cuda | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | DISP_DISTANCE | 0.03 | - | MC-rattle 目标位移 RMS(Å)。虚频可疑时用 0.01 / 0.05 各跑一次看是否敏感 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | FC2_SUPERCELL | (空) | - | 二阶专用大超胞（--dim-fc2），如 "6 6 6"。二阶的长程尾巴 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | KAPPA_MESH | 24 24 24 | - | BTE q 网格，写进 klmace_params 供 step4（2D 自动 kz=1） |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MAX_MULTIPLE | 8 | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MC_DMIN_SCALE | 0.85 | - | MC-rattle d_min = 最近邻 × 此系数 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MC_NITER | 10 | - | MC-rattle 迭代数 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | METHOD | random | - | random = hiphive MC-rattle 随机位移（帧数按 ALM nfree 反推，默认） |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MIN_SC_LEN | 12.0 | - | 超胞每个非真空方向的最小胞长(A)。与 kl-dft-cpu/CPU 版一致 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | N_RANDOM | auto | - | random 帧数：auto=按 ALM nfree 反推（N=max(10, ceil(Σnfree/DOF)*OVERSAMPLE)） |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | N_RANDOM_FC2 | auto | - | fc2 专用超胞随机帧数（auto=按 ALM 反推，只算二阶） |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | OVERSAMPLE | 3 | - | 随机位移过采样系数 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | RANDOM_SEED | 2025 | - | MC-rattle 随机种子 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | SUBTRACT_RESIDUAL | true | - | 扣掉未位移超胞的残余力。不扣会破坏声学求和规则 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "4 4 4"；空 = 按 MIN_SC_LEN 自动。2D 真空方向恒 1 |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | gres | gpu:1 | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | nodes | 1 | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | ntasks_per_node | 1 | - | - |
| kl-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | partition | gpu | - | - |
| kl-mace-gpu | templates/step3_fc/step.conf（步骤级） | BAND_POINTS | 101 | - | - |
| kl-mace-gpu | templates/step3_fc/step.conf（步骤级） | FIT | auto | - | auto（findiff -> --sym-fc，random -> --fc-calc symfc） |
| kl-mace-gpu | templates/step3_fc/step.conf（步骤级） | IMAG_THR | 0.10 | - | 虚频阈值(THz)：最小频率 < -IMAG_THR 判不稳定，step4 不启动 |
| kl-mace-gpu | templates/step3_fc/step.conf（步骤级） | NAC_BORN | (空) | - | 外部 BORN 文件路径。MACE 给不出 Born 电荷和 eps_inf，极性 |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | BTE | rta | - | rta = --br（默认）；lbte = --lbte 直接解，内存 O(N_mode^2) |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | EXTRA_ARGS | (空) | - | 原样附加给 phono3py-load，如 "--boundary-mfp 1e6" |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | ISOTOPE | true | - | 同位素散射（自然丰度） |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | MESH_OVERRIDE | (空) | - | 空 = 用 step2 写进 klmace_params 的 MESH |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | MESH_SCAN | (空) | - | 分号分隔的多套网格做收敛测试，如 "16 16 16; 20 20 20; 24 24 24" |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | SCALEBROAD | 1.0 | - | shengbte 高斯展宽系数（=ShengBTE 默认；0.1 更快但会漏过程） |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_CORES_PER_NUMA | 8 | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_EXE | ShengBTE | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | SHENGBTE_TOTAL_CORES | 96 | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | SOLVER | phono3py | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | T_MAX | 800 | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | T_MIN | 100 | - | - |
| kl-mace-gpu | templates/step4_kappa/step.conf（步骤级） | T_STEP | 100 | - | - |
| mlff-mace | templates/step.conf（技能级） | AUTO_U | auto | - | auto / true / false |
| mlff-mace | templates/step.conf（技能级） | BATCH_SIZE | 10 | - | autoplex MACE 默认；GPU OOM 自动降 1 |
| mlff-mace | templates/step.conf（技能级） | CELL_POLICY | primitive | - | step1 取原胞（超胞/声子都在原胞上展开） |
| mlff-mace | templates/step.conf（技能级） | CURVE_POINTS | 25,50,100,200,all | - | - |
| mlff-mace | templates/step.conf（技能级） | CURVE_TOL | 0.05 | - | 学习曲线判平阈值（相对改善） |
| mlff-mace | templates/step.conf（技能级） | DATA_MODE | scratch | - | scratch / extend |
| mlff-mace | templates/step.conf（技能级） | DEVICE | auto | - | auto / cpu / cuda（jzzn 无 GPU 分区，auto 恒落 cpu； |
| mlff-mace | templates/step.conf（技能级） | DIMENSION | auto | - | auto / 2d / 3d |
| mlff-mace | templates/step.conf（技能级） | DISP_WARN | 80 | - | - |
| mlff-mace | templates/step.conf（技能级） | DTYPE | float64 | - | ★ 不许改（float32 力误差会造出假虚频） |
| mlff-mace | templates/step.conf（技能级） | E0S_MODE | estimated | - | estimated / json（基座与 DFT 零点差大，用 estimated 对齐） |
| mlff-mace | templates/step.conf（技能级） | EDIFF | 1e-7 | - | 单点必须比常规静态更严 |
| mlff-mace | templates/step.conf（技能级） | ENCUT_OVERRIDE | (空) | - | 空 = ceil(1.5 × max ENMAX) |
| mlff-mace | templates/step.conf（技能级） | ENERGY_LIMIT | 0.005 | - | eV/atom，离群过滤（>10% WARN，>30% FAIL） |
| mlff-mace | templates/step.conf（技能级） | ENERGY_WEIGHT | auto | - | 3D→1.0，2D→10.0 |
| mlff-mace | templates/step.conf（技能级） | EPOCHS | 1500 | - | 实际轮数上限：mace 0.3.16 单卡 PATIENCE 只打日志不跳出 |
| mlff-mace | templates/step.conf（技能级） | FORCES_WEIGHT | 100.0 | - | - |
| mlff-mace | templates/step.conf（技能级） | FORCE_CONTINUE | false | - | 停机后强制继续（默认禁止） |
| mlff-mace | templates/step.conf（技能级） | FORCE_LIMIT | 40.0 | - | eV/Å，力分量上限（=autoplex force_max 默认，源码为准） |
| mlff-mace | templates/step.conf（技能级） | FUNC | pbesol | - | pbe / pbesol / pbe-d3（step1 与 step5 必须一致） |
| mlff-mace | templates/step.conf（技能级） | GENERATION | 0 | - | - |
| mlff-mace | templates/step.conf（技能级） | GEN_INCREMENT | 20 | - | 每代新增 rattle 帧数 |
| mlff-mace | templates/step.conf（技能级） | GRUNEISEN_STRAIN | 0.01 | - | Grüneisen 用 ±1% 晶格应变 |
| mlff-mace | templates/step.conf（技能级） | HUBER_DELTA | 0.05 | - | eV/Å。MACE 默认 0.01 对本数据太小，全程线性段=L1，会推高力 RMSE |
| mlff-mace | templates/step.conf（技能级） | IMPROVE_MIN | 0.02 | - | THz，判定「本代是否有实质改善」 |
| mlff-mace | templates/step.conf（技能级） | ISO_BOX | 15.0 | - | 孤立原子盒子边长（Å） |
| mlff-mace | templates/step.conf（技能级） | KPOINTS_GRID | (空) | - | 空 = Γ-only（超胞取力标准做法）；"2 2 2" 显式 |
| mlff-mace | templates/step.conf（技能级） | KSPACING_TOL | 0.20 | - | 指纹 k 点密度容差 |
| mlff-mace | templates/step.conf（技能级） | LOSS | huber | - | autoplex 同款（含应力项；换 stress/universal 需先确认你的 mace 版本） |
| mlff-mace | templates/step.conf（技能级） | LR | 0.001 | - | autoplex MACE 默认（multihead 模式需 --force_mh_ft_lr） |
| mlff-mace | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | 与 kl-mace-* / phonon-mace-cpu 对齐（原 2023-12-03-mace-128-L1_epoch-199.model） |
| mlff-mace | templates/step.conf（技能级） | MAX_ATOMS | 150 | - | - |
| mlff-mace | templates/step.conf（技能级） | MAX_GENERATION | 4 | - | - |
| mlff-mace | templates/step.conf（技能级） | MIN_ATOMS | 60 | - | - |
| mlff-mace | templates/step.conf（技能级） | MIN_DIST_RATIO | 0.75 | - | d_min = 比例 × (两原子共价半径之和) |
| mlff-mace | templates/step.conf（技能级） | MIN_VACUUM | 15.0 | - | 2D 真空厚度下限 |
| mlff-mace | templates/step.conf（技能级） | NCORE | 4 | - | 12 核：KPAR=1、NCORE=4 |
| mlff-mace | templates/step.conf（技能级） | N_COMMITTEE | 4 | - | seed 固定 1,2,3,4 |
| mlff-mace | templates/step.conf（技能级） | N_GPU | 0 | - | GPU 微调分卡数：0=auto（=N_COMMITTEE，seed s → GPU (s-1)%N_GPU）； |
| mlff-mace | templates/step.conf（技能级） | N_PER_CELL | 2 | - | 每个（应变, 幅度）格点的种子数 |
| mlff-mace | templates/step.conf（技能级） | PATIENCE | 100 | - | 验证损失连续无改善多少代即"早停"（单卡仅记录不跳出；模型取最佳 checkpoint） |
| mlff-mace | templates/step.conf（技能级） | PRE_XYZ_FILES | (空) | - | extend 用，逗号分隔（绝对/相对路径） |
| mlff-mace | templates/step.conf（技能级） | RATTLE_STD | auto | - | auto = step3 自校准 [0.5,1.0,1.6]×u_rms(300K) |
| mlff-mace | templates/step.conf（技能级） | RATTLE_STD_FALLBACK | 0.03,0.06,0.10 | - | - |
| mlff-mace | templates/step.conf（技能级） | REF_DISP | 0.1 | - | Å，=autoplex 默认（0.01 的位移力会被 rattle 帧在损失里淹没，谐波刚度学不到） |
| mlff-mace | templates/step.conf（技能级） | REF_FC2_PATH | (空) | - | 有则跳过 DFT 声子基准（displ 帧不生成） |
| mlff-mace | templates/step.conf（技能级） | REPLAY_XYZ | (空) | - | 多头微调 replay，必需（集群无外网） |
| mlff-mace | templates/step.conf（技能级） | RMS_MAX | 0.2 | - | THz，主收敛闸 |
| mlff-mace | templates/step.conf（技能级） | SEED_BASE | 2025 | - | - |
| mlff-mace | templates/step.conf（技能级） | START_SWA | 1200 | - | autoplex SWA 起点（USE_SWA=true 且 EPOCHS>START_SWA 才生效） |
| mlff-mace | templates/step.conf（技能级） | STRESS_WEIGHT | auto | - | 3D→10.0，2D→0（面外应力是垃圾不能训） |
| mlff-mace | templates/step.conf（技能级） | USE_SWA | false | - | [FIX-H3] 默认关：PATIENCE 不跳出时 EPOCHS=START_SWA 前的段才训到， |
| mlff-mace | templates/step.conf（技能级） | U_ANION_GATE | true | - | - |
| mlff-mace | templates/step.conf（技能级） | U_GATE_ANIONS | O F | - | - |
| mlff-mace | templates/step.conf（技能级） | U_OVERRIDE | (空) | - | 如  Mn:3.9 Fe:0 |
| mlff-mace | templates/step.conf（技能级） | VACUUM_AXIS_POLICY | error | - | 2D 真空不在 c 轴时报错（不会静默转走） |
| mlff-mace | templates/step.conf（技能级） | VALID_FRACTION | 0.10 | - | - |
| mlff-mace | templates/step.conf（技能级） | VOL_FACTORS | 0.94,0.96,0.98,1.00,1.02,1.04,1.06 | - | 7 档 ±6%：B0 是 required 闸，3 点(E0/B0/V0) |
| mlff-mace | templates/step1_relax/step.conf（步骤级） | nodes | 1 | - | - |
| mlff-mace | templates/step1_relax/step.conf（步骤级） | ntasks_per_node | 12 | - | - |
| mlff-mace | templates/step5_label/step.conf（步骤级） | nodes | 1 | - | - |
| mlff-mace | templates/step5_label/step.conf（步骤级） | ntasks_per_node | 12 | - | - |
| mlff-mace | mlff_common.py | ALGO | Normal | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | ALLOW_2D_FIXED_CELL | False | bool | SPEC 内建默认 |
| mlff-mace | relax_common.py | CALC_FORMATION | None | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | CALC_INTERCALATION | None | str | SPEC 内建默认 |
| mlff-mace | mlff_common.py | CONDA_ENV | - | str | SPEC 内建默认 |
| mlff-mace | mlff_common.py | CONDA_SH | - | str | SPEC 内建默认 |
| mlff-mace | mlff_common.py | FORCE_MH_FT_LR | True | bool | SPEC 内建默认 |
| mlff-mace | relax_common.py | GUEST_ELEMENT | None | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | HOST_DIR | None | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | HOST_ENERGY | None | float | SPEC 内建默认 |
| mlff-mace | relax_common.py | HOST_FORMULA | None | elemmap | SPEC 内建默认 |
| mlff-mace | mlff_common.py | MACE_MODEL_DIR | - | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_ALLOW_3D_TPL | false | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_DIPOL | auto | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_ENCUT_FLOOR | 0 | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_ISPIN | auto | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_KPOINTS | gamma | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MOL_MOMENT | 1.0 | str | SPEC 内建默认 |
| mlff-mace | relax_common.py | MU | None | elemmap | SPEC 内建默认 |
| mlff-mace | mlff_common.py | MULTIHEAD | False | bool | SPEC 内建默认 |
| mlff-mace | relax_common.py | MU_GUEST | None | float | SPEC 内建默认 |
| mlff-mace | mlff_common.py | NUM_SAMPLES_PT | 30000 | int | SPEC 内建默认 |
| mlff-mace | relax_common.py | STALL_MINUTES | 60 | int | SPEC 内建默认 |
| mlff-mace | relax_common.py | STD_CELL | None | str | SPEC 内建默认 |
| mlff-mace | gen_step2_supercell.py | STEP | STEP | str | SPEC 内建默认 |
| mlff-mace | mlff_common.py | STRESS_WEIGHT_3D | 1.0 | float | SPEC 内建默认 |
| mlff-mace | relax_common.py | U_UNKNOWN_POLICY | None | str | SPEC 内建默认 |
| opt-dft-cpu | templates/step1_opt/step.conf（步骤级） | FUNC | pbesol | - | - |
| opt-dft-cpu | templates/step3_energy/step.conf（步骤级） | CALC_FORMATION | true | - | - |
| opt-dft-cpu | templates/step3_energy/step.conf（步骤级） | CALC_INTERCALATION | true | - | - |
| opt-dft-cpu | gen_step3_energy.py | AUTO_U | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | CELL_POLICY | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | GUEST_ELEMENT | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | HOST_DIR | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | HOST_ENERGY | None | float | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | HOST_FORMULA | None | elemmap | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_ALLOW_3D_TPL | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_DIPOL | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_ENCUT_FLOOR | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_ISPIN | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_KPOINTS | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MOL_MOMENT | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MU | None | elemmap | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | MU_GUEST | None | float | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | STALL_MINUTES | None | int | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | STD_CELL | None | str | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | U_ANION_GATE | None | bool | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | U_GATE_ANIONS | None | words | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | U_OVERRIDE | None | elemmap | SPEC 内建默认 |
| opt-dft-cpu | gen_step3_energy.py | VACUUM_AXIS_POLICY | None | str | SPEC 内建默认 |
| opt-mace-cpu | templates/step.conf（技能级） | DEVICE | cpu | - | - |
| opt-mace-cpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| opt-mace-cpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | primitive | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | auto / 2d / 3d |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-3 | - | eV/A。结构优化+形成能 1e-3 足够（能量已收敛到 ~1e-5 eV）； |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | - |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | eV/A |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | STRESS_TOL | 0.05 | - | GPa。晶胞残余应力闸（RELAX_CELL=true 才有意义）； |
| opt-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| opt-mace-cpu | templates/step3_formation/step.conf（步骤级） | MU | Si:-5.386 | - | - |
| opt-mace-cpu | templates/step3_formation/step.conf（步骤级） | MU_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| opt-mace-cpu | gen_step1_mace_relax.py | CONDA_ENV | kc.DEFAULT_CONDA_ENV | str | SPEC 内建默认 |
| opt-mace-cpu | gen_step1_mace_relax.py | CONDA_SH | kc.DEFAULT_CONDA_SH | str | SPEC 内建默认 |
| opt-mace-cpu | gen_step1_mace_relax.py | MACE_MODEL_DIR | - | str | SPEC 内建默认 |
| opt-mace-gpu | templates/step.conf（技能级） | DEVICE | cuda | - | - |
| opt-mace-gpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| opt-mace-gpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | primitive | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | auto / 2d / 3d |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-3 | - | eV/A。结构优化+形成能 1e-3 足够（能量已收敛到 ~1e-5 eV）； |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | - |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | eV/A |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | STRESS_TOL | 0.05 | - | GPa。晶胞残余应力闸（RELAX_CELL=true 才有意义）； |
| opt-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| opt-mace-gpu | templates/step3_formation/step.conf（步骤级） | MU | Si:-5.386 | - | - |
| opt-mace-gpu | templates/step3_formation/step.conf（步骤级） | MU_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| opt-mace-gpu | gen_step1_mace_relax.py | CONDA_ENV | kc.DEFAULT_CONDA_ENV | str | SPEC 内建默认 |
| opt-mace-gpu | gen_step1_mace_relax.py | CONDA_SH | kc.DEFAULT_CONDA_SH | str | SPEC 内建默认 |
| opt-mace-gpu | gen_step1_mace_relax.py | MACE_MODEL_DIR | - | str | SPEC 内建默认 |
| phonon-dft-cpu | templates/step.conf（技能级） | FUNC | pbesol | - | - |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | FD_DISTANCE | 0.01 | - | phonopy 2 阶有限位移幅度(Å) |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | FUNC | auto | - | - |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | KSPACING | 0.03 | - | 超胞取力 K 点间距 |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | MAX_DISP | 200 | - | 位移帧数硬闸 |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | MAX_MULTIPLE | 6 | - | - |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | MIN_SC_LEN | 12.0 | - | 超胞每个非真空方向最小胞长(Å) |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | PHONON_MESH | 20 20 20 | - | q 网格（写进 kl_params 供 step3；2D 自动 kz=1） |
| phonon-dft-cpu | templates/step2_disp/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "3 3 3"；空=按 MIN_SC_LEN 自动。2D 真空方向压 1 |
| phonon-dft-cpu | templates/step3_phonon/step.conf（步骤级） | BAND_POINTS | 51 | - | - |
| phonon-dft-cpu | templates/step3_phonon/step.conf（步骤级） | IMAG_THR | 0.10 | - | 虚频阈值(THz)：min_freq >= -IMAG_THR 判稳定 |
| phonon-dft-cpu | gen_step2_disp.py | ENCUT | None | int | SPEC 内建默认 |
| phonon-dft-cpu | gen_step2_disp.py | ENCUT_FACTOR | 1.5 | float | SPEC 内建默认 |
| phonon-dft-cpu | gen_step2_disp.py | KSCHEME | 2 | str | SPEC 内建默认 |
| phonon-dft-cpu | gen_step2_disp.py | VASPKIT_EXE | vaspkit | str | SPEC 内建默认 |
| phonon-mace-cpu | templates/step.conf（技能级） | DEVICE | cpu | - | - |
| phonon-mace-cpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| phonon-mace-cpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | primitive | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-4 | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | - |
| phonon-mace-cpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | ALM_CUT2 | 7.0 | - | ALM 二阶力常数截断半径(Å)，避免搜索整个超胞 |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | DISP_DISTANCE | 0.01 | - | hiphive MC-rattle 目标位移 RMS(Å)；声子 2 阶用 0.01 |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MAX_MULTIPLE | 6 | - | - |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MC_DMIN_SCALE | 0.85 | - | d_min = 最近邻 × 此系数（防止原子靠太近） |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MC_NITER | 10 | - | MC 迭代数 |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | MIN_SC_LEN | 15.0 | - | 超胞每个非真空方向的最小胞长(Å) |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | N_DISP | auto | - | 随机位移帧数：auto=按 ALM 2 阶自由力常数反推（下限 10） |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | OVERSAMPLE | 3 | - | N = ceil(nfree_fc2/DOF) × OVERSAMPLE |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | RANDOM_SEED | 2025 | - | - |
| phonon-mace-cpu | templates/step2_disp_force/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "4 4 4"；空=按 MIN_SC_LEN 自动 |
| phonon-mace-cpu | templates/step3_phonon/step.conf（步骤级） | FC_CUTOFF | 6.0 | - | - |
| phonon-mace-cpu | templates/step3_phonon/step.conf（步骤级） | MAX_FREQ_MIN_THZ | 20.0 | - | - |
| phonon-mace-cpu | templates/step4_cutoff_conv/step.conf（步骤级） | CUTOFF_LIST | 4.0,6.0,8.0 | - | - |
| phonon-mace-cpu | gen_step2_disp_phonon.py | CONDA_ENV | kc.DEFAULT_CONDA_ENV | str | SPEC 内建默认 |
| phonon-mace-cpu | gen_step2_disp_phonon.py | CONDA_SH | kc.DEFAULT_CONDA_SH | str | SPEC 内建默认 |
| phonon-mace-cpu | gen_step2_disp_phonon.py | MACE_MODEL_DIR | - | str | SPEC 内建默认 |
| phonon-mace-gpu | templates/step.conf（技能级） | DEVICE | cuda | - | - |
| phonon-mace-gpu | templates/step.conf（技能级） | DTYPE | float64 | - | - |
| phonon-mace-gpu | templates/step.conf（技能级） | MACE_MODEL | MACE-matpes-pbe-omat-ft.model | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | CELL_POLICY | primitive | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | DIMENSION | auto | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FIX_SYMMETRY | true | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | FMAX | 1e-4 | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | MAX_STEPS | 2000 | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | OPTIMIZER | FIRE | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX | true | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RELAX_CELL | true | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | RESIDUAL_TOL | 2e-3 | - | - |
| phonon-mace-gpu | templates/step1_mace_relax/step.conf（步骤级） | SYMPREC | 1e-4 | - | - |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | DISP_DISTANCE | 0.01 | - | hiphive MC-rattle 目标位移 RMS(Å)；声子 2 阶用 0.01 |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MAX_MULTIPLE | 6 | - | - |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MC_DMIN_SCALE | 0.85 | - | d_min = 最近邻 × 此系数（防止原子靠太近） |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MC_NITER | 10 | - | MC 迭代数 |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | MIN_SC_LEN | 15.0 | - | 超胞每个非真空方向的最小胞长(Å) |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | N_DISP | auto | - | 随机位移帧数：auto=按 ALM 2 阶自由力常数反推（下限 10） |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | OVERSAMPLE | 3 | - | N = ceil(nfree_fc2/DOF) × OVERSAMPLE |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | RANDOM_SEED | 2025 | - | - |
| phonon-mace-gpu | templates/step2_disp_force/step.conf（步骤级） | SUPERCELL | (空) | - | 显式 "4 4 4"；空=按 MIN_SC_LEN 自动 |
| phonon-mace-gpu | gen_step2_disp_phonon.py | CONDA_ENV | kc.DEFAULT_CONDA_ENV | str | SPEC 内建默认 |
| phonon-mace-gpu | gen_step2_disp_phonon.py | CONDA_SH | kc.DEFAULT_CONDA_SH | str | SPEC 内建默认 |
| phonon-mace-gpu | gen_step2_disp_phonon.py | MACE_MODEL_DIR | - | str | SPEC 内建默认 |
| unihamgnn | templates/step.conf（技能级） | CONDA_ENV | ML | - | - |
| unihamgnn | templates/step.conf（技能级） | CONDA_SH | /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh | - | - |
| unihamgnn | templates/step.conf（技能级） | DEVICE | cpu | - | predict 设备 cpu/cuda（3090 纯 CPU 跑） |
| unihamgnn | templates/step.conf（技能级） | DFT_DATA | /home/wangchaoyue852/software/Uni-HamGNN/DFT_DATA19 | - | - |
| unihamgnn | templates/step.conf（技能级） | DIMENSION | auto | - | auto / 2d / 3d |
| unihamgnn | templates/step.conf（技能级） | ELECTRONIC_TEMP | 300.0 | - | K |
| unihamgnn | templates/step.conf（技能级） | ENERGY_CUTOFF | 200.0 | - | Ry |
| unihamgnn | templates/step.conf（技能级） | HAMGNN_DIR | /home/wangchaoyue852/software/Uni-HamGNN/HamGNN | - | - |
| unihamgnn | templates/step.conf（技能级） | KGRID | 5 5 5 | - | - |
| unihamgnn | templates/step.conf（技能级） | MAX_SCF_ITER | 300 | - | - |
| unihamgnn | templates/step.conf（技能级） | MKL_LIB | /home/wangchaoyue852/miniconda3/pkgs/mkl-2025.3.1-h0e700b2_10/lib | - | - |
| unihamgnn | templates/step.conf（技能级） | MPIRUN | /home/wangchaoyue852/miniconda3/envs/openmx_build/bin/mpirun | - | - |
| unihamgnn | templates/step.conf（技能级） | NAO_MAX | 26 | - | OpenMX 最大轨道数：14 / 19 / 26 |
| unihamgnn | templates/step.conf（技能级） | NK | 120 | - | band_cal 能带路径 k 点数 |
| unihamgnn | templates/step.conf（技能级） | NPROC | 1 | - | - |
| unihamgnn | templates/step.conf（技能级） | NTHREADS | 16 | - | predict 的 OMP 线程数 |
| unihamgnn | templates/step.conf（技能级） | OPENMX_POSTPROCESS | (空) | - | - |
| unihamgnn | templates/step.conf（技能级） | READ_OPENMX | (空) | - | - |
| unihamgnn | templates/step.conf（技能级） | SCF_CRITERION | 1.0e-7 | - | Hartree |
| unihamgnn | templates/step.conf（技能级） | SOC | true | - | 通用 SOC 模型：true = 生成 non-SOC + SOC 两份 graph_data |
| unihamgnn | templates/step.conf（技能级） | UNI_MODEL | /home/wangchaoyue852/software/Uni-HamGNN/uni-hamgnn_2_1.pkl | - | - |
| unihamgnn | templates/step.conf（技能级） | XC | GGA-PBE | - | - |
| zt-dft-cpu | step.conf（技能级） | BANDGAP | hse | - | - |
| zt-dft-cpu | step.conf（技能级） | FUNC | pbesol | - | - |
| zt-dft-cpu | step20_zt/step.conf（步骤级） | KL_CONST_T | 300 | - | - |
| zt-dft-cpu | step20_zt/step.conf（步骤级） | KTEMP_MODE | interp | - | - |

统计：20 个技能，493 条 step.conf/SPEC 默认值，45 条声明参数。
