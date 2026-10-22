# unihamgnn 技能状态：暂不可用（2026-10-22 决定放弃）

## 决定

**放弃 `unihamgnn` 技能的实际计算支持**，将其加入全局技能黑名单（`setting/tf.yaml` 的
`disabled_skills: [unihamgnn]`，该文件 gitignore 不入库）。已完成的环境部署（工具、模型、基组数据）
已从各集群移除，避免 `autozt check` 出现“就绪但不可用”的误导。

## 为什么放弃

1. **问题不在部署，在工具本身**：`openmx_postprocess`（HamGNN 自带的 OpenMX 修改版）在 Si 结构上
   始终卡在基组截断阶段：
   ```
   <truncation> CpyCell=14 ct_AN=   1 FNAN SNAN 13670   0
   <truncation> CpyCell=14 ct_AN=   2 FNAN SNAN 13670   0
   ```
   日志显示 `PAO/VPS of species Si were normally found`（数据读取正常），但 240 秒内不出现
   `Finish calculating S & H0`，即 truncation/基组处理死循环（NaN 计数 13670）。
2. **在参考主机上同样复现**：把同一份 `openmx.dat` 拿到最早部署、此前被认为可用的 3090 上，
   用其原版二进制、原版 `DFT_DATA19`、正确 `MKL_LIB` 运行，结果同样卡死超时；
   用新生成的结构数据复测也失败（rc=124）。
3. **已排除的假设**：
   - 不是“复制了错误二进制”——在 jzzn 上用本地工具链从源码重编了 `openmx_postprocess/read_openmx`
     （系统 gcc/gfortran 8.3 + OpenMPI 4.0.5 + MKL 2019 + 自编 GSL 2.7.1 + 系统 FFTW 3.3.8，
     修了 gfortran8 不支持 `-fallow-argument-mismatch`、ELPA 并行模块竞态、GSL include 三个构建问题），
     编出的二进制仍然同样卡死。
   - 不是“数据缺失/损坏”——两套 `DFT_DATA19` 内容逐文件一致（566 文件：PAO 341 + VPS 224 + vps_info），
     日志确认 PAO/VPS 均正常读取。
   - 不是“资源不够”——缩短 K 网格/截断能、单线程、MKL 路径显式指定后均无改善。

结论：现有 HamGNN/OpenMX 版本组合下，`openmx_postprocess` 对 Si 的基组截断处理存在兼容问题，
需要上游修复或改用作者的 `honpas_1.2_H0` 路线；继续盲调投入产出比过低。

## 已尝试的部署（供将来重启时参考）

- 分发方式：以 yue 的 `~/software/AutoZT/Uni-HamGNN` 打包为可移植 bundle（5 个可执行文件 +
  依赖动态库 + HamGNN 源码 `HamGNN-main`），上传到 jzzn / hfeshell / clus_xa / hanhai25 / 3090；
- Python 侧：`torch-geometric 2.8.1`、`natsort`、`opt_einsum_fx`、`torch-ema`、`lightning-utilities`
  等依赖已装进各集群环境；
- 数据：`uni-hamgnn_2_1.pkl`（885 MB，Zenodo 17239078）+ `DFT_DATA19`（382,244,449 字节 / 566 文件）；
- 3090 参考部署：conda `ML` 环境（torch 1.11）+ conda `openmx_build`（gcc/gfortran 14、OpenMPI、GSL、MKL）。

## 恢复步骤（若将来上游修复）

1. 复核上游 HamGNN/OpenMX 版本与 `DFT_DATA19` 的基组命名兼容性，或改用 `honpas` 路线；
2. 重新部署工具与数据（参考上面的分发方式，或按 `skill/unihamgnn/README.md` 的 3090 配方本地编译）；
3. 在 `setting/tf.yaml` 中移除 `disabled_skills` 里的 `unihamgnn`；
4. 用 Si 等小体系跑 S1→S2→S3 端到端验收，确认 `openmx_postprocess` 能写出
   `Finish calculating S & H0` 后再释放给用户。

## 影响

- 其余 20 个 CPU 技能不受影响（五台 CPU 超算的 VASP/MACE/声子/eph-QE/cohp 等均正常）；
- 3090 的 `eph-qe-cpu`（QE 未编译完成）与 `kl-gpumd`（gpumd/nep 未编译）仍是独立遗留项。
