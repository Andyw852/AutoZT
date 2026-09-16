# step8.4_amset2d -- AMSET 二维散射核（插件式）

> 只在 2D 项目里打开；三维项目走原版 `step8_amset`，完全不受影响。
> 开关：项目 `project_setting` 里写 `amset2d: true`（技能里该组 `default: false`）。

## 这个步骤解决什么

原版 AMSET 的散射核是三维的。把它直接用在真空 slab 上会有三个问题：

1. **弹性张量顺序**（三维也中招，已在 gen_step10 / gen_step14 一起修）
   VASP `TOTAL ELASTIC MODULI` 的行列序是 `XX YY ZZ XY YZ ZX`，而 AMSET 经 pymatgen
   按标准 Voigt `XX YY ZZ YZ XZ XY` 解析。不重排就把面内剪切 `C66(XY)` 塞进 `YZ` 槽位。
   实测 CrS2 单层：`C66 = 22.29 GPa` 被换成 `C_zxzx = -0.10 GPa`（面内 TA 支刚度近零）。
2. **ADP 不该乘 c/t**。AMSET 的散射积分在超胞里做，二维核的正确输入是**原始 slab
   弹性常数**（`= C_2D/c`）。乘了 c/t 会让 ADP 散射率偏小 `t/c` 倍（CrS2 约 1/3，
   对应迁移率偏大 3.06 倍）。
3. **散射核的 q 依赖**。二维 Fröhlich / 压电 / 带电杂质 / 声学形变势在 `q->0` 与
   大 `q` 的行为都与三维不同；真空方向的 `q_z` 在二维里必须积分掉（映射系数见下）。

## 映射关系（一句话）

AMSET 的三维积分里，真空方向 `q_z` 的积分贡献因子 `2*pi/c`，因此

```
G_s = c * G_2D
```

`G_s` 是插件塞进 AMSET `prefactor x factor` 的量。对 ADP 的等价说法：直接用原始
slab 弹性常数（`C_2D/c`）即可，不要再乘 c/t。推导与各机制的式子见
`../METHODOLOGY.md` 第 4 节；实测验证见 `VERIFICATION.md`。

## 目录里的文件

| 文件 | 作用 |
|---|---|
| `gen_step14_amset2d.py` | 登录节点 gen：写 settings.yaml + 2d_correction.json + submit.sh，复制插件 |
| `amset2d_plugin.py` | 运行期插件：把 ADP/POP/IMP/PIE 四类换成二维核（同名注册表替换，不改 AMSET 安装） |
| `amset2d_pop_freq.py` | 取 Gamma 点**面内极性模**有效频率（AMSET 自带口径会把面外 ZO 模算进去） |
| `test_kernels.py` | 四类核的自检（旋转不变、极限行为、护栏），`python test_kernels.py` 应全 PASS |
| `example_2d_correction.json` | 自检用的最小 2d_correction.json 示例 |
| `VERIFICATION.md` | 在真实 2D 数据上做的端到端验证记录（T1/T2 + 三材料对比） |

## 运行前自检

```bash
cd skill/ke-dft-cpu/step8.4_amset2d
python test_kernels.py          # 需要 amset + pymatgen 环境（如 conda amset_clean）
```

集群上确认 `python -c "import amset; print(amset.__version__)"`：插件在
**0.4.19 与 0.5.1** 上都核对过（两者的 `scattering/elastic.py`、`inelastic.py`、
`common.py` 逐字节相同）。其它版本会打印 WARN 但继续跑，请自行核对散射接口。

跑完后在 `amset.log` 里应看到：

```
[amset2d] amset 0.4.19 | c=20.000 A | r_inf=[[...]] | mechanisms=ADP,POP,IMP
[amset2d] worker: ADP=AcousticDeformation2D POP=PolarOptical2D IMP=IonizedImpurity2D PIE=Piezoelectric2D
```

没有这两行说明插件没生效（子进程尤其要看第二行 -- ADP 是以 reference 传给 spawn
子进程的，插件没被导入就会**静默**退回原版）。

## 已知局限（写论文要写进方法学）

- 偶极近似，没有四极矩；文献显示这会带来 20-50% 的迁移率差异。
- 极性散射只用**单一有效声子模**。
- 忽略 Janus 结构的面外偶极耦合。
- 自由载流子屏蔽用局域场 + 二维 Thomas-Fermi（文献对比偏低 10-20%）。
- 要求层法向平行于笛卡尔 z（插件会检查）。
- 仍受 AMSET 自身限制：无非极性光学形变势（ODP）、无谷间散射。

要作为定量结论发表，建议选 1-2 个材料用 EPW / Perturbo（开二维库仑截断）对照。
