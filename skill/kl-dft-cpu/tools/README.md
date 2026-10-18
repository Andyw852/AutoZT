# kl-dft-cpu/tools —— 已有结果的补出图/补制表工具

这些脚本**只读取已算好的结果**，不重算、不提交作业，`step6_kappa` 跑完后随时可跑。

## make_kappa_plots.py

kl-dft-cpu 的 `step6_kappa` 只产出 `kappa_summary.json` / `kappa-m*.hdf5`，
**技能内没有 κL 绘图的 gen 步骤**；本脚本补齐这部分，并顺带把 `step5_phonon_plot`
（S5.1_plot 自动出的声子谱）一并归集。

```bash
python make_kappa_plots.py <kl-dft-cpu 目录> [<输出目录>] [<材料名>]
```

- `<kl-dft-cpu 目录>`：含 `step6_kappa/kappa_summary.json` 的目录；
  也接受材料根目录（自动在 `./` 或 `./*/` 下找 `step6_kappa`）。
- `<输出目录>`：缺省 `<材料根>/kl_plots`。
- `<材料名>`：缺省由目录名推断。

### 产物

| 文件 | 内容 |
|---|---|
| `<M>_kl_vs_T.png` | κL(T)：面内（左轴）+ 2D 归一化 κ·h/d（右轴） |
| `<M>_kl_vs_T.csv` | 温度依赖长表（含 xx/yy/zz 分量） |
| `<M>_kl_mesh_convergence.png` | q 网格收敛 + κ∞ 外推 |
| `<M>_kl_cutoff_selection.png` | FC3 截断扫描（`plateau_ok` 绿圈标记，标题带 status） |
| `<M>_kl_anisotropy.png` | xx / yy / zz 分量（2D 的 zz 无物理意义，仅参考） |
| `<M>_kl_summary.csv` / `.json` | 300 K 标量汇总（文献格式） |
| `phonon_band/` | S5.1_plot 声子谱图（存在才复制） |

### 说明

- 只写 `<输出目录>`，绝不写进 `step6_kappa/`。
- `cutoff_selection.records` 的字段名是 `cut` / `kappa` / `kappa_err` / `plateau_ok`。
- 2D 的 `kappa_zz` 无物理意义（真空方向），仅供核对各向异性。
- 原始 κL 依赖 `Lz`；2D 归一化 `κ·h/d` 与 `Lz` 无关，跨算例可比。
