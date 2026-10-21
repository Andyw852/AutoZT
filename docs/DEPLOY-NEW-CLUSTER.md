# 新超算 / 新账号接入 AutoZT —— 通用部署手册

> 目标：拿到一个新集群（或同集群新账号）后，按本手册可以**不靠原作者**把 AutoZT 布通：
> 跑得了 `autozt check` → 提交得了作业 → 结果数值可对照。
>
> 本手册**不假设**集群长什么样：SLURM / PBS / 无调度器、Intel / GCC、MKL / OpenBLAS、
> NVIDIA / 海光 DCU / 无加速卡、有无 module 系统，都按“先侦察、再分支”的方式处理。
> 每节都给了**判据**（怎么算这一步完成）和**常见变体**（别人集群不一样时怎么办）。
>
> 配套实例：`docs/CLUS_XA-DEPLOY-20261019.md`（SCNet 西安 clus_xa 的实测记录，含所有坑的现场数据）。
> 配套工具：`tools/probe_cluster.sh`（只读侦察，见第 1 节）。

---

## 0. 总览：新集群要交付什么

| # | 交付物 | 位置 | 判据 |
|---|--------|------|------|
| 1 | SSH 别名可用（BatchMode 免密） | 本机 `~/.ssh/config` | `ssh <别名> true` 返回 0 |
| 2 | 编译/运行环境 | 远端 `~/software/AutoZT/` | `ldd` 无 not found，二进制能跑 |
| 3 | POTCAR 库 + VASPKIT 配置 | `~/software/AutoZT/vasp_pseudopotentials` + `~/.vaspkit` | `autozt check` 不报 potcar 缺失 |
| 4 | Python 环境（venv/conda） | 远端 `~/software/AutoZT/venvs/<name>` | 非交互 `python3 -c "import numpy,pymatgen,spglib"` 通过 |
| 5 | 集群配置 | 本机 `setting/<hpc>.yaml` | `autozt doctor` 认得出该集群 |
| 6 | 提交模板 | 本机 `setting/<hpc>/templates/` | 生成的 `submit.sh` 能提交 |
| 7 | 验收记录 | `docs/<HPC>-DEPLOY-<日期>.md` | `autozt check` 全绿 + 冒烟数据 |

**原则**：
1. **只读优先**：侦察阶段全部只读（`ssh` 跑 `sinfo/ldd/ls` 可以，`sbatch` 只用 `--test-only`）。
2. **最小可用再扩展**：先 VASP 全链路，再 QE/amset/MACE 等可选组件。
3. **每步留判据**：不要“看起来装好了”——用 `ldd`、`autozt check`、一个真实小作业说话。
4. **集群私有信息不进仓库**：`setting/*` 默认被 `.gitignore` 排除；文档里用占位符代替个人路径。

---

## 1. 第一步：侦察集群（不要凭经验猜）

```bash
# 本机
tools/probe_cluster.sh

# 远端（先把 SSH 别名配好，见第 2 节）
tools/probe_cluster.sh <ssh别名>

# 加 --test-only 时才会向调度器做“空提交预检”（不占资源）
tools/probe_cluster.sh <ssh别名> --test-only
```

侦察报告覆盖：基本信息、调度器、分区/GRES/时限、账户配额、存储、module/编译器/MPI、
数学库/FFT、加速卡、Python/conda、AutoZT 现状、POTCAR 候选、可提交性预检。

### 1.1 必须手工确认的 6 个数

| 问题 | 命令（只读） | 记到 yaml |
|------|--------------|-----------|
| 能提交哪个分区？ | 逐分区 `sbatch --test-only`（probe 的 --test-only） | `partition` |
| 每节点多少核/多少内存？ | `sinfo -o "%c %m"`、`scontrol show node` | `max_cpus`（每节点上限） |
| 账户并发上限？ | `sacctmgr show assoc` 的 GrpTRES；或超配提交报错反推 | 写进 yaml 注释 |
| 单作业最长多久？ | `sinfo -o "%l"` / `scontrol show partition` | 模板 `--time` |
| 有没有 GPU/DCU？几张？ | `nvidia-smi -L` / `hy-smi` / `sinfo GRES` | 是否用 GPU 技能 |
| 家目录/工作目录在哪？多大？ | `df -h ~`、`quota` | `work_dir` |

> 变体：集群没有 `sinfo`（PBS/SGE/LSF/裸机）→ 用 `qstat -Q` / `bqueues` / 问管理员；
> AutoZT 的 `start` 需要 `sbatch/squeue/scancel` 语义，非 SLURM 需适配层（第 9.5 节）。

---

## 2. SSH：别名、密钥、直连与跳板

### 2.1 最小可用

```sshconfig
Host <alias>
    HostName <host>
    Port <port>
    User <user>
    IdentityFile ~/.ssh/id_rsa_<name>
    IdentitiesOnly yes
    ServerAliveInterval 30
    ServerAliveCountMax 6
    ControlMaster auto
    ControlPath ~/.ssh/cm-%r@%h:%p
    ControlPersist 4h
```

判据：`ssh -o BatchMode=yes <alias> true && echo OK`。密钥权限 `600`、目录 `700`。

### 2.2 直连被墙/只能跳板

```sshconfig
Host <alias>
    HostName <real-host>
    ProxyCommand ssh -W %h:%p <jump-alias>
```

> 注意：跳板链上任何一跳的密钥/网络失效都会表现为“超时”，先分别验证每一跳。

### 2.3 ★ 非交互会话的两个坑（AutoZT 采集/提交都离不开）

1. **非交互 PATH 可能不含你的 venv**：AutoZT 远端采集用 `ssh host 'python3 ...'`，
   若登录 shell 才 source 的 `.bashrc` 没生效，就找不到 python3。
   - 方案 A（通用）：`mkdir -p ~/.ssh && echo 'PATH=...' > ~/.ssh/environment` 并让 sshd 允许（`PermitUserEnvironment yes`）；
   - 方案 B（不依赖 sshd 配置）：把 `python3` 软链到**本来就非交互在 PATH 里**且你可写的目录，例如 `~/perl5/bin/python3`（先看 `probe_cluster.sh` 第 8 节的“第一个可写 PATH 目录”）。
2. **module 只在交互 shell 生效**：作业脚本里必须显式 `module load ...`；`ssh host 'command'` 不会自动 load。

判据：`ssh <alias> 'python3 -V && which python3'` 在**非登录式**命令里也能成功（AutoZT 的 collect 就长这样）。

---

## 3. 远端目录布局（AGENTS 铁律 11）

```
~/software/AutoZT/
├── vasp/                    # 源码 + 编译产物（每个版本一个子目录）
│   ├── vasp.6.4.3-pristine/
│   ├── vasp.6.4.3-optcell/  # 若需要 2D 变胞
│   └── makefile.include     # 集群专用编译配置（留档）
├── vasp_pseudopotentials/   # potpaw_PBE_64/ 等
├── venvs/                   # 每个 Python 环境一个子目录（autozt / mace_cpu / cogito ...）
├── conda_envs/              # 若用 conda（amset 等）
├── mlff_models/             # 机器学习势权重
├── qe-*/ pseudo-*/          # 可选组件
└── work/                    # AutoZT 远端工作根（work_dir）
```

容量估算：VASP 编译 ~2 GB（源码+对象+3 个二进制×2 版本）；POTCAR ~300 MB；venv 0.6–2 GB；
MACE 模型 0.1–1 GB；QE ~0.6 GB；UniHamGNN 全套 ~3 GB。**先 `df -h` 确认家目录配额**。

判据：`du -sh ~/software/AutoZT/*` 与 `df -h ~` 都在预算内。

---

## 4. 编译/库选型决策表（不同集群差异最大的地方）

### 4.1 编译器

| 集群有啥 | 建议 |
|----------|------|
| 只有系统 gcc（常是 4.8/7.x 老版本） | 找 module（gcc/9+）或用 `scl enable devtoolset-*`；**老 gcc 会 ICE**（见 13 节） |
| Intel oneAPI / parallel_studio | 用 `mpiicc/mpiifort`，数值库直接 MKL，最省事 |
| 混合（gcc 编译 + Intel MKL） | 可行：`-mkl=sequential` 或手动链 MKL 的三件套 |

**ABI 原则**：MPI 库由谁编译，就用谁家的编译器和 Fortran runtime；`libgfortran` 版本混用会
“能编译、运行时报 undefined symbol / segfault”。

### 4.2 BLAS / LAPACK / FFT

| 来源 | 优点 | 注意 |
|------|------|------|
| MKL | 最快、最省心 | 与 gcc 混用要选 `intel_lp64+intel_thread` 或 `gnu_thread` 变体；`libiomp5` 在 `compiler/lib/intel64*` 下 |
| 系统 BLAS/LAPACK | `-l:liblapack.so.3 -l:libblas.so.3` | 老 LAPACK(<3.8) 有 `dgegv_`；新版本删了它 → 补源码（第 12.2 节） |
| OpenBLAS | 多线程好 | 注意 `OPENBLAS_NUM_THREADS` 与 MPI 线程超订 |
| FFTW3 | 通用 | 头文件/库分 double 版：`.../fftw3/double/gnu/{include,lib}` |

判据：`ldd <二进制> | grep "not found"` 为空；跑冒烟作业不出现 `undefined symbol`。

---

## 5. 源码构建通用流程（以 VASP 为例，其它同理）

```bash
module load <compiler> <mpi>          # 显式加载；作业脚本里也要写
cd ~/software/AutoZT/vasp
# 1) 解包源码
# 2) 写 makefile.include（本集群专用，留档在包里）
# 3) 先小规模试编，再并行
make -k -j 16 std > build_std.log 2>&1 || true
```

| 症状（build 日志） | 根因 | 处理 |
|---|---|---|
| `internal compiler error: Segmentation fault` | 编译器太老/与源码不兼容 | 换新版 gcc/Intel |
| `cannot find -lxxx` | 库路径/名字不对 | 用 `ldconfig -p`、module show 找真实路径，必要时 `-L... -Wl,-rpath,...` |
| `undefined reference to 'xxx_'`（Fortran 符号） | 缺库/链接顺序/ABI | 把库放最后；确认 MPI/Fortran runtime 与编译器同一套 |
| 大量 `.mod` 找不到、并发编译时随机失败 | 并行 make 的模块依赖竞态 | `make -k -j N` **迭代多轮**直到收敛（每轮补一批） |
| `invalid command 'bdist_wheel'`（Python 包） | 环境缺 wheel | `pip install wheel` 后重试 |

**安装后必做**：`ldd` 检查 + 用**最小输入**真跑一次（第 6.3 节），别只看“编过了”。

---

## 6. VASP 接入（第一个要通的东西）

### 6.1 makefile.include 关键项（通用骨架）

```make
FC = mpif90
FCL = mpif90
CPP_OPTIONS = -DHOST=\"LinuxGNU\" -DMPI -DMPI_BLOCK=8000 -Duse_collective \
              -DCACHE_SIZE=4000 -Davoidalloc -Dvasp6 -Duse_bse_te -Dtbdyn -Dfock_dblbuf
BLASPACK = -L<blas/lapack dir> -l:liblapack.so.3 -l:libblas.so.3
LLIBS   += $(BLASPACK) -L<fftw lib> -Wl,-rpath,<fftw lib> -l:libfftw3.so.3
INCS    += -I<fftw include>
OBJECTS_LIB = linpack_double.o dgegv.o dggsvd.o   # 仅当 LAPACK 太新、缺 dgegv/dggsvd 时
```

> 别照抄别人的 `makefile.include`——库名/路径/可选宏都跟集群相关；用 `make -n` 或编译日志确认。

### 6.2 2D 变胞（optcell）

若材料需要“面内固定、真空方向放开”的变胞，用打过补丁的源码编一份（`relax_2d` profile 指向它），
否则用原版。两份都保留在 `vasp/` 下，yaml 里用两个 profile 区分。

### 6.3 冒烟作业（每个新集群必做）

```bash
# 目录里放 POSCAR(2 原子 Si 立方)/INCAR/KPOINTS/POTCAR，用集群模板提交
# INCAR 最小集： ENCUT=250 ISMEAR=0 SIGMA=0.05 EDIFF=1e-5 NELM=20 NCORE=<核数>
# 判据： exit code 0、OUTCAR 含 "reached required accuracy"、General timing 有数
```

**数值基线**：同一材料在另一台机器上跑过的总能用同一个 POSCAR/INCAR/KPOINTS/POTCAR 比
（Si 立方 2 原子、ENCUT 250、2×2×2 Γ 中心常作基准）。差 >1e-4 eV 要查 POTCAR/ENCUT/版本。

### 6.4 并行参数（通用原则）

- `NCORE`（VASP 5+）≈ sqrt(每节点核数) 或 4–8；`KPAR` 只在 k 点足够多时开。
- 单节点核数 ≥ 64 时，整节点提交更划算；但**整节点要等完全空闲节点**，忙时排队更久。
- 账户并发上限 = 核数预算：`jobcores × max_jobs ≤ 配额`。例：320 核 → 32 核×10 个作业。
- 作业时长设账户/分区上限（如 3 天）；短作业（<1 天）通常排得更快。

---

## 7. POTCAR、VASPKIT 与 Python 环境

### 7.1 POTCAR 三种布局（`potcar_dir` 要指向**哪一层**看你的目录结构）

```
A) potcar_dir/potpaw_PBE_64/<El>/POTCAR     ← 推荐（AutoZT/defect 类按这个找）
B) potcar_dir/potpaw_PBE/<El>/POTCAR
C) potcar_dir/<El>/POTCAR                   ← 解包时没建父目录的情况（要搬一层）
```
解包官方 tar 包通常得到的是 C 形态（元素目录直接铺开），记得 `mkdir potpaw_PBE_64 && mv */ potpaw_PBE_64/`，
搬之前 `chmod -R u+w`（官方包目录常是只读 `r-x`）。

### 7.2 VASPKIT（生成 POTCAR 的工具，也是 `autozt check` 的检查点之一）

- VASPKIT 1.x 的 `~/.vaspkit` 用**空格分隔**：`PBE_PATH    /abs/path/potpaw_PBE_64`；
  写成 `KEY = VALUE` 时 vaspkit 会把值读成 `=`，**静默不生成 POTCAR**。
- 生成 POTCAR 的快速验证：任意目录放 POSCAR，`printf '1\n103\n' | vaspkit`，看有没有 `POTCAR` 且首行是 `PAW_PBE Si ...`。
- 若 `autozt check` 报 `vaspkit_potcar` 缺失但 vaspkit 明明能用：那是 check 的解析格式与 vaspkit 原生格式不一致
  （仓库已修：两种格式都认）；老版本可用 `KEY = VALUE` 临时通过 check，但会破坏 vaspkit，**不推荐**。

### 7.3 Python 环境（venv 还是 conda）

| 场景 | 建议 |
|------|------|
| 只用 AutoZT 内置技能 | venv 足够：`python -m venv`，pip 装 numpy/scipy/matplotlib/spglib/pymatgen/h5py/seekpath |
| 要用 amset / phonopy 链 / MACE | 组件各自独立环境（numpy 版本需求常冲突，见 12 节） |

**版本原则**：
- `numpy`：老 amset/pheasy/mace 生态多为 1.x；cogito-dft 0.3.13 要求 ≥2 → **不要塞同一个环境**；
- `pymatgen`：许多老工具不兼容太新的主版本，装完用 `import pymatgen`+实际脚本验证；
- 装完必须测**非交互**调用（第 2.3 节）。

判据：`ssh <alias> 'python3 -c "import numpy,scipy,pymatgen,spglib,h5py;print(1)"'`。

---

## 8. AutoZT 集群配置：`setting/<hpc>.yaml` 字段逐项

> 这是“这台机器怎么跑”的唯一事实来源。字段不要求全填，按需；**必填**见下表。
> 不同集群的例子（CPU/SLURM、GPU/SLURM、无调度器、DCU）见 8.2。

### 8.1 字段表

| 字段 | 必填 | 含义 / 取值 |
|------|------|-------------|
| `name` | ✅ | 集群名（`setting/<name>.yaml`、`-hpc`、`check --cluster` 都用它） |
| `ssh_host` | ✅ | `~/.ssh/config` 别名；**空字符串** = 本机模式（AutoZT 会接 fakeslurm） |
| `template_map` | | 逻辑模板名→文件映射；一般留 `{}`，模板按 `setting/<name>/templates/` 同名取 |
| `vasp.standard` | ✅* | 原版 VASP 的 `std/gam/ncl` 绝对路径 + `cell_constraint: none` |
| `vasp.relax_2d` | | 变胞版（`cell_constraint: ioptcell_tag`）；没有 2D 需求可不写 |
| `work_dir` | ✅ | 远端工作根（与 `setting` 里别的集群**不能相同**） |
| `partition` | ✅ | 默认提交分区 |
| `max_cpus` | | **每节点**核数上限（preflight 用来发现“模板申请超过单节点核数”） |
| `conda_sh` | | conda 的 `profile.d/conda.sh`；纯 venv 可留空 |
| `conda_env` | ✅* | gen 脚本用的 Python 环境；**含 `/` 按 venv 激活**，否则 conda env 名/前缀 |
| `pheasy_env` / `kl_env` | | 拟合/声子专用环境（缺省同上） |
| `amset_env` | | ke-dft-cpu S8 的 amset 环境（应为**独立**环境，见 12.2） |
| `potcar_dir` | ✅* | POTCAR 根（第 7.1 节三种布局） |
| `mlff_model_dir` | | MACE 权重目录（MLFF 技能用） |
| `remote_path_prefix` | | 远端命令前缀 PATH（如放 vaspkit 软链的 `pybin` + venv `bin`） |
| `qe_bin` / `qe_pseudo_dir` / `qe_setup` | | eph-qe-cpu 用：`pw.x/ph.x` 目录、UPF 目录、运行前 export/module 命令串 |
| 其它（`mlff_*`/`amset_*`/`abc_*` 等） | | 新技能可能需要；以 `autozt check` 报告“缺什么”为准 |

> ✅* = 该功能用到时才必填：只用 band → vasp+potcar+conda_env；用 ke → 再加 amset_env；用 eph → qe_*。

### 8.2 四类集群的最小 yaml 骨架

```yaml
# A) CPU + SLURM + module（最常见）
name: mycluster
ssh_host: mycluster
work_dir: /work/<user>/software/AutoZT/work
partition: cpu64
max_cpus: 64
vasp:
  standard: {cell_constraint: none, std: /work/<user>/software/AutoZT/vasp/.../vasp_std, gam: .../vasp_gam, ncl: .../vasp_ncl}
conda_sh: ""
conda_env: /work/<user>/software/AutoZT/venvs/autozt
potcar_dir: /work/<user>/software/AutoZT/vasp_pseudopotentials
```

```yaml
# B) GPU + SLURM（模板要带 --gres=gpu:N）
name: mygpu
ssh_host: mygpu
work_dir: /home/<user>/AutoZT/work
partition: a800
max_cpus: 128
conda_env: /home/<user>/AutoZT/venvs/mace_gpu   # 含 torch+cuda
mlff_model_dir: /home/<user>/AutoZT/mlff_models
```

```yaml
# C) 无调度器（裸服务器）：ssh_host 仍写别名；模板用 tools/fakeslurm 语义
name: myserver
ssh_host: myserver
work_dir: /home/<user>/AutoZT/work
conda_env: /home/<user>/AutoZT/venvs/autozt
# 提交模板用 setting/myserver/templates/ 下自备的（fakeslurm 兼容 #SBATCH 变量）
```

```yaml
# D) 本机（local 模式）
name: local
ssh_host: ""
work_dir: /tmp/autozt_work
conda_env: ""
```

### 8.3 模板目录

```
setting/<hpc>/templates/
├── submit_std_3d.tpl   submit_std_2d.tpl   submit_std_0d.tpl
├── submit_ncl_3d.tpl   submit_ncl_2d.tpl
└── step1_opt/submit_std_2d.tpl
```

模板里一般包含：`#SBATCH -p <partition>`、`-N/-n/--ntasks-per-node`、`--time`、`module load`、
`export OMP_NUM_THREADS/OMPI_MCA_*`、`mpirun ... vasp_std`、可选 `cd $SLURM_SUBMIT_DIR`。
**从别的集群复制模板时，逐项核对**：分区名、module 名、mpirun 路径/参数、VASP 路径（AutoZT 渲染时会按
`setting/<hpc>.yaml` 的 vasp 路径替换，但模板里写死的 `module` 不会）。

判据：`autozt -p <材料> -j <步骤> init` 后，检查远端 `submit.sh` 的分区/核数/路径 都对，再 `start`。

---

## 9. 提交与并行：跨集群差异最大的一节

### 9.1 核数/时长/并发

1. 先看分区上限（`sinfo`）、再看账户并发（`GrpTRES cpu=`）。
2. `每作业核数 × 同时作业数 ≤ 账户并发`；预留 1 个槽给调试。
3. 单作业跨节点：`-N n` + `-n n*ppn`，注意 `--ntasks-per-node` 与 MPI 进程数一致。
4. 长作业上限：模板 `--time` 不要超过分区 `MaxTime`，否则 sbatch 直接拒。

### 9.2 module / 环境

- 作业脚本里**显式** `module load`；不要依赖登录 shell。
- 需要 venv/conda 的步骤：脚本里 `source <env>/bin/activate` 或 `conda activate`，
  并在 `step.conf` 里用 `CONDA_ENV` 指定（AutoZT 会把集群 yaml 的 `conda_env` 注入）。
- 运行时库路径：优先 `-Wl,-rpath` 写进二进制；其次 `module load`；最后才是 `export LD_LIBRARY_PATH`。

### 9.3 ★ MPI 启动器与“绑核”坑（实测高频）

| 现象 | 原因 | 处理 |
|---|---|---|
| 作业 3 秒失败，stderr `binding more processes than cpus` | OpenMPI 默认绑核与节点/cgroup 不一致 | `mpirun --bind-to none -np $SLURM_NTASKS ...`（先在一台节点上冒烟确认） |
| Intel MPI 程序用 OpenMPI 的 `mpirun` 启动 | ABI 不匹配（能起进程、计算崩） | 用 Intel MPI 的 `mpiexec/mpirun`，或 `module load mpi/intelmpi` 后直接用 |
| `mpirun -np N` 启动 N 个完整程序 | 模板写成了“每个进程算全量” | 确认 VASP 模板是 `mpirun -np $SLURM_NTASKS vasp_std` 一类；不要双重并行 |
| 单节点共享内存 | `OMPI_MCA_btl=self,vader,tcp`、`pml=ob1` | 单节点作业常设（多节点作业不要设） |

### 9.4 GPU / DCU

- NVIDIA：模板 `#SBATCH --gres=gpu:N`；环境装 `torch+cuda`；检查 `nvidia-smi` 可见。
- 海光 DCU：`--gres=dcu:N`，DTK 工具链（`compiler/dtk/*`），PyTorch 要用对应 DTK 版本（不是官方 CUDA wheel）。
- 判据：作业里 `nvidia-smi -L`（或 `hy-smi`）看得到卡，且程序 `torch.cuda.is_available()` 为真。

### 9.5 非 SLURM 集群

AutoZT 通过 `sbatch/squeue/scancel` 语义管理作业：
- 裸服务器/无调度器：用仓库 `tools/fakeslurm`（AutoZT 自动接入本机模式；远端 `remote_path_prefix` 注入）；
- PBS/SGE/LSF：需写一层适配脚本（把 `sbatch` 翻译成 `qsub`），或在文档中标注“不支持自动推进”。

---

## 10. 步骤级/材料级切换（同一材料跨集群跑）

- 材料级：`autozt -p <材料> hpc <cluster>` 改全部步骤（会重跑/需产物）——大动作，先想清楚。
- 步骤级：项目 `tf_<项目>.yaml` 里写
  ```yaml
  step_hpc:
    S7_deform: clus_xa_yue
    S7.1_read: clus_xa_yue
  ```
  只动指定步骤，其余留原集群；依赖产物由 AutoZT 按 needs 跨集群回传。
- 项目级 step.conf 覆盖：`project_setting/templates/step.conf`（参数）与
  `project_setting/templates/<步骤>/step.conf`（步骤级，优先级最高）。
- 改完**先在 `autozt list` 看 HPC 列/任务图**，再 `retry + start`。

---

## 11. 验收清单（新集群布完必须逐条过）

| # | 检查 | 命令 | 期望 |
|---|------|------|------|
| 1 | 连通 | `ssh -o BatchMode=yes <alias> true` | exit 0 |
| 2 | 非交互 python | `ssh <alias> 'python3 -V'` | 打印版本 |
| 3 | 配置预检 | `autozt doctor` | 无报错，认得出新集群 |
| 4 | 技能就绪 | `autozt check --cluster <hpc>` | VASP 系技能 ✅（可选组件可先 ❌） |
| 5 | 提交模板 | `autozt -p <材料> -j <步骤> init` 后看远端 submit.sh | 分区/核数/路径正确 |
| 6 | 冒烟计算 | 最小 Si 作业（第 6.3 节） | OUTCAR reached required accuracy，能量与基线一致 |
| 7 | 并行上限 | `--test-only` / 实提 2 个作业 | 不被 QOS/分区拒；同提交者作业互不影响 |
| 8 | 清理 | `squeue -u <user>` | 0 个测试作业残留 |
| 9 | 记录 | 写 `docs/<HPC>-DEPLOY-<date>.md` | 含实测配额/坑/验收数字 |

---

## 12. 可选组件（用到再装，别一次全上）

### 12.1 QE（quantum-espresso，eph-qe-cpu 用）

1. 优先：集群已有 QE module（`module avail | grep -i quantum`）→ `qe_bin` 指向其 `bin`，
   `qe_setup` 里写 `module load` 或显式 `LD_LIBRARY_PATH`（Intel MKL/impi 常需要）。
2. 次选：conda-forge `quantum-espresso`（自包含，但要网络/镜像；注意 conda 求解内存）。
3. 判据：`pw.x -h`、`ph.x -h` 正常；跑一个 Si SCF（ONCV 赝势），`convergence has been achieved`。
4. 坑：集群 module 的 QE 常与已加载的 MPI 冲突（`module unload mpi/...` 后再 load），
   或者用 `qe_setup` 显式给 `LD_LIBRARY_PATH`（MKL + Intel MPI 的 lib 目录）。

### 12.2 amset（ke-dft-cpu 的 S8）

- **独立环境**（amset 0.5.x 通常要 numpy 1.x、pymatgen 老版本；与 MACE 的 numpy/torch 常常冲突）。
- pip 路线：`python -m venv` → `pip install amset==0.5.1 "pymatgen<2025" "numpy<2" scipy ...`；
  装完 `import amset` 与 `amset deform create --help` 验证。
- 坑：`BoltzTraP2` 从 PyPI sdist 编译常失败（自带 spglib 子模块缺 CMakeLists / 需要 netCDF4 头文件）。
  解决：用 conda-forge 预编译包解包后拷进 venv 的 site-packages（见 clus_xa 记录），或干脆不用 BT2 的插值功能。
- 集群 yaml：`amset_env: <env 绝对路径>`；`check` 会验证 `amset`/`BoltzTraP2` 可导入。

### 12.3 MACE（MLFF 系列）

- 环境：`torch`（CPU 集群装 `--index-url .../whl/cpu`）、`mace-torch`、`e3nn`、`ase`、`matscipy`；
  **注意 torch 版本与 e3nn 兼容性**（torch≥2.6 默认 `weights_only=True`，老 e3nn 会报
  `UnpicklingError: Unsupported global` → 降 torch 到 2.4.x 或升级 e3nn）。
- 权重：技能默认模型名见 `skill/*/templates/step.conf` 的 `MACE_MODEL`；从官方 release 下载到
  `mlff_model_dir`；`MACECalculator(model_paths=<file>)` 冒烟算一个 Si 单点。
- 集群 yaml：`conda_env`（含 torch/mace 的环境）或单独 `mlff_env`（若技能支持）。

### 12.4 COGITO（cohp-cogito）

- `pip install cogito-dft==<版本>`（纯 Python + numpy2）；版本与 jzzn 对齐（如 0.3.13）。
- **独立 venv**（要求 numpy≥2，别和 MACE/amset 混装）；把 `COGITO/COGITOanalyze/COGITOpost`
  软链进集群 `conda_env` 的 `bin`，或在项目 step.conf 里指定专用环境。
- 前置：VASP `LWAVE=.TRUE.` 的静态 WAVECAR、`ISYM` 用约化网格（不要 `ISYM=-1`）、NBANDS 足够。

### 12.5 UniHamGNN（unihamgnn 技能）

- 组件：HamGNN 仓库、通用模型 `.pkl`（~1 GB）、OpenMX 基组库 `DFT_DATA19`、自编译的
  `openmx_postprocess/read_openmx`、Python 环境（torch CPU + torch_geometric + lightning + e3nn）。
- 编译 `openmx_postprocess` 的坑：需要 GSL（静态链可省运行时 libcblas）、MPI、LAPACK、FFTW；
  若缺 ScaLAPACK/BLACS，可链 MKL 的 scalapack+blacs；ELPA 的 Fortran 模块并行编译会竞态，用
  `make -k -j N` 迭代多轮。
- 集群 yaml 里通常要写 `mlff_model_dir` 之外的项目级 step.conf（`HAMGNN_DIR/UNI_MODEL/DFT_DATA/
  OPENMX_POSTPROCESS/READ_OPENMX`），因为这些是技能参数，不是集群字段。

### 12.6 其它

- `POTCAR` 之外的赝势库（SSSP/ONCV/LDA）：按 `potcar_dir` 布局放；
- `VASPKIT` 新版本、`phonopy/phono3py/ALM`、`pheasy`：装进对应技能环境，用 `autozt check` 暴露缺项。

---

## 13. 故障速查表（症状 → 原因 → 处理）

| 症状 | 可能原因 | 处理 |
|------|----------|------|
| `sbatch: error: Invalid account or account/partition combination` | 账户无该分区权限 | 用 `--test-only` 逐分区找可用分区；问管理员加权限 |
| `CPU count per node can not be satisfied` | 单节点申请核数 > 节点核数（模板写死 192/256 等） | 改模板/step.conf 的每节点核数；或用 `max_cpus` 让 preflight 提前警告 |
| `group max tres limit ...` | 超账户 GPU/CPU 配额 | 降并发或降核数（需说明影响，AGENTS 铁律 9） |
| 作业 3 秒失败 + `binding more processes than cpus` | OpenMPI 绑核 | `mpirun --bind-to none` |
| `error while loading shared libraries` | 编译机与运行机库不同 / 没 rpath | `ldd` 找缺的库；加 `-Wl,-rpath` 或 module |
| `undefined symbol: dgegv_` | LAPACK≥3.8 删了老接口 | 补 `dgegv.f/dggsvd.f` 进源码（netlib） |
| VASP 一启动就 `I REFUSE TO CONTINUE` | POSCAR 未更新/结构不一致 | 检查 push 的 POSCAR 与本地一致（AutoZT 有 sha 校验） |
| vaspkit 不生成 POTCAR | `~/.vaspkit` 用了 `KEY = VALUE` | 改空格分隔；或检查 `PBE_PATH` 指向元素目录的父级 |
| `autozt check` 说缺 python 包但明明装了 | 检查的是**非交互** PATH/指定环境 | 用 venv 绝对路径/软链到 PATH；或检查 yaml 的 `conda_env` |
| 远端采集 `Connection closed` / `timeout` | 跳板/代理/ControlMaster 失效 | 重建连接（`ssh -O exit`）；检查跳板密钥有效期 |
| 两个会话互相杀作业 | 同一材料被两个 autozt 会话 retry/start | **一个材料同一时间只留一个会话操作**；用 `step_hpc`/worktree 隔离 |

---

## 14. 文档与记录

- 每个新集群写一份 `docs/<HPC>-DEPLOY-<YYYYMMDD>.md`：
  **实测配额表**（并发/单作业/时限/存储）、**构建踩坑**（编译器/库/MPI）、**验收数字**（冒烟能量/带隙）、
  **可选组件状态**（QE/amset/MACE/…）。直接复制 **`docs/DEPLOY-RECORD-TEMPLATE.md`** 按实测填
  （模板已把"本集群特殊项 / 通用项"分开）。
- `setting/*` 是本地文件（.gitignore）：可以在仓库里留一份**脱敏模板**（用 `<user>/<work_dir>` 占位），
  真实路径不提交。
- 与 AGENTS 的对应：安装位置=铁律 11；改核数流程=铁律 9；只读诊断=铁律 4；
  临时代码/脚本放 `tmp/`。

---

## 附录 A：最小 `setting/<hpc>.yaml` 模板（复制后填占位符）

```yaml
name: <HPC>
ssh_host: <ssh别名>
template_map: {}
vasp:
  standard:
    cell_constraint: none
    std: <work>/software/AutoZT/vasp/<ver>/bin/vasp_std
    gam: <work>/software/AutoZT/vasp/<ver>/bin/vasp_gam
    ncl: <work>/software/AutoZT/vasp/<ver>/bin/vasp_ncl
work_dir: <work>/software/AutoZT/work
partition: <分区>
max_cpus: <每节点核数>
conda_sh: ""
conda_env: <work>/software/AutoZT/venvs/autozt
potcar_dir: <work>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: <work>/software/AutoZT/pybin:<work>/software/AutoZT/venvs/autozt/bin
# 用 ke 技能再加：amset_env: ...
# 用 eph 技能再加：qe_bin/qe_pseudo_dir/qe_setup
# 用 MLFF 技能再加：mlff_model_dir: ...
```

## 附录 B：新集群部署顺序（一页纸）

1. `tools/probe_cluster.sh <别名>` → 抄下 6 个数（第 1.1 节）
2. 配 SSH + `~/.ssh/environment`/python3 软链（第 2 节）
3. `mkdir -p ~/software/AutoZT/{vasp,vasp_pseudopotentials,venvs,work}`（第 3 节）
4. 编 VASP（原版+可选 optcell），`ldd` 检查（第 5、6 节）
5. 解包 POTCAR、配 `~/.vaspkit`、建 venv（第 7 节）
6. 写 `setting/<HPC>.yaml` + 复制/改写提交模板（第 8、9 节）
7. `autozt doctor` → `autozt check --cluster <HPC>` → 冒烟作业（第 11 节）
8. 写部署记录，跑第一批真实材料（小体系 1 步验证吞吐/排队/IO）
9. 需要 QE/amset/MACE/COGITO/UniHamGNN 时按第 12 节逐个装，每装一个就 `autozt check` 一次
