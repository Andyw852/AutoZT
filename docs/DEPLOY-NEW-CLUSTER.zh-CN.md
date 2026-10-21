# 在新超算站点部署 AutoZT：移植与验收指南（中文版）

**文档版本：** 1.1（事实性修订：第 5.2、7.1、7.2、8.1、8.3、9.1、9.3.1、10.1、11.3、13、14.6 节；修订记录见仓库历史）  
**适用版本：** AutoZT v2.0 及以后  
**读者：** 负责在计算集群或独立节点上部署 AutoZT 的站点管理员与用户  
**配套文件：** `tools/probe_cluster.sh`（只读站点勘察脚本）、`docs/DEPLOY-RECORD-TEMPLATE.md`（站点记录模板）、`docs/CONFIGURING.md`（配置参考）  
**英文版：** `docs/DEPLOY-NEW-CLUSTER.md`（内容与本文对应）  
**许可：** MIT（见 `LICENSE`）

---

## 1. 适用范围

AutoZT 随软件分发的是一组模板，而不是可直接使用的站点配置。凡标识某一台机器的信息——
作业调度器、分区名称、模块系统、编译器、数学库与 FFT 库、赝势路径、暂存目录——均与站点相关，
必须由部署者提供。本指南给出一个可复现的部署流程：从首次只读勘察，到形成有据可查的验收记录。

本指南不对目标机器作任何假设，覆盖以下情形：

- 批处理调度器（SLURM、PBS/SGE/LSF，或无调度器）；
- 编译器栈（GNU、Intel，或二者混用）；
- 数学库（MKL、OpenBLAS、厂商 BLAS/LAPACK）；
- 加速卡（NVIDIA GPU、海光 DCU，或无加速卡）；
- 模块系统（Lmod/Tcl module、Spack，或仅使用环境变量）。

每一节都给出该步骤的**验收判据**，并列明站点之间的常见差异。实测示例见附录 C。

## 2. 术语

| 术语 | 本文中的含义 |
|---|---|
| 站点（site）/ 集群 | 通过 SSH 访问的计算资源（同一台机器可存在多个账号） |
| 站点配置文件（site profile） | 描述站点的 YAML 文件 `setting/<cluster>.yaml` |
| 技能（skill） | AutoZT 的工作流模块（`skill/<name>/skill.yaml`），定义步骤、判据与产物 |
| 运行时环境 | 供输入生成脚本使用的 Python 环境（virtualenv 或 conda 前缀） |
| 工作目录（work_dir） | 站点上用于创建材料目录的可写文件系统根目录 |
| 赝势库 | 电子结构计算所用的 POTCAR（或 UPF）数据集 |
| 提交模板 | `setting/<cluster>/templates/` 下的调度器脚本，在生成输入时渲染 |
| 验收测试 | 证明站点配置正确的证据（见第 13 节与第 16 节） |

## 3. 与既有实践的关系

AutoZT 采用科学软件中的常见做法：代码保持可移植，把与机器相关的参数集中在声明式配置文件中，
并以专门命令加以验证。代表性实例如下。

| 项目 | 站点相关配置 | 验证 / 诊断 |
|---|---|---|
| ASE | `siteconfig.py`（随包提供模板，每个站点生成一次） | `ase info`、导入测试 |
| AiiDA | `verdi computer setup`（主机名、调度器、工作目录、`mpirun` 命令、前后置命令） | `verdi computer test` |
| Quantum ESPRESSO | `make.inc` 与 `arch/` 机器文件 | 测试套编译 |
| CP2K | `arch/` 文件，每个平台一份 | `make -j` 与回归测试 |
| LAMMPS | `src/MAKE/MACHINES/Makefile.*` | `unittest/` 单元测试 |
| Snakemake | 执行器插件与站点 profile（`profiles/<site>/config.yaml`） | `--dry-run` |
| Nextflow | `nextflow.config` 中按站点命名的 `profiles`，以 `-profile` 选择 | `-preview` |
| jobflow-remote | `jobflow.yaml`（runner、workers、调度器、工作目录、前后置命令） | `jf runner check` |

上列项目共有、本文也遵循的约定为：

1. **分发模板，不分发可用配置。** 站点路径、主机名与凭据不进入版本库；`setting/*` 已被排除。
2. **每个站点一份声明式文件**，集中描述访问方式、调度器、工作目录、运行时环境与启动命令。
3. **通过站点自身机制发现环境**（模块系统、Spack 或环境变量），不做硬编码假设。
4. **提供预检命令。** 提交作业之前报告配置错误（`autozt check`，见第 13 节）。
5. **参考文档与操作指南分离。** 字段语义见 `docs/CONFIGURING.md` 与本文第 10 节；逐步操作见本指南。
6. **保留验收证据。** 版本、编译参数与冒烟测试输出记入站点记录（`docs/<SITE>-DEPLOY-<date>.md`），以保证可复现。

## 4. 站点勘察（只读）

在改动站点之前先完成勘察。`tools/probe_cluster.sh` 通过 SSH 收集以下信息，且不提交作业
（给出 `--test-only` 时除外）：

```bash
tools/probe_cluster.sh                       # 本机
tools/probe_cluster.sh <ssh别名>              # 远端站点
tools/probe_cluster.sh <ssh别名> --test-only   # 附加调度器可提交性检查
```

报告内容包括：主机与操作系统；可用的批处理调度器；SLURM 分区的节点数、核数、内存、GRES 与时限；
账户关联与 QOS 限额；文件系统与配额；模块系统及编译器/MPI 版本；数学库与 FFT 库；加速卡；
非交互 Python 环境；既有 AutoZT 安装；候选赝势目录。

### 4.1 需人工确认的六个量

| 问题 | 只读命令 | 记录位置 |
|---|---|---|
| 本账号可提交哪些分区？ | 逐分区 `sbatch --test-only` | `partition` |
| 每节点核数与内存？ | `sinfo -o "%c %m"`、`scontrol show node` | `max_cpus` |
| 账户级并发上限？ | `sacctmgr show assoc`（GrpTRES），或以超大 `--test-only` 请求反推 | 站点配置文件注释 |
| 单作业最长时限？ | `sinfo -o "%l"`、`scontrol show partition` | 模板中的 `--time` |
| 每节点加速卡数量？ | `nvidia-smi -L`、`hy-smi`、`sinfo` 的 GRES | 使用加速卡的技能 |
| 工作目录与家目录的位置、容量与配额？ | `df -h`、`quota` | `work_dir` |

批量探测多个分区时，两次请求之间应间隔一至两秒；部分调度器会限制提交频率
（`sbatch: request too frequently`），否则可能把可用分区误判为不可用。

## 5. 访问方式与非交互执行

### 5.1 SSH 配置

```sshconfig
Host <别名>
    HostName <主机>
    Port <端口>
    User <账号>
    IdentityFile ~/.ssh/id_rsa_<名称>
    IdentitiesOnly yes
    ServerAliveInterval 30
    ServerAliveCountMax 6
    ControlMaster auto
    ControlPath ~/.ssh/cm-%r@%h:%p
    ControlPersist 4h
```

验收判据：`ssh -o BatchMode=yes <别名> true` 返回 0。私钥仅属主可读（`chmod 600`），
所在目录权限为 `chmod 700`。

若集群只能经跳板机访问，使用 `ProxyCommand ssh -W %h:%p <跳板别名>`，并逐跳独立验证；
否则任一跳的故障都只表现为超时。

### 5.2 非交互会话的两项要求

AutoZT 以非交互方式调用远端主机（`ssh host 'python3 ...'`）完成作业采集与提交，因此：

1. **运行时 Python 必须位于非交互 `PATH` 中。** 仅在交互式启动文件中导出的变量不可见。
   常用处理方式：把 `PATH` 写入 `~/.ssh/environment`（需在 `sshd_config` 中启用
   `PermitUserEnvironment yes`，多数超算出于安全考虑不允许）；或在已位于非交互 `PATH` 且账号
   可写的目录中创建解释器软链（勘察脚本会报告第一个这样的目录）；或把运行时环境的 `bin`
   通过 `remote_path_prefix` 注入。
2. **模块命令不会被继承。** 每个作业脚本都必须显式加载模块；非交互 SSH 命令不执行登录 shell 的模块初始化。

验收判据：`ssh <别名> 'python3 -V && command -v python3'` 成功，且报告的解释器为预期运行时环境。

## 6. 目录布局与存储

约定布局如下：

```
~/software/AutoZT/
├── vasp/                    # 源代码与二进制，每个版本一个目录
│   ├── vasp.6.4.3-pristine/
│   ├── vasp.6.4.3-optcell/  # 仅在需要变胞弛豫时
│   └── makefile.include     # 站点专用编译配置，留存以备追溯
├── vasp_pseudopotentials/   # 如 potpaw_PBE_64/
├── venvs/                   # 每个组件一个 Python 环境
├── conda_envs/              # 如使用 conda 前缀
├── mlff_models/             # 机器学习势权重
├── qe-*/ pseudo-*/          # 可选组件
└── work/                    # AutoZT 的 work_dir
```

容量参考：VASP 构建树约 2 GB（源码、对象文件，以及原版与变胞版共六个二进制）；赝势库约 0.3 GB；
单个 Python 环境 0.6–2 GB；机器学习势 0.1–1 GB；可选组件最多约 3 GB。安装前用 `df -h`、`quota` 核对配额。

验收判据：`du -sh ~/software/AutoZT/*` 不超过可用配额。

## 7. 工具链选择

### 7.1 编译器栈

| 站点条件 | 建议选择 |
|---|---|
| 仅有系统 GCC（常为 4.8–7.x） | 加载较新模块（`gcc/9` 或更高）或使用 developer toolset；过旧的编译器可能报内部编译错误 |
| Intel oneAPI / Parallel Studio | 使用 Intel MPI 编译器包装（经典编译器为 `mpiicc`/`mpiifort`；oneAPI 2024 起 `ifort` 已移除，改用 `mpiifx`）并链接 MKL，配置量最小 |
| 混合（GNU 编译器 + Intel MKL） | 可行；显式链接相应 MKL 变体（`gnu_thread`） |

ABI 原则：MPI 库决定了应使用的编译器与 Fortran 运行时。混用不同版本的 `libgfortran` 常导致
链接成功但运行失败的二进制。

### 7.2 数学与 FFT 库

| 来源 | 优点 | 注意事项 |
|---|---|---|
| MKL | 性能好、自包含 | 与 GNU 编译器混用时须链接 `mkl_gf_lp64` 接口层（而非 `mkl_intel_lp64`），线程层选 `mkl_sequential` 或 `mkl_gnu_thread`，具体参数以 Intel MKL Link Line Advisor 为准；Intel 编译器使用 `-qmkl`（旧版为 `-mkl`）；`libiomp5` 位于编译器库目录 |
| 系统 BLAS/LAPACK | 无需安装 | 以 `-l:liblapack.so.3 -l:libblas.so.3` 链接。LAPACK 3.6.0 起，已弃用的 `DGEGV`、`DGGSVD` 默认不再编译；若链接时缺少这些符号，处理方法见第 8.1 节 |
| OpenBLAS | 多线程性能好 | 需注意与 MPI 进程的线程超订 |
| FFTW3 | 可移植 | 使用双精度安装（`.../fftw3/double/...`） |

验收判据：`ldd <二进制> | grep "not found"` 无输出，且冒烟运行不出现未定义符号错误。

## 8. 电子结构引擎的构建

以下流程以 VASP 为例，同样适用于其他以源码分发的引擎。

```bash
module load <编译器> <MPI>             # 作业脚本中同样需要
cd ~/software/AutoZT/vasp
# 1) 解包源码
# 2) 为本站点编写 makefile.include，并与源码树一同保留
# 3) 先编译单个目标，再提高并行度
make -k -j 16 std > build_std.log 2>&1 || true
```

### 8.1 构建失败诊断表

| 构建日志中的现象 | 可能原因 | 处理方式 |
|---|---|---|
| `internal compiler error: Segmentation fault` | 编译器过旧，与源码不兼容 | 加载较新的编译器 |
| `cannot find -lxxx` | 库路径或名称不对 | 用 `ldconfig -p`、`module show` 定位；补 `-L`，并优先加 `-Wl,-rpath` |
| `undefined reference to 'name_'` | 缺库、链接顺序错误或 ABI 不匹配 | 库放在最后；确认 MPI 与 Fortran 运行时和编译器匹配 |
| 并行编译时大量 `.mod` 文件缺失 | 并行构建的模块依赖竞争 | 重复执行 `make -k -j N` 直至收敛；每一轮补齐一批对象 |
| `invalid command 'bdist_wheel'`（Python 包） | 环境中缺少 `wheel` | 先 `pip install wheel`，再重试 |

若链接时报 `undefined reference to 'dgegv_'`（或 `dggsvd_`），说明所用 LAPACK 未编译这些已弃用例程：
从 netlib 取回 `dgegv.f` 与 `dggsvd.f`，加入引擎源码树，并在 `makefile.include` 的
`OBJECTS_LIB` 中加入相应对象。

安装完成后必须动态验证：用 `ldd` 检查未解析库，再运行最小算例（第 8.3 节）。链接成功不足以证明二进制可用。

### 8.2 变胞弛豫（可选）

需要在固定真空轴条件下做面内弛豫的材料，需使用打过补丁的源码另行构建。两个版本并存，
并在站点配置文件中通过两个 VASP profile（`standard`、`relax_2d`）分别引用。

### 8.3 冒烟测试与数值基线

使用 Si 金刚石结构的面心立方（fcc）原胞（2 个原子；常规立方胞含 8 个原子）与最小输入集
（`ENCUT=250`、`ISMEAR=0`、`SIGMA=0.05`、`EDIFF=1e-5`、`NELM=20`），经站点模板提交。验收判据：
退出码为 0；`OUTCAR` 中出现 `aborting loop because EDIFF is reached`（电子步收敛）；写出计时信息。
注：`reached required accuracy` 仅在离子弛豫（`NSW > 0`）收敛时输出，静态自洽计算不会出现该行。

记录总能，并与另一台机器在相同输入下得到的结果比较。若偏差超过 1e-4 eV，应检查赝势文件、
平面波截断能与引擎版本。

### 8.4 并行配置

- `NCORE`（VASP 5 及以后）常取 4–8，或约为每节点核数的平方根；仅在 k 点足够多时使用 `KPAR`。
- 整节点申请效率更高，但在繁忙系统上排队时间更长；两种策略应实测比较（第 11.2 节）。
- 账户并发上限约束“作业规模 × 同时作业数”；应预留一个槽位用于调试。

## 9. 赝势、后处理与 Python 运行时

### 9.1 赝势库布局

赝势有两种使用方式，两者要求的目录布局不同：

- **`defect-dft-cpu`** 直接读取站点配置中的 `potcar_dir`，按元素/变体依次查找：
  `potpaw_PBE_64/<变体>/POTCAR`、`potpaw_PBE_54/<变体>/POTCAR`、`potpaw_PBE/<变体>/POTCAR`、
  `potpaw_PBE/POTCAR_<变体>` 与 `potpaw_PBE.54/POTCAR_<变体>`。
- **使用 VASPKIT 的技能**（例如 `band-dft-cpu`）通过 VASPKIT 生成 `POTCAR`，VASPKIT 读取
  `~/.vaspkit`（见第 9.2 节）；对这些技能而言 `potcar_dir` 不被使用。

若元素目录直接位于 `potcar_dir` 之下（`potcar_dir/<元素>/POTCAR`，即官方压缩包解压后的常见形态），
`defect-dft-cpu` **不会识别**。应把压缩包解压到以库名命名的子目录（例如 `potpaw_PBE_64/`），
并预先执行 `chmod -R u+w`，因为压缩包保留了只读权限。

验收判据：`autozt check --cluster <站点>` 对在用技能不报告赝势数据缺失，且后处理工具生成的 `POTCAR`
首行给出预期的势名。
### 9.2 VASPKIT 配置

VASPKIT 1.x 以空白分隔的 `KEY VALUE` 格式读取 `~/.vaspkit`。写成 `KEY = VALUE` 时其值会被解析为
`=`，赝势将静默地不被生成。可用单个结构验证：

```bash
printf '1\n103\n' | vaspkit      # 应生成首行为 'PAW_PBE Si ...' 的 POTCAR
```

早期版本还有一个诊断问题：`autozt check` 只解析 `KEY = VALUE` 行，因此在 VASPKIT 本身可用时
仍报告 `vaspkit_potcar` 缺失。现已同时兼容两种格式；但 VASPKIT 本身要求空白分隔格式。

### 9.3 Python 运行时环境

| 用途 | 建议环境 |
|---|---|
| AutoZT 核心输入生成 | 含 numpy、scipy、matplotlib、spglib、pymatgen、h5py、seekpath 的 virtualenv |
| 组件专用包（amset、声子相关、MACE 等） | 每个组件一个独立环境；版本要求经常冲突 |

版本策略：

- `numpy`：传统科学计算包（amset 0.5.x、pheasy）一般要求 1.x，而 `cogito-dft` 0.3.13 要求 2.x，二者不可共用环境。
- `pymatgen`：许多较旧的工具与新版主版本不兼容；应通过组件自身入口验证，而不只测试 `import pymatgen`。
- 安装完成后一律以非交互方式验证（第 5.2 节）。

验收判据：`ssh <别名> 'python3 -c "import numpy,scipy,pymatgen,spglib,h5py; print(1)"'`。

#### 9.3.1 以 `conda create -p` 创建的环境

AutoZT 把含 `/` 的运行时环境路径视为 virtualenv，并以 `source <env>/bin/activate` 激活。
新版 conda 创建的环境不含 virtualenv 风格的 `bin/activate`，激活静默失败后会退回系统解释器，
表现为 `autozt check` 报告缺包、生成脚本导入失败。

处理方式：补充一个最小激活脚本：

```bash
# 环境路径在生成 activate 时直接写入：不要用 $0 推算，脚本被 source 时 $0 指向 shell 本身。
ENV=<环境绝对路径>
cat > "$ENV/bin/activate" <<EOF
export VIRTUAL_ENV="$ENV"
export CONDA_PREFIX="$ENV"
export PATH="$ENV/bin:\$PATH"
export LD_LIBRARY_PATH="$ENV/lib:\${LD_LIBRARY_PATH:-}"
EOF
chmod +x "$ENV/bin/activate"
```

验证：`bash -c 'source <env>/bin/activate; python -V; python -c "import <包名>"'` 必须使用该环境的解释器。

## 10. 站点配置参考（`setting/<cluster>.yaml`）

站点配置文件是 AutoZT 在某站点运行方式的唯一事实来源。除标注“必填”外各字段可选；
缺省值依次回退到项目配置与软件包默认值，优先级为：`project_setting/` > 用户配置 > 软件包默认值
（见 `docs/CONFIGURING.md`）。

### 10.1 字段参考

| 字段 | 必填 | 含义 |
|---|---|---|
| `name` | 是 | 集群标识，供 `setting/<name>.yaml`、`-hpc` 与 `check --cluster` 使用 |
| `ssh_host` | 是 | `~/.ssh/config` 中的别名；空字符串表示在本机执行 |
| `template_map` | 否 | 逻辑模板到文件的映射；通常为 `{}`，模板按名称在 `setting/<name>/templates/` 下解析 |
| `vasp.<profile>` | DFT 技能必填 | VASP profile：`standard`（`cell_constraint: none`）与可选的 `relax_2d`（`cell_constraint: ioptcell_tag`），各自给出 `std`、`gam`、`ncl` 的绝对路径 |
| `work_dir` | 是 | 站点可写根目录；不得与其他站点相同 |
| `partition` | 是 | 默认提交分区 |
| `max_cpus` | 建议 | 每节点核数（亦可写作 `max_cpus_per_task`）；预检据此发现模板申请超过单节点核数的情形 |
| `cores` | 否 | 作业默认核数；步骤级与技能级设置优先 |
| `references_dir` | 条件 | `defect-dft-cpu` 的凸包参考相共享目录 |
| `conda_sh` | 否 | `profile.d/conda.sh` 路径；仅用 virtualenv 的站点可留空 |
| `conda_env` | 输入生成必填 | 供输入生成脚本使用的 Python 环境；含 `/` 的路径按 virtualenv 激活 |
| `kl_env` | 随技能 | 基于 phono3py 的声子计算专用环境，会被注入步骤配置；`autozt check` **不**评估该环境，需手工验证（见第 13 节） |
| `amset_env` | 随技能 | AMSET 专用环境（应为独立环境）；`autozt check` 以 `py:<包名>@amset_env` 报告其缺项 |
| `potcar_dir` | 条件 | 赝势库根目录，仅 `defect-dft-cpu` 读取；使用 VASPKIT 的技能读取 `~/.vaspkit`（布局见第 9.1 节） |
| `mlff_model_dir` | MLFF 技能 | 机器学习势权重目录 |
| `remote_path_prefix` | 建议 | 远端命令的 `PATH` 前缀（例如存放工具软链的 `pybin`，其后接环境 `bin`） |
| `qe_bin`、`qe_pseudo_dir`、`qe_setup` | eph-qe 技能 | Quantum ESPRESSO 可执行文件、UPF 目录与运行前置环境语句 |
| 其他字段 | 视报告而定 | `autozt check` 会指出仍缺失的项，按技能需要补充 |

### 10.2 示例配置

```yaml
# (a) 带模块系统的 CPU + SLURM 集群
name: examplecpu
ssh_host: examplecpu
work_dir: /work/<账号>/software/AutoZT/work
partition: cpu64
max_cpus: 64
vasp:
  standard:
    cell_constraint: none
    std: /work/<账号>/software/AutoZT/vasp/vasp.6.4.3-pristine/bin/vasp_std
    gam: .../vasp_gam
    ncl: .../vasp_ncl
conda_sh: ""
conda_env: /work/<账号>/software/AutoZT/venvs/autozt
potcar_dir: /work/<账号>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: /work/<账号>/software/AutoZT/pybin:/work/<账号>/software/AutoZT/venvs/autozt/bin
```

```yaml
# (b) GPU 集群：模板需带 --gres=gpu:N
name: examplegpu
ssh_host: examplegpu
work_dir: /home/<账号>/AutoZT/work
partition: gpu
max_cpus: 128
conda_env: /home/<账号>/AutoZT/venvs/mace_gpu
mlff_model_dir: /home/<账号>/AutoZT/mlff_models
```

```yaml
# (c) 无调度器的服务器：提交由 tools/fakeslurm 模拟
name: exampleserver
ssh_host: exampleserver
work_dir: /home/<账号>/AutoZT/work
conda_env: /home/<账号>/AutoZT/venvs/autozt
```

```yaml
# (d) 本机执行
name: local
ssh_host: ""
work_dir: /tmp/autozt_work
conda_env: ""
```

## 11. 提交模板与资源策略

### 11.1 模板与渲染

模板位于 `setting/<cluster>/templates/`（例如 `submit_std_3d.tpl`、`submit_ncl_2d.tpl`、
`step1_opt/submit_std_2d.tpl`）。生成输入时，AutoZT 按当前 VASP profile 替换可执行文件路径，
但**不解析** `module load` 行与调度器指令。从其他站点复制模板时，至少应核对：分区名、模块名、
MPI 启动器及其参数、时限。

验收判据：执行 `autozt -t <技能> -p <材料> -j <步骤> init`（只生成、不提交）后，检查站点上生成的
`submit.sh`，确认分区、核数与可执行文件路径无误，再执行 `start`。

### 11.2 资源规模

1. 确认分区限额（`sinfo`）与账户并发上限（GrpTRES）。
2. 使“作业规模 × 同时作业数”不超过账户上限；预留一个槽位用于调试。
3. 多节点作业的 `-N`、`-n` 与 `--ntasks-per-node` 必须一致。
4. 时限应低于分区上限，否则提交会被拒绝。

部分节点与整节点的取舍：在繁忙系统上，整节点申请需要等待完全空闲的节点，排队时间更长；
“较多、较小”的作业组合通常总体吞吐更高，但最优策略与站点相关，应实测确定（第 13 节）。

### 11.3 MPI 启动

| 现象 | 可能原因 | 处理方式 |
|---|---|---|
| 作业在数秒内失败；stderr 出现 `binding more processes than cpus` | OpenMPI 默认绑核与节点或 cgroup 布局不兼容 | `mpirun --bind-to none -np $SLURM_NTASKS ...`，先在单节点确认 |
| 以 OpenMPI 的 `mpirun` 启动 Intel MPI 程序 | ABI 不匹配 | 使用 Intel MPI 的 `mpiexec`/`mpirun`，或加载对应的 `mpi/intelmpi` 模块 |
| 程序被启动了多份，或输出互相覆盖 | 启动器与程序所链接的 MPI 实现不是同一套，导致每个进程各自启动一份完整实例 | 使用与程序所链接 MPI 实现一致的 `mpirun`/`mpiexec` 启动 |
| 单节点通信异常 | InfiniBand/UCX 初始化问题 | 单节点作业可设 `OMPI_MCA_btl=self,vader,tcp` 与 `OMPI_MCA_pml=ob1` |

### 11.4 加速卡

- NVIDIA GPU：以 `#SBATCH --gres=gpu:N` 申请，安装 CUDA 版运行时，并在作业内用 `nvidia-smi` 确认可见性。
- 海光 DCU：以 `#SBATCH --gres=dcu:N` 申请，使用 DTK 工具链（`compiler/dtk/...`），并安装与 DTK 匹配的
  PyTorch 构建，而不是上游 CUDA wheel。
- 验收判据：作业内能列出加速卡，且运行时报告其可用（`torch.cuda.is_available()` 或厂商等价接口）。

### 11.5 无 SLURM 的站点

AutoZT 以 `sbatch`、`squeue`、`scancel` 的语义管理作业。对无调度器的服务器，仓库提供
`tools/fakeslurm`（本地模式下自动启用，远端节点可通过 `remote_path_prefix` 注入）。
对 PBS/SGE/LSF，需要一层把 `sbatch` 语义转换为本地调度器的适配层；在适配层可用之前，
不支持自动推进作业，这一限制应记入站点记录。

## 12. 多集群运行与项目级覆盖

### 12.1 步骤级路由

同一材料的不同步骤可在不同站点运行，在项目配置（`tf_<项目>.yaml`）中声明：

```yaml
step_hpc:
  S7_deform: examplegpu
  S7.1_read: examplegpu
```

只有列出的步骤被路由，其余步骤仍使用材料默认站点。下游步骤所需的产物由 AutoZT 按依赖关系传输。
对材料整体切换站点（`autozt -p <材料> hpc <集群>`）会影响所有步骤，可能需要重跑前置步骤，
该决定应记入站点记录。

### 12.2 项目级步骤配置

参数可在项目级覆盖，优先级由低到高为：`project_setting/templates/step.conf`（共用）与
`project_setting/templates/<步骤>/step.conf`（步骤专用）。修改单个参数可执行：

```bash
autozt -t <技能> -p <材料> -j <步骤> conf --set params.<KEY>=<VALUE>   # 写入项目配置
autozt -t <技能> -p <材料> -j <步骤> conf                              # 只读，显示合并后的值
```

写入项目配置在 AutoZT 中属于受控操作，须由材料的负责人提出。

### 12.3 操作纪律

作业提交命令彼此之间不具备幂等性。`retry`、`init`、`rerun` 会取消同一材料排队中的作业
（`kill_if_queued`）。若有两位操作者或两个会话同时作用于同一材料，后执行的一方会取消先提交的作业，
结果将无法归因。为此应遵循：

1. 同一材料同一时间只由一位操作者处理；执行前先查看任务图（`autozt graph`）。
2. 并行工作树之间隔离站点配置。`setting/*` 不纳入版本控制，因此每个 git worktree 拥有各自的站点配置副本。
3. 在站点或项目记录中注明该次运行的操作者。

## 13. 验收判据

下列各项全部通过，即认为站点可用；证据应粘贴进站点记录。

| # | 检查项 | 命令 | 期望结果 |
|---|---|---|---|
| 1 | 连通性 | `ssh -o BatchMode=yes <别名> true` | 退出码 0 |
| 2 | 非交互 Python | `ssh <别名> 'python3 -V'` | 输出预期版本 |
| 3 | 配置预检 | `autozt doctor` | 无错误，且列出新站点 |
| 4 | 技能就绪 | `autozt check --cluster <站点>` | DFT 相关技能通过 |
| 5 | 模板渲染 | `autozt -p <材料> -j <步骤> init` 后检查 `submit.sh` | 分区、核数与路径正确 |
| 6 | 冒烟计算 | 最小 Si 算例（第 8.3 节） | 电子步收敛（`aborting loop because EDIFF is reached`）；能量与基线一致 |
| 7 | 并发 | 提交两个探测作业 | 均被接受；无 QOS 或分区拒绝 |
| 8 | 清理 | `squeue -u <账号>` | 无测试作业残留 |
| 9 | 记录 | `docs/<SITE>-DEPLOY-<日期>.md` | 配额、问题与验收数字齐备 |

`autozt check` 输出的解读：

- **硬性缺失**：技能无法运行（缺少引擎二进制、赝势库、核心 Python 包，或该技能专用环境）。
- **可选项缺失**：次要功能不可用（`optional exe:`/`optional py:`、无 GPU 等），不影响主流程。
- 缺包信息是**针对该技能实际使用的环境**报告的：AMSET 依赖在 `amset_env` 中检查，并带 `@amset_env` 后缀；
  其余依赖在 `conda_env` 中检查（站点配置未设 `conda_env` 时，使用技能 `skill.yaml` 中 `requires.conda` 声明的环境）。
  `kl_env` 当前不在 `autozt check` 的检查范围内，需手工验证。不带 `--cluster` 时，`autozt check` 检查的是本机。

## 14. 可选组件

组件按需安装；每小节给出安装途径、验收测试与已知问题。

### 14.1 Quantum ESPRESSO（技能 `eph-qe-cpu`）

1. 优先使用站点模块（`module avail | grep -i quantum-espresso`）：将 `qe_bin` 指向其 `bin` 目录，
   并在 `qe_setup` 中表达运行环境（模块加载，或显式设置 MKL 与 Intel MPI 的 `LD_LIBRARY_PATH`）。
2. 否则使用 conda-forge 包（自包含；需要网络访问，并为 conda 求解器准备足够内存）。
3. 验收：`pw.x -h` 与 `ph.x -h` 可运行；小体系 Si 自洽计算达到 `convergence has been achieved`；
   `qe_pseudo_dir` 中存在模守恒赝势。
4. 已知问题：预先加载的 MPI 模块可能与 QE 模块冲突（先卸载前者）；MKL/Intel MPI 运行库有时需显式导出。

### 14.2 AMSET（技能 `ke-dft-cpu`）

应使用独立环境：AMSET 0.5.x 通常要求 numpy 1.x 与较旧的 pymatgen，与机器学习环境冲突。基于 pip 的安装示例：

```bash
python -m venv <env>
<env>/bin/pip install 'amset==0.5.1' 'pymatgen<2025' 'numpy<2' scipy matplotlib h5py numba interpolation sumo pyFFTW
```

已知问题：`BoltzTraP2` 从其 PyPI 源码分发构建时常失败（自带 spglib 子模块缺少 CMakeLists，且需要 netCDF4 头文件）。
可将 conda-forge 的二进制包解包后把模块复制进环境；输运计算本身并不需要其插值功能。
验收：`import amset`、`amset deform create --help` 正常，并在站点配置文件中登记 `amset_env`。

### 14.3 MACE（MLFF 技能）

先安装 CPU 或 CUDA 版 PyTorch，再安装 `mace-torch`、`e3nn`、`ase` 与 `matscipy`。注意 PyTorch 与 e3nn 的兼容性：
PyTorch 2.6 将 `torch.load` 的默认值改为 `weights_only=True`，会使旧版 e3nn 报 `UnpicklingError`；
应使用 PyTorch 2.4.x 或更新版 e3nn。技能期望的默认模型名见 `skill/*/templates/step.conf`（`MACE_MODEL`），
将对应权重下载到 `mlff_model_dir`，并用单点计算验证。

### 14.4 COGITO（技能 `cohp-cogito`）

将 `cogito-dft`（例如 0.3.13）安装到独立 virtualenv（其要求 numpy 2.x）。
把 `COGITO`、`COGITOanalyze`、`COGITOpost` 暴露给该技能使用的环境（软链即可）。
前置条件：开启 `LWAVE=.TRUE.` 的静态 VASP 计算、对称性约化的 k 点网格（`ISYM` 不为 `-1`）、足够的能带数。

### 14.5 UniHamGNN（技能 `unihamgnn`）

该组件包括：HamGNN 源码树、通用模型权重（约 1 GB）、OpenMX 基组库 `DFT_DATA19`、
辅助可执行文件 `openmx_postprocess` 与 `read_openmx`，以及含 PyTorch（CPU）、`torch_geometric`、
`lightning`、`e3nn` 的 Python 环境。

构建 `openmx_postprocess` 需要 GSL、MPI、LAPACK 与 FFTW。缺少 ScaLAPACK/BLACS 符号时，
可用 MKL 的 ScaLAPACK 与 BLACS 库；`libiomp5` 位于 Intel 编译器库目录；静态链接 GSL 可避免运行时依赖 `libcblas`。
ELPA 的 Fortran 模块存在并行依赖竞争，需要多轮 `make -k -j N`。
技能参数（`HAMGNN_DIR`、`UNI_MODEL`、`DFT_DATA`、`OPENMX_POSTPROCESS`、`READ_OPENMX`）属于项目级步骤参数，
应写入 `project_setting/templates/step.conf`，而不是站点配置文件。

### 14.6 声子与力常数拟合（phonopy、phono3py、ALM、pheasy）

- `phonopy`、`phono3py`：安装二进制 wheel（`pip install --only-binary=:all: phonopy phono3py`）；
  确认 `phono3py` 可执行文件存在且两个模块均可导入。
- `ALM`：PyPI 上同名包是 2.4 KB 的空包，不可安装。conda-forge 只提供较旧的构建（2.0.0_dev），且仅覆盖部分 Python 版本。
  当无法把匹配的构建安装进目标环境时，可解包对应 Python 版本的 conda-forge 包，把其 `site-packages/alm/` 复制进环境，
  并将环境中 `spglib` 自带的 `libsymspg.so.2` 软链为所需的 `libsymspg.so.1`
  （扩展按 `$ORIGIN/../../..` 解析库）。
- `pheasy`：从源码安装（先 `pip install wheel`，再 `pip install --no-deps --no-build-isolation .`）；
  确认 `pheasy` 可执行文件可用。

## 15. 故障排查表

| 现象 | 可能原因 | 处理方式 |
|---|---|---|
| `Invalid account or account/partition combination` | 账号未被授权使用该分区 | 用 `--test-only` 找出可用分区；申请权限 |
| `CPU count per node can not be satisfied` | 申请核数超过单节点核数（如模板硬编码 192） | 调整模板或步骤配置中的每节点核数；使 `max_cpus` 一致 |
| `group max tres limit ...` | 超过账户并发上限 | 减少同时作业数或作业规模；记录该变更 |
| 作业数秒内失败并报 `binding more processes than cpus` | OpenMPI 默认绑核 | `mpirun --bind-to none` |
| `error while loading shared libraries` | 运行时缺少库 | 用 `ldd` 定位；以 `-Wl,-rpath` 重新构建，或加载对应模块 |
| `undefined symbol: dgegv_` | LAPACK 3.8 及以后移除了该例程 | 从 netlib 取 `dgegv.f`、`dggsvd.f` 加入源码 |
| 启动时报 `I REFUSE TO CONTINUE` | 结构不一致或 POSCAR 损坏 | 核对传输的 POSCAR 与本地文件一致（AutoZT 会做校验） |
| VASPKIT 不生成 POTCAR | `~/.vaspkit` 使用了 `KEY = VALUE` | 改为 `KEY VALUE`；核对 `PBE_PATH` 根目录 |
| `autozt check` 报告缺包但实际已安装 | 非交互 `PATH`；技能专用环境；conda 前缀缺少 `bin/activate` | 见第 5.2、13、9.3.1 节 |
| `conda create` 被终止（exit 137）或极慢 | 求解器内存或索引访问 | 改用 micromamba 或 `--solver=libmamba`；直接解包 conda 包；使用本地镜像 |
| `sbatch: request too frequently` | 提交频率受限 | 每次提交间隔一至两秒，分批提交 |
| 远端采集报 `Connection closed` 或超时 | 中间跳、代理或 ControlMaster 失效 | 重建连接（`ssh -O exit` 后重连）；检查各跳密钥有效期 |
| 两个会话提交的作业相互取消 | 两个会话作用于同一材料（`kill_if_queued`） | 同一材料只保留一位操作者；用 worktree 隔离 `setting/`；操作前查看任务图 |

## 16. 站点记录与可复现性

每个站点应创建 `docs/<SITE>-DEPLOY-<YYYYMMDD>.md`（由 `docs/DEPLOY-RECORD-TEMPLATE.md` 复制）。
记录应包含：实测资源上限、软件版本与编译参数、验收输出、尚存的限制。该文档是使站点可被其他团队复现的依据。

可复现性要求：

1. 固定引擎、库、Python 包与模型权重的版本。
2. 站点配置不进入版本控制（`setting/*` 已被忽略）；分发模板与占位符。
3. 引用软件版本：见 `README.md` 与 `CITATION.cff` 中的版本化 DOI。
4. 报告结果时，同时给出站点名称与站点记录文档，以及输入文件。

---

## 附录 A：最小站点配置模板

```yaml
name: <SITE>
ssh_host: <ssh别名>
template_map: {}
vasp:
  standard:
    cell_constraint: none
    std: <根目录>/software/AutoZT/vasp/<版本>/bin/vasp_std
    gam: <根目录>/software/AutoZT/vasp/<版本>/bin/vasp_gam
    ncl: <根目录>/software/AutoZT/vasp/<版本>/bin/vasp_ncl
work_dir: <根目录>/software/AutoZT/work
partition: <分区>
max_cpus: <每节点核数>
conda_sh: ""
conda_env: <根目录>/software/AutoZT/venvs/autozt
potcar_dir: <根目录>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: <根目录>/software/AutoZT/pybin:<根目录>/software/AutoZT/venvs/autozt/bin
# ke-dft-cpu 需追加：      amset_env: ...
# eph-qe-cpu 需追加：      qe_bin / qe_pseudo_dir / qe_setup
# MLFF 技能需追加：        mlff_model_dir: ...
```

## 附录 B：一页式流程

1. 勘察站点：`tools/probe_cluster.sh <别名>`；记录第 4.1 节的六个量。
2. 配置 SSH 与非交互运行时（第 5 节）。
3. 创建 `~/software/AutoZT/{vasp,vasp_pseudopotentials,venvs,work}`（第 6 节）。
4. 构建引擎（原版，必要时加变胞版），以 `ldd` 验证（第 7、8 节）。
5. 安装赝势库、配置 VASPKIT、创建 Python 环境（第 9 节）。
6. 编写站点配置与提交模板（第 10、11 节）。
7. 依次运行 `autozt doctor`、`autozt check --cluster <站点>` 与冒烟计算（第 13 节）。
8. 建立站点记录（第 16 节）。
9. 以一条真实材料短链运行，测量吞吐与排队情况。
10. 按需安装可选组件，每装一项重新运行 `autozt check`（第 14 节）。

## 附录 C：实测示例（数值参考）

以下数值测自 SCNet 西安 `clus_xa`（2026-10-19），详细记录见 `docs/CLUS_XA-DEPLOY-20261019.md`。
这些数值仅作为量级参考，不应照搬为配置。

| 项目 | 实测值 | 影响 |
|---|---|---|
| 节点规格 | 64 核、255 GB | 决定 `max_cpus` 与模板规模 |
| 分区时限 | 3 天 | `--time=3-00:00:00` |
| 账户并发 | 320 核（单作业同样可达） | 例如 32 核作业 × 10 并发 |
| 提交模板 | 每作业 32 核 | 整节点作业需等待完全空闲节点 |
| MPI 问题 | 部分节点 OpenMPI 默认绑核失败 | `mpirun --bind-to none` |
| 编译器 | GCC 9.3.0（系统 7.3.1 会报内部编译错误） | 通过模块加载 |
| MPI 库 | HPC-X 2.11 | 与编译器匹配 |
| BLAS/LAPACK | 系统 3.4.2，并补入 `dgegv`/`dggsvd` 源码 | 新版 LAPACK 缺少这些例程 |
| 引擎 | VASP 6.4.3，原版与 IOPTCELL 补丁版 | `standard` 与 `relax_2d` 两个 profile |
| 赝势 | potpaw_PBE_64，346 个元素 | 解压后需调整权限 |
| VASPKIT | 1.2.4，`~/.vaspkit` 使用空白分隔 | 用于生成 POTCAR |
| Quantum ESPRESSO | 站点模块 6.8（Intel MPI）、MKL 运行时、`Si_ONCV_PBE-1.2.upf` | 自洽 6 步收敛 |
| AMSET | 0.5.1（独立环境，pymatgen 2024.11、numpy 1.26） | `ke-dft-cpu` 就绪 |
| COGITO | cogito-dft 0.3.13（独立环境，numpy 2.x） | `cohp-cogito` 就绪 |
| MACE | mace-torch 0.3.16、PyTorch 2.4.1 CPU、e3nn 0.4.4；`MACE-matpes-pbe-omat-ft.model` | `mlff`/`opt-mlff-cpu` 就绪 |
| UniHamGNN | HamGNN 2.3.0、885 MB 模型、`DFT_DATA19`、自编译 OpenMX 工具 | `unihamgnn` 就绪 |

## 附录 D：参考文献

- ASE，安装定制（`siteconfig.py`）：https://ase-lib.org/ase/installation.html
- AiiDA，计算机配置与测试：https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html
- Quantum ESPRESSO 用户手册（安装、`make.inc`、`arch/`）：https://www.quantum-espresso.org/Doc/user_guide/
- CP2K 手册（`arch/` 编译）：https://manual.cp2k.org/
- LAMMPS 安装与机器 Makefile：https://docs.lammps.org/Install.html
- Snakemake 执行器插件与 profile：https://snakemake.github.io/snakemake-plugin-catalog/
- Nextflow 配置与 profiles：https://www.nextflow.io/docs/latest/config.html
- jobflow-remote 配置：https://matgenix.github.io/jobflow-remote/
- AutoZT 仓库与版本化 DOI：见 `README.md` 与 `CITATION.cff`。
