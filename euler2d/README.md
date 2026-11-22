# euler2d

二维 Euler 方程七个经典算例的 WCNS 格式测试程序。
Seven classical two-dimensional Euler test cases for the WCNS schemes.

---

## 中文

### 内容

`euler2d` 把学位论文中的七个二维 Euler 算例统一到一个驱动程序中：
`double_mach`、`kelvin_helmholtz`、`rayleigh_taylor`、`riemann1`、
`riemann2`、`shock_bubble`、`shock_vortex`。算例只提供初值、边界、
源项和默认网格；空间插值格式与差分套子由 `input.nml` 在运行时选择。

默认配置是 `scheme = 'QR'`、`difference_operator = 'AUTO'`，
即 **`QR + AUTO(HCS4)`**。`AUTO` 按插值阶数配对：三阶 `QR` 配 `HCS4`。

### 算例表

域写成 `[x_min,x_max] × [z_min,z_max]`，网格比 `nx/nz` 由算例固定。
“论文配置”是 `scripts/run_results.sh` 为复现论文结果写入的插值与套子；
它们不是 `input.nml` 的默认值。驱动程序的默认始终是 `QR + AUTO(HCS4)`。

| 算例 | 域 | `nx/nz` | 默认 `nz` | 默认 `tend` | 默认 `gamma` | 论文配置（插值 + 套子） | 特殊项 |
|------|----|---------|-----------|-------------|--------------|--------------------------|--------|
| `double_mach` | `[0,4] × [0,1]` | 4 | 200 | 0.2 | 1.4 | `MOZ + HCS6` | 顶部激波边界随 `t` 移动 |
| `kelvin_helmholtz` | `[0,1] × [0,1]` | 1 | 400 | 0.5 | 1.4 | `Z5 + HCS6` | 两层密度、正弦 `w` 扰动 |
| `rayleigh_taylor` | `[0,0.25] × [0,1]` | 1/4 | 400 | 1.95 | 5/3 | `MOC53 + HCS6` | 重力源项；`nx = nz/4` |
| `riemann1` | `[0,1] × [0,1]` | 1 | 400 | 0.8 | 1.4 | `MOC53 + HCS6` | 四象限初值 |
| `riemann2` | `[0,1] × [0,1]` | 1 | 1600 | 0.3 | 1.4 | `MOC33 + HCS4` | 四象限初值 |
| `shock_bubble` | `[-0.1,1.9] × [-0.5,0.5]` | 2 | 800 | 0.4 | 1.4 | `MOC33 + HCS4` | 双输出时刻 `tend2 = 0.2` |
| `shock_vortex` | `[-0.25,1.75] × [-0.5,0.5]` | 2 | 180 | 0.6 | 1.4 | `Z3 + HCS4` | 静止涡与运动激波 |

六个算例只写一个输出时刻；`shock_bubble` 在 `tend2` 和 `tend` 各写一次。
坐标原点按初值定义：`shock_bubble` 为 `(-0.1,-0.5)`，`shock_vortex`
为 `(-0.25,-0.5)`，其余为 `(0,0)`。`shock_vortex` 的旧导出程序曾把
坐标写成 `(0,0)`，与初值不一致；本实现统一使用初值原点。

### 数值方法概要

控制方程为二维 Euler 方程组。网格是 `nx+1 × nz+1` 的单元中心网格，
`h = 1/nz`，物理域为 `[0,nx/nz] × [0,1]` 加算例原点平移。数组 halo
为 **`mt = 7`** 层：八点特征插值模板在 `j = -3` 要读 `u(-6)`，在
`j = n+2` 要读 `u(n+6)`，因此 `mt = 5` 会越界。

空间离散的流程是：

1. 对每个半节点，从八点模板投影到特征变量；
2. 调用 `../interpolation/` 的统一插值接口，把四个特征场插值到半节点；
3. 反投影回守恒变量，得到左右状态的半节点值；
4. 用修正 Steger–Warming 分裂得到半节点通量 `Eflux`、`Gflux`；
5. 用所选差分套子作用于半节点通量和整数节点通量，得到 `-dE/dx - dG/dz`。

时间推进固定为四步四阶 Runge–Kutta（RK4），
`dt = CFL·h / max(|u| + c)`，`c = sqrt(gama*p/rho)`（`gama` 为代码中的
比热比变量）；
末步截断到目标时刻。两层修正是：

- 特征插值后，若半点密度或压力低于模板最小值的 0.1%，回退到中心点；
- 每个 RK 子步后，若单元 `rho < 1e-8` 或 `p < 1e-8`，用旧状态四邻点平均替换。

完整公式、特征矩阵、边界表与差分系数见
[`docs/method_zh.md`](docs/method_zh.md) 和
[`docs/method_en.md`](docs/method_en.md)。

### 插值格式、套子与 AUTO 配对

`&scheme_nml` 的 `scheme` 复用 `interpolation/` 的格式名，可选：
`MUSCL`、`JS3`、`Z3`、`JS5`、`Z5`、`QR`、`MOC322`、`MOC33`、`MOC53`、
`MOZ`、`LINEAR3_M1P1`、`LINEAR3_M2P0`、`LINEAR3_0P2`、`LINEAR4_M2P1`、
`LINEAR4_M1P2`、`LINEAR5_M2P2`、`LINEAR6_M2P3`、`RBF_GA3`、`RBF_MQ3`、
`RBF_IQ3`、`RBF_GA5`、`RBF_MQ5`、`RBF_IQ5`。

`&difference_nml` 的 `difference_operator` 可选 `AUTO`、`WCNSE2`、
`WCNSE4`、`WCNSE6`、`HCS4`、`HCS6`。仓库的 `WCNSE2/4/6` 就是论文
语境中的 **E2/E4/E6**。`AUTO` 的配对规则是：

| 插值阶数 | `AUTO` 套子 | 典型插值 |
|----------|-------------|----------|
| 2 | `WCNSE2` | `MUSCL` |
| 3 | `HCS4` | `JS3`、`Z3`、`QR`、`MOC33`、`MOC322` |
| 4 | `WCNSE4` | `LINEAR4_M2P1`、`LINEAR4_M1P2` |
| 5 | `HCS6` | `JS5`、`Z5`、`MOC53`、`MOZ`、`LINEAR5_M2P2` |
| 6 | `WCNSE6` | `LINEAR6_M2P3` |

显式选择一个与插值阶数不匹配的套子时，程序打印
`warning: ... is not the default pairing` 但继续运行，便于做“同插值不同套子”
的对照。套子是运行级选择，不按算例锁死；RT 的重力源项在五个套子下都生效。

论文配置与本仓库标准件的差异（`MOZ` 参数名、`MOC53` 混合权重、
`MOC33` 的 THINC 非线性候选等）逐项记录在
[`docs/formula_audit.md`](docs/formula_audit.md)。

### 构建与运行

```sh
# 在 euler2d/ 目录中
make                 # 构建 bin/dp/euler2d 与 test_control
make test            # 单元测试 + 文档检查
make smoke           # 七个算例的小网格端到端检查
./bin/dp/euler2d     # 读取当前目录的 input.nml

# 从仓库根目录
make -C euler2d
make -C euler2d test
```

切换精度或编译器前先清理：

```sh
make -C euler2d clean
make -C euler2d test PRECISION=qp
make -C euler2d clean
make -C euler2d test FC=gfortran
```

`PRECISION` 可取 `dp`（默认）或 `qp`；`FC` 默认 `ifx`，也可用 `gfortran`。
构建目录按精度隔离：`build/$(PRECISION)`、`bin/$(PRECISION)`。

结果矩阵脚本（默认 `cfl = 0.1`）：

```sh
sh scripts/run_results.sh --smoke   # nz=16，写 results/smoke/
sh scripts/run_results.sh --quick   # nz=32，写 results/quick/
sh scripts/run_results.sh --full    # 算例默认 nz 与 tend，写 results/full/
python3 scripts/plot_fields.py results/smoke/<case>_<nz>x<nx>_t1.dat <case>.png
```

只提交 `results/smoke/` 的样本；`results/quick/` 与 `results/full/` 由
`.gitignore` 忽略。

### 配置（`input.nml`）

```fortran
&euler2d_nml
  case_2d = 'riemann1'      ! 七个算例之一
  mesh_num = 0              ! z 方向单元数；0 = 算例默认；>= 8
  cfl = 0.3d0               ! 必须为正
  tend = -1.0d0             ! <= 0 = 算例默认
  tend2 = -1.0d0            ! 仅 shock_bubble：<0 默认，0 关闭，>0 第二输出
  gamma = 0.0d0             ! <= 0 = 算例默认
  output_dir = 'results'    ! 输出目录，程序自动创建
  write_sample = .false.    ! true 时另写 t0 初值样本
/
&scheme_nml
  scheme = 'QR'             ! 插值格式；其余参数见 interpolation/README.md
/
&difference_nml
  difference_operator = 'AUTO'   ! AUTO | WCNSE2 | WCNSE4 | WCNSE6 | HCS4 | HCS6
/
```

`nx` 由 `mesh_num` 和算例比例推出。`rayleigh_taylor` 的比例为 1/4，
因此 `mesh_num` 必须能被 4 整除。`tend2` 只允许用于 `shock_bubble`，
且必须小于 `tend`。

输出为 Tecplot ASCII 点格式，变量 `x, z, rho`，文件名
`<case>_<nz>x<nx>_t1.dat`（第二个输出为 `_t2.dat`，可选 `_t0.dat`）。
日志 `<case>.log` 是键值对，包含 `case`、`scheme`、`operator`、`mesh`、
`cfl`、`tend`、`gama`、`steps`、`interp_corrections`、`cell_corrections`、
`nan_count`、`cpu_seconds`。

### 结果与插图

`results/smoke/` 是提交到仓库的小网格样本：每个算例 `nz = 16`、
`tend = 0.001`，`shock_bubble` 另写 `tend2 = 0.0005`，全部
`nan_count=0`。它们只用于展示输出格式与完整管线，不是论文结果。

`docs/figures/` 中的插图**复制自作者的学位论文 LaTeX 源**
（`figures/chap4/`），保留论文原图网格与格式，不是本仓库重新生成的：

| 算例 | 论文插图（PDF） | 原图网格 |
|------|------------------|----------|
| `double_mach` | [`Double_mach_QR_1920_480.pdf`](docs/figures/Double_mach_QR_1920_480.pdf) | 1920×480 |
| `double_mach`（对照） | [`Double_mach_Z3_1920_480.pdf`](docs/figures/Double_mach_Z3_1920_480.pdf) | 1920×480 |
| `kelvin_helmholtz` | [`KH_800_QR.pdf`](docs/figures/KH_800_QR.pdf) | 800×800 |
| `rayleigh_taylor` | [`RT_240_960_QR.pdf`](docs/figures/RT_240_960_QR.pdf) | 240×960 |
| `riemann1` | [`RM_type1_800_QR.pdf`](docs/figures/RM_type1_800_QR.pdf) | 800×800 |
| `riemann2` | [`RM_Type2_1600_QR.pdf`](docs/figures/RM_Type2_1600_QR.pdf) | 1600×1600 |
| `shock_bubble`（`t=0.2`） | [`Shock_bubble_1600_800_QR_t02.pdf`](docs/figures/Shock_bubble_1600_800_QR_t02.pdf) | 1600×800 |
| `shock_bubble`（`t=0.4`） | [`Shock_bubble_1600_800_QR_t04.pdf`](docs/figures/Shock_bubble_1600_800_QR_t04.pdf) | 1600×800 |
| `shock_vortex` | [`SV_360_180_all.pdf`](docs/figures/SV_360_180_all.pdf) | 360×180 |

每个 PDF 旁边有同名 `.png` 预览，由
`pdftoppm -png -r 150 -singlefile` 生成。论文插图使用论文自己的网格
（例如 DM 的 `1920×480` 与本仓库算例表的 `nz = 200` 不同），因此不要把
它们当作本仓库 `run_results.sh` 在默认 `nz` 下的输出。`shock_vortex`
在论文中没有单独的 `Z3` 图，只有 `SV_360_180_all.pdf`。

学位论文 PDF 本身不放入仓库。插图用于说明算例形态与格式效果；
定量、逐值的新旧对比口径见 `docs/formula_audit.md`。

### 精度与编译器

默认 `dp`（双精度），`qp`（四精度）可构建并通过同一组测试。默认编译器
`ifx`；`gfortran` 也可用。提交的 smoke 样本由 `dp` 生成。若要做
`-fcheck=all -ffpe-trap=invalid,zero,overflow` 边界检查，先 `make clean`，
再指定 `FC=gfortran FFLAGS="-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"`。

### 文件与参考

- 方法说明：[`docs/method_zh.md`](docs/method_zh.md)、[`docs/method_en.md`](docs/method_en.md)；
- 与论文程序的差异：[`docs/formula_audit.md`](docs/formula_audit.md)；
- 插值格式实现：`../interpolation/`；
- 参考程序为作者开发算例时使用的版本，未随本仓库发布；文档与代码中的
  公式均可独立核对。

---

## English

### Contents

`euler2d` collects seven two-dimensional Euler test cases from the degree
thesis into one driver: `double_mach`, `kelvin_helmholtz`,
`rayleigh_taylor`, `riemann1`, `riemann2`, `shock_bubble`, and
`shock_vortex`. Each case supplies only its initial state, boundary rule,
source term, and default mesh; the interpolation scheme and difference
operator are selected at run time through `input.nml`.

The default configuration is `scheme = 'QR'` with
`difference_operator = 'AUTO'`, i.e. **`QR + AUTO(HCS4)`**. `AUTO` pairs a
third-order interpolation such as `QR` with `HCS4`.

### Case table

Domains are written `[x_min,x_max] × [z_min,z_max]`; the grid ratio
`nx/nz` is fixed by the case. The “paper configuration” column is what
`scripts/run_results.sh` writes to reproduce the thesis runs; it is not the
`input.nml` default. The driver default is always `QR + AUTO(HCS4)`.

| Case | Domain | `nx/nz` | Default `nz` | Default `tend` | Default `gamma` | Paper config (scheme + operator) | Special item |
|------|--------|---------|--------------|----------------|-----------------|----------------------------------|--------------|
| `double_mach` | `[0,4] × [0,1]` | 4 | 200 | 0.2 | 1.4 | `MOZ + HCS6` | top shock boundary moves with `t` |
| `kelvin_helmholtz` | `[0,1] × [0,1]` | 1 | 400 | 0.5 | 1.4 | `Z5 + HCS6` | two density layers, sinusoidal `w` |
| `rayleigh_taylor` | `[0,0.25] × [0,1]` | 1/4 | 400 | 1.95 | 5/3 | `MOC53 + HCS6` | gravity source; `nx = nz/4` |
| `riemann1` | `[0,1] × [0,1]` | 1 | 400 | 0.8 | 1.4 | `MOC53 + HCS6` | four-quadrant initial state |
| `riemann2` | `[0,1] × [0,1]` | 1 | 1600 | 0.3 | 1.4 | `MOC33 + HCS4` | four-quadrant initial state |
| `shock_bubble` | `[-0.1,1.9] × [-0.5,0.5]` | 2 | 800 | 0.4 | 1.4 | `MOC33 + HCS4` | two output times, `tend2 = 0.2` |
| `shock_vortex` | `[-0.25,1.75] × [-0.5,0.5]` | 2 | 180 | 0.6 | 1.4 | `Z3 + HCS4` | stationary vortex, moving shock |

Six cases write one output time; `shock_bubble` writes once at `tend2` and
once at `tend`. Coordinate origins follow the initial states:
`shock_bubble` uses `(-0.1,-0.5)`, `shock_vortex` uses `(-0.25,-0.5)`, and
the others use `(0,0)`. The old `shock_vortex` export routine wrote
coordinates from `(0,0)`, inconsistent with its initial state; this
implementation uses the initial-state origin everywhere.

### Numerical method at a glance

The governing equations are the two-dimensional Euler equations. The grid
is cell-centred with `nx+1 × nz+1` points, spacing `h = 1/nz`, and physical
domain `[0,nx/nz] × [0,1]` shifted by the case origin. The halo is
**`mt = 7`** cells: the eight-point characteristic stencil at `j = -3`
reads `u(-6)`, and the stencil at `j = n+2` reads `u(n+6)`, so `mt = 5`
would read out of bounds.

The spatial pipeline is:

1. project the eight-point stencil to characteristic variables at each
   half node;
2. interpolate the four characteristic fields to the half node through the
   shared `../interpolation/` dispatcher;
3. project back to conserved variables to obtain left and right states;
4. apply the modified Steger–Warming split to build the half-node fluxes
   `Eflux` and `Gflux`;
5. apply the selected difference operator to the half-node fluxes and the
   integer-node fluxes to form `-dE/dx - dG/dz`.

Time integration is fixed four-stage, fourth-order Runge–Kutta (RK4) with
`dt = CFL·h / max(|u| + c)` and `c = sqrt(gama*p/rho)` (the code variable
`gama` holds the ratio of specific heats); the last step is
clipped to the target time. Two correction layers are active:

- after characteristic interpolation, a half-node density or pressure below
  0.1% of the stencil minimum falls back to the centre state;
- after every RK substage, an interior cell with `rho < 1e-8` or
  `p < 1e-8` is replaced by the average of its four neighbours from the
  saved old state.

Full formulas, characteristic matrices, the boundary table, and the
difference coefficients are in [`docs/method_zh.md`](docs/method_zh.md) and
[`docs/method_en.md`](docs/method_en.md).

### Interpolation schemes, operators, and AUTO pairing

`scheme` in `&scheme_nml` reuses the names from `interpolation/`:
`MUSCL`, `JS3`, `Z3`, `JS5`, `Z5`, `QR`, `MOC322`, `MOC33`, `MOC53`,
`MOZ`, `LINEAR3_M1P1`, `LINEAR3_M2P0`, `LINEAR3_0P2`, `LINEAR4_M2P1`,
`LINEAR4_M1P2`, `LINEAR5_M2P2`, `LINEAR6_M2P3`, `RBF_GA3`, `RBF_MQ3`,
`RBF_IQ3`, `RBF_GA5`, `RBF_MQ5`, and `RBF_IQ5`.

`difference_operator` in `&difference_nml` accepts `AUTO`, `WCNSE2`,
`WCNSE4`, `WCNSE6`, `HCS4`, and `HCS6`. The repository names
`WCNSE2/4/6` are the thesis-context **E2/E4/E6**. `AUTO` pairs:

| Interpolation order | `AUTO` operator | Typical schemes |
|---------------------|-----------------|-----------------|
| 2 | `WCNSE2` | `MUSCL` |
| 3 | `HCS4` | `JS3`, `Z3`, `QR`, `MOC33`, `MOC322` |
| 4 | `WCNSE4` | `LINEAR4_M2P1`, `LINEAR4_M1P2` |
| 5 | `HCS6` | `JS5`, `Z5`, `MOC53`, `MOZ`, `LINEAR5_M2P2` |
| 6 | `WCNSE6` | `LINEAR6_M2P3` |

An explicit operator that does not match the interpolation order prints
`warning: ... is not the default pairing` and continues, so “same
interpolation, different operator” studies remain possible. The operator
is a run-level choice rather than a per-case lock; the RT gravity source is
active with all five operators.

Differences between the thesis configurations and the repository standard
schemes (`MOZ` parameter names, `MOC53` mixing weights, the `MOC33` THINC
nonlinear candidate) are recorded item by item in
[`docs/formula_audit.md`](docs/formula_audit.md).

### Build and run

```sh
# inside euler2d/
make                 # build bin/dp/euler2d and test_control
make test            # unit tests + documentation checks
make smoke           # end-to-end check of the seven cases on a small mesh
./bin/dp/euler2d     # reads input.nml from the current directory

# from the repository root
make -C euler2d
make -C euler2d test
```

Clean before changing precision or compiler:

```sh
make -C euler2d clean
make -C euler2d test PRECISION=qp
make -C euler2d clean
make -C euler2d test FC=gfortran
```

`PRECISION` is `dp` (default) or `qp`; `FC` defaults to `ifx` and also
supports `gfortran`. Build directories are separated by precision:
`build/$(PRECISION)`, `bin/$(PRECISION)`.

Result-matrix scripts (they use `cfl = 0.1`):

```sh
sh scripts/run_results.sh --smoke   # nz=16, writes results/smoke/
sh scripts/run_results.sh --quick   # nz=32, writes results/quick/
sh scripts/run_results.sh --full    # case-default nz and tend, results/full/
python3 scripts/plot_fields.py results/smoke/<case>_<nz>x<nx>_t1.dat <case>.png
```

Only `results/smoke/` samples are committed; `results/quick/` and
`results/full/` are ignored by `.gitignore`.

### Configuration (`input.nml`)

```fortran
&euler2d_nml
  case_2d = 'riemann1'      ! one of the seven cases
  mesh_num = 0              ! z cells; 0 = case default; >= 8
  cfl = 0.3d0               ! must be positive
  tend = -1.0d0             ! <= 0 = case default
  tend2 = -1.0d0            ! shock_bubble only: <0 default, 0 off, >0 second output
  gamma = 0.0d0             ! <= 0 = case default
  output_dir = 'results'    ! created automatically
  write_sample = .false.    ! true also writes a t0 initial-state sample
/
&scheme_nml
  scheme = 'QR'             ! interpolation; other parameters in interpolation/README.md
/
&difference_nml
  difference_operator = 'AUTO'   ! AUTO | WCNSE2 | WCNSE4 | WCNSE6 | HCS4 | HCS6
/
```

`nx` follows from `mesh_num` and the case ratio. `rayleigh_taylor` has ratio
1/4, so `mesh_num` must be divisible by 4. `tend2` is accepted only for
`shock_bubble` and must be smaller than `tend`.

Output is Tecplot ASCII point format with variables `x, z, rho`, in files
named `<case>_<nz>x<nx>_t1.dat` (`_t2.dat` for the second output, and an
optional `_t0.dat`). The log `<case>.log` uses key=value lines and records
`case`, `scheme`, `operator`, `mesh`, `cfl`, `tend`, `gama`, `steps`,
`interp_corrections`, `cell_corrections`, `nan_count`, and `cpu_seconds`.

### Results and figures

`results/smoke/` holds the committed small-mesh samples: `nz = 16` and
`tend = 0.001` for every case, with `tend2 = 0.0005` for `shock_bubble`,
and `nan_count=0` throughout. They demonstrate the output format and the
full pipeline; they are not thesis results.

The figures in `docs/figures/` are **copied from the author's degree-thesis
LaTeX source** (`figures/chap4/`) and keep the thesis grids and formats.
They are not regenerated by this repository:

| Case | Thesis figure (PDF) | Original grid |
|------|---------------------|---------------|
| `double_mach` | [`Double_mach_QR_1920_480.pdf`](docs/figures/Double_mach_QR_1920_480.pdf) | 1920×480 |
| `double_mach` (comparison) | [`Double_mach_Z3_1920_480.pdf`](docs/figures/Double_mach_Z3_1920_480.pdf) | 1920×480 |
| `kelvin_helmholtz` | [`KH_800_QR.pdf`](docs/figures/KH_800_QR.pdf) | 800×800 |
| `rayleigh_taylor` | [`RT_240_960_QR.pdf`](docs/figures/RT_240_960_QR.pdf) | 240×960 |
| `riemann1` | [`RM_type1_800_QR.pdf`](docs/figures/RM_type1_800_QR.pdf) | 800×800 |
| `riemann2` | [`RM_Type2_1600_QR.pdf`](docs/figures/RM_Type2_1600_QR.pdf) | 1600×1600 |
| `shock_bubble` (`t=0.2`) | [`Shock_bubble_1600_800_QR_t02.pdf`](docs/figures/Shock_bubble_1600_800_QR_t02.pdf) | 1600×800 |
| `shock_bubble` (`t=0.4`) | [`Shock_bubble_1600_800_QR_t04.pdf`](docs/figures/Shock_bubble_1600_800_QR_t04.pdf) | 1600×800 |
| `shock_vortex` | [`SV_360_180_all.pdf`](docs/figures/SV_360_180_all.pdf) | 360×180 |

Each PDF has a same-named `.png` preview generated by
`pdftoppm -png -r 150 -singlefile`. The thesis figures use the thesis grids
(for example, DM uses `1920×480`, whereas the case table default is
`nz = 200`), so they are not the output of `run_results.sh` at the default
`nz`. The thesis has no separate `Z3` figure for `shock_vortex`; only
`SV_360_180_all.pdf` is available.

The thesis PDF itself is not placed in the repository. The figures document
the case morphology and scheme behaviour; the quantitative old/new
comparison protocol is in `docs/formula_audit.md`.

### Precision and compilers

The default precision is `dp` (double); `qp` (quad) builds and passes the
same tests. The default compiler is `ifx`; `gfortran` is also supported.
The committed smoke samples were produced with `dp`. For a bounds and
floating-point check, clean first and then build with
`FC=gfortran FFLAGS="-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"`.

### Files and references

- Method notes: [`docs/method_zh.md`](docs/method_zh.md),
  [`docs/method_en.md`](docs/method_en.md);
- Differences from the thesis programs: [`docs/formula_audit.md`](docs/formula_audit.md);
- Interpolation implementations: `../interpolation/`;
- The reference programs are the versions used by the author when the cases
  were developed and are not distributed with this repository. Every formula
  in the documentation and code can be checked independently.
