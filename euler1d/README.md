# euler1d

一维 Euler 方程经典 Riemann 算例：Sod、Lax、Shu-Osher。
Classical one-dimensional Euler Riemann problems: Sod, Lax, and Shu-Osher.

---

## 中文

### 内容

| 程序 | 算例 | 域 | 终止时间 | 初值 |
|------|------|----|----------|------|
| `euler1d` | `sod` | `[0,1]` | 0.15 | `(1,0,1)` / `(0.125,0,0.1)` |
| `euler1d` | `lax` | `[0,1]` | 0.15 | `(0.445,0.698,3.528)` / `(0.5,0,0.571)` |
| `euler1d` | `shu_osher` | `[0,10]` | 1.8 | 左 `(3.857143,2.629369,10.3333)`，右 `ρ=1+0.2 sin(5x)` |
| `euler1d_exact` | `sod`、`lax` | `[0,1]` | 0.15 | 精确 Riemann 解 |

`euler1d_exact` 只支持有闭式解的 `sod` 与 `lax`。三个算例共用同一套
Steger-Warming 分裂、特征变量插值、RK4 与零阶外推边界。

### 控制方程

```
∂U/∂t + ∂E/∂x = 0
U = (ρ, ρu, ρE)^T
E = (ρu, ρu² + p, (ρE + p)u)^T
p = (γ - 1)(ρE - 0.5 ρu²),  γ = 1.4
```

网格为单元中心 `x_j = j h`，`h = L/N`，默认 `N = 200`。

### 空间离散

半节点通量由特征变量插值给出：

1. 用中心点状态构造 3×3 左右特征矩阵；
2. 对三个特征场分别调用
   `../interpolation/src/interpolation/` 的 `scheme_dispatch`；
3. 反投影回守恒变量；
4. 用 Steger-Warming 分裂得到 `E_{j+1/2}`。

节点右端项 `K_j = -dE/dx` 由差分套子给出：

| 名称 | 阶数 | 说明 |
|------|------|------|
| `WCNSE2` | 2 | 二阶通量差分 |
| `WCNSE4` | 4 | 四阶通量差分 |
| `WCNSE6` | 6 | 六阶通量差分 |
| `HCS4` | 4 | 一对半点 + 一对整点 |
| `HCS6` | 6 | 一对半点 + 两对整点 |

`difference_operator = 'AUTO'` 时，三阶插值配 `HCS4`，五阶插值配 `HCS6`。
默认 `scheme = 'QR'`、`difference_operator = 'AUTO'`，即 `QR + HCS4`。

### 时间推进与边界

RK4，`dt = CFL h / max_j(|u_j| + sqrt(p_j/ρ_j))`，默认 `CFL = 0.1`。
外边界采用 `halo = 7` 的零阶外推。每个 RK 阶段后执行正性保护：
`ρ < 1e-8` 或 `p < 1e-8` 的单元用旧状态邻点平均替换。

### 结果矩阵

`scripts/run_results.sh` 默认跑三组代表格式：

| 组 | 插值 | 差分 |
|----|------|------|
| 二阶 | `MUSCL` | `WCNSE2` |
| 三阶（默认） | `QR` | `HCS4` |
| 五阶 | `Z5` | `HCS6` |

`MOC53` 不参与结果矩阵。结果只做数值解与精确解/参考解的曲线对比，
不计算 L1/L2/L∞。

- `results/csv/<case>_<scheme>_<operator>.csv`：`x,rho,u,p`；
- `results/csv/exact_<case>.csv`：Sod/Lax 精确解；
- `results/csv/reference_shu_osher.csv`：Shu-Osher 的 `N=1600` 数值参考解；
- `results/logs/`：网格、格式、CFL、步数、保护计数、极值与 `nan_count`；
- `results/figures/`：每格式三面板 `ρ,u,p` 图，以及每算例的密度对比图。

### 构建与运行

```sh
make                 # 构建 euler1d 与 euler1d_exact（默认 dp）
make test            # 单元测试
make smoke           # 三算例小网格端到端检查
make results         # 生成结果矩阵与图
./bin/dp/euler1d     # 运行 input.nml 中的配置
```

默认编译器 `FC=ifx`，默认精度 `dp`；`make FC=gfortran` 与
`make PRECISION=qp` 可切换。构建目录按精度隔离：`build/$(PRECISION)`、
`bin/$(PRECISION)`。切换编译器前先 `make clean`。

### 配置（`input.nml`）

```
&euler1d_nml
  case_1d = 'sod'          ! sod | lax | shu_osher
  mesh_num = 200
  cfl = 0.1d0
  tend = -1.0d0            ! <=0 用算例默认值
  output_dir = 'results'
  write_sample = .false.
/
&scheme_nml
  scheme = 'QR'
/
&difference_nml
  difference_operator = 'AUTO'
/
```

`euler1d_exact` 复用同一组 `&euler1d_nml`，只接受 `sod` 与 `lax`，
忽略 `cfl` 与 `write_sample`。

### 精度

默认 `dp`；`qp` 可构建与测试。提交的结果由 `dp` 生成。

### 文件与参考

- 方法说明：`docs/method_zh.md`、`docs/method_en.md`；
- 与旧程序的差异：`docs/formula_audit.md`；
- 参考程序为作者开发算例时使用的版本，未随本仓库发布；文档与代码中的公式
  均可独立核对。

---

## English

### Contents

| Program | Case | Domain | Final time | Initial state |
|---------|------|--------|------------|---------------|
| `euler1d` | `sod` | `[0,1]` | 0.15 | `(1,0,1)` / `(0.125,0,0.1)` |
| `euler1d` | `lax` | `[0,1]` | 0.15 | `(0.445,0.698,3.528)` / `(0.5,0,0.571)` |
| `euler1d` | `shu_osher` | `[0,10]` | 1.8 | left `(3.857143,2.629369,10.3333)`, right `ρ=1+0.2 sin(5x)` |
| `euler1d_exact` | `sod`, `lax` | `[0,1]` | 0.15 | exact Riemann solution |

`euler1d_exact` supports only the closed-form Sod and Lax cases. All three
cases share the same Steger-Warming split, characteristic interpolation, RK4,
and zero-order extrapolation boundary.

### Governing equations

```
∂U/∂t + ∂E/∂x = 0
U = (ρ, ρu, ρE)^T
E = (ρu, ρu² + p, (ρE + p)u)^T
p = (γ - 1)(ρE - 0.5 ρu²),  γ = 1.4
```

The grid is cell-centred: `x_j = j h`, `h = L/N`, with `N = 200` by default.

### Spatial discretization

Half-node fluxes come from characteristic interpolation:

1. build the 3×3 left/right characteristic matrices from the centre state;
2. interpolate each characteristic field with
   `../interpolation/src/interpolation/` through `scheme_dispatch`;
3. project back to conserved variables;
4. apply the Steger-Warming split to obtain `E_{j+1/2}`.

The nodal residual `K_j = -dE/dx` uses one of `WCNSE2`, `WCNSE4`, `WCNSE6`,
`HCS4`, or `HCS6`. With `difference_operator = 'AUTO'`, third-order
interpolation pairs with `HCS4` and fifth-order with `HCS6`. The default is
`scheme = 'QR'` and `difference_operator = 'AUTO'`, i.e. `QR + HCS4`.

### Time advance and boundary

RK4 with `dt = CFL h / max_j(|u_j| + sqrt(p_j/ρ_j))` and `CFL = 0.1` by
default. The halo is `halo = 7` with zero-order extrapolation. After every
RK stage, cells with `ρ < 1e-8` or `p < 1e-8` are replaced by the average of
their old neighbours.

### Result matrix

The default run uses three representative pairs: `MUSCL + WCNSE2`,
`QR + HCS4`, and `Z5 + HCS6`. `MOC53` is not part of the result matrix.
Only numerical-versus-exact/reference curves are produced; no L1/L2/Linf
norms are computed.

### Build and run

```sh
make                 # build both drivers (dp by default)
make test            # unit tests
make smoke           # small end-to-end checks for all three cases
make results         # generate the result matrix and figures
./bin/dp/euler1d     # run the configuration in input.nml
```

The default compiler is `ifx` and the default precision is `dp`;
`make FC=gfortran` and `make PRECISION=qp` switch them. Run `make clean`
before switching compilers.

### Precision

The committed results use `dp`. `qp` is supported for builds and tests.
