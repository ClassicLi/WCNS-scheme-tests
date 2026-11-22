# 复杂波时间推进测试：WCNS 插值与可选差分套子
# Complex-Wave Time-Advancement Test: WCNS Interpolation with Selectable Difference Operators

本子项目求解一维周期线性对流方程，用复杂初值检验 WCNS 格式在长时间推进下对
光滑峰、方波、三角波和半圆的综合表现。差分套子可选，推进周期数可自定义，
每个周期输出一次波形快照并与初值对比。

This subproject solves the one-dimensional periodic linear advection equation with a
complex initial profile containing a smooth bump, a square block, a triangle, and a
semicircle. The difference operator is selectable, the number of periods is
user-defined, and one waveform snapshot is written per period for comparison with the
initial profile.

---

## 中文

### 控制方程与网格

```
u_t + u_x = 0,   x ∈ [-1, 1],   周期边界,   波速 1
x_j = -1 + j h,  h = 2 / N,     j = 0..N
```

一个波形周期为 `T = 2`。`u(N)` 与 `u(0)` 相同，halo 由周期映射填充。

### 初值

初值由四段组成：

- `x ≈ -0.7` 附近的平滑高斯峰；
- `[-0.4, -0.2]` 上高度为 1 的方波；
- 顶点在 `x = 0.1` 的三角波；
- `x = 0.5` 附近的平滑半圆。

高斯峰与半圆使用三点模板 `(f(x-δ) + 4 f(x) + f(x+δ)) / 6` 平滑，`δ = 0.005`
（`N = 200` 时为 `h/2`），与参考程序一致。

### 空间离散

半节点通量 `E_{j+1/2}` 由 WCNS 插值给出，插值格式直接复用
`../interpolation/src/interpolation/` 中的 23 个标准件。节点空间右端项
`K_j = -dE/dx` 由下列差分套子之一给出：

| 名称 | 公式 | 阶数 |
|------|------|------|
| `WCNSE2` | `(E_{j+1/2} - E_{j-1/2}) / h` | 2 |
| `WCNSE4` | `[(9/8)(E_{j+1/2} - E_{j-1/2}) - (1/24)(E_{j+3/2} - E_{j-3/2})] / h` | 4 |
| `WCNSE6` | `[(75/64)(E_{j+1/2} - E_{j-1/2}) - (25/384)(E_{j+3/2} - E_{j-3/2}) + (3/640)(E_{j+5/2} - E_{j-5/2})] / h` | 6 |
| `HCS4` | `[(4/3)(E_{j+1/2} - E_{j-1/2}) - (1/6)(u_{j+1} - u_{j-1})] / h` | 4 |
| `HCS6` | `[(64/45)(E_{j+1/2} - E_{j-1/2}) - (2/9)(u_{j+1} - u_{j-1}) + (1/180)(u_{j+2} - u_{j-2})] / h` | 6 |

所有公式前的整体负号已经包含在 `K_j = -dE/dx` 中。`HCS6` 是“一对半点 +
两对整点”的混合紧致套子，系数由 Taylor 展开匹配到六阶得到。

### 时间推进与周期数

时间推进为四步四阶 Runge–Kutta（RK4）。时间步长

```
dt = CFL * h / max_j |u_j|
```

默认 `CFL = 0.3`。`periods`（周期数）在 `input.nml` 中设置，默认 4。程序在
`t = 0, T, 2T, ..., periods*T` 各输出一次快照；最后一步会截断到整周期，
保证快照恰好落在周期边界上。本子项目只比较波形，不计算误差范数。

### 默认结果

结果由 `scripts/run_results.sh` 生成，网格 `N = 200`，推进 4 个周期：

| 组 | 格式 | 差分套子 |
|----|------|----------|
| 五阶 | `JS5`、`Z5` | `HCS6` |
| 三阶 | `JS3`、`Z3`、`QR` | `HCS4` |

每个算例输出：

- `results/csv/<scheme>_<operator>_t<k>.csv`：`x, u`，`k = 0..periods`；
- `results/dat/<scheme>_<operator>.dat`：多 zone Tecplot 文件；
- `results/logs/<scheme>_<operator>.log`：配置、步数、极值与 `nan_count`；
- `results/figures/*.png`：每格式的演化图，以及三阶/五阶的四周期对比图。

### 构建、运行与测试

```bash
make                      # 默认双精度 dp
make FC=gfortran          # 可选的 gfortran 构建
make PRECISION=qp         # 四精度构建
make test                 # 单元测试 + 参数/文档检查
make smoke                # 端到端烟测
make results              # 生成 5 个算例与图
./bin/dp/complex_wave     # 运行 input.nml 中的配置
```

切换编译器前请先 `make clean`。

### 精度

本子项目默认使用双精度 `dp`；`qp` 仍可构建和测试，但提交的结果由 `dp` 生成。
构建目录按精度隔离：`build/$(PRECISION)`、`bin/$(PRECISION)`。

### 文件与参考

- 设计说明：`docs/method_zh.md`、`docs/method_en.md`；
- 与参考程序的差异：`docs/formula_audit.md`；
- 参考程序为作者开发本章算例时使用的版本，未随本仓库发布；本文档与代码中的公式
  均可独立核对。

---

## English

### Governing equation and grid

```
u_t + u_x = 0,   x ∈ [-1, 1],   periodic boundary,   wave speed 1
x_j = -1 + j h,  h = 2 / N,     j = 0..N
```

One wave period is `T = 2`. The endpoint satisfies `u(N) = u(0)`, and halo points
are filled by periodic mapping.

### Initial condition

The profile has four features: a smoothed Gaussian bump near `x = -0.7`, a unit
square block on `[-0.4, -0.2]`, a triangle with apex at `x = 0.1`, and a smoothed
semicircle near `x = 0.5`. The Gaussian and semicircle use the three-point stencil
`(f(x-δ) + 4 f(x) + f(x+δ)) / 6` with `δ = 0.005`, matching the reference program.

### Spatial discretization

Half-node fluxes `E_{j+1/2}` come from the WCNS interpolation standards shared with
`../interpolation/src/interpolation/`. The spatial right-hand side `K_j = -dE/dx`
uses one of five difference operators: `WCNSE2`, `WCNSE4`, `WCNSE6`, `HCS4`, or
`HCS6` (formulas in the Chinese table above). `HCS6` combines one pair of half
nodes with two pairs of integer nodes and is sixth-order accurate.

### Time advance and period count

RK4 is used with

```
dt = CFL * h / max_j |u_j|
```

and default `CFL = 0.3`. The number of periods (`periods`, default 4) is set in
`input.nml`. Snapshots are written at `t = 0, T, 2T, ..., periods*T`; the last
step is clipped to land exactly on each period boundary. No error norms are
computed in this subproject.

### Default result set

`scripts/run_results.sh` uses `N = 200` and four periods:

| Group | Schemes | Difference operator |
|-------|---------|---------------------|
| Fifth order | `JS5`, `Z5` | `HCS6` |
| Third order | `JS3`, `Z3`, `QR` | `HCS4` |

Each run writes per-period CSV snapshots, a multi-zone Tecplot file, a log with
`nan_count`, and waveform figures.

### Build, run, and test

```bash
make
make FC=gfortran
make PRECISION=qp
make test
make smoke
make results
./bin/dp/complex_wave
```

Run `make clean` before switching compilers.

### Precision

This subproject defaults to double precision `dp`; `qp` remains buildable and
testable, but committed results are generated in `dp`. Build directories are
separated by precision: `build/$(PRECISION)` and `bin/$(PRECISION)`.
