# accuracy

格式精度测试：一维周期对流与二维 Euler（等熵涡 / 正弦密度波）。
Scheme accuracy tests: 1D periodic convection and 2D Euler (isentropic vortex /
sine density wave).

---

## 中文

### 内容

| 程序 | 算例 | 方程 |
|------|------|------|
| `accuracy_1d` | 一维周期对流 | `u_t + u_x = 0`，`x ∈ [-1,1]`，`Tend = 2` |
| `accuracy_2d` | `vortex` | 二维 Euler，`[0,10]²`，等熵涡，`Tend = 1` |
| `accuracy_2d` | `sine` | 二维 Euler，`[0,2π]²`，`ρ = ρ0 + A sin(x+y)`，`Tend = 1` |

两个程序都用仓库的 `../interpolation/src/interpolation/` 插值标准件，并在
`src/wcns/difference_operators.f90` 中提供可选差分套子。

### 目录

```
accuracy/
├── app/                 两个驱动
├── src/common/          配置与误差/阶数模块
├── src/convection1d/    一维初值、精确解、RHS/RK4
├── src/euler2d/         二维场、两套算例、通量、求解器
├── src/wcns/            差分套子（与 adr/complex_wave 同源）
├── tests/               单元测试、smoke、文档检查
├── scripts/             结果矩阵脚本
├── results/             CSV 与日志
└── docs/                方法与公式说明
```

### 构建与运行

```sh
make                 # 构建 accuracy_1d 与 accuracy_2d
make test            # 单元测试 + smoke + 文档检查
make smoke           # 只跑小网格端到端检查
make results         # 跑完整结果矩阵，写入 results/
sh scripts/run_results.sh --quick   # 小网格结果，写入 results/quick/
```

默认编译器 `FC=ifx`，默认精度 `dp`；`make FC=gfortran` 与
`make PRECISION=qp` 可切换。构建目录按精度隔离：`build/$(PRECISION)`、
`bin/$(PRECISION)`。切换编译器前先 `make clean`。

### 配置（`input.nml`）

```
&accuracy_nml
  write_sample = .false.        ! 是否输出一个 Tecplot 样本
  output_dir   = 'results'
/
&accuracy_1d_nml
  initial_1d = 'sin2'           ! sin2|sin|sin3|sin_phase|gauss|sqrt|rational
  nbase = 10
  levels = 8
  cfl = 0.1d0
  tend = 2.0d0
/
&accuracy_2d_nml
  case_2d = 'vortex'            ! 'vortex' | 'sine'
  nbase = 40
  levels = 4
  cfl = 0.1d0
  tend = 1.0d0
  sine_rho0 = 1.0d0
  sine_amp  = 0.2d0
/
&scheme_nml
  scheme = 'Z5'
/
&difference_nml
  difference_operator = 'AUTO'
/
```

`difference_operator` 支持 `AUTO`、`WCNSE2`、`WCNSE4`、`WCNSE6`、`HCS4`、
`HCS6`。`AUTO` 按阶数配对：三阶→`HCS4`，五阶→`HCS6`。显式选择与阶数不匹配时
只打印警告，仍按用户选择运行，便于做“同插值不同套子”的对照。

### 一维初值

所有初值都是 `[-1,1]` 上的连续光滑周期函数，精确解为 `u0(x - t)`：

| ID | 公式 |
|----|------|
| `sin2` | `sin²(πx)` |
| `sin` | `sin(πx)` |
| `sin3` | `sin³(πx)` |
| `sin_phase` | `sin(πx - sin(πx)/π)` |
| `gauss` | `exp(-20x²)` |
| `sqrt` | `sqrt(1 + 20x²)` |
| `rational` | `1/(1 + 20x²)` |

`INITIALS="sin2 gauss" ./scripts/run_results.sh` 可只跑指定初值。

### 方法矩阵

| 算例 | 插值 | 默认套子 |
|------|------|----------|
| 1D 对流 / 2D `vortex` / 2D `sine` | `Z3` | `HCS4` |
| 同上 | `QR` | `HCS4` |
| 同上 | `Z5` | `HCS6` |

`Z5` 在光滑正弦算例上显示约 5 阶；`Z3`/`QR` 在临界点附近按已知规律降阶
（L1 约 2、L∞ 约 1–2），与 `interpolation/` 子项目的 Z3/QR 结果一致。
二维强等熵涡上 `Z5` 出现非线性权重降阶，L1 阶数从约 3 逐步升到 `n=320`
的约 4.3，与旧 `js5.dat` 的 3–4 阶趋势一致；二维管线用
`LINEAR5_M2P2 + HCS6` 对旧 `linear5.dat` 逐值验证，相对差约 `1e-15`。
边界处理：一维 `halo = 7`、二维 `mt = 7`，覆盖插值与套子的最远读取；
详见 `docs/method_zh.md`。

### 结果

`results/csv/`：每个算例/初值/格式/套子一个 CSV，列为
`n,L1,order_L1,L2,order_L2,Linf,order_Linf,cpu_time`。
`results/logs/`：`scheme`、`operator`、`precision`、网格层级、CFL、Tend 等
键值对。一维文件名包含 `initial`，二维包含 `case`。

### 精度

默认 `dp`；`qp` 可构建。二维 `qp` 内存约为 `dp` 的四倍，完整矩阵建议用 `dp`。
一维最细层若在 `dp` 下触底，可用 `make results PRECISION=qp` 重跑。

---

## English

### Contents

| Program | Case | Equation |
|---------|------|----------|
| `accuracy_1d` | 1D periodic convection | `u_t + u_x = 0`, `x ∈ [-1,1]`, `Tend = 2` |
| `accuracy_2d` | `vortex` | 2D Euler, `[0,10]²`, isentropic vortex, `Tend = 1` |
| `accuracy_2d` | `sine` | 2D Euler, `[0,2π]²`, `ρ = ρ0 + A sin(x+y)`, `Tend = 1` |

Both programs use the shared interpolation standards in
`../interpolation/src/interpolation/`. Difference operators live in
`src/wcns/difference_operators.f90`.

### Build and run

```sh
make
make test
make smoke
make results
sh scripts/run_results.sh --quick
```

The default compiler is `ifx` and the default precision is `dp`;
`make FC=gfortran` and `make PRECISION=qp` switch them. Builds are isolated by
precision in `build/$(PRECISION)` and `bin/$(PRECISION)`. Run `make clean`
before switching compilers.

### Configuration

The five namelist groups are shown above. `difference_operator` accepts
`AUTO`, `WCNSE2`, `WCNSE4`, `WCNSE6`, `HCS4`, `HCS6`. `AUTO` pairs third-order
interpolation with `HCS4` and fifth-order interpolation with `HCS6`; an
explicit mismatch prints a warning but still runs.

### 1D initial profiles

All seven profiles (`sin2`, `sin`, `sin3`, `sin_phase`, `gauss`, `sqrt`,
`rational`) are continuous and periodic on `[-1,1]`; the exact solution is
`u0(x - t)`.

### Method matrix

| Case | Interpolation | Default operator |
|------|---------------|------------------|
| 1D / 2D `vortex` / 2D `sine` | `Z3` | `HCS4` |
| same | `QR` | `HCS4` |
| same | `Z5` | `HCS6` |

`Z5` reaches about fifth order on the smooth sine case. `Z3`/`QR` show the known
critical-point order reduction (L1 about 2, Linf about 1–2), consistent with the
`interpolation/` Z3/QR results. On the strong isentropic vortex, `Z5` shows
the nonlinear-weight order reduction: the L1 order rises from about 3 to about
4.3 at `n=320`, consistent with the old `js5.dat` (3–4). The 2D pipeline itself
is validated by `LINEAR5_M2P2 + HCS6` against the old `linear5.dat` to about
`1e-15` relative. The 1D halo is 7 and the 2D `mt` is 7, covering the farthest
interpolation/operator read; see `docs/method_en.md`.

### Results

`results/csv/` holds one CSV per case/profile/scheme/operator with columns
`n,L1,order_L1,L2,order_L2,Linf,order_Linf,cpu_time`; `results/logs/` holds
key=value run metadata. 1D filenames include the profile and 2D filenames
include the case.

### Precision

`dp` is the default; `qp` is supported. A full 2D `qp` matrix needs roughly
four times the memory of `dp`, so the committed results use `dp`.
