# WCNS-scheme-tests

WCNS 数值格式测试程序集合 / Numerical scheme test programs for WCNS research

本仓库包含六个自包含的子项目，共用一套可复现的 Fortran 构建与结果脚本。
This repository contains six self-contained subprojects with a shared, reproducible
Fortran build and result workflow.

---

## 中文

### 子项目

| 目录 | 内容 | 文档 |
|------|------|------|
| `interpolation/` | WCNS 非线性插值精度测试：收敛阶、非线性误差、权/光滑因子诊断 | [`interpolation/README.md`](interpolation/README.md) |
| `adr/` | ADR 频谱测试：`adr_weno`（WENO 重构）与 `adr_wcns`（WCNS 插值 + 可选差分套子） | [`adr/README.md`](adr/README.md) |
| `complex_wave/` | 复杂波时间推进测试：一维周期对流，WCNS 插值 + 可选差分套子，按周期输出波形 | [`complex_wave/README.md`](complex_wave/README.md) |
| `accuracy/` | 格式精度测试：一维周期对流、二维 Euler（等熵涡 / 正弦密度波），输出误差与收敛阶 | [`accuracy/README.md`](accuracy/README.md) |
| `euler1d/` | 一维 Euler 经典 Riemann 算例：Sod、Lax、Shu-Osher，附精确解程序与曲线对比 | [`euler1d/README.md`](euler1d/README.md) |
| `euler2d/` | 二维 Euler 七个经典算例：双马赫反射、Kelvin–Helmholtz、Rayleigh–Taylor、两个 Riemann 问题、激波–气泡、激波–涡，含 Tecplot 场输出 | [`euler2d/README.md`](euler2d/README.md) |

### 构建与测试

```bash
make                 # 构建六个子项目
make test            # 运行六个子项目的测试
make all-precisions  # qp 与 dp 各构建一遍
make FC=gfortran     # 可选的 gfortran 构建
make clean
```

### 精度选择

`interpolation/` 与 `adr/` 默认 `PRECISION=qp`（四精度）；`complex_wave/`、
`accuracy/`、`euler1d/` 与 `euler2d/` 默认 `PRECISION=dp`（双精度）。六个子项目都支持 `qp` 与 `dp`：

```bash
make PRECISION=dp
make -C interpolation test PRECISION=dp
make -C adr test PRECISION=dp
make -C complex_wave test PRECISION=qp
```

构建目录按精度隔离：`build/$(PRECISION)`、`bin/$(PRECISION)`。`interpolation/` 与
`adr/` 的提交结果由 qp 生成；`complex_wave/`、`accuracy/`、`euler1d/` 与 `euler2d/` 的提交结果由 dp 生成。
qp/dp 都可用于快速验证与对照。

### 结果

- `interpolation/results/`：收敛阶 CSV、误差曲线与诊断图。
- `adr/results/`：ADR 色散/耗散 CSV 与图。
- `complex_wave/results/`：每周期波形快照、日志与图。
- `accuracy/results/`：一维/二维误差与收敛阶 CSV、运行日志。
- `euler1d/results/`：Sod/Lax/Shu-Osher 的数值解、精确/参考解 CSV、日志与对比图。
- `euler2d/results/`：七个二维算例的 Tecplot 场数据、运行日志与提交的 smoke 样本。

### 引用与许可

各子项目的公式来源见其 README 与 `docs/`。本仓库使用 MIT License，见 `LICENSE`。
联系方式：haomin.li@buaa.edu.cn。

---

## English

### Subprojects

| Directory | Content | Documentation |
|-----------|---------|---------------|
| `interpolation/` | WCNS nonlinear interpolation accuracy test: convergence order, nonlinear error, weight/smoothness diagnostics | [`interpolation/README.md`](interpolation/README.md) |
| `adr/` | ADR spectral test: `adr_weno` (WENO reconstruction) and `adr_wcns` (WCNS interpolation + selectable difference operator) | [`adr/README.md`](adr/README.md) |
| `complex_wave/` | Complex-wave time-advancement test: 1D periodic advection with WCNS interpolation, selectable difference operators, per-period waveforms | [`complex_wave/README.md`](complex_wave/README.md) |
| `accuracy/` | Scheme accuracy tests: 1D periodic convection and 2D Euler (isentropic vortex / sine density wave), with error and convergence-order tables | [`accuracy/README.md`](accuracy/README.md) |
| `euler1d/` | Classical 1D Euler Riemann problems: Sod, Lax, and Shu-Osher, with an exact-solution driver and profile comparisons | [`euler1d/README.md`](euler1d/README.md) |
| `euler2d/` | Seven classical 2D Euler cases: double Mach reflection, Kelvin–Helmholtz, Rayleigh–Taylor, two Riemann problems, shock-bubble, and shock-vortex, with Tecplot field output | [`euler2d/README.md`](euler2d/README.md) |

### Build and test

```bash
make                 # build all six subprojects
make test            # run all six test suites
make all-precisions  # build qp and dp
make FC=gfortran     # optional gfortran build
make clean
```

### Precision

`interpolation/` and `adr/` default to `PRECISION=qp` (quad); `complex_wave/`,
`accuracy/`, `euler1d/`, and `euler2d/` default to `PRECISION=dp` (double). All six subprojects support both:

```bash
make PRECISION=dp
make -C interpolation test PRECISION=dp
make -C adr test PRECISION=dp
make -C complex_wave test PRECISION=qp
```

Build directories are separated by precision: `build/$(PRECISION)`, `bin/$(PRECISION)`.
Committed results come from qp in `interpolation/` and `adr/`, and from dp in
`complex_wave/`, `accuracy/`, `euler1d/`, and `euler2d/`. Both precisions remain available for quick
checks and comparison.

### Results

- `interpolation/results/`: convergence CSVs, error curves and diagnostics.
- `adr/results/`: ADR dispersion/dissipation CSVs and figures.
- `complex_wave/results/`: per-period waveform snapshots, logs and figures.
- `accuracy/results/`: 1D/2D error and convergence-order CSVs, run logs.
- `euler1d/results/`: Sod/Lax/Shu-Osher numerical, exact/reference CSVs, logs and comparison figures.
- `euler2d/results/`: Tecplot fields, run logs, and committed smoke samples for the seven 2D cases.

### Citation and license

Formula sources are documented in each subproject's README and `docs/`. This repository
uses the MIT License; see `LICENSE`. Contact: haomin.li@buaa.edu.cn.
