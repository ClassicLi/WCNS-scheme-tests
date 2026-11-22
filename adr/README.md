# ADR 频谱测试：WENO 重构与 WCNS 格式 / ADR Spectral Test: WENO Reconstruction and WCNS Schemes

本目录提供两套 ADR（Approximate Dispersion Relation，近似色散关系）测试程序：
`adr_weno` 与 `adr_wcns`。两套程序共用同一套周期边界、DFT、时间推进与输出骨架，
只替换空间格式。

This directory provides two ADR (Approximate Dispersion Relation) test programs:
`adr_weno` and `adr_wcns`. They share the same periodic boundary, DFT, time stepping
and output harness; only the spatial scheme differs.

---

## 中文

### 项目简介

ADR 方法由 Pirozzoli (2006) 提出，用于近似非线性格式的修正波数。程序在周期网格上
取单模态初值，推进一个很小的时间步，用离散 Fourier 变换提取增长因子，得到
色散（实部）与耗散（虚部）曲线。

### 测试方法

- 网格：`N = mesh_num`，默认 `827`。**`mesh_num` 必须为奇数**；偶数网格可能出现振荡。
- 初值：`u_0(x_j) = cos(n · 2π j / N)`，`j = 0..N-1`。
- 时间推进：一步三步三阶 Runge–Kutta；`dt` 取足够小以消除时间离散误差。
- DFT：`û_n = (1/N) Σ u_j exp(-i n 2π j/N)`。
- 修正波数：`Phi_n = 2π i ln(û_t/û_0) / (dt · N)`。
- 横轴：`varphi_n = 2π n / N`；`Re(Phi)` 为色散，`Im(Phi)` 为耗散。

### 两套程序

| 程序 | 空间格式 | 差分/导数 |
|------|----------|-----------|
| `adr_weno` | WENO 重构半节点通量 | 固定二阶通量差分 |
| `adr_wcns` | WCNS 半节点插值 | 可选 `WCNSE2` / `WCNSE4` / `WCNSE6` / `HCS4` / `HCS6` |

`adr_wcns` 的插值部分直接使用 `../interpolation/src/interpolation/` 中的标准件。

### 差分套子（重要）

**同一插值格式、不同差分套子，ADR 色散/耗散曲线不同。**
做格式对比时应选择实际计算中该格式对应的差分套子；跨套子比较没有意义。

### 保留的 WENO 格式

`WENO5_LINEAR`、`WENO5_JS`、`WENO5_M`、`WENO5_Z`、`WENO_C`、
`WENO3_LINEAR`、`WENO3_JS`、`WENO3_TRIG`、`WENO3_LLR`、
`WENO3_GA_RBF_DEFINE`、`WENO3_GA_RBF_UNDEFINE`、`BAEZA_WENO3`、`THINC`、
`TENO_333`、`TENO_333_RBF`、`TENO_MOD_SEMI_ENO`、`TENO_MOD_SEMI_ENO_RBF`、
`RBF5_LINEAR`、`WENO5_DLR`。

### 保留的 WCNS 格式

`MUSCL`、`JS3`、`Z3`、`JS5`、`Z5`、`QR`、`MOC322`、`MOC33`、`MOC53`、`MOZ`、
`LINEAR3_M1P1`、`LINEAR3_M2P0`、`LINEAR3_0P2`、`LINEAR4_M2P1`、`LINEAR4_M1P2`、
`LINEAR5_M2P2`、`LINEAR6_M2P3`、`RBF_GA3`、`RBF_MQ3`、`RBF_IQ3`、
`RBF_GA5`、`RBF_MQ5`、`RBF_IQ5`。

### 未验证格式

以下子程序保留、可编译、可选，但标注为“实验/未验证”，不进入默认结果矩阵：

- `TENO_333_RBF`：原始实现中 `eps2h2` 未初始化。
- `THINC`：原始表达式的量纲存疑，在对称模板上可能产生非有限值。
- `WENO5_DLR`：原始文件已注明五阶 DLR 有问题；在 ADR 数据上可能产生 `NaN`。

### 精度与 `dt`

- 构建：`make PRECISION=qp`（默认，四精度）或 `make PRECISION=dp`（双精度）。
- `dt = 0` 表示按精度自动选择：qp → `1e-20`，dp → `1e-10`。
- 构建目录按精度隔离：`build/$(PRECISION)`、`bin/$(PRECISION)`。

### 构建与运行

```bash
make                    # 构建 adr_weno 与 adr_wcns（默认 qp）
make test               # Fortran 单元测试
make smoke              # 两个程序的端到端烟测
make PRECISION=dp       # 双精度构建
```

`input.nml` 分为四组：`&adr_nml`、`&weno_nml`、`&scheme_nml`、`&difference_nml`。
`adr_weno` 读取前两组，`adr_wcns` 读取 `&adr_nml`、`&scheme_nml`、`&difference_nml`。

### 结果

- `results/csv/`：`varphi, phi_re, phi_im`。
- `results/figures/`：`weno3`、`weno5`、`wcns3`、`wcns5`、`operators` 的色散/耗散图。
- `results/quick/`：`--quick` 的本地快速结果（不随仓库提交）。

```bash
sh scripts/run_results.sh --quick      # 33 点快速矩阵
sh scripts/run_results.sh              # 默认精选矩阵（827 点）
sh scripts/run_results.sh --all-schemes
```

### 引用

- S. Pirozzoli. On the spectral properties of shock-capturing schemes.
  *Journal of Computational Physics*, 2006.

---

## English

### Overview

The ADR method of Pirozzoli (2006) estimates the modified wavenumber of a nonlinear
scheme. A single Fourier mode is advanced by one small time step on a periodic grid;
a discrete Fourier transform gives the growth factor, whose real and imaginary parts
describe dispersion and dissipation.

### Method

- Grid: `N = mesh_num`, default `827`. **`mesh_num` must be odd**; even grids may oscillate.
- Initial data: `u_0(x_j) = cos(n · 2π j / N)`.
- Time stepping: one three-stage third-order Runge–Kutta step with a very small `dt`.
- DFT: `û_n = (1/N) Σ u_j exp(-i n 2π j/N)`.
- Modified wavenumber: `Phi_n = 2π i ln(û_t/û_0) / (dt · N)`.
- Horizontal axis: `varphi_n = 2π n / N`; `Re(Phi)` is dispersion, `Im(Phi)` is dissipation.

### Two programs

| Program | Spatial scheme | Difference operator |
|---------|----------------|---------------------|
| `adr_weno` | WENO reconstruction of half-node fluxes | fixed second-order flux difference |
| `adr_wcns` | WCNS half-node interpolation | selectable `WCNSE2` / `WCNSE4` / `WCNSE6` / `HCS4` / `HCS6` |

### Difference operators (important)

**With the same interpolation, a different difference operator gives a different ADR
curve.** When comparing schemes, choose the operator that the production scheme uses;
comparing across operators is not meaningful.

### Retained WENO schemes

`WENO5_LINEAR`, `WENO5_JS`, `WENO5_M`, `WENO5_Z`, `WENO_C`,
`WENO3_LINEAR`, `WENO3_JS`, `WENO3_TRIG`, `WENO3_LLR`,
`WENO3_GA_RBF_DEFINE`, `WENO3_GA_RBF_UNDEFINE`, `BAEZA_WENO3`, `THINC`,
`TENO_333`, `TENO_333_RBF`, `TENO_MOD_SEMI_ENO`, `TENO_MOD_SEMI_ENO_RBF`,
`RBF5_LINEAR`, `WENO5_DLR`.

### Retained WCNS schemes

`MUSCL`, `JS3`, `Z3`, `JS5`, `Z5`, `QR`, `MOC322`, `MOC33`, `MOC53`, `MOZ`,
`LINEAR3_M1P1`, `LINEAR3_M2P0`, `LINEAR3_0P2`, `LINEAR4_M2P1`, `LINEAR4_M1P2`,
`LINEAR5_M2P2`, `LINEAR6_M2P3`, `RBF_GA3`, `RBF_MQ3`, `RBF_IQ3`,
`RBF_GA5`, `RBF_MQ5`, `RBF_IQ5`.

### Unvalidated schemes

The following subroutines are kept and selectable, but are marked experimental and are
not part of the default result matrix:

- `TENO_333_RBF`: the original implementation leaves `eps2h2` uninitialized.
- `THINC`: the original expression is dimensionally suspect and can produce non-finite
  values on symmetric stencils.
- `WENO5_DLR`: the original source marks the fifth-order DLR as problematic; it can
  produce `NaN` on ADR data.

### Precision and `dt`

- Build with `make PRECISION=qp` (default, quad) or `make PRECISION=dp` (double).
- `dt = 0` selects the precision default: qp → `1e-20`, dp → `1e-10`.
- Build directories are separated by precision: `build/$(PRECISION)`, `bin/$(PRECISION)`.

### Build and run

```bash
make                 # adr_weno and adr_wcns (default qp)
make test            # Fortran unit tests
make smoke           # end-to-end smoke tests for both programs
make PRECISION=dp    # double-precision build
```

### Results

- `results/csv/`: `varphi, phi_re, phi_im`.
- `results/figures/`: dispersion and dissipation figures for `weno3`, `weno5`,
  `wcns3`, `wcns5`, `operators`.
- `results/quick/`: local `--quick` output (not committed).

```bash
sh scripts/run_results.sh --quick
sh scripts/run_results.sh
sh scripts/run_results.sh --all-schemes
```

### Reference

- S. Pirozzoli. On the spectral properties of shock-capturing schemes.
  *Journal of Computational Physics*, 2006.
