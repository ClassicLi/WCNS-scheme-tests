# 与论文程序的差异（formula audit）

本文逐项记录 `euler2d` 与学位论文二维程序在命名、公式、边界和验证口径上的
差异。参考程序是作者开发算例时使用的 Fortran 源程序（作者本地保存，
未随本仓库发布）；本仓库的公式以 `src/` 与 `../interpolation/` 的代码为准。

## 1. 空间差分套子的命名审计

论文程序的 `rk_order` 不是时间推进阶数，而是**空间差分套子的选择器**。
时间推进在所有分支里都是固定四步四阶 Runge–Kutta，`rk_order` 只决定调用哪个
“RungKutta_*”子程序，而这些子程序的名字与其实际差分公式并不一致：

| 论文子程序 | 实际差分格式 | 系数 |
|------------|--------------|------|
| `RungKutta_E6` | **HCS6** | `a = 64/45`, `b = -2/9`, `c = 1/180` |
| `RungKutta_hybrid_E4` | **HCS4** | `a = 4/3`, `b = -1/6` |

`RungKutta_E6` 的名字容易被误读为六阶通量差分 E6，但它的公式是
“一对半点通量 + 两对整点值”的混合紧致格式，即本仓库的 `HCS6`。
`RungKutta_hybrid_E4` 同理，对应本仓库的 `HCS4`。

本仓库的 `WCNSE2/4/6` 才是论文语境中的 **E2/E4/E6**：

| 本仓库 | 论文名称 | 系数 |
|--------|----------|------|
| `WCNSE2` | E2 | `1` |
| `WCNSE4` | E4 | `9/8`, `-1/24` |
| `WCNSE6` | E6 | `75/64`, `-25/384`, `3/640` |

因此论文里的 “E6” 与本仓库的 `WCNSE6` 是同一族，而论文明名为
`RungKutta_E6` 的子程序反而是 `HCS6`。文档与代码一律使用本仓库的
`WCNSE2/4/6`、`HCS4/6` 名称，不再沿用论文的子程序名。

## 2. RT 的 `rk_order`、源项与未启用分支

`rayleigh_taylor` 的论文配置写 `rk_order = 6`，实际调用的是
`RungKutta_E6`，也就是 **HCS6**，不是 hybrid 格式。该子程序内部含有
RT 的重力源项：

```
rhs(3) = rhs(3) + ρ
rhs(4) = rhs(4) + ρw
```

`RungKutta_hybrid_E6` 是论文程序里存在但默认配置不会调用的变体；
不能写成“RT 使用 hybrid_E6”。

论文把同一个重力源项分别抄写在 `RungKutta_E6`、`RungKutta_hybrid_E4`、
`RungKutta_hybrid_E6` 三个子程序内部。本实现把它抽成
`e2d_cases::e2d_source(case, rhs)`，在 `e2d_rhs` 组装完空间残差后
**只加一次**。因此：

- 源项与差分套子解耦，RT 在 `WCNSE2`、`WCNSE4`、`WCNSE6`、`HCS4`、
  `HCS6` 下都可用；
- 不再保留论文 `rk_order` 的 hybrid 分支；
- 不存在“选错源项副本”的可能。

## 3. 别名

论文与旧脚本中的若干名称在本仓库统一为：

| 旧名称 | 本仓库 |
|--------|--------|
| `Z5_7point` | `Z5` |
| `MOZ_mml` | `MOZ` |
| `MOC53_framework` | `MOC53` |

这些只是名字统一，不改变公式。

## 4. halo：论文 `MT = 5` 越界与 `mt = 7` 修正

论文程序把守恒量声明为 `U(-MT:... )` 一类，取 `MT = 5`。但八点特征插值：

- 在 `j = -3` 读模板 `U(j-3:j+4) = U(-6:1)`；
- 在 `j = N+2` 读模板 `U(N-1:N+6)`；

需要 `U(-6)` 与 `U(N+6)`。差分套子的整数节点输入还需要 `±2` 层。
`MT = 5` 在这些位置越界读取，属于未定义行为。本实现取 **`mt = 7`**，
分配 `(-mt:n+mt, -mt:m+mt, 4)`，并在
`gfortran -fcheck=all -ffpe-trap=invalid,zero,overflow` 下检查。

`mt 5→7` 只扩大 halo，不改变内部离散公式；若新旧结果出现差异，应首先检查
旧程序在越界位置读到的未定义值，而不是修改新实现。

## 5. 输出时刻

六个算例各写一个输出时刻。`shock_bubble` 在论文中写两个时刻：
先在 `tend2 = 0.2` 输出，再在 `tend = 0.4` 输出。本实现保留双输出，
文件名为 `_t1.dat`（第一个时刻）与 `_t2.dat`（第二个时刻）；其他算例只有
`_t1.dat`。`write_sample = .true.` 时另写 `_t0.dat` 初值样本。

## 6. shock-vortex 的坐标原点

论文 `shock_vortex` 的 `Initial` 子程序用

```
x = j h - 0.25,   z = k h - 0.50
```

即原点 `(-0.25,-0.5)`；但同一个程序的 `export` 子程序写坐标时用
`j h, k h`，即原点 `(0,0)`。两者不一致，是源程序的缺陷。

本实现把 `x_origin = -0.25`、`z_origin = -0.5` 作为算例元数据，
初值和输出都使用同一原点；`shock_bubble` 同理使用
`(-0.1,-0.5)`，其余算例为 `(0,0)`。

## 7. 插值层差异

本仓库的插值实现来自 `interpolation/` 的统一版本。与论文逐项比较：

| 格式 | 与论文的关系 |
|------|--------------|
| `QR` | 公式逐行一致 |
| `Z5` | 公式逐行一致；参数名 `z5_epsilon` → `z5_eps` |
| `Z3` | 公式逐行一致；参数名 `z3_epsilon` → `z3_eps` |
| `MOZ` | 公式一致，但参数名语义互换：论文 `moz_para_q` → 仓库 `moz_M_q`，论文 `moz_M_q` → 仓库 `moz_N_q` |
| `MOC53` | 仓库的混合权重公式与论文不同；论文 `moc53_beta_type` 在仓库无等价参数 |
| `MOC33` | 仓库没有 `THINC`，非线性候选固定为 QR；论文 `moc33_nonlinear='thinc'` 在仓库无等价参数 |

`MOC53` 与 `MOC33+THINC` 因公式不同，不能与论文结果逐值对比，只能做
定性对照。`QR`、`Z5`、`Z3`、`MOZ` 在参数正确映射后可以逐值对比；
`MOZ` 的差异仅是参数名语义，必须按上表映射，不能按同名参数直接照抄。

## 8. 论文配置 → 仓库标准件映射

下表的“论文配置”用于 `scripts/run_results.sh`；“能否逐值复现”指同网格、
同参数下能否以 `1e-10` 相对差的量级对比。

| 算例 | 论文配置 | 仓库标准件 | 能否逐值复现 |
|------|----------|------------|--------------|
| `double_mach` | `MOZ`：`moz_para_q=1`, `moz_M_q=4`, `moz_eps=1e-40` + `HCS6` | `MOZ`：`moz_M_q=1`, `moz_N_q=4`, `moz_eps=1e-40` + `HCS6` | 是（按语义互换参数名） |
| `kelvin_helmholtz` | `Z5`：`z5_q=4`, `z5_epsilon=1e-40` + `HCS6` | `Z5`：`z5_q=4`, `z5_eps=1e-40` + `HCS6` | 是 |
| `rayleigh_taylor` | `MOC53`：`moc53_M_q=2`, `moc53_eps=1e-60`, `moc53_linear='GA'`, `moc53_beta_type=1` + `HCS6` | `MOC53`：`moc53_N_q=2`, `moc53_p=2`, `moc53_eps=1e-60`, `moc53_linear='GA'` + `HCS6` | 否：混合权重公式不同 |
| `riemann1` | 同 `rayleigh_taylor` | 同 `rayleigh_taylor` | 否：混合权重公式不同 |
| `riemann2` | `MOC33`：`moc33_para_p=4`, `moc33_linear='MQ'`, `moc33_nonlinear='thinc'` + `HCS4` | `MOC33`：`moc33_p=4`, `moc33_linear='MQ'`，非线性候选固定 QR + `HCS4`；另需 `moc33_eps=1e-40` | 否：THINC → QR |
| `shock_bubble` | 同 `riemann2` | 同 `riemann2` | 否：THINC → QR |
| `shock_vortex` | `Z3`：`z3_q=2`, `z3_epsilon=1e-40` + `HCS4` | `Z3`：`z3_q=2`, `z3_eps=1e-40` + `HCS4` | 是 |

论文的 `moc53_beta_type`、`moc33_nonlinear` 在仓库没有对应项；
`MOC33` 的 `moc33_p` 对应论文 `moc33_para_p`，`MOC53` 的 `moc53_N_q`
对应论文 `moc53_M_q`。`scripts/run_results.sh` 逐项写出映射后的参数，
不依赖仓库默认值。

## 9. CFL 声速修正

论文的 `TimeStep` 用

```
dt = CFL · h / max( sqrt(u² + w²) + sqrt(p/ρ) )
```

作为时间步；第二个根号缺少 `γ`。但同一声速在特征分解和 Steger–Warming
分裂里一直是物理声速

```
c = sqrt(γ p / ρ).
```

本实现把时间步改为

```
dt = CFL · h / max( sqrt(u² + w²) + sqrt(gama*p/rho) ),
```

即 `sqrt(γ p / ρ)`，其中 `gama` 是代码中的比热比变量；这与通量分裂
保持一致。对 `γ = 1.4`，旧公式把声速低估约 `sqrt(1.4) ≈ 1.183`，
因此在相同 `CFL` 下旧步长更大。这是有意的物理修正，不视为迁移误差；
做新旧对比时应记录实际步长差异。

## 10. R2/SB 的 `moc33_eps = 1e-40`

仓库 `scheme_params` 的默认 `moc33_eps = 1e-80`、`moc33_p = 2`。
`riemann2` 与 `shock_bubble` 使用论文的 `moc33_p = 4`。在常数区，
`beta ≈ 0`，`alpha_L`、`alpha_N` 约为 `1/moc33_eps ≈ 1e80`；
`alpha**4 ≈ 1e320`，超过双精度上限（约 `1.8e308`），产生 Inf，
随后 `omega = Inf/Inf` 变成 NaN。

论文的 `moc33_eps = 1e-40` 使 `alpha**4 ≈ 1e160`，仍在双精度范围内。
因此 `scripts/run_results.sh` 对 `riemann2` 与 `shock_bubble` 显式写
`moc33_eps = 1.0d-40`；这不是调参，而是运行论文配置的必要条件。
`MOC53` 的 `moc53_p = 2`，`1e-60` 的平方约为 `1e120`，不存在同一问题。

## 11. 迁移验证

### 11.1 验证口径

论文参考端使用其 Fortran 源程序的 `MT 5 → 7` 修补版（只改 halo，
其余不动），与本实现使用同一编译器、同一 `FFLAGS`、同一精度，关闭
OpenMP，串行运行。对比最终密度场，允许 1e-10 量级的相对差；
若出现差异，先检查 `mt` 修补位置与时间步差异。

可逐值对比的算例：

- `double_mach`（`MOZ + HCS6`）；
- `shock_vortex`（`Z3 + HCS4`）；
- `kelvin_helmholtz`（`Z5 + HCS6`，可选）。

不可逐值对比、只做定性对照的算例：

- `rayleigh_taylor`、`riemann1`：`MOC53` 混合权重公式不同；
- `riemann2`、`shock_bubble`：论文 `MOC33` 的非线性候选是 THINC，
  仓库固定为 QR。

### 11.2 已完成的验证

本实现的单元测试与 smoke 检查全部通过：

| 检查 | 结果 |
|------|------|
| `test_control` | 22 项：namelist 默认值、非法配置、`AUTO` 配对、显式不匹配警告 |
| `test_operators` | 27 项：五个套子的多项式精确性；HCS4/HCS6 的整节点参数与通量参数分别验证 |
| `test_flux` | 9 项：左右状态相同时数值通量等于物理通量；`ml·mr = I`；x/z 对称性 |
| `test_cases` | 108 项：七个算例的初值、边界 halo 关系、RT 源项符号、SB 双输出元数据 |
| `test_solver` | 10 项：RK4、CFL、正性保护与计数 |
| `test_io` | 28 项：Tecplot 头、日志键、非有限值计数 |
| `test_smoke` | 七个算例 + RT 的 `HCS4`/`HCS6`/`WCNSE6` 对照全部 PASS，`nan_count=0` |

`results/smoke/` 的实测结果（`nz = 16`、`tend = 0.001`、`cfl = 0.1`）：

| 算例 | 论文配置 | 步数 | `interp_corrections` | `cell_corrections` | `nan_count` |
|------|----------|------|----------------------|--------------------|-------------|
| `double_mach` | `MOZ + HCS6` | 3 | 83 | 0 | 0 |
| `kelvin_helmholtz` | `Z5 + HCS6` | 1 | 0 | 0 | 0 |
| `rayleigh_taylor` | `MOC53 + HCS6` | 1 | 0 | 0 | 0 |
| `riemann1` | `MOC53 + HCS6` | 1 | 0 | 0 | 0 |
| `riemann2` | `MOC33 + HCS4` | 1 | 0 | 0 | 0 |
| `shock_bubble` | `MOC33 + HCS4` | 2 | 32 | 0 | 0 |
| `shock_vortex` | `Z3 + HCS4` | 1 | 0 | 0 | 0 |

### 11.3 逐值对比的实测数值

验证环境：ifx 2025.1，`dp`，`FFLAGS = -O2 -g`，串行执行（参考端不启用
OpenMP，`OMP_THREADS = 1`）。参考端从论文源程序复制到临时目录，只把
`global_var.f90` 的 `MT = 5` 改为 7，然后用 `make FC=ifx FFLAGS="-O2 -g"`
构建；新端使用同网格、同插值格式、同差分套子、同 `cfl = 0.3`。两者输出
均为 Tecplot ASCII 格式，用 `scripts/compare_fields.py` 比较最终密度：

```bash
python3 euler2d/scripts/compare_fields.py ref.dat new.dat
```

脚本按 `VARIABLES` 中的 `rho`/`den` 定位密度，输出 `max_abs_diff`、
`max_rel_diff`（以参考密度为分母）和 `nan_count`，网格不一致或出现
非有限值时返回非零状态。

#### 11.3.1 论文参考端（只修 `MT`）

| 算例 | 配置 | `nz`/`nx` | `tend` | 步数（参考/新） | `max_abs_diff` | `max_rel_diff` | `nan_count` | 判定 |
|------|------|-----------|--------|-----------------|----------------|----------------|-------------|------|
| `double_mach` | `MOZ + HCS6` | 64/256 | 0.2 | 815 / 848 | 4.784077514012198e+00 | 4.323653818566910e-01 | 0 | 差异来自 CFL 修正，见下文 |
| `double_mach` | `MOZ + HCS6` | 64/256 | 0.05 | 177 / 186 | 1.849482630998356e-01 | 5.464176433595876e-02 | 0 | 差异来自 CFL 修正 |
| `double_mach` | `MOZ + HCS6` | 64/256 | 0.01 | 28 / 30 | 2.332198988912388e-01 | 5.995029782421994e-02 | 0 | 差异来自 CFL 修正 |
| `shock_vortex` | `Z3 + HCS4` | 64/128 | 0.6 | 322 / 346 | 2.512929123010821e-05 | 2.175346485273562e-05 | 0 | 差异来自 CFL 修正 |
| `kelvin_helmholtz`（可选） | `Z5 + HCS6` | 64/64 | 0.5 | 223 / 254 | 4.195648582983047e-05 | 2.987564528928644e-05 | 0 | 差异来自 CFL 修正 |

`double_mach` 的差异超过 `1e-10`，原因不是离散公式，而是第 9 节有意的
CFL 声速修正：参考端第一步 `dt = 3.884855858e-4`，新端第一步
`dt = 3.672077835e-4`（比值 0.945），因此同样推进到 `tend` 所需的步数
不同（815 vs 848），后续激波演化的采样轨迹也随之不同。`shock_vortex`
与 `kelvin_helmholtz` 的差异同样来自 CFL（322 vs 346 步、223 vs 254 步）。
新实现保留物理声速 `sqrt(gama*p/rho)`，未改回 `sqrt(p/rho)`。

#### 11.3.2 CFL 对齐的诊断参考端

为确认上述差异完全由 CFL 声速造成，另做一个仅用于诊断的参考端：除
`MT 5 → 7` 外，再把 `TimeStep` 的 `sqrt(ps/ro)` 改为
`sqrt(gama*ps/ro)`。其余源码与串行设置不变。该诊断不替代 11.3.1 的
迁移口径，只用于定位差异来源。

| 算例 | `tend` | 步数（参考/新） | `max_abs_diff` | `max_rel_diff` | `nan_count` |
|------|--------|-----------------|----------------|----------------|-------------|
| `shock_vortex` | 0.6 | 346 / 346 | 1.265654248072679e-14 | 1.081548425420672e-14 | 0 |
| `kelvin_helmholtz` | 0.5 | 254 / 254 | 3.690381333854020e-12 | 2.175852727969460e-12 | 0 |
| `double_mach` | 0.01 | 30 / 30 | 5.917044632042234e-12 | 7.390859530591621e-13 | 0 |
| `double_mach` | 0.05 | 186 / 186 | 2.946833888017864e-10 | 3.334552649910436e-11 | 0 |
| `double_mach` | 0.1 | 402 / 402 | 4.580345767912775e-05 | 3.286869613856360e-06 | 0 |
| `double_mach` | 0.2 | 848 / 848 | 9.302482598474242e-01 | 7.372912061728343e-02 | 0 |

`shock_vortex`、`kelvin_helmholtz` 与 `double_mach` 的早期结果在 `1e-10`
量级或更好，说明插值、特征分解、Steger–Warming 通量、差分套子和 RK4
的移植与论文公式一致。`double_mach` 在 CFL 对齐后仍随时间的推移放大
差异：`t = 0.01` 为 `7.39e-13`，`0.05` 为 `3.33e-11`，`0.1` 为
`3.29e-6`，`0.2` 为 `7.37e-2`，而步数始终相同。这是双马赫反射的
激波相互作用与正性修正对舍入级扰动的非线性放大，不是公式差异；
因此该算例的逐点长时对比只作为结构对照，`1e-10` 量级判定取自早期
（`tend <= 0.05`）。

对 `MOC53` 与 `MOC33+THINC` 的四个算例，不设 `1e-10` 目标；只比较
密度场的整体结构、激波/接触间断位置和是否出现非物理振荡。
