# 格式精度测试：方法与实现

## 1. 一维周期对流

求解

```
u_t + u_x = 0,   x ∈ [-1, 1],
```

周期边界，`Tend = 2`。网格 `x_j = -1 + j h`，`h = 2/N`，`N = nbase·2^(k-1)`，
`k = 1..levels`，默认 `nbase = 10`、`levels = 8`。

### 1.1 初值

七个连续光滑的周期初值，均由 `c1d_profile(x)` 提供：

| ID | 公式 |
|----|------|
| `sin2` | `sin²(πx)` |
| `sin` | `sin(πx)` |
| `sin3` | `sin³(πx)` |
| `sin_phase` | `sin(πx - sin(πx)/π)` |
| `gauss` | `exp(-20x²)` |
| `sqrt` | `sqrt(1 + 20x²)` |
| `rational` | `1/(1 + 20x²)` |

精确解为 `u0(x - t)`，按周期回绕到 `[-1,1]`。`c1d_initial` 与 `c1d_exact`
都调用同一个 `c1d_profile`，因此初值与精确解不会走样。

### 1.2 空间离散

半节点插值：

```
E_{j+1/2} = Interp(u_{j-3..j+4}),   j = -3..N+2.
```

`Interp` 由 `scheme_dispatch` 按 `scheme` 选择，使用 `interpolation/`
标准件。差分套子：

```
K_j = difference_apply(op, h, E_{j-3/2..j+5/2}, u_{j-2..j+2}).
```

`HCS4`/`HCS6` 的整节点参数在标量情形就是 `u_j`。

### 1.3 时间推进

四步四阶 Runge–Kutta：

```
k1 = R(u^n)
k2 = R(u^n + dt/2 k1)
k3 = R(u^n + dt/2 k2)
k4 = R(u^n + dt k3)
u^{n+1} = u^n + dt (k1 + 2 k2 + 2 k3 + k4)/6
```

`dt = CFL·h/max|u|`，默认 `CFL = 0.1`；末步截断到 `Tend`。

### 1.4 误差与阶数

```
L1  = Σ|u_num - u_ex|/(N+1)
L2  = sqrt(Σ(u_num - u_ex)²/(N+1))
L∞  = max|u_num - u_ex|
p_k = log(e_{k-1}/e_k)/log(2)
```

### 1.5 边界与 halo

插值在 `j = -3` 需要 `u(-6)`，在 `j = N+2` 需要 `u(N+6)`，因此 `halo = 7`。
周期 halo：`u(-i) = u(N-i)`、`u(N+i) = u(i)`，`i = 1..7`。旧程序的
`U(-5:N+5)` 会越界，本实现不再沿用。

## 2. 二维 Euler

二维 Euler 方程，周期边界，`Tend = 1`，`CFL = 0.1`，
`N = M = nbase·2^(k-1)`，默认 `nbase = 40`、`levels = 4`。

### 2.1 `vortex`

域 `[0,10]²`，`h = 10/N`。等熵涡参数
`EKC = 0.3`、`RC = 0.5`、`RFA = 0.35`、`γ = 1.4`，中心 `(5,5)`。
初值与精确解由同一套公式给出：

```
x̄ = j h - 5 - t,   ȳ = k h - 5
r = sqrt(x̄² + ȳ²)
T = 1 - (γ-1) EKC² exp(2 RFA (1 - (r/RC)²)) / (4 RFA γ)
ρ = T^(1/(γ-1))
u = 1 + (EKC/RC) exp(RFA (1 - (r/RC)²)) ȳ
w = -(EKC/RC) exp(RFA (1 - (r/RC)²)) x̄
p = ρ^γ
```

密度误差在 `(N+1)×(M+1)` 个节点上统计。

### 2.2 `sine`

域 `[0,2π]²`，`h = 2π/N`，`u = w = 1`，`p = 1`，
`ρ = ρ0 + A sin(x + y - 2t)`，默认 `ρ0 = 1`、`A = 0.2`。
初值取 `t = 0`；旧 AIA 数据对应 `ρ0 = 2`、`A = 1`。

### 2.3 通量

保留参考程序的流程：

1. `E5_X`/`E5_Z` 用中心点的特征矩阵把 7 点守恒量模板投影到特征变量，
   逐特征场调用 `interpolate_dispatch`，再投影回半节点守恒量；
2. `InterFlux` 对左、右两侧各重构一次，用修改的 Steger–Warming 分裂
   （`eps2 = 0.12²`、`sig = 0.5/γ`）得到 `Eflux/Gflux`；
3. `e2d_rhs` 把套子作用于半节点通量与整节点通量：

```
K_j = difference_apply(op, h, Eflux_{j-3/2..j+5/2}, Eint_{j-2..j+2})
    + difference_apply(op, h, Gflux_{k-3/2..k+5/2}, Gint_{k-2..k+2})
```

注意：二维 Euler 中 `HCS4`/`HCS6` 的整节点参数是**物理通量**
`Eint/Gint`，不是守恒变量 `Ubar`。标量一维测试中二者相同，二维中不同。

限幅：若重构出的密度或压力低于模板最小值的 0.1%，回退到中心点状态。
该逻辑与参考程序一致。

### 2.4 边界与 halo

`mt = 7`，覆盖插值与套子在 `j = -3..N+2` 上的最远读取
（`j-3 = -6`）。周期 halo 在 `e2d_boundary` 中填充。

### 2.5 误差与阶数

与一维相同的 L1/L2/L∞ 定义，分母为 `(N+1)(M+1)`；对密度 `ρ` 统计。

## 3. 实测阶数

- `Z5 + HCS6` 在光滑正弦算例上约 5 阶（例如 5.11、5.00）。
- `Z3`/`QR + HCS4` 在临界点附近降到 L1 约 2、L∞ 约 1–2；这与
  `interpolation/` 子项目 Z3/QR 的插值结果一致，是格式性质。
- 二维强等熵涡上 `Z5` 的 L1 阶数从约 3 逐步升到 `n=320` 的约 4.3，
  与旧 `js5.dat` 的 3–4 阶趋势一致；旧 `linear5.dat` 的干净 5 阶由本实现的
  `LINEAR5_M2P2 + HCS6` 逐值复现，相对差约 `1e-15`，说明二维管线移植正确。

## 4. 选择建议

- 想看设计阶数，用 `sine` + `Z5 + HCS6`。
- 想看强梯度/非线性权重行为，用 `vortex` + `Z5`，并对照 `LINEAR5_M2P2`。
- `AUTO` 只按阶数配对；显式套子不匹配时程序会警告，但仍运行。
