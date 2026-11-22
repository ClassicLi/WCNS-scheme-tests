# euler2d 方法说明

本文说明 `euler2d` 实际使用的离散方法。公式与常数均取自
`src/euler2d/`、`src/common/`、`src/wcns/` 与 `../interpolation/` 的代码。

## 1. 控制方程

二维 Euler 方程组写作

```
∂U/∂t + ∂E/∂x + ∂G/∂z = S,
U = (ρ, ρu, ρw, E)^T,
Eflux = (ρu, ρu² + p, ρuw, (E + p)u)^T,
Gflux = (ρw, ρuw, ρw² + p, (E + p)w)^T,
p = β (E - 0.5 ρ (u² + w²)),   β = γ - 1.
```

默认 `γ = 1.4`；`rayleigh_taylor` 用 `γ = 5/3`。`S` 只有 RT 非零：
`S = (0, 0, ρ, ρw)^T`，由 `e2d_source` 在空间残差组装完成后加入。

## 2. 网格、halo 与坐标原点

设 `n = nx`、`m = nz`、`h = 1/m`。单元中心坐标为

```
x_j = x_origin + j h,   j = 0..n,
z_k = z_origin + k h,   k = 0..m.
```

物理域是 `[x_origin, x_origin + nx/nz] × [z_origin, z_origin + 1]`。
`nx = nz · ratio_num / ratio_den` 由算例比例推出。

守恒量 `ubar`、整数节点通量 `eint/gint` 分配为
`(-mt:n+mt, -mt:m+mt, 4)`，半节点通量 `eflux`、`gflux` 分配为
`(-3:n+2, 0:m, 4)` 与 `(0:n, -3:m+2, 4)`。halo 取 **`mt = 7`**：

- 八点特征插值模板在 `j = -3` 读 `ubar(j-3:j+4) = ubar(-6:1)`；
  在 `j = n+2` 读 `ubar(n-1:n+6)`，因此需要 `-6` 与 `n+6`；
- 差分套子的整数节点输入在 `j = 0` 读 `eint(-2:2)`，在 `j = n` 读
  `eint(n-2:n+2)`。

论文原程序用 `MT = 5`，在上述两个位置会越界；本实现统一改为
`mt = 7`，并在 `gfortran -fcheck=all` 下检查。

坐标原点按初值定义：

| 算例 | `x_origin` | `z_origin` |
|------|------------|------------|
| `double_mach` | 0 | 0 |
| `kelvin_helmholtz` | 0 | 0 |
| `rayleigh_taylor` | 0 | 0 |
| `riemann1` | 0 | 0 |
| `riemann2` | 0 | 0 |
| `shock_bubble` | -0.1 | -0.5 |
| `shock_vortex` | -0.25 | -0.5 |

输出坐标统一使用初值原点。`shock_vortex` 的旧 `export` 子程序曾用
`(0,0)` 输出，与 `Initial` 的 `(-0.25,-0.5)` 不一致；这里以初值为准。

## 3. 特征变量插值

对每个半节点 `j+1/2`（或 `k+1/2`），从八点模板 `U_{-3}, ..., U_4`
投影到特征变量。记中心点（模板下标 0）的

```
u = (ρu)/ρ,   w = (ρw)/ρ,
α = 0.5 (u² + w²),
p = β (E - α ρ),
c² = γ p / ρ,   c = sqrt(c²).
```

x 方向左特征矩阵与右特征矩阵为

```
ml_x =
[ αβ - c²,   -βu,      -βw,   β ]
[ -w,         0,         1,    0 ]
[ αβ - c u,   c - βu,  -βw,   β ]
[ αβ + c u,  -c - βu,  -βw,   β ]

mr_x = (1/(2c²)) *
[ -2,     0,            1,                    1                   ]
[ -2u,    0,            u + c,                u - c               ]
[ -2w,    2c²,          w,                    w                   ]
[ -2α,    2w c²,        α + c u + c²/β,       α - c u + c²/β      ]
```

z 方向为

```
ml_z =
[ αβ - c²,   -βu,      -βw,   β ]
[ -u,         1,         0,    0 ]
[ αβ - c w,  -βu,      c - βw, β ]
[ αβ + c w,  -βu,     -c - βw, β ]

mr_z = (1/(2c²)) *
[ -2,     0,            1,                    1                   ]
[ -2u,    2c²,          u,                    u                   ]
[ -2w,    0,            w + c,                w - c               ]
[ -2α,    2u c²,        α + c w + c²/β,       α - c w + c²/β      ]
```

代码中矩阵按列存储为 `mr(i, p)`；`ml_x`、`mr_x`、`ml_z`、`mr_z` 的
逐项公式见 `e2d_flux.f90` 的 `eigen_matrix_x`、`eigen_matrix_z`。
四个特征场

```
q_p = ml(p,1) U_1 + ml(p,2) U_2 + ml(p,3) U_3 + ml(p,4) U_4
```

分别调用 `../interpolation/` 的 `interpolate_dispatch` 得到 `q_half(p)`，
再反投影

```
U_half(1:4) = mr(:,1) q_half(1) + ... + mr(:,4) q_half(4).
```

左右状态各做一次：左状态模板为 `ubar(j-3:j+4)`，右状态模板为
`ubar(j+4:j-3:-1)`。这样得到 `U_L` 与 `U_R` 后进入通量分裂。

### 第一层修正：插值正性回退

若反投影得到的半点 `ρ_half` 或 `p_half` 不大于模板内
`j-2..j+2` 五点最小值的 `0.1%`，即

```
ρ_half <= 0.001 min(ρ_{j-2..j+2})   或
p_half <= 0.001 min(p_{j-2..j+2}),
```

则令 `U_half = U_0`（中心点），并把 `interp_corrections` 加一。
x、z 两个方向都按此执行。

## 4. 修正 Steger–Warming 分裂

记 `eps² = 0.12²`、`sig = 1/(2γ)`。对 x 方向，

```
λ1 = u,   λ4 = u + c,   λ5 = u - c,
λ± = (λ ± sqrt(λ² + eps²))/2,
diss = -sig (2λ1 - λ4 - λ5),
disp =  sig (λ4 - λ5) c.
```

左状态取 `λ+`，得到 `fp`：

```
fp(1) = diss + λ1
fp(2) = u (diss + λ1) + disp
fp(3) = w (diss + λ1)
fp(4) = (α + c²/β) diss + (α + p/(ρβ)) λ1 + u disp
```

右状态取 `λ-`，用同样的式子得到 `fn`。半节点通量为

```
Eflux = ρ_L fp + ρ_R fn.
```

z 方向把 `λ1 = w`，并把 `disp` 放在第三分量：

```
fp(1) = diss + λ1
fp(2) = u (diss + λ1)
fp(3) = w (diss + λ1) + disp
fp(4) = (α + c²/β) diss + (α + p/(ρβ)) λ1 + w disp
```

`Gflux = ρ_L fp + ρ_R fn`。这与原始 Steger–Warming 的
`E^+ + E^-` 写法不同：代码把特征速度的贡献除以密度后重组，
`fp/fn` 是“单位密度”的通量因子。

## 5. 差分套子与 AUTO 配对

半节点通量数组按 `f_{j-5/2} ... f_{j+5/2}` 传入差分套子，整数节点值按
`u_{j-2} ... u_{j+2}` 传入。五个套子的节点残差为：

```
WCNSE2:
  K_j = -(f_{j+1/2} - f_{j-1/2}) / h

WCNSE4:
  K_j = -[ (9/8)(f_{j+1/2} - f_{j-1/2})
           - (1/24)(f_{j+3/2} - f_{j-3/2}) ] / h

WCNSE6:
  K_j = -[ (75/64)(f_{j+1/2} - f_{j-1/2})
           - (25/384)(f_{j+3/2} - f_{j-3/2})
           + (3/640)(f_{j+5/2} - f_{j-5/2}) ] / h

HCS4:
  K_j = -[ a (f_{j+1/2} - f_{j-1/2})
           + b (u_{j+1} - u_{j-1}) ] / h,
  a = 4/3,  b = -1/6

HCS6:
  K_j = -[ a (f_{j+1/2} - f_{j-1/2})
           + b (u_{j+1} - u_{j-1})
           + c (u_{j+2} - u_{j-2}) ] / h,
  a = 64/45,  b = -2/9,  c = 1/180
```

`difference_apply` 返回的是 `-dE/dx`。二维残差为

```
rhs(j,k,:) = difference_apply(op, h, eflux(j-3:j+2,k,:),
                              eint(j-2:j+2,k,:))
           + difference_apply(op, h, gflux(j,k-3:k+2,:),
                              gint(j,k-2:k+2,:)).
```

`AUTO` 按 `scheme_order` 配对：

| 插值阶数 | 套子 |
|----------|------|
| 2 | `WCNSE2` |
| 3 | `HCS4` |
| 4 | `WCNSE4` |
| 5 | `HCS6` |
| 6 | `WCNSE6` |

默认 `QR` 的阶数为 3，因此 `AUTO` 给出 `HCS4`，即
`QR + AUTO(HCS4)`。显式套子与阶数不匹配时打印警告并继续。

## 6. 时间推进、CFL 与源项

时间推进固定为四步四阶 Runge–Kutta：

```
u^(1) = u^n + 0.5 dt k1
u^(2) = u^n + 0.5 dt k2
u^(3) = u^n + dt k3
u^(n+1) = u^n + dt (k1 + 2k2 + 2k3 + k4)/6
k_r = RHS(u^(r), t_n + dt)
```

每个 RK 子步前重新填充边界，残差中包含 `e2d_source` 的源项。RT 的
重力源项加在 z 动量与能量方程上，与差分套子无关。

时间步为

```
dt = CFL · h / max_{j,k} ( sqrt(u² + w²) + sqrt(γ p / ρ) ).
```

第二个根号是物理声速 `c = sqrt(gama*p/rho)`（即 `sqrt(γ p / ρ)`；
`gama` 为代码中的比热比变量），与特征分解和 Steger–Warming 分裂使用的
声速一致。末步若越过目标时刻，则截断为 `t_end - t`。
默认 `CFL = 0.3`；结果脚本使用 `0.1`。

## 7. 第二层修正：单元正性保护

每个 RK 子步结束后扫描内部单元 `0..n, 0..m`。按“先密度、后压力”的顺序：

```
若 ρ < 1e-8:
    U(j,k,:) = (U0(j-1,k,:) + U0(j+1,k,:)
              + U0(j,k-1,:) + U0(j,k+1,:)) / 4
    计数加一
否则:
    计算 u = ρu/ρ, w = ρw/ρ, p = β(E - 0.5ρ(u²+w²))
    若 p < 1e-8:
        U(j,k,:) = 四邻点平均（取自本步开始时的 U0）
        计数加一
```

四邻点来自本步开始时保存的 `u0`，而不是刚更新过的状态。
`cell_corrections` 记录触发次数。先判断密度可避免用零或负密度作除数。

## 8. 七个算例的边界类型

边界填充由 `e2d_boundary(case, t)` 分发。所有算例都填充 `mt = 7` 层 halo。

| 算例 | 左边界 | 右边界 | 下边界 | 上边界 |
|------|--------|--------|--------|--------|
| `double_mach` | 固定入流 `ρ=8, u=7.14471, w=-4.125, p=116.5` | 零阶外推 | `x < 1/6` 固定入流；否则镜像反射（`w` 反号） | 随 `t` 移动的激波：`x_s(t) = 1/6 + (1+20t)/√3`；`x > x_s` 用 `(1.4,0,0,1)`，`x <= x_s` 用 `(8,7.14471,-4.125,116.5)` |
| `kelvin_helmholtz` | 镜像复制 | 镜像复制 | 镜像复制 | 镜像复制 |
| `rayleigh_taylor` | 反射（`u` 反号） | 反射（`u` 反号） | 固定 `ρ=2, p=1, u=w=0` | 固定 `ρ=1, p=2.5, u=w=0` |
| `riemann1` | 镜像复制 | 镜像复制 | 镜像复制 | 镜像复制 |
| `riemann2` | 镜像复制 | 镜像复制 | 镜像复制 | 镜像复制 |
| `shock_bubble` | 固定入流 `ρ=3.81, u=2.85, w=0, p=10` | 零阶外推 | 反射（`w` 反号） | 反射（`w` 反号） |
| `shock_vortex` | 镜像复制 | 镜像复制 | 镜像复制 | 镜像复制 |

“镜像复制”指源程序的 `U(:, -p) = U(:, p-1)`、
`U(:, m+p) = U(:, m-p+1)`（x 方向同理），不是带平移的严格周期边界。
`double_mach` 底部在 `x - 1/6 < 1e-8` 时用入流，否则把下方 halo
写成 `U(j,-p) = U(j,p-1)` 并把 `w` 反号。

## 9. 输出与日志

Tecplot ASCII 点格式，变量 `x, z, rho`。主文件
`<case>_<nz>x<nx>_t1.dat`；`shock_bubble` 在 `tend2` 写 `_t1.dat`、
在 `tend` 写 `_t2.dat`。`write_sample = .true.` 时另写 `_t0.dat`。

日志 `<case>.log` 记录算例、插值、套子、网格、CFL、终止时间、
`gamma`、步数、`interp_corrections`、`cell_corrections`、`nan_count`
与 CPU 时间。非有限值计数只扫描物理单元，不扫描 halo。
