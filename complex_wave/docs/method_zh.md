# 复杂波时间推进测试：方法与实现

## 1. 控制方程与网格

求解一维周期线性对流方程

```
u_t + u_x = 0,   x ∈ [-1, 1]
```

波速为 1，因此一个波形周期为

```
T = 2 / 1 = 2
```

网格为

```
x_j = -1 + j h,   h = 2 / N,   j = 0..N
```

周期边界通过端点重复实现：`u(N) = u(0)`。半节点插值需要左右各 6 个 halo 点，
本实现统一使用 `halo = 7`：

```
u(-i) = u(N-i),   u(N+i) = u(i),   i = 1..halo
```

## 2. 初值

初值由四段组成，参数与参考程序一致：`δ = 0.005`、`z = -0.7`、`a = 0.5`。

```
G(x, z) = exp(-log(2) (x - z)^2 / 0.03^2)
F(x, a) = sqrt(max(0, 1 - (10 (x - a))^2))

k = N/10 .. 2N/10   :  u_k = (G(x_k, z-δ) + 4 G(x_k, z) + G(x_k, z+δ)) / 6
k = 3N/10 .. 4N/10  :  u_k = 1
k = 5N/10 .. 6N/10  :  u_k = max(0, 1 - |10 x_k - 1|)
k = 7N/10 .. 8N/10  :  u_k = (F(x_k, a-δ) + 4 F(x_k, a) + F(x_k, a+δ)) / 6
其余                 :  u_k = 0
```

三角波的 `max(0, ...)` 只在 `N` 不是 10 的倍数时起作用：这时分段端点可能
略微落在 `[0, 0.2]` 之外，参考公式会给出很小的负值。

## 3. 空间离散

### 3.1 半节点插值

```
E_{j+1/2} = Interp(u_{j-3..j+4}),   j = -3 .. N+2
```

`Interp` 为 `../interpolation/src/interpolation/` 中的 WCNS 标准件，由
`&scheme_nml` 的 `scheme` 选择。本算例使用 `JS3`、`JS5`、`Z3`、`Z5`、`QR`。

### 3.2 差分套子

节点空间右端项为 `K_j = -dE/dx`。设

```
ΔE_m = E_{j+m+1/2} - E_{j-m-1/2},   m = 0, 1, 2
Δu_m = u_{j+m} - u_{j-m},           m = 1, 2
```

五个可选套子为：

```
WCNSE2:  K_j = -ΔE_0 / h
WCNSE4:  K_j = -[(9/8) ΔE_0 - (1/24) ΔE_1] / h
WCNSE6:  K_j = -[(75/64) ΔE_0 - (25/384) ΔE_1 + (3/640) ΔE_2] / h
HCS4:    K_j = -[(4/3) ΔE_0 - (1/6) Δu_1] / h
HCS6:    K_j = -[(64/45) ΔE_0 - (2/9) Δu_1 + (1/180) Δu_2] / h
```

`HCS6` 的系数由 Taylor 展开匹配到六阶得到。把

```
ΔE_0/h = u' + (h^2/24) u''' + (h^4/1920) u^(5) + ...
Δu_1/h = 2u' + (h^2/3) u''' + (h^4/60) u^(5) + ...
Δu_2/h = 4u' + (8h^2/3) u''' + (8h^4/15) u^(5) + ...
```

代入 `a ΔE_0/h + b Δu_1/h + c Δu_2/h`，令 `u'` 的系数为 1、`u'''` 与 `u^(5)`
的系数为 0，得到

```
a + 2b + 4c = 1
a/24 + b/3 + 8c/3 = 0
a/1920 + b/60 + 8c/15 = 0
```

解得 `a = 64/45`、`b = -2/9`、`c = 1/180`。

## 4. 时间推进

采用四步四阶 Runge–Kutta：

```
k1 = R(u^n)
k2 = R(u^n + dt/2 k1)
k3 = R(u^n + dt/2 k2)
k4 = R(u^n + dt   k3)
u^{n+1} = u^n + dt (k1 + 2 k2 + 2 k3 + k4) / 6
```

其中 `R` 为第 3 节的空间右端项。时间步长按参考程序取

```
dt = CFL h / max_j |u_j|,   CFL = 0.3（默认）
```

每个周期末端的最后一步会截断为 `target - t`，使快照恰好落在 `t = kT`。
若最大幅值为 0 或非有限，程序报错退出，避免 `dt = 0` 或 NaN 造成死循环。

## 5. 输出

对 `k = 0..periods`：`t = kT`，输出：

- `results/csv/<scheme>_<operator>_t<k>.csv`：`x, u`；
- `results/dat/<scheme>_<operator>.dat`：每个快照一个 Tecplot zone；
- `results/logs/<scheme>_<operator>.log`：`version`、`precision`、`scheme`、
  `operator`、`mesh_num`、`periods`、`cfl`、`dt_min`、`dt_max`、`steps`、
  `umin`、`umax`、`nan_count`。

本算例不输出误差范数；整数周期后解析解等于初值，因此把快照与初值画在一起
即可看出累积误差。

## 6. 参数与精度

`input.nml` 分为三组：

- `&complex_wave_nml`：`mesh_num`、`periods`、`cfl`、`output_dir`；
- `&scheme_nml`：`scheme` 与 23 个插值标准件的参数；
- `&difference_nml`：`difference_operator`。

默认 `mesh_num = 200`、`periods = 4`、`cfl = 0.3`、`scheme = 'Z5'`、
`difference_operator = 'HCS6'`。本子项目默认双精度 `dp`，`qp` 可构建但结果
以 `dp` 为准。
