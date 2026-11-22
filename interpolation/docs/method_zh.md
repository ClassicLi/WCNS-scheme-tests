# 插值精度测试方法

本文档说明 `interpolation_test` 的测试口径、误差定义、非线性误差和诊断量。
数值公式依据作者已发表的工作（见 README 的“公式来源与验收”）；逐格式的
公式与实现选择见 `docs/formula_audit.md`。

## 1. 半节点插值

在均匀网格 `x_j = x_ini + j h` 上，插值格式由整节点值
`u_{j-3}, ..., u_{j+4}` 计算半节点值

```text
u_{j+1/2} = I(u_j-3, ..., u_j+4).
```

测试区间为 `[-1, 1]`。每个格式只负责半节点插值；完整 WCNS 求解器的通量分裂、
紧致差分和边界处理不在本仓库范围内。

## 2. 误差与收敛阶

对每个网格数 `N`，在半节点上统计

```text
L1     = (1/N) * sum_i |u_exact(x_i) - u_interp(x_i)|
L2     = sqrt( (1/N) * sum_i (u_exact(x_i) - u_interp(x_i))^2 )
Linf   = max_i |u_exact(x_i) - u_interp(x_i)|
```

网格逐级加倍：`N_i = ini_mesh_num * 2^(i-1)`，`i = 1..iter_times`。收敛阶为

```text
p_i = log2( e_{i-1} / e_i ),
```

其中 `e` 为对应的 L1/L2/L∞ 误差。程序同时输出 `En` 和 `En_inf` 的阶数。

## 3. 非线性误差

按文献 [3] 的定义，非线性误差是

```text
En(x) = |u_interp(x) - u_linear^p(x)|,
```

其中 `u_linear^p` 是**与被测格式同阶**的纯线性多项式插值：

- 2 阶：`(u_j + u_{j+1}) / 2`；
- 3 阶：`poly3`；
- 5 阶：`linear5th`；
- 6 阶：`linear6th`。

`En` 是半节点上的平均，`En_inf` 是最大绝对值。默认
`nonlinear_baseline = 'same_order'`；`'poly3'` 仅用于复现旧版固定三阶背景的结果。
线性格式自身的 `En` 记为 0。

## 4. 诊断量

格式例程通过可选的 `diag` 参数返回以下量（不适用的字段保持 0）：

| 字段 | 含义 |
|------|------|
| `beta(1:3)` | 子模板光滑因子 |
| `omega(1:3)` | 非线性权 |
| `tau` | 全局光滑因子（Z5/Z3/MOZ 等） |
| `beta_L`, `beta_N` | 混合框架的主/非线性候选光滑因子 |
| `omega_L`, `omega_N` | 混合框架的主/非线性候选权重 |
| `lambda_h2` | RBF 扰动形状参数 |
| `h_noscale` | MOZ/MOC53 的无量纲量阶调节器 `H_A` |

`delta_omega_mml` 和 `delta_omega_bjw` 按格式实际权数计算：二权格式用
`(1/4, 3/4)`，三权格式用 `(1/16, 10/16, 5/16)`；无权重格式写 0。

## 5. 测试函数

`test_func_type` 选择测试函数：

1. 临界点函数：`exp(0.75 (x-1)) * x^(n+1)`，`n = 0/1/2`。
2. 标准波函数：分段光滑/间断波（`wave_func`）。
3. UDF：`udf_id = 1..5` 对应 `udf_1..udf_5`，`udf_id = 6` 对应无编号 `udf`。
   默认 `udf_id = 3`（`cos(x)`，光滑）。
4. omega 偏移函数：`x^2 exp(0.75(x-1))`，在 `x=0.525` 处加偏移。

所有函数位于 `src/test/standard_function.f90`。

## 6. 格式与权重

### JS3/Z3、JS5/Z5

三点或五点模板拆分为子模板，子模板光滑因子 `beta_k` 和线性权 `C_k`：

- JS 型：`alpha_k = C_k / (beta_k + eps)^p`；
- Z 型：`alpha_k = C_k [1 + (tau / (beta_k + eps))^q]`；
- `omega_k = alpha_k / sum(alpha)`。

三阶 `C=(1/4,3/4)`；五阶 `C=(1/16,10/16,5/16)`。`tau` 为全局光滑因子
（三阶 `|beta_0-beta_1|`，五阶 `|beta_0-beta_2|`）。

### QR

类有理分式插值在半节点上使用

```text
u_{j+1/2} = u_j
  + [3 |u_j-u_{j-1}| (u_{j+1}-u_j) + |u_{j+1}-u_j| (u_j-u_{j-1})]
    / (2 [|u_{j+1}-u_j| + 3 |u_j-u_{j-1}|] + qr_smv).
```

### RBF 扰动

RBF 扰动格式先取同阶线性插值，再用理想形状参数的差分近似加扰动项。
GA/MQ/IQ 的系数和截断阈值见 `src/interpolation/interp_rbf.f90`；
`|lambda_h2| >= rbf_floor_value` 时退化为线性插值。

### MOC322

文献 [3] 的三阶 MOC 方法：主候选为 `poly3`，两个二阶候选为左右两点插值。
光滑因子为

```text
eta_k = |u_{j+k+1} - 2u_{j+k} + u_{j+k-1}| + eps,  k=-1,0,1
r = ((eta_1 - eta_-1) / (eta_1 + eta_-1))^2
beta_00 = r eta_0
beta_L  = |u_j - u_{j-1}| + (1-r)/(4r+eps) eta_0
beta_R  = |u_{j+1} - u_j| + (1-r)/(4r+eps) eta_0
alpha_00 = 1/(beta_00+eps)^Mq, alpha_L = (1/4)/(beta_L+eps)^Mq,
alpha_R = (3/4)/(beta_R+eps)^Mq
u = (alpha_00 poly3 + alpha_L u_L + alpha_R u_R) / sum(alpha)
```

### MOC33

文献 [4] 的三阶混合框架：线性/RBF 候选 + 三点 QR 候选，非线性候选固定 QR。
`beta_00 = kappa eta_0`，`beta_0^H = |u_j-u_{j-1}| + (1-kappa)/2 eta_0`，
`beta_1^H = |u_{j+1}-u_j| + (1-kappa)/2 eta_0`，
`beta_01 = beta_0^H beta_1^H / (beta_0^H + beta_1^H)`，
`omega = alpha^p / sum(alpha^p)`。

### MOC53

MOC53 五阶混合框架（暂未发表，代码中一并提供）：五点线性/RBF 候选 + 三点 QR 候选。
`D_L = u_{j+1}-3u_j+3u_{j-1}-u_{j-2}`，
`D_R = u_{j+2}-3u_{j+1}+3u_j-u_{j-1}`，
`H_A = |(D_L^2-D_R^2)/(D_L^2+D_R^2)|`，
`tau = H_A^Nq |D_L^2-D_R^2|`，
`beta_00 = tau`，
`beta_10 = beta_0 beta_1 beta_2 / (beta_0 beta_1 + beta_1 beta_2 + beta_0 beta_2)`，
权重用 `alpha=1/(beta+eps)`、`omega=alpha^p/sum(alpha^p)`。

### MOZ

文献 [1] 的五阶改进 Z 型加权：用 `tau_5^I` 替换 Z5 的 `tau_5^Z`，
`alpha_k = C_k [1 + (tau_5^I/(beta_k+eps))^Mq]`，
`C=(1/16,10/16,5/16)`。`Nq` 在 `tau_5^I` 内，`Mq` 是 Z 型权指数。

### 线性与 MUSCL

线性模板 3/4/5/6 阶直接使用多项式插值系数；2 阶 MUSCL 使用 minmod 限制器。
它们作为低阶基线和同阶非线性误差背景。

## 7. 参数与默认值

全部格式参数在 `src/interpolation/scheme_params.f90` 中定义，并由
`&scheme_nml` 覆盖。默认值与文献见 `docs/formula_audit.md`。常用参数：

- `scheme`：格式 ID；
- `moc322_M_q`、`moc33_p`、`moc53_p`、`moc53_N_q`、`moz_N_q`、`moz_M_q`；
- `moc33_linear`、`moc53_linear`：`GA/MQ/IQ/LINEAR`；
- `rbf_num_deriv`、`rbf_denom_deriv`、`rbf_floor_value`；
- `js3_eps/js3_p`、`z3_eps/z3_q`、`js5_eps/js5_p`、`z5_eps/z5_q`、`qr_smv`。

## 8. 复现与验收

```bash
make
make test
make smoke
sh scripts/run_results.sh --quick
python3 scripts/plot_convergence.py --quick
```

验收以逻辑正确为准：能编译、能跑完网格加密序列、收敛阶符合设计、
临界点行为合理、无 NaN/崩溃、结果可复现；不要求逐位复现文献表格。

精度：默认四精度（qp）；`make PRECISION=dp` 可构建双精度。提交结果由 qp 生成；
dp 在高分辨率算例（如 160→320 的 GA RBF）上受舍入限制，相关测试阈值按精度放宽。
