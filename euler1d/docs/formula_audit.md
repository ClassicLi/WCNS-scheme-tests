# 与旧程序的差异（formula audit）

参考程序：

- `D:\numerical_scheme\5_chap234_detail\code_simplify\1d`（WSL
  `/mnt/d/numerical_scheme/5_chap234_detail/code_simplify/1d`）；
- `D:\numerical_scheme\3_chap3\code\shock_tube_exact_solution`（WSL
  `/mnt/d/numerical_scheme/3_chap3/code/shock_tube_exact_solution`）。

## 1. 算例与格式

| 项目 | 旧程序 | 本子项目 |
|------|--------|----------|
| Sod | `Interpolation_character_modi` + `muscl_van_albada` + `Difference_E2` | 统一特征插值；默认 `QR + HCS4` |
| Lax | `Interpolation_character` + `MOC53` + `Difference_E4_Hybrid` | 默认 `QR + HCS4`；`MOC53` 不参与结果矩阵 |
| Shu-Osher | `Interpolation_character` + `MOC53` + `Difference_hcs6` | 默认 `QR + HCS4`；参考解 `QR + HCS4`，`N = 1600` |
| Sedov、TT | 有独立目录 | 不纳入本子项目 |

旧程序的格式是逐例硬编码；本子项目统一为 `input.nml` 可选，默认 `QR`。

## 2. 算子等价关系

| 旧名称 | 本子项目 |
|--------|----------|
| `Difference_E2`（`a = 1`） | `WCNSE2` |
| `Difference_E4_Hybrid`（`a = 4/3, b = -1/6`） | `HCS4` |
| `Difference_hcs6`（`a = 64/45, b = -2/9, c = 1/180`） | `HCS6` |
| `Difference_E4` | 无直接对应；可用 `WCNSE4` 对照 |
| `Difference_E6`（`75/64, -25/384, 3/640`） | `WCNSE6` |

本子项目的 `difference_operators.f90` 与 `adr/`、`complex_wave/`、`accuracy/`
同源。

## 3. 特征投影

旧 Sod 目录的 `Interpolation_character_modi` 与 Lax/Shu-Osher 的
`Interpolation_character` 是两套不同的投影/反投影系数。本子项目统一采用
`accuracy/src/euler2d/e2d_flux.f90` 的一维退化版（去掉横向速度后的
3×3 矩阵），并对 `mr·ml = I` 做单元测试。

因此 Sod 的新结果不会与旧图逐点一致；这是设计中选择“统一可选”的预期取舍。

## 4. 边界与 halo

旧程序 `U(-5:N+5)`，但插值在 `j = -3` 需要 `u(-6)`、在 `j = N+2` 需要
`u(N+6)`。本子项目使用 `halo = 7`，并在 `-fcheck=all` 下检查越界。

## 5. 时间步

旧 `TimeStep` 与仓库 `e2d_timestep` 都用 `sqrt(p/ρ)` 作为时间步速度：

```
dt = CFL h / max_j (|u_j| + sqrt(p_j/ρ_j))
```

特征分解与 Steger-Warming 仍用物理声速 `c = sqrt(γp/ρ)`。两者不混用；
文档中明确写出，避免把“有效 CFL 略大”误认为实现错误。

## 6. 精确解程序

旧程序的问题与本子项目的修正：

| 旧行为 | 修正 |
|--------|------|
| `global_var` 中写死 `TYPE`、`MESH_NUM`，输出 `Exact.dat` | `input.nml` 选择 `case_1d`、`mesh_num`、`tend`、`output_dir` |
| 双膨胀非真空分支使用未赋值的 `r_star_L/r_star_R` | 改为 `ρ_L(p*/p_L)^(1/γ)`、`ρ_R(p*/p_R)^(1/γ)` |
| 真空分支使用未赋值的 `c_star` | 用 `u_L + 2c_L/(γ-1)`、`u_R - 2c_R/(γ-1)` 作为波尾 |
| `ZONE` 标题硬编码 `"Sod"` | 由算例名生成 |
| list-directed 输出折行 | 固定 `es24.16` CSV，四列 `x,rho,u,p` |
| 波坐标 `k = x - 0.5` 只隐含在代码里 | 明确 `ξ = (x - 0.5)/t`，并有端点/波位置测试 |
| 牛顿迭代靠 `1e-40` 容差退出 | 精度相关相对容差 + 最多 200 次 + 二分兜底 |

## 7. 结果

本子项目不计算误差范数，只输出数值解与精确解/参考解的 CSV、日志和图。
Shu-Osher 的参考解是同求解器在 `N = 1600` 的数值解，不是解析解；
参考解的网格、格式与算子记录在 `results/logs/reference_shu_osher.log`。
