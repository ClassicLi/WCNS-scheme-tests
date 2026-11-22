# 公式与参数审计 / Formula and Parameter Audit

本文件记录 ADR 两套程序的格式来源、默认参数与已知改动。代码注释统一使用英文。

## 1. WENO 格式

| Scheme ID | 类别 | 阶数 | 关键参数 | 实现 |
|-----------|------|------|----------|------|
| `WENO5_LINEAR` | 线性最优重构 | 5 | — | `weno5_linear_reconstruct` |
| `WENO5_JS` | JS 型权 | 5 | `js5_eps=1e-6`, `js5_p=2` | `weno5_js` |
| `WENO5_M` | M 型权 | 5 | `js5_eps`, `js5_p` | `weno5_m` |
| `WENO5_Z` | Z 型权 | 5 | `z5_eps=1e-40`, `z5_q=2` | `weno5_z` |
| `WENO_C` | WENO-C | 5 | `weno_c_c4=0.2`, `weno_c_c5=0.1`, `weno_c_eps=1e-6` | `weno_c_deheilan` |
| `WENO3_LINEAR` | 线性 | 3 | — | `weno3_linear` |
| `WENO3_JS` | JS 型权 | 3 | `js3_eps=1e-6`, `js3_p=2` | `weno3_js` |
| `WENO3_TRIG` | 三角重构 | 3 | — | `weno3_trig` |
| `WENO3_LLR` | LLR | 3 | `llr_eps=1e-20` | `weno3_llr` |
| `WENO3_GA_RBF_DEFINE` | RBF-GA（ε 由数据定义） | 3 | — | `weno3_ga_rbf_eps_define` |
| `WENO3_GA_RBF_UNDEFINE` | RBF-GA（固定 ε） | 3 | `rbf_eps=75` | `weno3_ga_rbf_eps_undefine` |
| `BAEZA_WENO3` | Baeza 型 | 3 | 内部 `vareps=1e-40` | `baeza_weno3` |
| `THINC` | THINC | 3 | `thinc_beta=2` | `get_thinc` |
| `TENO_333` | TENO | 3 | `teno_q=6`, `teno_ct=1e-4` | `teno_original_333` |
| `TENO_333_RBF` | TENO + RBF | 3 | `teno_q`, `teno_ct` | `teno_original_333_rbf` |
| `TENO_MOD_SEMI_ENO` | 改进 TENO | 3 | `teno_q`, `teno_ct` | `teno_mody_semi_eno` |
| `TENO_MOD_SEMI_ENO_RBF` | 改进 TENO + RBF | 3 | `teno_q`, `teno_ct` | `teno_mody_semi_eno_rbf_ga` |
| `RBF5_LINEAR` | RBF 扰动线性 | 5 | 内部 ε 由数据定义 | `rbf_5_linear` |
| `WENO5_DLR` | DLR | 5 | 内部 `tol=1.0`, `vareps=1e-8` | `weno5_dlr` |

说明：

- WENO 的 JS/Z 默认参数与插值测试标准件一致：JS `eps=1e-6, p=2`；Z `eps=1e-40, q=2`。
- `WENO5_LINEAR` 是 WENO 最优线性重构（单元平均意义），不是点值五点插值。
- `WENO5_DLR` 为实验子程序，见第 5 节。

## 2. WCNS 格式

WCNS 侧直接复用 `../interpolation/src/interpolation/` 的 23 个标准件：

`MUSCL`、`JS3`、`Z3`、`JS5`、`Z5`、`QR`、`MOC322`、`MOC33`、`MOC53`、`MOZ`、
7 个线性格式、6 个 RBF 格式。

公式、默认参数与参数审计见 `../interpolation/docs/formula_audit.md`。
`MOC53` 使用论文定稿形式（`N_q=4`、`p=2`、默认线性候选 `LINEAR`）。

## 3. 差分套子

`flux_half(1:6)` 映射 `f_{j-5/2}..f_{j+5/2}`；`u_node(1:5)` 映射 `u_{j-2}..u_{j+2}`。
所有算子返回 `R_j = -(...)/h`。

旧 `adr_wcns` 使用 halo=5，实际会读到 `u(-6:nx+6)`；本实现改用 halo=7，修正了该越界。

```
WCNSE2: -(f_{j+1/2} - f_{j-1/2}) / h
WCNSE4: -[ 9/8 (f_{j+1/2}-f_{j-1/2}) - 1/24 (f_{j+3/2}-f_{j-3/2}) ] / h
WCNSE6: -[ 75/64 (f_{j+1/2}-f_{j-1/2}) - 25/384 (f_{j+3/2}-f_{j-3/2})
          + 3/640 (f_{j+5/2}-f_{j-5/2}) ] / h
HCS4:   -[ 4/3 (f_{j+1/2}-f_{j-1/2}) - 1/6 (u_{j+1}-u_{j-1}) ] / h
HCS6:   -[ 64/45 (f_{j+1/2}-f_{j-1/2}) - 2/9 (u_{j+1}-u_{j-1})
          + 1/180 (u_{j+2}-u_{j-2}) ] / h
```

## 4. 已修正的笔误

- `TENO_MOD_SEMI_ENO_RBF`：原代码 `4.0d0*(-1)` 改为 `4.0d0*u(-1)`。
- `WENO_C`：原代码在 `do loop=1,2` 中写 `alpha_4 = ...`（缺下标），
  改为 `alpha_4(loop) = ...`；并移除调试用 `write "error!!!"` 检查。
- `adr_io.write_csv`：输出改为逗号分隔的 CSV。

## 5. 未验证格式

以下子程序保留、可编译、可选，但不进入默认结果矩阵：

- `TENO_333_RBF`：`eps2h2` 未初始化，结果依赖未定义内存。
- `THINC`：原始表达式量纲存疑；在对称模板上可能除零并产生非有限值。
- `WENO5_DLR`：原始文件已注明五阶 DLR 有问题；在 ADR 数据上可能产生 `NaN`。

## 6. 旧版参数 vs 标准件

旧 `adr_wcns` 目录与 `interpolation/` 标准件的默认参数不同。本仓库以标准件为准：

| 项目 | 旧 adr_wcns | 本仓库标准件 |
|------|-------------|--------------|
| MOC53 N_q / M_q | 8 | 4 |
| MOC53 beta_N | 16 b1 b2 b3 / (b1b2 + 10b1b3 + 5b2b3 + eps) | b1b2b3 / (b1b2 + b1b3 + b2b3 + eps) |
| MOC53 p | — | 2 |
| MOC53 线性候选 | linear5th | LINEAR |
| RBF floor value | 1000 | 30 |
| JS5 p | 1 | 2 |
| WCNS halo | 5 | 7 |

旧 `.dat` 不作为复现目标；由参数差异导致的曲线变化不是 bug。

## 7. 旧 `.dat` 输出

旧版程序目录中的 `.dat` 文件不作为本程序的复现或回归目标。旧文件的网格、
活跃调用和边界处理与当前实现不同；验收以单元测试、解析校验和无 `NaN` 为准。
