# Interpolation accuracy test method

This document describes the test conventions, error definitions, nonlinear
error, and diagnostics used by `interpolation_test`. Formulas follow the
author's published work (see "Formula source and acceptance" in the README).
Per-scheme formula and implementation choices are recorded in
`docs/formula_audit.md`.

## 1. Half-point interpolation

On a uniform grid `x_j = x_ini + j h`, a scheme computes the half-point value

```text
u_{j+1/2} = I(u_{j-3}, ..., u_{j+4})
```

from integer-node values. The test domain is `[-1, 1]`. The repository tests
the interpolation operator only; the flux splitting, compact differencing,
and boundary treatment of a full WCNS solver are out of scope.

## 2. Errors and convergence order

For each grid size `N`, the half-point errors are

```text
L1   = (1/N) * sum_i |u_exact(x_i) - u_interp(x_i)|
L2   = sqrt( (1/N) * sum_i (u_exact(x_i) - u_interp(x_i))^2 )
Linf = max_i |u_exact(x_i) - u_interp(x_i)|
```

The grid is doubled: `N_i = ini_mesh_num * 2^(i-1)`, `i = 1..iter_times`.
The convergence order is

```text
p_i = log2( e_{i-1} / e_i ),
```

where `e` is the corresponding L1/L2/L∞ error. The program also reports
orders for `En` and `En_inf`.

## 3. Nonlinear error

Following the published definition [3], the nonlinear error is

```text
En(x) = |u_interp(x) - u_linear^p(x)|,
```

where `u_linear^p` is the pure linear polynomial interpolation of the
**same order** as the tested scheme:

- order 2: `(u_j + u_{j+1}) / 2`;
- order 3: `poly3`;
- order 5: `linear5th`;
- order 6: `linear6th`.

`En` is the half-point mean and `En_inf` is the maximum absolute value.
The default is `nonlinear_baseline = 'same_order'`; `'poly3'` is a
compatibility mode that reproduces the old fixed third-order background.
Linear schemes have `En = 0`.

## 4. Diagnostics

Scheme routines return the following quantities through the optional `diag`
argument (unused fields stay zero):

| Field | Meaning |
|-------|---------|
| `beta(1:3)` | sub-template smoothness indicators |
| `omega(1:3)` | nonlinear weights |
| `tau` | global smoothness indicator |
| `beta_L`, `beta_N` | main/nonlinear candidate indicators |
| `omega_L`, `omega_N` | main/nonlinear candidate weights |
| `lambda_h2` | RBF perturbation shape parameter |
| `h_noscale` | MOZ/MOC53 normalization factor `H_A` |

`delta_omega_mml` and `delta_omega_bjw` use the actual number of weights:
two-weight schemes use `(1/4, 3/4)`, three-weight schemes use
`(1/16, 10/16, 5/16)`, and schemes without weights write zero.

## 5. Test functions

`test_func_type` selects the test function:

1. Critical-point function: `exp(0.75 (x-1)) * x^(n+1)`, `n = 0/1/2`.
2. Standard wave function: piecewise smooth/discontinuous wave (`wave_func`).
3. UDF: `udf_id = 1..5` selects `udf_1..udf_5`; `udf_id = 6` selects the
   unnumbered `udf`. Default `udf_id = 3` (`cos(x)`, smooth).
4. Omega-offset function: `x^2 exp(0.75(x-1))` with an offset at `x=0.525`.

All functions live in `src/test/standard_function.f90`.

## 6. Schemes and weights

### JS3/Z3 and JS5/Z5

The three- or five-point stencil is split into sub-stencils. With indicators
`beta_k` and linear weights `C_k`:

- JS type: `alpha_k = C_k / (beta_k + eps)^p`;
- Z type: `alpha_k = C_k [1 + (tau / (beta_k + eps))^q]`;
- `omega_k = alpha_k / sum(alpha)`.

For the third order `C=(1/4,3/4)`; for the fifth order
`C=(1/16,10/16,5/16)`. `tau` is the global indicator (`|beta_0-beta_1|`
for the third order, `|beta_0-beta_2|` for the fifth order).

### QR

Quasi-rational interpolation uses

```text
u_{j+1/2} = u_j
  + [3 |u_j-u_{j-1}| (u_{j+1}-u_j) + |u_{j+1}-u_j| (u_j-u_{j-1})]
    / (2 [|u_{j+1}-u_j| + 3 |u_j-u_{j-1}|] + qr_smv).
```

### RBF perturbation

RBF perturbed schemes start from the same-order linear interpolation and add
a perturbation based on a finite-difference approximation of the ideal shape
parameter. The GA/MQ/IQ coefficients and cutoff are in
`src/interpolation/interp_rbf.f90`; when
`|lambda_h2| >= rbf_floor_value` the scheme falls back to linear
interpolation.

### MOC322

The third-order MOC method [3] uses `poly3` as the main
candidate and the left/right two-point interpolations as the two second-order
candidates:

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

The third-order hybrid framework [4] combines a
linear/RBF candidate with a three-point QR candidate; the nonlinear
candidate is fixed to QR. With `kappa` the normalized indicator ratio,

```text
beta_00 = kappa eta_0
beta_0^H = |u_j-u_{j-1}| + (1-kappa)/2 eta_0
beta_1^H = |u_{j+1}-u_j| + (1-kappa)/2 eta_0
beta_01 = beta_0^H beta_1^H / (beta_0^H + beta_1^H)
omega = alpha^p / sum(alpha^p)
```

### MOC53

The MOC53 fifth-order hybrid framework (unpublished; included for
completeness) combines a five-point linear/RBF candidate with a three-point
QR candidate. With

```text
D_L = u_{j+1}-3u_j+3u_{j-1}-u_{j-2}
D_R = u_{j+2}-3u_{j+1}+3u_j-u_{j-1}
H_A = |(D_L^2-D_R^2)/(D_L^2+D_R^2)|
tau = H_A^Nq |D_L^2-D_R^2|
beta_00 = tau
beta_10 = beta_0 beta_1 beta_2 / (beta_0 beta_1 + beta_1 beta_2 + beta_0 beta_2)
```

the weights use `alpha=1/(beta+eps)` and `omega=alpha^p/sum(alpha^p)`.

### MOZ

The improved fifth-order Z-type weighting [1] replaces the
Z5 global indicator with `tau_5^I`:

```text
alpha_k = C_k [1 + (tau_5^I/(beta_k+eps))^Mq],
C = (1/16, 10/16, 5/16).
```

`Nq` is the exponent inside `tau_5^I`; `Mq` is the Z-type weight exponent.

### Linear and MUSCL

The 3/4/5/6-order linear stencils use polynomial interpolation coefficients;
the second-order MUSCL scheme uses the minmod limiter. They serve as
low-order baselines and as the same-order nonlinear-error backgrounds.

## 7. Parameters and defaults

All scheme parameters are defined in
`src/interpolation/scheme_params.f90` and can be overridden through
`&scheme_nml`. Defaults and references are in
`docs/formula_audit.md`. Common parameters:

- `scheme`: scheme ID;
- `moc322_M_q`, `moc33_p`, `moc53_p`, `moc53_N_q`, `moz_N_q`, `moz_M_q`;
- `moc33_linear`, `moc53_linear`: `GA/MQ/IQ/LINEAR`;
- `rbf_num_deriv`, `rbf_denom_deriv`, `rbf_floor_value`;
- `js3_eps/js3_p`, `z3_eps/z3_q`, `js5_eps/js5_p`, `z5_eps/z5_q`, `qr_smv`.

## 8. Reproducibility and acceptance

```bash
make
make test
make smoke
sh scripts/run_results.sh --quick
python3 scripts/plot_convergence.py --quick
```

Acceptance is logical correctness: the code builds, completes the refinement
sequence, shows the designed orders, behaves sensibly at critical points,
produces no NaN/crash, and its results are reproducible. Bitwise reproduction
of the published tables is not required.

Precision: the default is quad precision (qp); `make PRECISION=dp` builds double
precision. The committed results use qp; dp is roundoff-limited on the finest meshes
(e.g. GA RBF at 160→320), so those test thresholds are relaxed for dp.
