# Formula audit

This file records, per retained scheme, the published reference, the formula
actually implemented, and the parameter defaults. The published work [1]-[4]
is the formula source; the legacy chapter-3 program and `double_mach` are
implementation references only.

## References

[1] Li H, Mao M, Yan Z-G, Zhu H, Bai J, Dong W, Deng X. A novel high-order
nonlinear weighting with quasi-rational interpolation for capturing
discontinuities. *Journal of Computational Physics*, 2026, 565: 115207.
DOI: 10.1016/j.jcp.2026.115207.

[2] 黎昊旻, 毛枚良, 白进维, 董维中, 朱华君. 非多项式高阶非线性加权格式研究综述.
*空气动力学学报*, 2026, 44(1): 6-18. DOI: 10.7638/kqdlxxb-2025.0070.

[3] 黎昊旻, 毛枚良, 燕振国, 朱华君, 白进维, 董维中. 基于多精度加权方法改进的
三阶精度HWCNS格式. *空气动力学学报*, 2026, 网络首发.
DOI: 10.7638/kqdlxxb-2025.0243.

[4] 黎昊旻, 毛枚良, 燕振国, 刘庆宗, 白进维, 董维中. 基于类有理分式的强鲁棒性
三阶精度激波捕捉格式研究. *空气动力学学报*, 2026, 网络首发（2026-11-19）.
DOI: 10.7638/kqdlxxb-2026.0127.

## MOC322

- Reference: [3] (multi-precision weighting / third-order HWCNS); the general
  MOC weight form is also reviewed in [2].
- Implemented:
  - `eta(-1)=|u(0)-2u(-1)+u(-2)|+eps`, `eta(0)=|u(1)-2u(0)+u(-1)|+eps`,
    `eta(1)=|u(2)-2u(1)+u(0)|+eps`.
  - `r=((eta(1)-eta(-1))/(eta(1)+eta(-1)))^2`.
  - `beta_00=r*eta(0)`;
    `beta_L=|u(0)-u(-1)|+(1-r)/(4r+eps)*eta(0)`;
    `beta_R=|u(1)-u(0)|+(1-r)/(4r+eps)*eta(0)`.
  - `alpha_00=1/(beta_00+eps)^M_q`,
    `alpha_L=C_L/(beta_L+eps)^M_q`,
    `alpha_R=C_R/(beta_R+eps)^M_q`, `C_L=1/4`, `C_R=3/4`.
  - `omega=alpha/sum(alpha)`;
    `u_half=omega_00*poly3+omega_L*u_L+omega_R*u_R`.
- Defaults: `moc322_eps=1e-80`, `moc322_M_q=2`.
- Indexing note: the source derivation writes `beta_{1,k}=|u_{j+k-1}-u_{j+k}|, k=1,2`.
  The two 2nd-order candidates are the left pair `(u_{j-1},u_j)` and the
  right pair `(u_j,u_{j+1})`; the implemented indicators use
  `|u_j-u_{j-1}|` and `|u_{j+1}-u_j|`, matching the legacy
  `MOC_322_original`.  The `moc322_eps` term in every `alpha` denominator
  keeps flat data finite.

## MOC33

- Reference: [4] (quasi-rational robust third-order shock-capturing); see also [2].
- Implemented:
  - `eta(-1)`, `eta(0)`, `eta(1)` as in MOC322.
  - `kappa=|（eta(1)-eta(-1)）/（eta(1)+eta(-1)）|`.
  - `beta_00=kappa*eta(0)`;
    `beta_H(k)=|u(k)-u(k-1)|+(1-kappa)/2*eta(0)`, `k=0,1`;
    `beta_N=beta_H(0)*beta_H(1)/(beta_H(0)+beta_H(1)+eps)`.
  - `alpha_L=1/(beta_00+eps)`, `alpha_N=1/(beta_N+eps)`;
    `omega_L=alpha_L^p/(alpha_L^p+alpha_N^p)`, `omega_N=1-omega_L`.
  - Linear candidate selected by `moc33_linear` (`GA/MQ/IQ/LINEAR`);
    nonlinear candidate fixed to `quasi_rational`.
- Defaults: `moc33_eps=1e-80`, `moc33_p=2`, `moc33_linear='GA'`.
- The `eps` terms in `beta_N` and the `alpha` denominators keep flat
  data finite.

## MOC53

- Reference: the MOC53 framework is not yet published and is included for
  completeness. Its improved global indicator and five-point sub-template
  indicators follow [1] and [2].
- Implemented:
  - `D_L=u(1)-3u(0)+3u(-1)-u(-2)`,
    `D_R=u(2)-3u(1)+3u(0)-u(-1)`.
  - `num=D_L^2-D_R^2`, `den=D_L^2+D_R^2`,
    `H_A=|num/(den+eps)|`, `tau=H_A^{N_q}|num|`.
  - `beta_00=tau`;
    `beta_10=beta_0 beta_1 beta_2/(beta_0 beta_1+beta_1 beta_2+beta_0 beta_2+eps)`
    (no 16/10/5 factor).
  - `alpha_L=1/(beta_00+eps)`, `alpha_N=1/(beta_10+eps)`;
    `omega_L=alpha_L^p/(alpha_L^p+alpha_N^p)`, `omega_N=1-omega_L`.
  - Five-point candidate selected by `moc53_linear` (`GA/MQ/IQ/LINEAR`);
    nonlinear candidate fixed to `quasi_rational`.
- Defaults: `moc53_eps=1e-80`, `moc53_N_q=4`, `moc53_p=2`,
  `moc53_linear='LINEAR'`.
- `double_mach`'s `beta_N=16 beta_0 beta_1 beta_2/(beta_0 beta_1+
  10 beta_0 beta_2+5 beta_1 beta_2+eps)` is not used; the published work and the
  legacy chapter-3 program use the product form above.
- The `eps` terms in `H_A`, `beta_10`, and the `alpha` denominators keep
  flat data finite.

## MOZ

- Reference: [1]; see also [2].
- Implemented:
  - `tau=tau_5^I` and `h_noscale=H_A` from the same helper as MOC53.
  - `beta_k` are the three 3rd-order sub-template indicators.
  - `alpha_k=C_k[1+(tau/(beta_k+eps))^{M_q}]`,
    `C=(1/16,10/16,5/16)`, `omega_k=alpha_k/sum(alpha)`.
  - `u_half=sum(omega_k*u_k)` with the three 3rd-order sub-stencils.
- Defaults: `moz_eps=1e-80`, `moz_N_q=4` (inside `tau`), `moz_M_q=2`
  (Z-type weight exponent).
- The `eps` term in the beta denominator keeps flat data finite.

## JS3 / Z3

- Reference: [2] (review); the JS/Z weights are standard.
- Implemented:
  - `beta_0=(u_j-u_{j-1})^2`, `beta_1=(u_{j+1}-u_j)^2`,
    `tau=|beta_0-beta_1|`, `C=(1/4,3/4)`.
  - JS3: `alpha_k=C_k/(beta_k+js3_eps)^{js3_p}`.
  - Z3: `alpha_k=C_k[1+(tau/(beta_k+z3_eps))^{z3_q}]`.
  - `omega_k=alpha_k/sum(alpha)`; `u_half=sum(omega_k u_k)`.
- Defaults: `js3_eps=1e-6`, `js3_p=2`, `z3_eps=1e-40`, `z3_q=2`.

## JS5 / Z5

- Reference: [2] (review).
- Implemented:
  - Three-point sub-stencils; `beta_0..beta_2` as in the published formulation;
    `tau=|beta_0-beta_2|`, `C=(1/16,10/16,5/16)`.
  - JS5: `alpha_k=C_k/(beta_k+js5_eps)^{js5_p}`.
  - Z5: `alpha_k=C_k[1+(tau/(beta_k+z5_eps))^{z5_q}]`.
- Defaults: `js5_eps=1e-6`, `js5_p=2`, `z5_eps=1e-40`, `z5_q=2`.

## QR

- Reference: [2], [4]; the quasi-rational form is standard.
- Implemented:
  - `u_half=u(0)+[3|u0-u-1|(u1-u0)+|u1-u0|(u0-u-1)]/
    (2[|u1-u0|+3|u0-u-1|]+qr_smv)`.
- Default: `qr_smv=1e-80`.
- Exact for constant/linear data; third-order on smooth data.

## RBF perturbed schemes

- Reference: [2] (review); implementation coefficients follow
  the simplified `double_mach/interp_rbf.f90`.
- Implemented:
  - Three-point variants start from `linear_3rd_minus1_to_plus1`; five-point
    variants start from `linear_5th_minus2_to_plus2`.
  - `lambda_h2` uses 3rd/1st derivatives (3-point) or 5th/3rd derivatives
    (5-point), with `left/mid/right` selectable.
  - If `|lambda_h2| >= rbf_floor_value`, `lambda_h2=0` and the scheme falls
    back to linear interpolation.
- Defaults: `rbf_floor_value=30`, `rbf_num_deriv='right'`,
  `rbf_denom_deriv='right'`.
- Note: the simplified GA5 variant keeps only the `lambda_h2` perturbation;
  the legacy `lambda_h4` term is not used.

## MUSCL

- Standard second-order minmod interpolation, used as a low-order baseline;
  not a contribution of this work.
- `u_half=u(0)+0.5*minmod(u(1)-u(0), u(0)-u(-1))`.

## Linear stencils

- Reference: standard polynomial interpolation; see [2].
- Implemented: 3/4/5/6-order stencils with the coefficients in
  `src/interpolation/interp_linear.f90`; `poly3` is the symmetric third-order
  stencil and `linear_5th_minus2_to_plus2` / `linear_6th_minus2_to_plus3`
  are the fifth/sixth-order backgrounds.
- Used as low-order baselines and as same-order nonlinear-error backgrounds.

## Test functions

- Reference: [2], [3] for the critical-point and standard
  wave functions.
- Implemented in `src/test/standard_function.f90`:
  - critical-point function `exp(0.75(x-1)) x^{n+1}`, `n=0/1/2`;
  - standard wave function;
  - UDFs: `udf_id=1..5` = `udf_1..udf_5`, `udf_id=6` = the unnumbered `udf`;
    default `udf_id=3` (`cos`);
  - omega-offset function.
- The UDF numbering is an implementation choice made to keep `udf_id=3` as
  the smooth `cos` used by the quick tests and `udf_id=4` as the current
  discontinuous `udf_4`.

## Nonlinear error and delta_omega

- Reference: [3] for the nonlinear-error definition; see also [2].
- Implemented: `En=|u_interp-u_linear^p|`, where `u_linear^p` is the pure
  linear interpolation of the same order. Default `nonlinear_baseline` is
  `'same_order'`; `'poly3'` is a compatibility mode for the old fixed
  third-order background.
- `delta_omega_mml/bjw` use the actual number of weights: two-weight schemes
  use `(1/4,3/4)`, three-weight schemes use `(1/16,10/16,5/16)`.
