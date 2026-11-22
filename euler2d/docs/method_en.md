# euler2d method notes

This document describes the discretization used by `euler2d`. Every formula
and constant is taken from the code in `src/euler2d/`, `src/common/`,
`src/wcns/`, and `../interpolation/`.

## 1. Governing equations

The two-dimensional Euler equations are

```
∂U/∂t + ∂E/∂x + ∂G/∂z = S,
U = (ρ, ρu, ρw, E)^T,
Eflux = (ρu, ρu² + p, ρuw, (E + p)u)^T,
Gflux = (ρw, ρuw, ρw² + p, (E + p)w)^T,
p = β (E - 0.5 ρ (u² + w²)),   β = γ - 1.
```

The default is `γ = 1.4`; `rayleigh_taylor` uses `γ = 5/3`. Only RT has a
non-zero source, `S = (0, 0, ρ, ρw)^T`, added by `e2d_source` after the
spatial residual is assembled.

## 2. Grid, halo, and coordinate origins

Let `n = nx`, `m = nz`, and `h = 1/m`. The cell centres are

```
x_j = x_origin + j h,   j = 0..n,
z_k = z_origin + k h,   k = 0..m.
```

The physical domain is
`[x_origin, x_origin + nx/nz] × [z_origin, z_origin + 1]`, with
`nx = nz · ratio_num / ratio_den` fixed by the case.

The conservative state `ubar` and the integer-node fluxes `eint/gint` are
allocated as `(-mt:n+mt, -mt:m+mt, 4)`; the half-node fluxes `eflux` and
`gflux` are allocated as `(-3:n+2, 0:m, 4)` and `(0:n, -3:m+2, 4)`.
The halo is **`mt = 7`**:

- the eight-point characteristic stencil at `j = -3` reads
  `ubar(j-3:j+4) = ubar(-6:1)`, and at `j = n+2` it reads
  `ubar(n-1:n+6)`, so the state must extend to `-6` and `n+6`;
- the difference operators read integer-node values
  `eint(-2:2)` at `j = 0` and `eint(n-2:n+2)` at `j = n`.

The original thesis programs use `MT = 5`, which reads out of bounds at
those two locations. This implementation uses `mt = 7` and checks it under
`gfortran -fcheck=all`.

Coordinate origins follow the initial states:

| Case | `x_origin` | `z_origin` |
|------|------------|------------|
| `double_mach` | 0 | 0 |
| `kelvin_helmholtz` | 0 | 0 |
| `rayleigh_taylor` | 0 | 0 |
| `riemann1` | 0 | 0 |
| `riemann2` | 0 | 0 |
| `shock_bubble` | -0.1 | -0.5 |
| `shock_vortex` | -0.25 | -0.5 |

Output coordinates always use the initial-state origin. The old
`shock_vortex` `export` routine wrote coordinates from `(0,0)`, which is
inconsistent with the `Initial` routine's `(-0.25,-0.5)`; the initial state
is authoritative here.

## 3. Characteristic interpolation

At each half node `j+1/2` (or `k+1/2`), the eight-point stencil
`U_{-3}, ..., U_4` is projected to characteristic variables. With the
centre state (stencil index 0),

```
u = (ρu)/ρ,   w = (ρw)/ρ,
α = 0.5 (u² + w²),
p = β (E - α ρ),
c² = γ p / ρ,   c = sqrt(c²).
```

The x-direction left and right characteristic matrices are

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

The z-direction matrices are

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

The code stores the right matrices column-wise as `mr(i, p)`; the
entry-by-entry formulas are in `eigen_matrix_x` and `eigen_matrix_z` in
`e2d_flux.f90`. The four characteristic fields

```
q_p = ml(p,1) U_1 + ml(p,2) U_2 + ml(p,3) U_3 + ml(p,4) U_4
```

are interpolated independently to `q_half(p)` by `interpolate_dispatch`
from `../interpolation/`, then projected back:

```
U_half(1:4) = mr(:,1) q_half(1) + ... + mr(:,4) q_half(4).
```

This is done twice: the left state uses the stencil `ubar(j-3:j+4)` and the
right state uses `ubar(j+4:j-3:-1)`. The resulting `U_L` and `U_R` enter the
flux split.

### First correction layer: interpolation positivity fallback

If the reconstructed half-node `ρ_half` or `p_half` is no larger than 0.1%
of the five-point stencil minimum over `j-2..j+2`,

```
ρ_half <= 0.001 min(ρ_{j-2..j+2})   or
p_half <= 0.001 min(p_{j-2..j+2}),
```

then `U_half` is replaced by the centre state `U_0` and
`interp_corrections` is incremented. This is applied in both x and z.

## 4. Modified Steger–Warming split

Let `eps² = 0.12²` and `sig = 1/(2γ)`. In the x direction,

```
λ1 = u,   λ4 = u + c,   λ5 = u - c,
λ± = (λ ± sqrt(λ² + eps²))/2,
diss = -sig (2λ1 - λ4 - λ5),
disp =  sig (λ4 - λ5) c.
```

The left state uses `λ+` and gives `fp`:

```
fp(1) = diss + λ1
fp(2) = u (diss + λ1) + disp
fp(3) = w (diss + λ1)
fp(4) = (α + c²/β) diss + (α + p/(ρβ)) λ1 + u disp
```

The right state uses `λ-` with the same expressions to give `fn`. The
half-node flux is

```
Eflux = ρ_L fp + ρ_R fn.
```

In the z direction, `λ1 = w` and `disp` moves to the third component:

```
fp(1) = diss + λ1
fp(2) = u (diss + λ1)
fp(3) = w (diss + λ1) + disp
fp(4) = (α + c²/β) diss + (α + p/(ρβ)) λ1 + w disp
```

and `Gflux = ρ_L fp + ρ_R fn`. This differs from the original
Steger–Warming `E^+ + E^-` bookkeeping: the code divides the characteristic
velocity contributions by density and regroups them, so `fp/fn` are
per-unit-density flux factors.

## 5. Difference operators and AUTO pairing

The half-node fluxes are passed as `f_{j-5/2} ... f_{j+5/2}`, and the
integer-node values as `u_{j-2} ... u_{j+2}`. The five nodal residuals are

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

`difference_apply` returns `-dE/dx`. The two-dimensional residual is

```
rhs(j,k,:) = difference_apply(op, h, eflux(j-3:j+2,k,:),
                              eint(j-2:j+2,k,:))
           + difference_apply(op, h, gflux(j,k-3:k+2,:),
                              gint(j,k-2:k+2,:)).
```

`AUTO` pairs by interpolation order:

| Interpolation order | Operator |
|---------------------|----------|
| 2 | `WCNSE2` |
| 3 | `HCS4` |
| 4 | `WCNSE4` |
| 5 | `HCS6` |
| 6 | `WCNSE6` |

The default `QR` has order 3, so `AUTO` selects `HCS4`, i.e.
`QR + AUTO(HCS4)`. An explicit operator that disagrees with the order prints
a warning and continues.

## 6. Time advance, CFL, and source terms

Time integration is fixed four-stage, fourth-order Runge–Kutta:

```
u^(1) = u^n + 0.5 dt k1
u^(2) = u^n + 0.5 dt k2
u^(3) = u^n + dt k3
u^(n+1) = u^n + dt (k1 + 2k2 + 2k3 + k4)/6
k_r = RHS(u^(r), t_n + dt)
```

The boundary halo is refreshed before every RK substage, and the residual
includes the `e2d_source` term. The RT gravity source acts on the z
momentum and energy equations and is independent of the difference
operator.

The time step is

```
dt = CFL · h / max_{j,k} ( sqrt(u² + w²) + sqrt(γ p / ρ) ).
```

The second square root is the physical sound speed
`c = sqrt(gama*p/rho)` (that is, `sqrt(γ p / ρ)`, where the code variable
`gama` holds the ratio of specific heats), consistent with the
characteristic decomposition and the Steger–Warming split. If the next step would pass the target time, it is clipped to
`t_end - t`. The default `CFL = 0.3`; the result scripts use `0.1`.

## 7. Second correction layer: cell positivity guard

After every RK substage, the interior cells `0..n, 0..m` are scanned.
Density is tested first, then pressure:

```
if ρ < 1e-8:
    U(j,k,:) = (U0(j-1,k,:) + U0(j+1,k,:)
              + U0(j,k-1,:) + U0(j,k+1,:)) / 4
    increment the count
else:
    compute u = ρu/ρ, w = ρw/ρ, p = β(E - 0.5ρ(u²+w²))
    if p < 1e-8:
        U(j,k,:) = four-neighbour average from U0
        increment the count
```

The neighbours come from `U0`, the state saved at the start of the time
step, not from the just-updated state. `cell_corrections` records the
number of replacements. Testing density first avoids dividing by a zero or
negative density.

## 8. Boundary types for the seven cases

`e2d_boundary(case, t)` dispatches the halo fill. Every case fills
`mt = 7` layers.

| Case | Left | Right | Bottom | Top |
|------|------|-------|--------|-----|
| `double_mach` | fixed inflow `ρ=8, u=7.14471, w=-4.125, p=116.5` | zero-order extrapolation | fixed inflow for `x < 1/6`; otherwise mirror reflection with `w` sign flipped | moving shock `x_s(t) = 1/6 + (1+20t)/√3`; `x > x_s` uses `(1.4,0,0,1)`, `x <= x_s` uses `(8,7.14471,-4.125,116.5)` |
| `kelvin_helmholtz` | mirror copy | mirror copy | mirror copy | mirror copy |
| `rayleigh_taylor` | reflection (`u` sign flipped) | reflection (`u` sign flipped) | fixed `ρ=2, p=1, u=w=0` | fixed `ρ=1, p=2.5, u=w=0` |
| `riemann1` | mirror copy | mirror copy | mirror copy | mirror copy |
| `riemann2` | mirror copy | mirror copy | mirror copy | mirror copy |
| `shock_bubble` | fixed inflow `ρ=3.81, u=2.85, w=0, p=10` | zero-order extrapolation | reflection (`w` sign flipped) | reflection (`w` sign flipped) |
| `shock_vortex` | mirror copy | mirror copy | mirror copy | mirror copy |

“Mirror copy” means the source program's `U(:, -p) = U(:, p-1)` and
`U(:, m+p) = U(:, m-p+1)` (with the analogous x relations), not a strict
periodic shift. At the `double_mach` bottom, `x - 1/6 < 1e-8` uses the
inflow state; otherwise the lower halo is set to
`U(j,-p) = U(j,p-1)` with the `w` sign flipped.

## 9. Output and logs

Tecplot ASCII point format with variables `x, z, rho`. The primary file is
`<case>_<nz>x<nx>_t1.dat`; `shock_bubble` writes `_t1.dat` at `tend2` and
`_t2.dat` at `tend`. `write_sample = .true.` additionally writes `_t0.dat`.

The log `<case>.log` records the case, interpolation scheme, operator, mesh,
CFL, final time, `gamma`, step count, `interp_corrections`,
`cell_corrections`, `nan_count`, and CPU time. The non-finite count scans
physical cells only, not the halo.
