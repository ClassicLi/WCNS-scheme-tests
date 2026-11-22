# euler1d method notes

## 1. Governing equations

```
∂U/∂t + ∂E/∂x = 0
U = (ρ, ρu, ρE)^T
E = (ρu, ρu² + p, (ρE + p)u)^T
p = (γ - 1)(ρE - 0.5 ρu²),  γ = 1.4
```

## 2. Grid and initial states

Cell-centred grid `x_j = j h`, `j = 0..N`, `h = L/N`. Sod/Lax use `L = 1`;
Shu-Osher uses `L = 10`. The halo is filled by zero-order extrapolation.
The Sod/Lax discontinuity lies between `j = N/2` and `N/2+1`; the
Shu-Osher left segment is `j < N/10` and the right segment starts at
`j = N/10`.

## 3. Characteristic interpolation

For the centre state `U_0 = (ρ, ρu, ρE)`, with `u = ρu/ρ`, `c² = γp/ρ`,
`α = u²/2`, and `β = γ - 1`, the left and right matrices are

```
ml = [ αβ - c²,      -βu,      β ]
     [ αβ - c u,   c - βu,     β ]
     [ αβ + c u,  -c - βu,     β ]

mr = 0.5/c² *
     [ -2,          1,                 1                ]
     [ -2u,         u + c,             u - c            ]
     [ -2α,         α + cu + c²/β,     α - cu + c²/β    ]
```

They are inverses. Each characteristic field `q_p = ml(p,:)·U` is
interpolated to the half node by the shared `interpolate_dispatch`, then
`U_{j+1/2} = mr·q_half`.

If the reconstructed `ρ` or `p` is below 0.1% of the stencil minimum, the
centre state is used instead.

## 4. Steger-Warming split

```
eps2 = 0.12²,  sig = 0.5/γ
λ1 = u,  λ4 = u + c,  λ5 = u - c
λk± = (λk ± sqrt(λk² + eps2))/2
diss = -sig (2λ1 - λ4 - λ5)
disp = sig (λ4 - λ5) c
fp(1) = diss + λ1
fp(2) = u diss + u λ1 + disp
fp(3) = (α + c²/β) diss + (α + p/(ρβ)) λ1 + u disp
E_{j+1/2} = ρ_L fp + ρ_R fn
```

`fp` uses the left state with positive eigenvalues and `fn` the right state
with negative eigenvalues.

## 5. Difference operators

The nodal residual `K_j = -dE/dx` uses `WCNSE2`, `WCNSE4`, `WCNSE6`,
`HCS4`, or `HCS6` (formulas in the Chinese section, identical coefficients).
`AUTO` pairs third-order interpolation with `HCS4` and fifth-order with
`HCS6`.

## 6. Time advance, boundary, and positivity guard

Four-stage Runge–Kutta (RK4) with

```
dt = CFL h / max_j (|u_j| + sqrt(p_j/ρ_j))
```

and default `CFL = 0.1`. The last step is clipped to `Tend`. The boundary is
seven layers of zero-order extrapolation. After every RK stage, a cell with
`ρ < 1e-8` or `p < 1e-8` is replaced by `(U_{j-1} + U_{j+1})/2` from the old
state. The log records how many corrections fired.

## 7. Exact solution

The exact solver is used for Sod and Lax only. Newton iteration finds the
star pressure:

```
f(p*, p_L, ρ_L) + f(p*, p_R, ρ_R) = u_L - u_R
```

with the shock-side and rarefaction-side forms of `f`. The iteration uses a
precision-dependent tolerance and falls back to bisection if needed. It then
gives `u*`, `ρ*_L`, `ρ*_R`, and the wave speeds.

Samples must use `ξ = (x - 0.5)/t`, because the initial discontinuity is at
`x = 0.5`. The four wave patterns (two shocks, two rarefactions including
vacuum, left shock/right rarefaction, left rarefaction/right shock) are
selected from `p*` relative to `p_L` and `p_R`. At `t = 0` the solver returns
the left or right initial state directly.
