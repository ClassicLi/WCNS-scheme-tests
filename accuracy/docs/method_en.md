# Accuracy tests: methods and implementation

## 1. One-dimensional periodic convection

Solve

```
u_t + u_x = 0,   x ∈ [-1, 1],
```

with periodic boundaries and `Tend = 2`. The grid is `x_j = -1 + j h`,
`h = 2/N`, `N = nbase·2^(k-1)`, `k = 1..levels`; defaults are `nbase = 10`
and `levels = 8`.

### 1.1 Initial profiles

Seven continuous periodic profiles are provided by `c1d_profile(x)`:
`sin2 = sin²(πx)`, `sin = sin(πx)`, `sin3 = sin³(πx)`,
`sin_phase = sin(πx - sin(πx)/π)`, `gauss = exp(-20x²)`,
`sqrt = sqrt(1 + 20x²)`, and `rational = 1/(1 + 20x²)`.
The exact solution is `u0(x - t)` with periodic wrapping. `c1d_initial` and
`c1d_exact` both call `c1d_profile`, so the two formulas cannot drift apart.

### 1.2 Spatial discretisation

Half-node interpolation:

```
E_{j+1/2} = Interp(u_{j-3..j+4}),   j = -3..N+2.
```

`Interp` is selected by `scheme_dispatch` from the shared `interpolation/`
standards. The difference operator is

```
K_j = difference_apply(op, h, E_{j-3/2..j+5/2}, u_{j-2..j+2}).
```

For `HCS4`/`HCS6`, the integer-node argument is the flux; in the scalar case
it is simply `u_j`.

### 1.3 Time integration

Four-stage fourth-order Runge–Kutta:

```
k1 = R(u^n)
k2 = R(u^n + dt/2 k1)
k3 = R(u^n + dt/2 k2)
k4 = R(u^n + dt k3)
u^{n+1} = u^n + dt (k1 + 2 k2 + 2 k3 + k4)/6
```

with `dt = CFL·h/max|u|`, default `CFL = 0.1`, and the last step truncated to
`Tend`.

### 1.4 Errors and orders

```
L1  = Σ|u_num - u_ex|/(N+1)
L2  = sqrt(Σ(u_num - u_ex)²/(N+1))
Linf = max|u_num - u_ex|
p_k = log(e_{k-1}/e_k)/log(2)
```

### 1.5 Halo and boundary

Interpolation at `j = -3` needs `u(-6)`, and at `j = N+2` it needs `u(N+6)`,
so `halo = 7`. The periodic halo is `u(-i) = u(N-i)`, `u(N+i) = u(i)` for
`i = 1..7`. The old program's `U(-5:N+5)` read out of bounds; this
implementation does not.

## 2. Two-dimensional Euler

Periodic Euler on a square, `Tend = 1`, `CFL = 0.1`,
`N = M = nbase·2^(k-1)`, defaults `nbase = 40` and `levels = 4`.

### 2.1 `vortex`

Domain `[0,10]²`, `h = 10/N`. Isentropic vortex with
`EKC = 0.3`, `RC = 0.5`, `RFA = 0.35`, `γ = 1.4`, centred at `(5,5)`.
The initial and exact fields use the same formulas:

```
x̄ = j h - 5 - t,   ȳ = k h - 5
r = sqrt(x̄² + ȳ²)
T = 1 - (γ-1) EKC² exp(2 RFA (1 - (r/RC)²)) / (4 RFA γ)
ρ = T^(1/(γ-1))
u = 1 + (EKC/RC) exp(RFA (1 - (r/RC)²)) ȳ
w = -(EKC/RC) exp(RFA (1 - (r/RC)²)) x̄
p = ρ^γ
```

Density errors are taken over the `(N+1)×(M+1)` nodes.

### 2.2 `sine`

Domain `[0,2π]²`, `h = 2π/N`, `u = w = 1`, `p = 1`,
`ρ = ρ0 + A sin(x + y - 2t)`, defaults `ρ0 = 1`, `A = 0.2`.
The initial state is at `t = 0`; the old AIA data corresponds to
`ρ0 = 2`, `A = 1`.

### 2.3 Flux assembly

The reference pipeline is kept:

1. `E5_X`/`E5_Z` project a seven-point conservative stencil to characteristic
   variables using the centre state, interpolate every characteristic field
   with `interpolate_dispatch`, and project back to a half-node state;
2. `InterFlux` reconstructs left and right states and applies the modified
   Steger–Warming split (`eps2 = 0.12²`, `sig = 0.5/γ`) to obtain
   `Eflux/Gflux`;
3. `e2d_rhs` applies the difference operator to the half-node and integer-node
   fluxes:

```
K_j = difference_apply(op, h, Eflux_{j-3/2..j+5/2}, Eint_{j-2..j+2})
    + difference_apply(op, h, Gflux_{k-3/2..k+5/2}, Gint_{k-2..k+2})
```

In the 2D Euler case the integer-node argument of `HCS4`/`HCS6` is the
**physical flux** `Eint/Gint`, not the conserved state `Ubar`. The two are
identical only in the scalar 1D test.

The positivity fallback is unchanged: if the reconstructed density or pressure
is below 0.1% of the stencil minimum, the centre state is used.

### 2.4 Halo and boundary

`mt = 7`, covering the farthest read of the interpolation and operators at
`j = -3..N+2` (where `j-3 = -6`). `e2d_boundary` fills the periodic halo.

### 2.5 Errors and orders

The same L1/L2/Linf definitions as in 1D, with denominator `(N+1)(M+1)`, on
the density `ρ`.

## 3. Observed orders

- `Z5 + HCS6` reaches about fifth order on the smooth sine case (e.g. 5.11,
  5.00).
- `Z3`/`QR + HCS4` show the known critical-point reduction (L1 about 2,
  Linf about 1–2), consistent with the `interpolation/` Z3/QR results.
- On the strong isentropic vortex, the `Z5` L1 order rises from about 3 to
  about 4.3 at `n=320`, consistent with the old `js5.dat` (3–4). The clean
  fifth-order old `linear5.dat` is reproduced value-for-value by
  `LINEAR5_M2P2 + HCS6` to about `1e-15` relative, confirming the 2D port.

## 4. Choosing a case

- To see the designed order, use `sine` with `Z5 + HCS6`.
- To see strong-gradient/nonlinear-weight behaviour, use `vortex` with `Z5`
  and compare against `LINEAR5_M2P2`.
- `AUTO` pairs only by order; an explicit mismatch prints a warning but still
  runs.
