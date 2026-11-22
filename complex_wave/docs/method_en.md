# Complex-Wave Time-Advancement Test: Method and Implementation

## 1. Governing equation and grid

We solve the one-dimensional periodic linear advection equation

```
u_t + u_x = 0,   x ∈ [-1, 1]
```

with unit wave speed, so one wave period is

```
T = 2 / 1 = 2
```

The grid is

```
x_j = -1 + j h,   h = 2 / N,   j = 0..N
```

The periodic boundary repeats the endpoint: `u(N) = u(0)`. Half-node interpolation
needs six halo points on each side; this implementation uses `halo = 7`:

```
u(-i) = u(N-i),   u(N+i) = u(i),   i = 1..halo
```

## 2. Initial condition

The profile has four features, with the same parameters as the reference program:
`δ = 0.005`, `z = -0.7`, `a = 0.5`.

```
G(x, z) = exp(-log(2) (x - z)^2 / 0.03^2)
F(x, a) = sqrt(max(0, 1 - (10 (x - a))^2))

k = N/10 .. 2N/10   :  u_k = (G(x_k, z-δ) + 4 G(x_k, z) + G(x_k, z+δ)) / 6
k = 3N/10 .. 4N/10  :  u_k = 1
k = 5N/10 .. 6N/10  :  u_k = max(0, 1 - |10 x_k - 1|)
k = 7N/10 .. 8N/10  :  u_k = (F(x_k, a-δ) + 4 F(x_k, a) + F(x_k, a+δ)) / 6
otherwise           :  u_k = 0
```

The `max(0, ...)` in the triangle acts only when `N` is not a multiple of 10: the
segment endpoints can then fall slightly outside `[0, 0.2]`, where the reference
formula would return a small negative value.

## 3. Spatial discretization

### 3.1 Half-node interpolation

```
E_{j+1/2} = Interp(u_{j-3..j+4}),   j = -3 .. N+2
```

`Interp` is a WCNS standard from `../interpolation/src/interpolation/`, selected by
`scheme` in `&scheme_nml`. This test uses `JS3`, `JS5`, `Z3`, `Z5`, and `QR`.

### 3.2 Difference operators

The spatial right-hand side is `K_j = -dE/dx`. Define

```
ΔE_m = E_{j+m+1/2} - E_{j-m-1/2},   m = 0, 1, 2
Δu_m = u_{j+m} - u_{j-m},           m = 1, 2
```

The five selectable operators are

```
WCNSE2:  K_j = -ΔE_0 / h
WCNSE4:  K_j = -[(9/8) ΔE_0 - (1/24) ΔE_1] / h
WCNSE6:  K_j = -[(75/64) ΔE_0 - (25/384) ΔE_1 + (3/640) ΔE_2] / h
HCS4:    K_j = -[(4/3) ΔE_0 - (1/6) Δu_1] / h
HCS6:    K_j = -[(64/45) ΔE_0 - (2/9) Δu_1 + (1/180) Δu_2] / h
```

The `HCS6` coefficients follow from matching the Taylor series through sixth order.
Using

```
ΔE_0/h = u' + (h^2/24) u''' + (h^4/1920) u^(5) + ...
Δu_1/h = 2u' + (h^2/3) u''' + (h^4/60) u^(5) + ...
Δu_2/h = 4u' + (8h^2/3) u''' + (8h^4/15) u^(5) + ...
```

and requiring the coefficients of `u'`, `u'''`, and `u^(5)` to be `1`, `0`, and `0`
gives

```
a + 2b + 4c = 1
a/24 + b/3 + 8c/3 = 0
a/1920 + b/60 + 8c/15 = 0
```

with the solution `a = 64/45`, `b = -2/9`, `c = 1/180`.

## 4. Time advancement

The classical four-stage fourth-order Runge–Kutta scheme is

```
k1 = R(u^n)
k2 = R(u^n + dt/2 k1)
k3 = R(u^n + dt/2 k2)
k4 = R(u^n + dt   k3)
u^{n+1} = u^n + dt (k1 + 2 k2 + 2 k3 + k4) / 6
```

where `R` is the spatial right-hand side from Section 3. As in the reference
program, the step is

```
dt = CFL h / max_j |u_j|,   CFL = 0.3 by default
```

The last step of each period is clipped to `target - t`, so every snapshot lands
exactly on `t = kT`. If the maximum amplitude is zero or non-finite, the program
stops with an error instead of looping with `dt = 0` or NaN.

## 5. Output

For `k = 0..periods` at `t = kT` the driver writes

- `results/csv/<scheme>_<operator>_t<k>.csv`: `x, u`;
- `results/dat/<scheme>_<operator>.dat`: one Tecplot zone per snapshot;
- `results/logs/<scheme>_<operator>.log`: `version`, `precision`, `scheme`,
  `operator`, `mesh_num`, `periods`, `cfl`, `dt_min`, `dt_max`, `steps`, `umin`,
  `umax`, and `nan_count`.

No error norms are written. After an integer number of periods the exact solution
equals the initial profile, so plotting each snapshot next to the initial profile
shows the accumulated error directly.

## 6. Parameters and precision

`input.nml` has three groups:

- `&complex_wave_nml`: `mesh_num`, `periods`, `cfl`, `output_dir`;
- `&scheme_nml`: `scheme` and the 23 interpolation-standard parameters;
- `&difference_nml`: `difference_operator`.

Defaults are `mesh_num = 200`, `periods = 4`, `cfl = 0.3`, `scheme = 'Z5'`, and
`difference_operator = 'HCS6'`. This subproject defaults to double precision `dp`;
`qp` remains buildable, but results are generated in `dp`.
