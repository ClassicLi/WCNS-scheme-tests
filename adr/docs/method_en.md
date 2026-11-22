# ADR Method

## 1. Approximate dispersion relation

The ADR method of Pirozzoli (2006) estimates the modified wavenumber of a nonlinear
scheme. A single Fourier mode on a periodic domain,

```
u_0(x_j) = cos(n · 2π j / N),   j = 0..N-1,
```

is advanced by a small time step `dt`, giving `u_t`. Fourier transforms of the initial
and final fields are

```
û_n = (1/N) Σ_{j=0}^{N-1} u_j exp(-i n 2π j / N).
```

The modified wavenumber is

```
Phi_n = 2π i ln(û_t / û_0) / (dt · N).
```

The horizontal axis is `varphi_n = 2π n / N`. `Re(Phi)` is dispersion and `Im(Phi)` is
dissipation; for linear advection the exact relation is `Phi = varphi`.

## 2. Grid and modes

- `mesh_num` must be odd; the default is `827`. Even grids may oscillate.
- Mode numbers `n = 0..mesh_num/2` (integer division); `mesh_num/2 + 1` output points.
- An odd grid has no Nyquist mode; the last `varphi` is slightly below `π`.

## 3. Time stepping

A three-stage third-order Runge–Kutta step is used:

```
u^(1) = u^n + dt R(u^n)
u^(2) = 3/4 u^n + 1/4 u^(1) + 1/4 dt R(u^(1))
u^(n+1) = 1/3 u^n + 2/3 u^(2) + 2/3 dt R(u^(2))
```

The periodic halo is refreshed before each stage. `dt` is small enough that the time
discretization error does not contaminate `Phi`.

## 4. Precision and dt

Since `û_t/û_0` differs from 1 by order `dt`, extracting its logarithm loses digits:

- Quad precision (`PRECISION=qp`, default): `dt = 1e-20`, about 13 significant digits.
- Double precision (`PRECISION=dp`): `dt = 1e-10`, about 6 significant digits; noise is
  more visible at high wavenumbers.
- `dt = 0` in `input.nml` selects the precision default.

## 5. Spatial discretizations

### 5.1 WENO reconstruction (`adr_weno`)

The numerical flux `f_{j+1/2}` is reconstructed, then a fixed second-order flux
difference is applied:

```
R_j = -(f_{j+1/2} - f_{j-1/2}) / h.
```

The WENO side has no operator selection.

### 5.2 WCNS interpolation (`adr_wcns`)

Standard interpolation modules produce `u_{j+1/2}`; then a selectable difference
operator is applied:

```
WCNSE2: R_j = -(f_{j+1/2} - f_{j-1/2}) / h
WCNSE4: R_j = -[ 9/8 (f_{j+1/2}-f_{j-1/2}) - 1/24 (f_{j+3/2}-f_{j-3/2}) ] / h
WCNSE6: R_j = -[ 75/64 (f_{j+1/2}-f_{j-1/2}) - 25/384 (f_{j+3/2}-f_{j-3/2})
                + 3/640 (f_{j+5/2}-f_{j-5/2}) ] / h
HCS4:   R_j = -[ 4/3 (f_{j+1/2}-f_{j-1/2}) - 1/6 (u_{j+1}-u_{j-1}) ] / h
HCS6:   R_j = -[ 64/45 (f_{j+1/2}-f_{j-1/2}) - 2/9 (u_{j+1}-u_{j-1})
                + 1/180 (u_{j+2}-u_{j-2}) ] / h
```

`HCS4` and `HCS6` use both half-node fluxes and integer-node values; `HCS6` combines
one pair of half nodes with two pairs of integer nodes.

## 6. Choosing a difference operator

**With the same interpolation, a different difference operator gives a different ADR
curve.** `Phi` describes the combined interpolation-plus-difference operator, not the
interpolation alone. For scheme comparisons, use the operator that the production scheme
actually uses; comparing across operators is not meaningful.

## 7. Boundary and halo

- `adr_weno`: halo = 9, covering the widest `WENO5_DLR` and TENO stencils.
- `adr_wcns`: halo = 7, covering `interpolate_dispatch`'s `u(-3:4)` and `HCS4/HCS6`'s
  integer-node reads.

## 8. Limitations

- Unvalidated variants (`THINC`, `TENO_333_RBF`, `WENO5_DLR`) are kept as experimental
  subroutines; their results are not guaranteed.
- The WENO fifth-order linear reconstruction is the WENO optimal reconstruction in the
  cell-average sense, not a pointwise five-point interpolation; do not assert pointwise
  polynomial exactness for it.
- Legacy `.dat` output is not a reproduction or regression target.

## 9. Reference

- S. Pirozzoli. On the spectral properties of shock-capturing schemes.
  *Journal of Computational Physics*, 2006.
