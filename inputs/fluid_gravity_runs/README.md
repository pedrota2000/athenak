# 2D conformal relativistic turbulence for fluid/gravity

Data-generation runs whose output is meant to be fed into the fluid/gravity metric
coefficients: the four-velocity `u^mu` and the fluid-frame energy density `e` (the
eigenvalue of the stress tensor, `T = diag(-e, p, p, p)` in the rest frame), stored
pointwise, plus the pointwise residuals of the inviscid relativistic conservation law
`d_mu T^{mu nu} = 0`.

Everything below is implemented with what AthenaK already has, plus one new problem
generator and two small (backwards-compatible) additions to the turbulence driver; see
[Code changes](#code-changes).

## Files

| file | what it is |
| --- | --- |
| `2d_conformal_turb_256.athinput` | 256^2 shakedown run (minutes; use it to calibrate `dedt`) |
| `2d_conformal_turb_1024.athinput` | 1024^2 production run |
| `analysis_window_decay.athinput` | restart override: dense dumps, **all sources off — this is where the data set comes from** |
| `analysis_window_forced.athinput` | restart override: dense dumps, driving still on (diagnostic use only, see the EFE section) |
| `../../src/pgen/fg_turb2d.cpp` | problem generator (hot uniform conformal fluid + diagnostics) |
| `../../scripts/fg_residuals.py` | post-processor: `u^mu`, `e`, and the residuals, cell by cell |

## The conformal equation of state

AthenaK's relativistic hydro carries an ideal gas, `p = (Gamma-1) e_int`, with primitives
`(rho, u^x, u^y, u^z, e_int)` — note the velocity slots hold the **four**-velocity — and
`T^{mu nu} = (rho + e_int + p) u^mu u^nu + p eta^{mu nu}`.  The eigenvalue of that tensor
is therefore

```
e = rho + e_int ,      p = (Gamma - 1) e_int
T^mu_mu = -e + (d-1) p = -rho + [(d-1)(Gamma-1) - 1] e_int
```

Choosing `Gamma = d/(d-1)` kills the `e_int` term **identically**, so

```
Gamma = 3/2  (d = 3, i.e. 2+1)  ->  p = e_int/2 ,  T^mu_mu = -rho  exactly
p/e = 1 / [(d-1) (1 + rho/e_int)]
```

There is no exactly conformal EOS in the code, and there cannot be one while a rest-mass
density is evolved: the trace anomaly *is* `rho`.  So the runs are set up "hot",
`e_int/rho = 10^4`, which puts `p/e` within 0.01% of `1/(d-1)`.  Because `rho` is just an
advected scalar in this limit (`d_mu (rho u^mu) = 0`) and is written out with every dump,
the anomaly is known pointwise: it can be monitored (`anom` column of the `.hst` file,
`anomaly` field in the post-processor) or removed altogether (`fg_residuals.py
--traceless`, which uses the strictly traceless `e = (d-1) p`).

**Use the full `T^{mu nu}`, with `e = dens + eint`.**  That is the tensor the code conserves
exactly, so `d_mu T^{mu nu} = 0` holds in the continuum and the residual measures only
numerical error; `e` is then the fluid/gravity energy density to 1 part in `10^4`.  The
`--traceless` alternative gives a strictly traceless tensor but injects `-rho a^nu` into the
residual (see below), which would compromise the main diagnostic — keep it as an optional
cross-check, not the default.

Raising `temp0` shrinks the anomaly for free: `rho/e_int = 1/(2 temp0)`, and the conformal
sound speed is `sqrt(Gamma-1) = 0.707` whatever `temp0` is, so there is no timestep penalty.
The shipped value `temp0 = 5000` gives `10^-4`; `5e4` would give `10^-5`.  The limit is
conditioning of the c2p inversion, which works with `q = tau/D ~ 2 temp0` — comfortable to
`~1e4-1e5` in double precision.  Note `temp0` is the *ideal-gas* temperature `p/rho`, not
the CFT temperature `T = e^{1/d} ~ 1` that sets the horizon radius; the two differ by a
factor of thousands here, and confusing them would give a badly wrong horizon field.

**`d = 3` vs `d = 4`.** The input files use `d = 3` (`gamma = 1.5`, `p = e/2`): the fluid
lives in the same 2+1 space-time as the simulation, which is what "conformal EOS in 2D,
`p = e/(d-1)` with `d` the space-time dimensionality" says.  If instead the intended
object is a 3+1 conformal fluid in a translationally invariant (2D) flow — which is what
`T = diag(-e, p, p, p)` with four entries suggests — then set `<hydro>/gamma = 1.3333333`
and `<problem>/d_spacetime = 4.0`, and pass `--gamma 1.3333333 --d 4` to the
post-processor.  Nothing else changes.

## Why the driving is set up the way it is

The cascade in 2D is inverse, so energy has to be injected at **small** scales and removed
at **large** ones:

* **Injection at high k.**  `<turb_driving>/nlow,nhigh` band-limit the force.  The
  production run drives `22 <= |k| L/2pi <= 26`, i.e. at 1/24 of the box, leaving `k < 22`
  for the inverse cascade and `k > 26` for the direct enstrophy cascade.

* **Use `driving_type = 2` (added for these runs).**  It builds the force from a random
  stream function, `F = (-d_y psi, d_x psi, 0)`, so it is exactly divergence-free on the
  2D slice, has no out-of-plane component (`u^z` stays *identically* zero, verified), is
  isotropic, and only visits `k_z = 0` modes.
  Neither pre-existing driver is usable here:
  - `driving_type = 0` builds its solenoidal projection in 3D and gives every mode a
    random `F_z`; on a 2D mesh (`nx3 = 1`) the `k_z != 0` modes collapse onto the `z = 0`
    slice, the projection stops cancelling, and the force acquires both a compressive part
    and an out-of-plane component that drives `u^z` (turning the run into 2.5D).
  - `driving_type = 1` *is* in-plane and divergence-free on a 2D slice, but it solves for
    `F_y` from `F_x` by dividing by `k_y`, so `|F_y| ~ (k_x/k_y) |F_x|` and the small-`k_y`
    modes are wildly over-driven.  Measured on a 128^2 test: `rms(u^y)/rms(u^x) = 3.6`.
    With `driving_type = 2` the same test gives 0.94.  It is also `nhigh-nlow+1` times more
    expensive, because its mode selection requires `nlow <= n_z <= nhigh` and so realises
    every in-plane mode once per `n_z` (the 128^2 test ran 2.7x faster after the switch).

* **No large-scale sink — deliberately.**  `rel_cooling` is left off.  This is the one
  point where the literature is explicit: Westernacher-Schneider (arXiv:1710.04264) writes
  that *"the regime of validity of the fluid-gravity duality is that of an arbitrarily high
  Reynolds number and no large-scale friction"*, and chooses the forcing configuration that
  gives `E(k) ~ k^-2` on exactly that ground.  A drag term is a source at **every** scale,
  including the large-scale band the metric is reconstructed from, so it is far more
  damaging than the band-limited drive (see the EFE section below).  Instead this is a
  *transient* quasi-steady state, as in that paper: no sink, run until the inverse cascade
  approaches the box scale, analyse snapshots before the energy piles up there.
  In 2D this also means the heating worry is much smaller than a 3D intuition suggests —
  most of the injected energy travels up-scale into large-scale *kinetic* energy rather
  than being dissipated into heat — so `T` drifts slowly.  Check `p/e` and `anom` in the
  `.hst` file rather than assuming it, and watch `Ekin` for the box-scale pileup.

* **`dedt` is not the physical injection rate — it overshoots by `3*temp0`.**  The driver
  normalizes its amplitude against the *rest-mass* density `D = rho W` (`m0 = 0.5<D F^2>dt`)
  but does its work against the true relativistic momentum `S_i = (e+p)W^2 v_i`
  (`m1 = <S_i F^i>`), so the actual injection is
  `eps = dedt * rho/(e+p) = dedt/(3 temp0)`.  With no large-scale sink the energy
  *accumulates* rather than balancing a flux, `(e+p) v^2 ~ eps t`, so

  ```
  dedt = 3 * temp0 * (e+p) * v_target^2 / t_sample
  ```

  Measured at 128^2: `dedt = 1440`, `temp0 = 5000` gives `eps = 0.096`, predicting
  `E_kin = 0.384` at `t = 4` against a measured `(e+p)*usq = 0.32` — 83%, the rest
  dissipated.  Note `dedt` scales with `temp0`: raising the temperature to shrink the trace
  anomaly weakens the forcing at fixed `dedt`, and both files' values assume
  `temp0 = 5000`.  Read `usq` from the `.hst` file and rescale as `v^2`.

* **No physical viscosity.**  Do *not* set `<hydro>/viscosity` here:
  `src/diffusion/viscosity.cpp` builds the Newtonian stress tensor out of `w0`, whose
  velocity slots are the four-velocity in SR — it is not a relativistic Navier–Stokes or
  Israel–Stewart term.  These runs are inviscid; dissipation is numerical (wenoz + hlle +
  FOFC), which is also exactly what makes the measured residual meaningful.

## Build and run

The pgen is selected at build time, so use a build directory of its own — this leaves the
existing 3D `PROBLEM=turb` build untouched:

```bash
module load modules/2.4-20250724 cuda/12.5.1
cmake -S . -B build_fg -D PROBLEM=fg_turb2d -D Athena_ENABLE_MPI=ON \
      -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_HOPPER90=ON
cmake --build build_fg -j
```

Single-node GPU run (same environment fixes as the 3D runs — without them UCX spins on
InfiniBand registration and never reaches the task list):

```bash
ulimit -l unlimited
export UCX_TLS=self,sm,cuda_copy,cuda_ipc
export OMPI_MCA_btl=^openib
./build_fg/src/athena -i inputs/fluid_gravity_runs/2d_conformal_turb_256.athinput -d run256
```

Then the production run, and finally the dense-cadence analysis window from a restart:

```bash
./build_fg/src/athena -r run1024/rst/fgturb2d.00010.rst \
    -i inputs/fluid_gravity_runs/analysis_window_decay.athinput -d run1024_window
```

## Output and the four stored variables

`variable = hydro_w` writes, per cell, `dens, velx, vely, velz, eint` =
`rho, u^x, u^y, u^z, e_int`.  From these:

```
u^t = sqrt(1 + u^x^2 + u^y^2 + u^z^2)      (u^z = 0 in 2D)
p   = (Gamma-1) e_int
e   = rho + e_int                          (eigenvalue of T; = (d-1)p up to the anomaly)
```

so the four fluid/gravity variables are exactly recoverable, with no derived-variable
support needed in the code.  `scripts/fg_residuals.py` does this and writes them out
alongside the residuals; `--fields-only` writes just the fields (no time stencil needed).

The `*.user.hst` file carries the run-level diagnostics enrolled by the pgen (the usual
conserved-variable history still goes to `*.hydro.hst`):
`rho, e, pres, p/e, anom, usq, W, Ekin, Etot, omsq` — all true volume averages.  `p/e`
should sit at `1/(d-1)`, `anom` at the trace anomaly, `usq` at the target `v^2`, and
`omsq` (enstrophy) should saturate while the large-scale energy grows: that growth is the
inverse cascade.

## Residuals of the conservation law

```bash
python3 scripts/fg_residuals.py 'run1024_window/bin/fgturb2d.hydro_w.*.bin' \
        -o run1024_window/residuals --gamma 1.5 --d 3 --order 4 --report
```

The script forms `T^{mu nu} = (e+p) u^mu u^nu + p eta^{mu nu}`, differences it, and writes
per-snapshot `.npz` files containing `ut, ux, uy, e, p, rho, anomaly`, the residuals
`res_t, res_x, res_y` and the baryon residual `res_N`, plus `norm_*` (the sum of the
absolute values of the individual terms in each equation) so that `|res|/norm` is a
dimensionless measure of how well the ideal law holds.  Verified against an analytically
differentiated smooth field to 3e-7 relative.  It also writes the first-order geometry
(`theta`, `acc_*`, `sigma_*`, `knudsen`), the holographic quantities (`temp`, `r_plus`,
`eta`, `visc_*`) and shell spectra (`spec_k`, `spec_E`, `spec_res_*`) — see the EFE section
for what each is for.

Three things to keep in mind:

1. **The residual is only zero if no source terms act.**  With the driver on,
   `d_mu T^{mu nu}` equals the driving force pointwise, and that field is not written to
   disk.  Take the data from `analysis_window_decay.athinput`, which switches every source
   off (`dedt = 0` now short-circuits the whole driver, so the fluid is not touched at all);
   the residual is then purely the numerical/dissipative part — the quantity that has to be
   small for the gradient expansion to be trustworthy.  Keep each window well under one
   forcing-scale eddy time and take an ensemble of them — see
   [the power law](#does-switching-the-force-off-destroy-the-power-law) below.
2. **Time derivatives need dense dumps.**  AthenaK writes a dump on the first step at or
   after each target time, so dump times are *not* exactly equispaced; the script uses a
   3-point non-uniform stencil and reads the exact times from the headers.  `dt = 0.001` in
   the window files is ~3–4 hydro timesteps at 1024^2, giving 200 snapshots
   (~8 GB) over `Delta t = 0.2`.
3. **Use `--order 4`.**  With `wenoz` reconstruction the scheme's own truncation error is
   high order; a 2nd-order analysis stencil would dominate the answer and you would be
   measuring the post-processor rather than the simulation.

## Are these profiles solutions of the Einstein equations?

**Forced ones are not — not at any scale.**  A forced conformal fluid is simply not in the
image of the fluid/gravity map into pure Einstein gravity with a cosmological constant.
Use the *unforced* window for the data set; the forced run is a device for manufacturing a
turbulent state, not the data set itself.

The reason is structural.  In fluid/gravity the bulk metric is a derivative expansion
`g = sum_n eps^n g^(n)[u, T]`, and the *constraint* equations of the bulk Einstein system
are identically the conservation of the boundary stress tensor: at order `n` the metric
solves EFE to `O(eps^{n+1})` if and only if
`d_mu (T_(0) + ... + T_(n))^{mu nu} = 0`.  arXiv:1710.04264 states the same condition
itself — its metric "will solve Einstein's equations with arbitrary accuracy in the perfect
fluid limit **if `u^a` and `T` evolve according to conformal hydrodynamics on the
boundary**".  A driven fluid does not: it evolves according to conformal hydrodynamics
*plus an external force*, so the stated sufficient condition fails.  Put the other way
round: **the residual `R^nu` this post-processor computes IS the constraint violation of the
metric you would build from `(u^mu, e)`** — that is the useful reframing, and it is not
merely a code-quality diagnostic.

**A band-limited force does not fix this, and here is why not.**  It is true, and worth
measuring, that the violation carries no spectral power outside the forcing band: the
constraint is linear in the transform, so `i k_mu (T^{mu nu})^(k) = f^nu(k) = 0` for `k`
outside the support of `f`.  On the 128^2 test, 100.0% of the residual power sits in the
forcing band and 0.00% below it.  But that is a statement about the *spectral distribution*
of the violation, not about the equations holding anywhere:

* the constraint is a **pointwise** condition and Einstein's equations are **nonlinear**, so
  they do not decompose mode by mode.  "No violation power at `k < k_f`" does not give a
  solution on any subset of the domain — pointwise, `d_mu T^{mu nu}(x) = f^nu(x) != 0` in
  every cell, because the force is band-limited in `k`, not in `x`;
* and you cannot repair it by low-pass filtering, because `T^{mu nu}` is **quadratic** in
  `u`.  Filtering to `k < k_f` gives `d_mu <T^{mu nu}> = -d_mu (subfilter stress)`, not
  zero: you have traded an external force for an unclosed Reynolds-stress term.  This is
  the LES closure problem, not a loophole.

To have a genuine bulk dual of a *driven* fluid you must add whatever sources the force —
a slowly-varying boundary metric, or an external boundary gauge field acting on a charged
fluid, `d_mu T^{mu nu} = F^{nu la} J_la` (which is the form the AthenaK drive happens to
take, since it is exactly orthogonal to `u`: it adds `rho a^i dt` to the momentum and
`rho (a_i u^i/u^t) dt` to the energy, so `f_mu u^mu = 0` identically).  But that is
Einstein–**Maxwell** with a charged brane, a different bulk theory with its own first-order
transport, and the required `F^{mu nu}` would have to be reverse-engineered from the
stochastic force.  It is not the map anyone is actually using.

**What does survive the objection.**  Quantities that depend only on an inertial-range
scaling exponent are largely insensitive to how the turbulence was maintained, because
forced and freely-evolving turbulence share their inertial-range statistics.  The fractal
dimension in arXiv:1710.04264 is of that kind — it is read off the roughness of
`r_+ = 4 pi T/3`, and the two values it reports (2.584 and 2.645) track the two spectral
exponents (`k^-2` and `k^-5/3`), with the forcing serving only to select which.  So its
*numbers* plausibly stand as statements about the geometry of a `k^-beta` random field.
What its construction does not establish is the stronger reading its language invites, that
these snapshots are horizons of solutions of the Einstein equations.  That gap is real and
the paper does not address it.

**So the protocol here is: force to build the state, then switch the force off and take the
data from the unforced evolution** (`analysis_window_decay.athinput`, where `dedt = 0` now
short-circuits the driver entirely so the fluid is not touched at all).  In 2D this is
nearly free, which is the happy accident that makes the whole thing practical: freely
evolving 2D turbulence conserves energy almost exactly in the inviscid limit — enstrophy is
bounded and decaying, so `eps = nu <omega^2> -> 0` — and only enstrophy decays.  An unforced
window therefore stays statistically steady for many eddy times, unlike in 3D, where you
would be fighting a rapidly decaying flow.  During that window the data satisfies
`d_mu T^{mu nu} = 0` and `d_mu (rho u^mu) = 0` up to numerical error, the snapshots are in
the image of the map at perfect-fluid order, and the residual measures what it is supposed
to measure.

There are two further obstructions, which apply even then.

**1. Numerical dissipation is still a source.**  In the unforced window the ideal residual
is no longer an external force, but it is not zero either: it is the scheme's dissipation,
`R^nu ≈ -d_mu T^{mu nu}_num`.  That is a violation of the perfect-fluid constraints too —
just a small, broadband, and *measurable* one rather than an `O(1)` imposed one.  Measured
on the 128^2 test with all sources off: `rms(R^t) = 1.5e-3` and `<|R^t|/N_t> = 0.022`,
against `0.98` and `0.73` with the thermostat on.  This is the error budget of the data set,
and it is exactly what `fg_residuals.py` is for.

**2. Order of hydrodynamics vs order of the metric.**  Perfect-fluid data yields a metric
that solves EFE only to `O(eps)`.  This is a legitimate published choice rather than a
defect — arXiv:1710.04264 works at perfect-fluid order throughout, with the metric
`ds^2 = -2 u_a dx^a dr - r^2 (1 - r_+^3/r^3) u_a u_b dx^a dx^b + r^2 (eta_ab + u_a u_b) dx^a dx^b`
(`R_AdS = 1`), and cites agreement with full GR simulations at the 1% level for particular
boundary fluid data.  But it does bound what the data set can support.  To be a *first-order* solution the data must satisfy
relativistic Navier–Stokes with the **holographic** transport coefficients — `eta = s/4pi`
and `zeta = 0` (forced by conformality).  These runs are inviscid: the dissipation is
numerical, of an unknown magnitude and not of Navier–Stokes form.  So strictly this is a
*zeroth-order* fluid/gravity data set.  The ideal residual is measuring that numerical
dissipation, `R^nu ≈ -d_mu T^{mu nu}_num`, whereas first-order consistency would require
`R^nu ≈ d_mu (2 eta sigma^{mu nu})` with `eta = s/4pi`.

That is now directly testable: regress the decay-window `R^nu` against
`d_mu (2 sigma^{mu nu})` to get the scheme's *effective* `eta`, and compare it to `s/4pi`.
It will not match.  The honest responses are (a) restrict the analysis to the band of
scales where both terms are negligible compared with the ideal ones, or (b) give AthenaK a
genuine relativistic (Israel–Stewart) viscous sector with `eta = s/4pi` — which it does not
have, and which is a real project, not a patch.  `<hydro>/viscosity` is *not* that thing.

**3. The trace anomaly.**  `T^mu_mu = -rho != 0` means the boundary theory is not exactly a
CFT (it carries a relevant deformation), so the bulk is not exactly asymptotically AdS.
It is `10^-4` here and known pointwise.  Removing it with `--traceless` is *not* free: the
conformal tensor differs from the conserved one by a dust term, `T_full - T_conf = rho
u^mu u^nu`, so `d_mu T_conf^{mu nu} = -rho a^nu` (the force that was accelerating the dust
you deleted; `u^nu d_mu(rho u^mu)` vanishes by baryon conservation).  That is `~7e-5` of the
`(e+p) a^nu` term at `temp0 = 5000`, but it is a source where there was none.  Better to
shrink `rho/e_int` with `temp0` and keep the exactly conserved tensor.

**And the validity of the expansion itself.**  Fluid/gravity requires
`Kn = |gradient| / T << 1`.  In turbulence `Kn` grows toward the dissipation scale and is
`O(1)` there by construction, so the reconstructed metric is meaningful only over the band
of scales where `Kn` is small — which is a statement to be *checked*, not assumed.  The
post-processor therefore also writes the first-order data the metric is built from:

| field | meaning |
| --- | --- |
| `theta` | `d_mu u^mu`, the expansion |
| `acc_t, acc_x, acc_y` | `a^mu = u^nu d_nu u^mu`, the four-acceleration |
| `sigma_tt ... sigma_yy` (`sigma_zz` for `d = 4`) | shear tensor `sigma^{mu nu}`, the object `g^(1)` is built from |
| `sigma_norm` | `sqrt(sigma^{mu nu} sigma_{mu nu})` |
| `kn_reduced` | `sigma_norm / e^{1/d}` |
| `res_dot_u`, `res_perp_norm` | the residual split into fluid-frame heating and pure force |

`sigma^{mu nu}` satisfies `u_mu sigma^{mu nu} = 0` (to 1e-16) and `sigma^mu_mu = 0` (to
8e-9, i.e. to stencil accuracy) on a smooth test field.

**The temperature normalization.**  A conformal fluid has `e ~ T^d`, but the constant is
not fixed by the hydrodynamics — it encodes the central charge (equivalently the horizon
scale / `G_{d+1}`) in code units.  arXiv:1710.04264 fixes it by fiat for `d = 3`:
`T = rho^(1/3)` with `rho` the energy density, `R_AdS = 1`, and the horizon at
`r_+ = 4 pi T / 3`.  The post-processor now uses that convention by default
(`--tnorm C`, `T = C e^(1/d)`, default `C = 1`), which makes `knudsen`, `r_plus`, `eta`
and the first-order stress `visc_*` quantitative rather than reduced.  If your bulk-side
conventions differ, set `--tnorm` accordingly — everything downstream scales with it.

## Does switching the force off destroy the power law?

Partly, on a timescale you have to respect — this is a real constraint, not a detail.

What is safe: **energy** is nearly conserved in freely evolving 2D turbulence, because
enstrophy is bounded and decaying so `eps = nu <omega^2> -> 0` in the inviscid limit.  The
flow does not die.

What is not safe: **energy conservation does not freeze the spectrum.**  Nonlinear transfer
continues without any flux to sustain it, so the shape drifts — enstrophy decays and the
`k > k_f` range steepens quickly, coherent vortices form and reorganize the field, and the
`k^-5/3` / `k^-2` range at `k < k_f` becomes a *relic* rather than a maintained constant-flux
range.  A constant-flux inertial range requires injection, full stop.

The drift time is the eddy time of the energy-containing scale.  With `v_l ~ l^{(beta-1)/2}`
for `E(k) ~ k^-beta`:

| | `L/l_f = 24` (these runs) | `L/l_f = 85` (arXiv:1710.04264) |
| --- | --- | --- |
| `E ~ k^-2` | `tau_box/tau_f = 4.9` | 9.2 |
| `E ~ k^-5/3` | 8.3 | 19.3 |

So the large-scale spectrum survives a *few* forcing-scale eddy times — and the
intermediate scales that carry the scaling range have shorter eddy times than that, so they
degrade first.  My earlier "many eddy times" was wrong; what is nearly steady is the total
energy, not the spectral shape.

**The way out is an ensemble of short windows, not one long window.**  The budget is
comfortable: at 1024^2 one forcing-scale eddy time `t_f = l_f/v ≈ 0.14` is ~474 hydro
timesteps, and the analysis cadence is `dt_out = 0.001` ≈ 4 timesteps.  So

* one window of `0.10 = 0.7 t_f` yields **100 snapshots** whose spectra have barely moved;
* the spin-up writes a restart every `2.0 ≈ 14 t_f`, so 10 of them are effectively
  independent realizations;
* releasing each gives **~1000 constraint-satisfying snapshots** carrying the forced-steady
  spectrum.

This is the same ensemble structure arXiv:1710.04264 uses (20 flows), with the difference
that each member is released before being sampled.  And it is checkable rather than assumed:
`spec_E` is written for every snapshot, so plot it across a window and set the length from
the measurement.

**A prior question, though: do you actually need the power law?**  The fluid/gravity map
requires only that (i) the data solve the hydro equations and (ii) `Kn << 1`.  It does not
care whether the spectrum is a power law.  The exponent matters only if the *question* is
about inertial-range scaling — as in arXiv:1710.04264, where the fractal dimension is
essentially `D(beta)`.  If the goal is instead "pointwise metric coefficients from turbulent
data", freely decaying 2D turbulence is fully turbulent and perfectly adequate, and the
window can be as long as you like.  Worth settling before spending GPU hours.

**Why you cannot split the difference by weakening the force.**  In a steady state the
injection must balance the flux, `eps = w v_f^3 / l_f`, and the force needed to supply it is
`f ~ eps/v_f = w v_f^2 / l_f` — the same size as the nonlinear term at the forcing scale.
That is what makes it the injection scale.  So the constraint violation is `O(1)` relative
to the ideal terms there, for any forcing amplitude; the measured `<|R^x|/N_x> ≈ 0.98`
is not a symptom of a badly chosen `dedt`.  There is no weak-forcing limit in which a driven
steady state approximately satisfies the constraints.

## Relation to arXiv:1710.04264

That paper reconstructs a turbulent AAdS black brane horizon from forced 2D conformal fluid
data at perfect-fluid order and measures its fractal dimension.  Its setup and this one line
up closely, which is deliberate:

| | arXiv:1710.04264 | these runs |
| --- | --- | --- |
| fluid | (2+1)-d conformal, `P = rho/2`, `T^ab = (3/2) rho u^a u^b + (1/2) rho eta^ab` | same, `Gamma = 3/2`; `rho_rest` is an extra `10^-4` trace anomaly |
| grid | 2048^2, 2pi-periodic | 1024^2, unit-periodic |
| forcing | white-noise-in-time, homogeneous isotropic, narrow band at `k_f` | same (`driving_type = 2`, `tcorr = 0`) |
| `k_f` | 85 and 170 (`k_max = N/3 = 682`) | 24 (`k_max = 341`), i.e. `k_max/k_f = 14` |
| spectrum | `k^-2` (`k_f = 85`) and `k^-5/3` (`k_f = 170`) | `k^-2` expected at this ratio — check `spec_E` |
| large-scale sink | none, by argument | none |
| protocol | transient quasi-steady state, ensemble of 20, snapshots before box-scale pileup | same; use several seeds |
| dissipation | explicit 4th-order hyperdissipation | numerical (wenoz + FOFC) |
| metric order | perfect fluid | perfect fluid (`r_plus`, and `visc_*` if you want to go further) |

Three differences to keep in mind.  Their fluid carries *no* rest-mass density at all, so
their `T^ab` is exactly traceless while ours has a `10^-4` anomaly (which is why `temp0`
is set high rather than reaching for `--traceless`).
Their dissipation is an explicit hyperviscosity, so its magnitude is known; ours is the
scheme's, so it has to be measured.  And — the substantive one — they analyse *forced*
snapshots, whereas here the data set is taken from the unforced window, for the reason given
in the EFE section: forced data does not satisfy the constraint their own metric is
conditioned on.

## What has actually been tested

A 128^2 serial CPU build (`-D PROBLEM=fg_turb2d`, Kokkos serial) was run for 8 cycles with
these settings and post-processed:

* Initial state as designed at the shipped `temp0 = 5000`: `p/e = 0.499950`, trace anomaly
  `rho/e = 9.999e-5`, CFT temperature `e^(1/3) = 1.00003` (against the ideal-gas `p/rho =
  5000` — the two temperatures differ by that factor, do not confuse them).  The dust term
  `--traceless` would inject is `|rho a|/|(e+p) a| = 6.7e-5`.
* **Raising `temp0` costs no timestep.**  At `temp0 = 5000` the CFL step is `3.3147e-3` at
  128^2, matching `cfl*dx/c_s` to five digits and identical to the `temp0 = 500` run, since
  `c_s = sqrt(Gamma-1)` is independent of temperature.  (`rms(R^t)` also improved from
  `1.5e-3` to `3.8e-4`, and the baryon residual from `4.4e-7` to `1.4e-8`.)
* `max|u^z| = 0.0` exactly — the run is strictly 2D, not 2.5D.
* `rms(u^y)/rms(u^x) = 0.94` with `driving_type = 2` (3.64 with `driving_type = 1`).
* **The residual pipeline recovers a known source term.**  With `rel_cooling` on, the
  analytic energy sink is `-crate_rel * p * u^t`; the measured `mean(d_mu T^{mu t})` was
  `-0.9772` against an analytic `-0.9836`, i.e. 0.6% — which is simultaneously a check on
  the stress tensor, the non-uniform time stencil, the spatial stencil, the meshblock
  stitching and the primitive-variable conventions.  `mean(R^x)`, `mean(R^y)` came out at
  `~1e-7`, consistent with the driver's zero-net-force constraint.
* **The projections isolate the individual sources.**  `u_nu R^nu` came out at `+0.98419`
  against the analytic `crate_rel * p = +0.98359` (0.06%), confirming that the fluid-frame
  part of the residual is exactly `rel_cooling` and the orthogonal part
  (`res_perp_norm = 0.96`) is exactly the drive.
* **The constraint violation is band-limited, as the argument above requires.**  Of the
  residual power in `R^x`, 100.0% falls in the forcing band `10 <= k <= 14` and 0.00% at
  `k < 10` — so the inverse-cascade range the metric would be reconstructed from is clean.
* **Measured pointwise residual in a developed unforced state (128^2, `v_rms = 0.45`,
  `W = 1.11`):** `<|R^x|/N_x> = 5.7%` of the typical term magnitude
  (`rms(R^x) = 1.02` against `rms(norm_x) = 16.5`), and the same for `R^y`; baryon residual
  4.7%.  Three checks on what that number means:
  - **It is the scheme's error, not the analysis stencil's.**  Doubling the dump spacing
    (stride 1 -> 2) leaves it flat (5.7% -> 5.2%); only at stride 4-8 does it start growing
    (x2.2 per doubling) as the 2nd-order time stencil takes over.  So `dt_out ~ 1-2 dt_CFL`
    is the right cadence.  Spatial order 2 gives 7.6% against order 4's 5.7%.
  - **The box integral of the residual vanishes to 9 digits** —
    `mean(R^t)/rms(R^t) = 9.1e-9` — as required, since with the driver off the code conserves
    `Etot` to all 8 printed digits.  That is a strong end-to-end check of the stress tensor,
    both stencils, the meshblock stitching and the variable conventions at once.
  - **It lives at the grid scale.**  94% of the `R^x` power sits at `k > 21` (top third of
    wavenumbers) and 0.2% at `k <= 8`, while 87% of the *energy* is at `k <= 8`.  The
    truncation error is concentrated where the energy is not.
* **`Etot` is exactly conserved with the driver off, and *not* with it on.**  Over the forced
  phase `Etot` fell 2.4% while the drive should have added ~0.30; in the unforced window the
  drift is `0.000e+00`.  The loss is therefore an artifact of the driver's relativistic
  branch (its c2p round trip and global boost), not of the hydro scheme.  One more reason to
  take data from the unforced window.
* **Removing the thermostat drops the energy-constraint violation by ~700x.**  With
  `rel_cooling` on, `rms(R^t) = 0.98` and `<u.R> = 0.98`; with it off (the shipped
  configuration), `rms(R^t) = 1.5e-3` and `<u.R> = 8e-4`, and the relative violation
  `<|R^t|/N_t>` falls from 0.73 to 0.022.  What remains is broadband, i.e. truncation error
  rather than a source.  Meanwhile `e` stays put (1.00096) and `p/e` holds at 0.49950.
* `fg_residuals.py` also agrees with an analytically differentiated smooth field to 3e-7
  relative, and its shear tensor satisfies `u_mu sigma^{mu nu} = 0` and `sigma^mu_mu = 0`
  (both independent of AthenaK).

Not yet tested: the 1024^2 GPU build, and whether `dedt = 1.0` / `crate_rel = 2.0` land on
`v_rms ~ 0.3` in saturation.  Expect to iterate on those two numbers with the 256^2 case.

One transient to be aware of: at `t = 0` the flow is at rest, so the driver does no work
while `rel_cooling` is already removing energy at `crate_rel * p ~ 1` per unit time.  The
box therefore cools by tens of percent during spin-up before injection ramps up.  This is
self-correcting — `p = dedt/(crate_rel <u^t>)` is a stable fixed point, approached from
above — but do not read the first ~1 time unit as a steady state.

## Code changes

* **new** `src/pgen/fg_turb2d.cpp` — hot uniform conformal initial state (built through
  `SingleP2C_IdealSRHyd`, so the SR "evolve `E-D`" convention cannot be gotten wrong),
  guards for `special_rel` / 2D / `gamma = d/(d-1)`, and the history diagnostics above.
* **`src/srcterms/turb_driver.{cpp,hpp}`** — three additions, all backwards compatible:
  - `driving_type = 2`: 2D isotropic, exactly solenoidal, purely in-plane driving from a
    random stream function (see above for why 0 and 1 do not work).  It refuses to run on
    a 3D mesh, and `driving_type` outside `[0,2]` is now an error instead of silently
    producing a zero force.
  - `turb_driving/eint_over_dens_max` (default `40.0`) replaces a hard-coded
    `e_int <= 40 rho` clip that the relativistic branch of `AddForcing()` applied every
    stage.  That clip caps `T` at 20 and with it the trace anomaly at ~2.4%, which is fatal
    for a conformal run; the input files disable it (`1.0e6`).  The default reproduces the
    old behaviour exactly.
  - `InitializeModes()`/`AddForcing()` return immediately when `dedt <= 0`.  Previously a
    `dedt = 0` run still paid for the mode sum and, in SR, still pushed the fluid through a
    c2p round trip, the `e_int` clip and a global boost every stage.  This is what makes the
    decay window exact.

  Neither change affects the 3D Newtonian `PROBLEM=turb` runs (default value preserved; the
  relativistic branch is not reached).
