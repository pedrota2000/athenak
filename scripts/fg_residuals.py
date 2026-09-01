#!/usr/bin/env python3
"""
Pointwise fluid/gravity variables and inviscid conservation-law residuals from AthenaK
2D special-relativistic ("conformal") turbulence dumps.

For each snapshot this writes, cell by cell, the four quantities that enter the
fluid/gravity metric,

    u^t, u^x, u^y      (four-velocity; u^z = 0 in a 2D run)
    e                  (fluid-frame energy density: the eigenvalue of T^{mu nu},
                        T = diag(-e, p, p, p) in the rest frame)

together with the residuals of the inviscid relativistic conservation laws

    R^nu  = d_mu T^{mu nu},        T^{mu nu} = (e + p) u^mu u^nu + p eta^{mu nu}
    R_N   = d_mu (rho u^mu)        (baryon number, for reference)

with eta = diag(-1, +1, +1, +1) and c = 1.

In the fluid/gravity construction the bulk Einstein constraint equations ARE the
conservation of the boundary stress tensor, so R^nu is the constraint violation of the
metric built from (u^mu, e) -- not merely a code-quality diagnostic.  It is split into the
part along u (fluid-frame heating: this is what `rel_cooling` contributes) and the
orthogonal part (a pure force: this is what the turbulence driver contributes, whose source
four-vector is exactly orthogonal to u), reported as `res_dot_u` and `res_perp_norm`.

Also written, because they are the first-order data the fluid/gravity metric is built from
and they set the validity of the gradient expansion: the expansion `theta = d_mu u^mu`, the
four-acceleration `a^mu`, the shear tensor `sigma^{mu nu}` and its norm, the Knudsen number
`|sigma|/T`, the temperature `T = C e^(1/d)` and the perfect-fluid-order horizon radius
`r_+ = 4 pi T / d`, and the holographic shear viscosity `eta = s/4pi` with the first-order
stress `-2 eta sigma^{mu nu}` it implies.  The normalization `C` (--tnorm) defaults to 1,
the convention of Westernacher-Schneider, arXiv:1710.04264, who takes T = rho^(1/3) with
rho the energy density of a (2+1)-dimensional conformal fluid and R_AdS = 1.

AthenaK's SR-hydro primitives are (dens, velx, vely, velz, eint) = (rho, u^x, u^y, u^z,
e_int), i.e. the velocity slots hold the *four*-velocity.  The energy density that
diagonalizes T is therefore

    e = rho + e_int,        p = (gamma - 1) e_int,

and for the conformal choice gamma = d/(d-1) one has exactly  T^mu_mu = -rho,  so the
rest-mass density that AthenaK carries along IS the trace anomaly.  It is written out so
it can be monitored (--report).  Use the FULL tensor, e = rho + e_int: it is the one that is
exactly conserved, so the residual measures only numerical error.  --traceless replaces e by
the strictly traceless (d-1) p, but that tensor differs from the conserved one by a dust term
and satisfies d_mu T^{mu nu} = -rho a^nu instead of zero -- an optional cross-check, not the
default.  Shrink the anomaly with <problem>/temp0 instead (rho/e_int = 1/(2 temp0)).

Time derivatives use a 3-point non-uniform stencil (AthenaK writes a dump on the first
step at or after each target time, so the dump times are NOT exactly equispaced), and
residuals are reported at the middle time of each consecutive triple.  Spatial
derivatives are periodic central differences; use --order 4 so the analysis stencil is
not the dominant error.

Usage
-----
  python3 scripts/fg_residuals.py 'run/bin/fgturb2d.hydro_w.*.bin' -o run/residuals \
      --gamma 1.5 --d 3 --order 4 --report

Each output file  <outdir>/resid.<NNNNN>.npz  contains the fields listed above plus the
per-equation normalizer  norm_{t,x,y} = sum of |individual terms|, so that
|R^nu| / norm^nu is a dimensionless measure of how well the ideal law is satisfied.
"""

import argparse
import glob
import os
import sys
import types

import numpy as np

# bin_convert.py lives both in scripts/ and vis/python/ in this tree; prefer scripts/
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(os.path.dirname(_here), "vis", "python"), _here):
    if _p not in sys.path:
        sys.path.insert(0, _p)
try:
    import bin_convert  # noqa: E402
except ImportError as _exc:  # pragma: no cover
    # bin_convert imports h5py at module scope, but only its athdf *writer* needs it;
    # stub it out so the .bin readers work in environments without HDF5.
    if "h5py" not in str(_exc):
        raise
    sys.modules["h5py"] = types.ModuleType("h5py")
    import bin_convert  # noqa: E402


def load_snapshot(fname, per_rank=False):
    """Read one .bin dump and stitch the meshblocks into 2D (nx2, nx1) arrays."""
    if per_rank:
        fd = bin_convert.read_all_ranks_binary(fname)
    else:
        fd = bin_convert.read_binary(fname)

    if fd["Nx3"] != 1:
        raise ValueError(f"{fname}: Nx3 = {fd['Nx3']}, this script is 2D only")
    if fd["mb_logical"][:, 3].max() != 0:
        raise ValueError(f"{fname}: mesh refinement is not supported")

    nx1, nx2 = fd["nx1_out_mb"], fd["nx2_out_mb"]
    nx1g, nx2g = fd["Nx1"], fd["Nx2"]
    out = {}
    for var, data in fd["mb_data"].items():
        # read_binary returns a list of per-meshblock arrays; read_all_ranks_binary an
        # ndarray.  Index with [m][0] so both work (block m, k = 0).
        g = np.empty((nx2g, nx1g), dtype=np.float64)
        for m, (lx1, lx2, _lx3, _lev) in enumerate(fd["mb_logical"]):
            g[lx2 * nx2:(lx2 + 1) * nx2, lx1 * nx1:(lx1 + 1) * nx1] = data[m][0]
        out[var] = g

    for need in ("dens", "velx", "vely", "eint"):
        if need not in out:
            raise ValueError(f"{fname}: missing '{need}'; dump <output>/variable = hydro_w")

    out["time"] = fd["time"]
    out["dx1"] = (fd["x1max"] - fd["x1min"]) / nx1g
    out["dx2"] = (fd["x2max"] - fd["x2min"]) / nx2g
    out["x1v"] = fd["x1min"] + out["dx1"] * (0.5 + np.arange(nx1g))
    out["x2v"] = fd["x2min"] + out["dx2"] * (0.5 + np.arange(nx2g))
    return out


def state(snap, gamma, dst, traceless):
    """Primitives -> (u^mu, e, p) and the stress tensor pieces needed here."""
    rho = snap["dens"]
    ux, uy = snap["velx"], snap["vely"]
    uz = snap.get("velz")
    eint = snap["eint"]

    usq = ux * ux + uy * uy
    if uz is not None:
        usq = usq + uz * uz
    ut = np.sqrt(1.0 + usq)

    p = (gamma - 1.0) * eint
    e = (dst - 1.0) * p if traceless else rho + eint
    w = e + p                                    # enthalpy density
    return dict(rho=rho, ux=ux, uy=uy, ut=ut, eint=eint, p=p, e=e, w=w, usq=usq)


def ddx(f, dx, axis, order):
    """Periodic centered difference along `axis`."""
    if order == 2:
        return (np.roll(f, -1, axis) - np.roll(f, 1, axis)) / (2.0 * dx)
    if order == 4:
        return (8.0 * (np.roll(f, -1, axis) - np.roll(f, 1, axis))
                - (np.roll(f, -2, axis) - np.roll(f, 2, axis))) / (12.0 * dx)
    raise ValueError("order must be 2 or 4")


def ddt3(f0, f1, f2, t0, t1, t2):
    """Derivative at t1 from a 3-point stencil with unequal spacing."""
    h0, h1 = t1 - t0, t2 - t1
    if h0 <= 0.0 or h1 <= 0.0:
        raise ValueError("snapshot times must be strictly increasing")
    return (-h1 / (h0 * (h0 + h1)) * f0
            + (h1 - h0) / (h0 * h1) * f1
            + h0 / (h1 * (h0 + h1)) * f2)


def tt_fluxes(s):
    """The T^{t nu} components (and rho u^t), i.e. everything d/dt acts on."""
    return dict(
        Ttt=s["w"] * s["ut"] * s["ut"] - s["p"],
        Ttx=s["w"] * s["ut"] * s["ux"],
        Tty=s["w"] * s["ut"] * s["uy"],
        N=s["rho"] * s["ut"],
    )


ETA = (-1.0, 1.0, 1.0)          # diag of eta_{mu nu} on the (t,x,y) block


def horizon(e, dst, tnorm):
    """CFT temperature, horizon radius, entropy density and shear viscosity in code units.

    For a conformal fluid e ~ T^d, so T = tnorm * e^(1/d).  Westernacher-Schneider
    (arXiv:1710.04264) uses tnorm = 1 for the (2+1)-dimensional case -- there T = rho^(1/3)
    with rho the energy density -- together with an AdS radius R = 1, and locates the
    horizon of the perfect-fluid-order metric at r_+ = 4 pi T / d.  Entropy density follows
    from e + p = T s for a conformal fluid, and Einstein-gravity duals have eta/s = 1/4pi
    (the KSS value), so eta = (e + p) / (4 pi T).
    """
    temp = tnorm * np.power(np.maximum(e, 1e-300), 1.0 / dst)
    r_plus = 4.0 * np.pi * temp / dst
    p_of_e = e / (dst - 1.0)                 # conformal
    s_dens = (e + p_of_e) / temp
    eta = s_dens / (4.0 * np.pi)
    return temp, r_plus, s_dens, eta


def uvec(s):
    """Contravariant u^mu on the (t,x,y) block, from a state dict or a raw snapshot."""
    if "ut" in s:
        return [s["ut"], s["ux"], s["uy"]]
    ux, uy = s["velx"], s["vely"]
    uz = s.get("velz")
    usq = ux * ux + uy * uy
    if uz is not None:
        usq = usq + uz * uz
    return [np.sqrt(1.0 + usq), ux, uy]


def gradients(prev, cur, nxt, times, order):
    """du[a][b] = d_a u_b  on the (t,x,y) block (lower index on u)."""
    u0v, u1v, u2v = uvec(prev), uvec(cur), uvec(nxt)
    d1, d2 = cur["dx1"], cur["dx2"]
    t0, t1, t2 = times
    du = [[None] * 3 for _ in range(3)]
    for b in range(3):
        sgn = ETA[b]                        # u_b = eta_{bb} u^b
        du[0][b] = ddt3(sgn * u0v[b], sgn * u1v[b], sgn * u2v[b], t0, t1, t2)
        du[1][b] = ddx(sgn * u1v[b], d1, 1, order)
        du[2][b] = ddx(sgn * u1v[b], d2, 0, order)
    return du


def first_order(s, du, dst):
    """Expansion, four-acceleration and shear tensor: the first-order fluid/gravity data.

    theta      = d_mu u^mu
    a^mu       = u^nu d_nu u^mu
    sigma^{mu nu} = Delta^{mu al} Delta^{nu be} d_{(al} u_{be)} - Delta^{mu nu} theta/(d-1)

    with Delta^{mu nu} = eta^{mu nu} + u^mu u^nu.  sigma is the object the first-order term
    of the fluid/gravity metric is built from; it obeys u_mu sigma^{mu nu} = 0 and
    sigma^mu_mu = 0 (checked in the unit test).  For d = 4 with a 2D flow the only extra
    component is sigma^{zz} = -theta/(d-1).
    """
    u = uvec(s)
    theta = -du[0][0] + du[1][1] + du[2][2]

    # a^mu = eta^{mu mu} sum_nu u^nu d_nu u_mu
    acc = [ETA[m] * sum(u[n] * du[n][m] for n in range(3)) for m in range(3)]

    delta = [[(ETA[m] if m == n else 0.0) + u[m] * u[n] for n in range(3)]
             for m in range(3)]
    sym = [[0.5 * (du[a][b] + du[b][a]) for b in range(3)] for a in range(3)]

    sigma = [[None] * 3 for _ in range(3)]
    for m in range(3):
        for n in range(m, 3):
            acc_mn = 0.0
            for a in range(3):
                for b in range(3):
                    acc_mn = acc_mn + delta[m][a] * sym[a][b] * delta[b][n]
            sigma[m][n] = acc_mn - delta[m][n] * theta / (dst - 1.0)
    for m in range(3):
        for n in range(m):
            sigma[m][n] = sigma[n][m]

    # sigma^{mu nu} sigma_{mu nu} = sum_{mn} eta_{mm} eta_{nn} (sigma^{mn})^2
    s2 = sum(ETA[m] * ETA[n] * sigma[m][n] * sigma[m][n]
             for m in range(3) for n in range(3))
    sigma_zz = None
    if dst > 3.5:                       # d = 4: the ignorable z direction still shears
        sigma_zz = -theta / (dst - 1.0)
        s2 = s2 + sigma_zz * sigma_zz
    return theta, acc, sigma, np.sqrt(np.maximum(s2, 0.0)), sigma_zz


def shell_spectrum(fields, nbins=None):
    """Shell-averaged spectra on integer wavenumber bins n = |k| L / 2pi.

    Used to check *where in k* the constraint violation lives.  A force with support in a
    narrow band around k_f violates d_mu T^{mu nu} = 0 only in that band, so conclusions
    drawn from the source-free bands (for a 2D inverse cascade, k < k_f) survive -- this is
    the argument that makes a *forced* fluid/gravity data set defensible, and this function
    is how you verify it instead of assuming it.
    """
    ny, nx = fields[0].shape
    ky = np.fft.fftfreq(ny, d=1.0 / ny)
    kx = np.fft.fftfreq(nx, d=1.0 / nx)
    kmag = np.sqrt(kx[None, :] ** 2 + ky[:, None] ** 2)
    if nbins is None:
        nbins = min(nx, ny) // 2
    idx = np.minimum(np.rint(kmag).astype(int), nbins)
    out = []
    for f in fields:
        power = np.abs(np.fft.fft2(f) / (nx * ny)) ** 2
        out.append(np.bincount(idx.ravel(), weights=power.ravel(),
                               minlength=nbins + 1)[:nbins + 1])
    return np.arange(nbins + 1, dtype=np.float64), out


def project_residual(s, res):
    """Split R^mu into its fluid-frame heating part u_mu R^mu and the orthogonal rest.

    A pure force (f_mu u^mu = 0, which is what AthenaK's relativistic driver applies)
    contributes only to the orthogonal part; a fluid-frame energy sink (rel_cooling,
    f^mu ~ -u^mu) contributes only to u_mu R^mu.
    """
    u = uvec(s)
    ul = [ETA[m] * u[m] for m in range(3)]
    r = [res["res_t"], res["res_x"], res["res_y"]]
    dot = sum(ul[m] * r[m] for m in range(3))
    perp = [r[m] + u[m] * dot for m in range(3)]
    p2 = sum(ETA[m] * perp[m] * perp[m] for m in range(3))
    return dot, np.sqrt(np.maximum(p2, 0.0))


def residuals(prev, cur, nxt, times, gamma, dst, order, traceless, tnorm=1.0):
    """Residuals of d_mu T^{mu nu} = 0 and d_mu (rho u^mu) = 0 at the middle time."""
    s = state(cur, gamma, dst, traceless)
    d1, d2 = cur["dx1"], cur["dx2"]

    f0 = tt_fluxes(state(prev, gamma, dst, traceless))
    f1 = tt_fluxes(s)
    f2 = tt_fluxes(state(nxt, gamma, dst, traceless))
    t0, t1, t2 = times

    dt_Ttt = ddt3(f0["Ttt"], f1["Ttt"], f2["Ttt"], t0, t1, t2)
    dt_Ttx = ddt3(f0["Ttx"], f1["Ttx"], f2["Ttx"], t0, t1, t2)
    dt_Tty = ddt3(f0["Tty"], f1["Tty"], f2["Tty"], t0, t1, t2)
    dt_N = ddt3(f0["N"], f1["N"], f2["N"], t0, t1, t2)

    w, p, ut, ux, uy = s["w"], s["p"], s["ut"], s["ux"], s["uy"]
    # axis 1 is x1, axis 0 is x2
    dx_Txt = ddx(w * ux * ut, d1, 1, order)
    dy_Tyt = ddx(w * uy * ut, d2, 0, order)
    dx_Txx = ddx(w * ux * ux + p, d1, 1, order)
    dy_Tyx = ddx(w * uy * ux, d2, 0, order)
    dx_Txy = ddx(w * ux * uy, d1, 1, order)
    dy_Tyy = ddx(w * uy * uy + p, d2, 0, order)
    dx_Nx = ddx(s["rho"] * ux, d1, 1, order)
    dy_Ny = ddx(s["rho"] * uy, d2, 0, order)

    out = dict(
        time=np.float64(t1),
        # the four fluid/gravity variables, pointwise
        ut=ut, ux=ux, uy=uy, e=s["e"],
        # supporting thermodynamics / anomaly
        p=p, rho=s["rho"], anomaly=s["rho"] / (s["rho"] + s["eint"]),
        # residuals
        res_t=dt_Ttt + dx_Txt + dy_Tyt,
        res_x=dt_Ttx + dx_Txx + dy_Tyx,
        res_y=dt_Tty + dx_Txy + dy_Tyy,
        res_N=dt_N + dx_Nx + dy_Ny,
        # normalizers: sum of |terms| in each equation
        norm_t=np.abs(dt_Ttt) + np.abs(dx_Txt) + np.abs(dy_Tyt),
        norm_x=np.abs(dt_Ttx) + np.abs(dx_Txx) + np.abs(dy_Tyx),
        norm_y=np.abs(dt_Tty) + np.abs(dx_Txy) + np.abs(dy_Tyy),
        norm_N=np.abs(dt_N) + np.abs(dx_Nx) + np.abs(dy_Ny),
    )

    # first-order (gradient) data: what the O(gradient) fluid/gravity metric needs, and
    # what the ideal residual has to be compared against to decide whether the profile is
    # a first-order solution of the Einstein equations
    du = gradients(prev, cur, nxt, times, order)
    theta, acc, sigma, sigma_norm, sigma_zz = first_order(s, du, dst)
    out["theta"] = theta
    out["acc_t"], out["acc_x"], out["acc_y"] = acc
    for m, mn in enumerate("txy"):
        for n, nn in enumerate("txy"):
            if n >= m:
                out[f"sigma_{mn}{nn}"] = sigma[m][n]
    if sigma_zz is not None:
        out["sigma_zz"] = sigma_zz
    out["sigma_norm"] = sigma_norm
    # temperature, horizon radius and transport, in the normalization of arXiv:1710.04264
    temp, r_plus, s_dens, eta = horizon(s["e"], dst, tnorm)
    out["temp"] = temp
    out["r_plus"] = r_plus            # the perfect-fluid-order horizon: r_+ = 4 pi T / d
    out["eta"] = eta                  # = s/4pi, the holographic shear viscosity
    # Knudsen number: the gradient expansion the metric is built on needs this << 1
    out["knudsen"] = sigma_norm / temp
    # the first-order stress the data would have to carry to be a first-order solution
    out["visc_xx"] = -2.0 * eta * sigma[1][1]
    out["visc_xy"] = -2.0 * eta * sigma[1][2]
    out["visc_yy"] = -2.0 * eta * sigma[2][2]
    # residual split into fluid-frame heating and pure-force parts
    out["res_dot_u"], out["res_perp_norm"] = project_residual(s, out)

    # spectra: E(k) of the 3-velocity (the quantity compared with k^-2 / k^-5/3 in the
    # turbulence literature) and of each residual component, so that the band in which the
    # Einstein constraints are actually violated can be read off directly
    vx, vy = s["ux"] / s["ut"], s["uy"] / s["ut"]
    kbins, spec = shell_spectrum([vx, vy, out["res_t"], out["res_x"], out["res_y"]])
    out["spec_k"] = kbins
    out["spec_E"] = 0.5 * (spec[0] + spec[1])
    out["spec_res_t"], out["spec_res_x"], out["spec_res_y"] = spec[2], spec[3], spec[4]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="glob(s) of hydro_w .bin dumps, in time order")
    ap.add_argument("-o", "--outdir", default="fg_residuals", help="output directory")
    ap.add_argument("--gamma", type=float, default=1.5, help="<hydro>/gamma (default 1.5)")
    ap.add_argument("--d", dest="dst", type=float, default=3.0,
                    help="space-time dimension d in p = e/(d-1) (default 3)")
    ap.add_argument("--order", type=int, default=4, choices=(2, 4),
                    help="spatial difference order (default 4)")
    ap.add_argument("--tnorm", type=float, default=1.0,
                    help="CFT temperature normalization C in T = C e^(1/d) "
                         "(default 1.0, the convention of arXiv:1710.04264)")
    ap.add_argument("--traceless", action="store_true",
                    help="cross-check with the strictly traceless e = (d-1) p; note this tensor is "
                         "not the exactly conserved one (see the module docstring)")
    ap.add_argument("--per-rank", action="store_true",
                    help="dumps were written with single_file_per_rank")
    ap.add_argument("--report", action="store_true", help="print per-snapshot RMS summary")
    ap.add_argument("--fields-only", action="store_true",
                    help="write u^mu and e only, skip the residuals (no time stencil)")
    args = ap.parse_args()

    names = []
    for pat in args.files:
        names.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])
    if not names:
        sys.exit("no input files matched")
    os.makedirs(args.outdir, exist_ok=True)

    if abs(args.gamma - args.dst / (args.dst - 1.0)) > 1e-8:
        print(f"# WARNING: gamma = {args.gamma} is not the conformal value "
              f"d/(d-1) = {args.dst/(args.dst-1.0)} for d = {args.dst}", file=sys.stderr)

    if args.fields_only:
        for n, name in enumerate(names):
            snap = load_snapshot(name, args.per_rank)
            s = state(snap, args.gamma, args.dst, args.traceless)
            np.savez_compressed(os.path.join(args.outdir, f"fields.{n:05d}.npz"),
                                time=np.float64(snap["time"]), x1v=snap["x1v"],
                                x2v=snap["x2v"], ut=s["ut"], ux=s["ux"], uy=s["uy"],
                                e=s["e"], p=s["p"], rho=s["rho"])
            print(f"{name} -> fields.{n:05d}.npz  t = {snap['time']:.6g}")
        return

    if len(names) < 3:
        sys.exit("need at least 3 consecutive dumps for the time derivative")

    if args.report:
        print("#    t         rms(R^t)    rms(R^x)    rms(R^y)    rms(R^N)   "
              "<|R^t|/N_t>  <|R^x|/N_x>  <|R^y|/N_y>    <u.R>      <|R_perp|>  "
              "rms|sigma|   <Kn>       <anomaly>")

    window = [load_snapshot(names[0], args.per_rank),
              load_snapshot(names[1], args.per_rank)]
    for n in range(1, len(names) - 1):
        window.append(load_snapshot(names[n + 1], args.per_rank))
        prev, cur, nxt = window
        times = (prev["time"], cur["time"], nxt["time"])
        r = residuals(prev, cur, nxt, times, args.gamma, args.dst, args.order,
                      args.traceless, args.tnorm)
        r["x1v"], r["x2v"] = cur["x1v"], cur["x2v"]
        np.savez_compressed(os.path.join(args.outdir, f"resid.{n:05d}.npz"), **r)

        if args.report:
            def rms(a):
                return float(np.sqrt(np.mean(a * a)))

            def rel(a, b):
                return float(np.mean(np.abs(a) / np.maximum(b, 1e-300)))
            print(f"{r['time']:11.5g} {rms(r['res_t']):11.4e} {rms(r['res_x']):11.4e} "
                  f"{rms(r['res_y']):11.4e} {rms(r['res_N']):11.4e} "
                  f"{rel(r['res_t'], r['norm_t']):12.4e} "
                  f"{rel(r['res_x'], r['norm_x']):12.4e} "
                  f"{rel(r['res_y'], r['norm_y']):12.4e} "
                  f"{float(np.mean(r['res_dot_u'])):11.4e} "
                  f"{float(np.mean(r['res_perp_norm'])):11.4e} "
                  f"{rms(r['sigma_norm']):11.4e} "
                  f"{float(np.mean(r['knudsen'])):10.3e} "
                  f"{float(np.mean(r['anomaly'])):11.4e}")
        window.pop(0)


if __name__ == "__main__":
    main()
