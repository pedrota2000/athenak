from tqdm import tqdm
import numpy as np
import glob
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
from matplotlib import patheffects
from read_data import athdf
import tempfile
import os
import bin_convert


# ── config ────────────────────────────────────────────────────────────────────
data_path  = "/mnt/home/ptarancon/ceph/runs/turb_1024_viscous_compressible_run_v2_20260622_113949/bin/"
output_mp4 = "turbulence_6panel_1024_v2.mp4"
FPS        = 15
BOX_LO, BOX_HI = -0.5, 0.5

# Isothermal EOS: cs = const = 1.0  (from config: iso_sound_speed = 1.0)
# Expected Mach ~ 0.5  (dedt=0.25 → urms ~ (2*dedt*L_drive)^(1/3) ≈ 0.63,
# sol_fraction=1.0 → purely solenoidal driving, mild density variance σ_lnρ ~ b·M ~ 0.17)
CS = 1.0 #1.0

# ── panel metadata ────────────────────────────────────────────────────────────
# (title, colormap, diverging, log_scale)
#
# Colormaps chosen for WHITE backgrounds:
#   Sequential  →  single-hue matplotlib cmaps (white at vmin, saturated at vmax)
#   Diverging   →  bwr  (blue–white–red, pure white at zero)
#   Schlieren   →  hot_r (white→yellow→orange→red→black, classic shock visualisation)
#
# Temperature removed: isothermal EOS keeps P/ρ = cs² = const everywhere.
# Replaced with |∇ρ|  (numerical schlieren) which reveals the mild shocks / interfaces
# that still form at Mach ~ 0.5.

PANEL_META = [
    ("Density  ρ  [log]",      "cividis",   False, True),
    ("Enstrophy  ½|ω|²",       "magma",     False, False),
    ("Dilatation ∇·v",         "coolwarm",  True,  False),
    ("Kin. energy  ½ρ|v|²",    "inferno",   False, False),
    ("|∇ρ| (schlieren)",       "turbo",     False, False),
    ("Mach number |v|/cₛ",     "viridis",   False, False),
]
N_PANELS = len(PANEL_META)

# ── periodic finite-difference gradient ──────────────────────────────────────
def grad_periodic(f, dx, axis):
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2.0 * dx)


# ── I/O ───────────────────────────────────────────────────────────────────────
# Get all the files inside data_path, not only vtk #
all_files = sorted(
    glob.glob(data_path + "/*.vtk") +
    glob.glob(data_path + "/*.bin")
)

print(f'Found {len(all_files)} files')



def read_fields(vtk_file):
    """Return (dens, vx, vy, vz, dx, dy, dz).
    For isothermal runs pressure is not needed (cs = const)."""
    mesh = pv.read(vtk_file)
    Nx, Ny, Nz = mesh.dimensions - np.ones(3, dtype=int)
    dx, dy, dz = mesh.spacing

    def field(name):
        return mesh.cell_data[name].reshape(Nx, Ny, Nz).transpose(2, 1, 0)

    return (field('dens'), field('velx'), field('vely'), field('velz'),
            dx, dy, dz)

def read_fields_bin(bin_file):
    """
    Convert .bin → temporary .athdf → read → delete.
    Returns (dens, vx, vy, vz, dx, dy, dz)
    """

    # create temporary directory
    tmp_dir = tempfile.gettempdir()

    base = os.path.basename(bin_file).replace(".bin", "")
    athdf_file = os.path.join(tmp_dir, base + ".athdf")

    try:
        # 1) read binary
        filedata = bin_convert.read_binary(bin_file)

        # 2) write temporary athdf
        bin_convert.write_athdf(athdf_file, filedata)

        # 3) read it using your existing ATHDF pipeline
        data = athdf(athdf_file)

        dens = np.asarray(data["dens"])
        vx   = np.asarray(data["velx"])
        vy   = np.asarray(data["vely"])
        vz   = np.asarray(data["velz"])

        x = np.asarray(data["x1v"])
        y = np.asarray(data["x2v"])
        z = np.asarray(data["x3v"])

        dx = float(np.mean(np.diff(x)))
        dy = float(np.mean(np.diff(y)))
        dz = float(np.mean(np.diff(z)))

        # fix orientation if needed
        if dens.shape[0] != len(x):
            dens = dens.transpose(2,1,0)
            vx   = vx.transpose(2,1,0)
            vy   = vy.transpose(2,1,0)
            vz   = vz.transpose(2,1,0)

        return dens, vx, vy, vz, dx, dy, dz

    finally:
        # 4) cleanup (VERY IMPORTANT)
        if os.path.exists(athdf_file):
            os.remove(athdf_file)


# ── physics ───────────────────────────────────────────────────────────────────
def compute_quantities(dens, vx, vy, vz, dx, dy, dz):
    """Return 6-tuple matching PANEL_META order."""

    # ① Vorticity → enstrophy
    ox = grad_periodic(vz, dy, 1) - grad_periodic(vy, dz, 2)
    oy = grad_periodic(vx, dz, 2) - grad_periodic(vz, dx, 0)
    oz = grad_periodic(vy, dx, 0) - grad_periodic(vx, dy, 1)
    enstrophy = 0.5 * (ox**2 + oy**2 + oz**2)

    # ② Dilatation (divergence of velocity)
    # negative → compression / shocks; positive → rarefaction
    div_v = (grad_periodic(vx, dx, 0) +
             grad_periodic(vy, dy, 1) +
             grad_periodic(vz, dz, 2))

    # ③ Kinetic energy density
    v2 = vx**2 + vy**2 + vz**2
    kin_energy = 0.5 * dens * v2

    # ④ Density-gradient magnitude  (numerical schlieren)
    # highlights shocks and contact discontinuities
    gx = grad_periodic(dens, dx, 0)
    gy = grad_periodic(dens, dy, 1)
    gz = grad_periodic(dens, dz, 2)
    grad_dens = np.sqrt(gx**2 + gy**2 + gz**2)

    # ⑤ Mach number  (isothermal: cs = CS = const)
    mach = np.sqrt(v2) / CS

    # Compute mean mach number, not pointwise #
    mach_mean = np.sqrt(np.mean(v2)) / CS
    print(f"  Mean Mach number: {mach_mean:.3f}")

    return dens, enstrophy, div_v, kin_energy, grad_dens, mach


# ── pre-read all frames ───────────────────────────────────────────────────────
print("Reading VTK files and computing quantities...")
frames = []

for f in tqdm(all_files, desc="Loading"):
    try:
        # Check if the files are vtk or athdf and read accordingly #
        if f.endswith('.vtk'):
            dens, vx, vy, vz, dx, dy, dz = read_fields(f)
        else:
            dens, vx, vy, vz, dx, dy, dz = read_fields_bin(f)
        frames.append(compute_quantities(dens, vx, vy, vz, dx, dy, dz))
    except Exception as e:
        print(f"\n  Skipping {f}: {e}")

N                   = len(frames)
Nx, Ny, Nz          = frames[0][0].shape
mid_z               = Nz // 2
print(f"  Grid: {Nx}×{Ny}×{Nz}  |  {N} frames  |  showing z-midplane")

# ── colour limits (global across all frames) ──────────────────────────────────
norms = []
for qi, (_, _, diverging, log_scale) in enumerate(PANEL_META):
    all_vals = np.concatenate([f[qi].ravel() for f in frames])
    if diverging:
        vabs = np.percentile(np.abs(all_vals), 99)
        norms.append(Normalize(-vabs, vabs))
    elif log_scale:
        vmin = max(np.percentile(all_vals,  1), 1e-10)
        vmax = np.percentile(all_vals, 99)
        norms.append(LogNorm(vmin=vmin, vmax=vmax))
    else:
        norms.append(Normalize(np.percentile(all_vals,  1),
                               np.percentile(all_vals, 99)))

# ── figure (white background) ─────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.patch.set_facecolor('black')
fig.subplots_adjust(wspace=0.40, hspace=0.46, top=0.92)

BG = "black"
FG = "white"
TICK = "#e6e6e6"
SPINE = "#aaaaaa"
CBAR = "#ffffff"

ims = []
for qi, ((title, cmap, _, _), norm) in enumerate(zip(PANEL_META, norms)):
    row, col = divmod(qi, 3)
    ax = axes[row, col]
    ax.set_facecolor(BG)

    ax.tick_params(colors=TICK, labelsize=10, width=1.2, length=5)

    for spine in ax.spines.values():
        spine.set_edgecolor(SPINE)
        spine.set_linewidth(1.2)

    ax.xaxis.label.set_color(FG)
    ax.yaxis.label.set_color(FG)
    ax.title.set_color(FG)

    sl = frames[0][qi][:, :, mid_z]
    im = ax.imshow(sl.T, origin='lower',
                   extent=[BOX_LO, BOX_HI, BOX_LO, BOX_HI],
                   cmap=cmap, norm=norm, aspect='equal')

    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    cb.ax.tick_params(color=FG, labelcolor=FG, labelsize=9, width=1.2, length=4)
    cb.outline.set_edgecolor(SPINE)
    cb.outline.set_linewidth(1.2)

    ims.append(im)
    ax.grid(False)
    ax.set_title(
        title,
        color='white',
        fontsize=14,
        pad=12,
        fontweight='bold'
    )
    ax.title.set_path_effects([
        matplotlib.patheffects.withStroke(linewidth=3, foreground='black')
    ])

frame_title = fig.suptitle(
    "Isothermal compressible turbulence  (M ≈ 0.5, solenoidal)  |  "
    "z-midplane  |  frame 0",
    color='white', fontsize=12, y=0.97
)


# ── animation ─────────────────────────────────────────────────────────────────
def update(frame_idx):
    for qi in range(N_PANELS):
        ims[qi].set_data(frames[frame_idx][qi][:, :, mid_z].T)
    frame_title.set_text(
        f"Isothermal compressible turbulence  (M ≈ 0.5, solenoidal)  |  "
        f"z-midplane  |  frame {frame_idx + 1}/{N}"
    )
    return ims + [frame_title]


ani = animation.FuncAnimation(
    fig, update, frames=N, interval=1000 // FPS, blit=True
)
writer = FFMpegWriter(
    fps=FPS,
    metadata={"title": "Isothermal Turbulence – 6 panels"},
    bitrate=6000,
    extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p']
)

print("Saving mp4…")
ani.save(
    output_mp4, writer=writer,
    progress_callback=lambda i, n: print(f"  Frame {i+1}/{n}", end="\r")
)
plt.close(fig)
print(f"\nSaved → {output_mp4}")