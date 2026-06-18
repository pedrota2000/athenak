from tqdm import tqdm
import numpy as np
import glob
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import cmasher as cmr

# ── config ────────────────────────────────────────────────────────────────────
data_path  = "/mnt/ceph/users/ptarancon/runs/turb_128_run_20260223_130640/vtk/"
output_mp4 = "enstrophy_animation.mp4"
FPS        = 60
BOX_LO, BOX_HI = -0.5, 0.5


# ── blue-black-orange colormap (zero = black) ─────────────────────────────────
cmap_ebo = cmr.fusion


# ── periodic gradient ─────────────────────────────────────────────────────────
def grad_periodic(f, dx, axis):
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2.0 * dx)


# ── load files ────────────────────────────────────────────────────────────────
all_vtk_files = sorted(glob.glob(data_path + 'Turb.hydro_w.*.vtk'))
print(f'Found {len(all_vtk_files)} VTK files')


def read_velocity(vtk_file):
    mesh = pv.read(vtk_file)
    Nx, Ny, Nz = mesh.dimensions - np.ones(3, dtype=int)
    dx, dy, dz = mesh.spacing

    def field(name):
        return mesh.cell_data[name].reshape(Nx, Ny, Nz).transpose(2, 1, 0)

    return field('velx'), field('vely'), field('velz'), dx, dy, dz


def compute_enstrophy(vx, vy, vz, dx, dy, dz):
    ox = grad_periodic(vz, dy, axis=1) - grad_periodic(vy, dz, axis=2)
    oy = grad_periodic(vx, dz, axis=2) - grad_periodic(vz, dx, axis=0)
    oz = grad_periodic(vy, dx, axis=0) - grad_periodic(vx, dy, axis=1)
    return 0.5 * (ox**2 + oy**2 + oz**2)


# ── pre-read all frames ───────────────────────────────────────────────────────
print("Reading and computing enstrophy for all frames...")
frames = []

for f in tqdm(all_vtk_files, desc="Reading VTK files"):
    try:
        vx, vy, vz, dx, dy, dz = read_velocity(f)
        frames.append(compute_enstrophy(vx, vy, vz, dx, dy, dz))
    except Exception as e:
        print(f"\n  Skipping {f}: {e}")

N                   = len(frames)
Nx, Ny, Nz          = frames[0].shape
mid_x, mid_y, mid_z = Nx // 2, Ny // 2, Nz // 2

# ── global colour limits (99th percentile) ────────────────────────────────────
vmax = np.percentile(np.concatenate([e.ravel() for e in frames]), 99)
vmin = 0.0

# ── build figure ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.45)

slice_labels = [('y', 'z'), ('x', 'z'), ('x', 'y')]


def get_slices(e):
    return [e[mid_x, :, :], e[:, mid_y, :], e[:, :, mid_z]]


slices0 = get_slices(frames[0])
ims = []
for col, (sl, (xl, yl)) in enumerate(zip(slices0, slice_labels)):
    ax = axes[col]
    im = ax.imshow(sl.T, origin='lower',
                   extent=[BOX_LO, BOX_HI, BOX_LO, BOX_HI],
                   cmap=cmap_ebo, vmin=vmin, vmax=vmax, aspect='auto')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.08)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ims.append(im)


def update(frame_idx):
    slices = get_slices(frames[frame_idx])
    for col, sl in enumerate(slices):
        ims[col].set_data(sl.T)
    return ims


ani    = animation.FuncAnimation(fig, update, frames=N, interval=1000 // FPS, blit=True)
writer = FFMpegWriter(fps=FPS, metadata={"title": "Enstrophy"}, bitrate=2000)
print("Saving mp4...")
ani.save(output_mp4, writer=writer,
         progress_callback=lambda i, n: print(f"  Frame {i+1}/{n}", end="\r"))
plt.close(fig)
print(f"\nSaved: {output_mp4}")