"""
Star Formation Simulation
=========================
Simulates gravitational collapse of a gas cloud into a protostar.

Requirements:
    pip install astropy matplotlib numpy

Run:
    python star_formation_sim.py

Output:
    star_formation.gif  (saved in the same directory)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from astropy import units as u
from astropy import constants as const

# ── Reproducibility ──────────────────────────────────────────────────────────
np.random.seed(42)

# ── Physical setup (using astropy units for clarity) ─────────────────────────
G_si      = const.G.to(u.m**3 / (u.kg * u.s**2)).value   # m³ kg⁻¹ s⁻²
M_cloud   = (1.0 * u.M_sun).to(u.kg).value                # total cloud mass [kg]
R_cloud   = (0.1 * u.pc).to(u.m).value                    # initial cloud radius [m]
c_sound   = 200.0                                          # sound speed [m/s] (~0.2 km/s for cold clouds)
t_ff_yr   = (np.sqrt(3 * np.pi / (32 * G_si * (M_cloud /
             (4/3 * np.pi * R_cloud**3))))
             / (u.yr).to(u.s))                             # free-fall time [yr]

print(f"Cloud mass      : 1.0 M☉")
print(f"Cloud radius    : 0.1 pc  =  {R_cloud/u.au.to(u.m):.0f} AU")
print(f"Sound speed     : {c_sound:.0f} m/s")
print(f"Free-fall time  : {t_ff_yr:.0f} yr")

# ── Simulation parameters ─────────────────────────────────────────────────────
N_PARTICLES = 300          # number of gas particles
N_FRAMES    = 1080         # animation frames (60 sec at 18 fps)
SOFTENING   = R_cloud * 0.05   # gravitational softening length [m]
DT_FRAC     = 0.012        # timestep as fraction of free-fall time
dt          = DT_FRAC * t_ff_yr * (u.yr).to(u.s)   # timestep [s]

# Work in normalised units: length in R_cloud, mass in M_cloud, time in t_ff
# This keeps numbers friendly for numpy
L  = R_cloud
M  = M_cloud
T  = t_ff_yr * (u.yr).to(u.s)

# ── Initial conditions ────────────────────────────────────────────────────────
# Bonnor-Ebert-like: uniform sphere with small random velocities
r_sphere   = R_cloud * np.cbrt(np.random.rand(N_PARTICLES))
theta      = np.arccos(1 - 2*np.random.rand(N_PARTICLES))
phi        = 2 * np.pi * np.random.rand(N_PARTICLES)

pos = np.column_stack([
    r_sphere * np.sin(theta) * np.cos(phi),
    r_sphere * np.sin(theta) * np.sin(phi),
    r_sphere * np.cos(theta),
])

# Small turbulent velocities (~5% of virial velocity)
v_virial = np.sqrt(G_si * M_cloud / R_cloud)
vel = (np.random.randn(N_PARTICLES, 3) * 0.05 * v_virial +
       # Add a tiny solid-body rotation for realism
       np.column_stack([
           -pos[:, 1] * 0.03 * v_virial / R_cloud,
            pos[:, 0] * 0.03 * v_virial / R_cloud,
            np.zeros(N_PARTICLES)
       ]))

mass = np.full(N_PARTICLES, M_cloud / N_PARTICLES)   # equal masses

# ── Physics functions ─────────────────────────────────────────────────────────

def compute_acceleration(pos, mass, softening):
    """
    Direct N-body gravitational acceleration.
    Uses softened gravity: a_i = -G * Σ_j m_j * (r_i - r_j) / (|r_ij|² + ε²)^(3/2)
    Returns array (N, 3) of accelerations [m/s²].
    """
    N   = len(pos)
    acc = np.zeros_like(pos)
    for i in range(N):
        dr  = pos - pos[i]                          # (N, 3) displacement vectors
        r2  = np.sum(dr**2, axis=1) + softening**2  # softened distance²
        r2[i] = np.inf                              # exclude self
        inv_r3 = r2**(-1.5)
        acc[i] = G_si * np.sum(
            mass[:, None] * dr * inv_r3[:, None], axis=0
        )
    return acc


def leapfrog_step(pos, vel, acc, dt):
    """Leapfrog (Störmer–Verlet) integrator — symplectic, energy-conserving."""
    vel_half = vel + 0.5 * dt * acc
    pos_new  = pos + dt * vel_half
    acc_new  = compute_acceleration(pos_new, mass, SOFTENING)
    vel_new  = vel_half + 0.5 * dt * acc_new
    return pos_new, vel_new, acc_new


# ── Jeans length calculation ─────────────────────────────────────────────────
def compute_jeans_length(pos, mass):
    """
    Calculate Jeans length: λ_J = √(π c_s² / (G ρ))
    where ρ is estimated from particle volume distribution.
    Returns Jeans length in meters.
    """
    # Estimate effective radius of cloud from particle positions
    r_particles = np.linalg.norm(pos, axis=1)
    r_eff = np.percentile(r_particles, 68)  # ~1 sigma radius
    
    if r_eff < SOFTENING:
        r_eff = R_cloud
    
    # Estimate density
    volume = (4/3) * np.pi * r_eff**3
    density = np.sum(mass) / volume
    
    # Jeans length
    jeans_length = np.sqrt(np.pi * c_sound**2 / (G_si * density))
    return jeans_length


# ── Staging: pre-compute all frames ──────────────────────────────────────────
print(f"\nPre-computing {N_FRAMES} frames …")
frames_pos = []
frames_vel = []
frames_jeans = []

acc = compute_acceleration(pos, mass, SOFTENING)

for f in range(N_FRAMES):
    if f % 20 == 0:
        print(f"  frame {f}/{N_FRAMES}")
    frames_pos.append(pos.copy())
    frames_vel.append(vel.copy())
    jeans_len = compute_jeans_length(pos, mass)
    frames_jeans.append(jeans_len)
    pos, vel, acc = leapfrog_step(pos, vel, acc, dt)

print("Done. Rendering animation …\n")

# ── Visualisation setup ───────────────────────────────────────────────────────
AU   = u.au.to(u.m)
PLOT_RANGE = R_cloud / AU * 1.4    # axis range in AU

fig, axes = plt.subplots(1, 2, figsize=(12, 6),
                          facecolor='#05060f')
ax_xy, ax_xz = axes
for ax in axes:
    ax.set_facecolor('#05060f')
    ax.tick_params(colors='#8899bb', labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor('#223344')

ax_xy.set_title('XY plane  (face-on)', color='#aabbcc', fontsize=10, pad=6)
ax_xz.set_title('XZ plane  (edge-on)', color='#aabbcc', fontsize=10, pad=6)

for ax in axes:
    ax.set_xlabel('x  [AU]', color='#6688aa', fontsize=8)
ax_xy.set_ylabel('y  [AU]', color='#6688aa', fontsize=8)
ax_xz.set_ylabel('z  [AU]', color='#6688aa', fontsize=8)

# Time annotation
time_text = fig.text(0.5, 0.96,
    f'Elapsed: 0 yr  |  Free-fall time: {t_ff_yr:.0f} yr',
    ha='center', color='#ccddeeff', fontsize=9,
    fontfamily='monospace')

# Stage colour legend
stage_text = fig.text(0.5, 0.01,
    'Stage: Molecular cloud', ha='center',
    color='#88aacc', fontsize=9)

# Jeans length annotation
jeans_text = fig.text(0.05, 0.92,
    'Jeans length: --', ha='left', color='#ffaa66', fontsize=8,
    fontfamily='monospace')

plt.tight_layout(rect=[0, 0.04, 1, 0.95])


def get_particle_style(pos_f, vel_f):
    """
    Map each particle to a size and colour based on local density proxy
    (distance from centre): closer = hotter/brighter = larger dot.
    Returns sizes, colours (RGBA).
    """
    r    = np.linalg.norm(pos_f, axis=1)
    r_au = r / AU
    v    = np.linalg.norm(vel_f, axis=1)

    # Normalise: density proxy = 1/r (clamped)
    density = 1.0 / (r_au + 0.1 * PLOT_RANGE)
    density /= density.max()

    # Colour ramp: cold blue → warm orange → white hot core
    cmap = plt.cm.inferno
    colors = cmap(density**0.6)
    colors[:, 3] = np.clip(0.3 + 0.7 * density, 0, 1)   # alpha

    # Size proportional to density
    sizes = 2 + 18 * density**1.5

    return sizes, colors


def get_stage(t_elapsed_yr):
    """Return a human-readable stage label."""
    frac = t_elapsed_yr / t_ff_yr
    if frac < 0.15:
        return "Stage: Molecular cloud — slow contraction"
    elif frac < 0.45:
        return "Stage: Dense core forming — collapse accelerating"
    elif frac < 0.70:
        return "Stage: Protostellar collapse — disk hint emerging"
    elif frac < 0.90:
        return "Stage: Protostar forming — accretion ongoing"
    else:
        return "Stage: ★  Protostar born  ★"


# ── Animation function ────────────────────────────────────────────────────────
scatters_xy = [None]
scatters_xz = [None]

def init():
    for ax in axes:
        ax.set_xlim(-PLOT_RANGE, PLOT_RANGE)
        ax.set_ylim(-PLOT_RANGE, PLOT_RANGE)
    return []


def update(frame):
    pos_f = frames_pos[frame] / AU
    vel_f = frames_vel[frame]
    sizes, colors = get_particle_style(frames_pos[frame], vel_f)

    # Clear previous scatter
    for sc in scatters_xy:
        if sc:
            sc.remove()
    for sc in scatters_xz:
        if sc:
            sc.remove()

    scatters_xy[0] = ax_xy.scatter(pos_f[:, 0], pos_f[:, 1],
                                    s=sizes, c=colors, linewidths=0)
    scatters_xz[0] = ax_xz.scatter(pos_f[:, 0], pos_f[:, 2],
                                    s=sizes, c=colors, linewidths=0)

    t_yr = frame * dt / (u.yr).to(u.s)
    jeans_len = frames_jeans[frame]
    jeans_au = jeans_len / AU
    r_cloud_au = R_cloud / AU
    
    time_text.set_text(
        f'Elapsed: {t_yr:.0f} yr  |  Free-fall time: {t_ff_yr:.0f} yr  '
        f'({100*t_yr/t_ff_yr:.0f}% collapse)'
    )
    stage_text.set_text(get_stage(t_yr))
    jeans_text.set_text(
        f'Jeans λ_J = {jeans_au:.1f} AU  |  Cloud R = {r_cloud_au:.1f} AU'
    )
    return [scatters_xy[0], scatters_xz[0]]


ani = animation.FuncAnimation(
    fig, update, frames=N_FRAMES, init_func=init,
    blit=False, interval=60
)

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = 'star_formation.mp4'
ani.save(OUT, writer='ffmpeg', fps=18, dpi=120)
print(f"Saved: {OUT}")
plt.close()