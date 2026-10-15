"""Plot saved Phase 2 measurements; never execute the simulation on import.

CLI: python scientific_visualization.py DATA_DIRECTORY --output NEW_DIRECTORY
     [--view com|fixed] [--animation]
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

AU_M = 149597870700.0
YEAR_S = 31557600.0


def load_dataset(directory):
    directory = Path(directory)
    report = json.loads((directory/'diagnostics.json').read_text())
    path = directory/'trajectory.npz'
    if hashlib.sha256(path.read_bytes()).hexdigest() != report['trajectory_sha256']:
        raise ValueError('Trajectory hash does not match diagnostics')
    with np.load(path, allow_pickle=False) as stored:
        data = {key: stored[key] for key in stored.files}
    validate_dataset(data, report)
    return data, report


def validate_dataset(data, report):
    p, v, m, t, rho = (np.asarray(data[k]) for k in ('pos', 'vel', 'mass', 'time_s', 'density_kg_m3'))
    if p.ndim != 3 or p.shape[2] != 3 or p.shape[0] == 0 or p.shape[1] == 0 or v.shape != p.shape:
        raise ValueError('Expected nonempty (frames,N,3) positions and velocities')
    if m.shape != (p.shape[1],) or t.shape != (p.shape[0],) or rho.shape != p.shape[:2]:
        raise ValueError('Trajectory shapes do not match')
    if any(not np.isfinite(x).all() for x in (p, v, m, t)) or (m <= 0).any() or (np.diff(t) <= 0).any():
        raise ValueError('Invalid masses, states or times')
    if np.isinf(rho).any() or np.any(rho[np.isfinite(rho)] <= 0):
        raise ValueError('Densities must be positive or explicitly undefined (NaN)')
    if len(report['snapshots']) != len(t):
        raise ValueError('Snapshot count mismatch')
    if report['units'] != 'SI':
        raise ValueError('Only SI datasets are supported')
    accuracy = report['accuracy']
    budget = accuracy['energy_budget']
    if not np.isfinite(budget) or budget <= 0:
        raise ValueError('A finite positive energy budget is required')
    if not np.isclose(accuracy['sample_start_s'], t[0], rtol=1e-12, atol=0) or not np.isclose(
            accuracy['sample_end_s'], t[-1], rtol=1e-12, atol=0):
        raise ValueError('Accuracy interval does not match saved times')
    if len(t) > 1 and not np.allclose(np.diff(t), report['configuration']['dt_s'], rtol=1e-12, atol=0):
        raise ValueError('Recorded timestep does not match saved times')
    expected_com = np.sum(m[None, :, None]*p, axis=1)/m.sum()
    stored_com = np.asarray([s['com'] for s in report['snapshots']])
    if stored_com.shape != expected_com.shape or not np.allclose(
            stored_com, expected_com, rtol=1e-12, atol=1e-12*max(1., abs(p).max())):
        raise ValueError('Reported COM does not match particle masses and positions')
    energy = np.asarray([s['energy'] for s in report['snapshots']])
    if not np.isfinite(energy).all() or energy[0] == 0:
        raise ValueError('Energy series must be finite with nonzero reference')
    error = (energy-energy[0])/abs(energy[0])
    if not np.allclose(error, accuracy['energy_error'], rtol=1e-12, atol=1e-15):
        raise ValueError('Energy history does not match reported errors')
    if not np.isclose(max(abs(error)), accuracy['peak_energy_error'], rtol=1e-12, atol=0):
        raise ValueError('Peak error mismatch')
    status = 'FAIL' if max(abs(error)) > accuracy['energy_budget'] else 'PASS (sampled interval only)'
    if accuracy['energy_status'] != status:
        raise ValueError('Inconsistent numerical accuracy label')


class ProjectionView:
    """Fixed temporal log-density colour scale; equal AU axes and explicit view origin."""

    def __init__(self, data, report, *, view='com', half_width_au=None):
        validate_dataset(data, report)
        if view not in ('com', 'fixed'):
            raise ValueError('view must be com or fixed')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        self.data, self.report, self.view = data, report, view
        com = np.asarray([s['com'] for s in report['snapshots']])
        self.projected = (data['pos']-com[:, None, :] if view == 'com' else data['pos'])/AU_M
        if half_width_au is None:
            self.half_width = float(1.05*np.max(abs(self.projected))) or 1.0
        elif not np.isfinite(half_width_au) or half_width_au <= 0:
            raise ValueError('half_width_au must be positive')
        else:
            self.half_width = float(half_width_au)
        density = data['density_kg_m3']
        finite = density[np.isfinite(density)]
        self.has_density = finite.size > 0
        if self.has_density:
            lo, hi = float(finite.min()), float(finite.max())
            if lo == hi:
                lo, hi = lo/10, hi*10  # display range only; data are untouched
            self.norm = LogNorm(lo, hi)
        else:
            self.norm = None  # grey points, no invented density scale
        self.figure, self.axes = plt.subplots(1, 2, figsize=(12, 6))
        self.figure.subplots_adjust(left=.07, right=.86, bottom=.30, top=.81, wspace=.28)
        cmap = plt.colormaps['viridis'].copy()
        cmap.set_bad('0.6')
        self.scatters = []
        for ax, vertical in zip(self.axes, ('y', 'z')):
            ax.set(xlim=(-self.half_width, self.half_width), ylim=(-self.half_width, self.half_width),
                   xlabel=f'x {"- COM_x " if view == "com" else ""}[AU]',
                   ylabel=f'{vertical} {"- COM_"+vertical+" " if view == "com" else ""}[AU]',
                   title=f'X{vertical.upper()} projection')
            ax.set_aspect('equal')
            ax.grid(alpha=.2)
            scatter = ax.scatter([], [], s=12, c=[] if self.has_density else '0.6',
                                 **({'norm': self.norm, 'cmap': cmap} if self.has_density else {}))
            self.scatters.append(scatter)
        if self.has_density:
            cax = self.figure.add_axes([.89, .28, .018, .43])
            bar = self.figure.colorbar(self.scatters[0], cax=cax)
            bar.set_label('Neighbour particle-density proxy [kg/m³]\nFixed log scale over saved interval')
        self.heading = self.figure.text(.5, .95, '', ha='center', fontsize=11)
        self.accuracy_text = self.figure.text(.5, .895, '', ha='center', fontsize=9)
        self.measurements = self.figure.text(.07, .15, '', fontsize=9)
        self.clipping = self.figure.text(.07, .08, '', fontsize=9)
        self.figure.text(.07, .025, 'Collisionless gravity only; no pressure forces. Jeans scales are hypothetical gas indicators.', fontsize=9)
        self.update(0)

    def update(self, frame):
        if not isinstance(frame, (int, np.integer)) or not 0 <= frame < len(self.projected):
            raise ValueError('Invalid frame index')
        p = self.projected[frame]
        for scatter, j in zip(self.scatters, (1, 2)):
            scatter.set_offsets(p[:, [0, j]])
            if self.has_density:
                scatter.set_array(np.ma.masked_invalid(self.data['density_kg_m3'][frame]))
        s, a = self.report['snapshots'][frame], self.report['accuracy']
        config = self.report['configuration']
        def show(value, scale=1):
            return 'undefined' if value is None else f'{value/scale:.4g}'
        self.heading.set_text(f'Collisionless N-body cloud | elapsed {self.data["time_s"][frame]/YEAR_S:,.0f} yr | '
                              f'{"COM-following" if self.view == "com" else "fixed-origin"} view')
        self.accuracy_text.set_text(
            f'Energy criterion: {a["energy_status"]} | current ΔE/|E₀|={a["energy_error"][frame]:+.3%} | '
            f'run peak={a["peak_energy_error"]:.3%} / budget={a["energy_budget"]:.1%}\n'
            f'dt={config["dt_s"]/YEAR_S:.1f} yr | ε={config["softening_m"]/AU_M:.1f} AU | '
            f'scope: {a["sample_start_s"]/YEAR_S:,.0f}–{a["sample_end_s"]/YEAR_S:,.0f} yr only')
        self.accuracy_text.set_color('firebrick' if a['energy_status'] == 'FAIL' else 'black')
        self.measurements.set_text(
            f'r₅₀={show(s["r50_m"], AU_M)} AU | RMS radius={show(s["rms_radius_m"], AU_M)} AU | '
            f'extent={show(s["extent_m"], AU_M)} AU (not a boundary)\n'
            f'Hypothetical λ_J={show(s["jeans_length_m"], AU_M)} AU; L=2r₅₀, L/λ_J={show(s["diameter_half_over_jeans"])} | '
            f'geometric 1−r₅₀/r₅₀(0)={show(s["contraction_half"])}')
        counts = [int(np.count_nonzero((abs(p[:, 0]) > self.half_width) | (abs(p[:, j]) > self.half_width))) for j in (1, 2)]
        undefined = int(np.count_nonzero(~np.isfinite(self.data['density_kg_m3'][frame])))
        self.clipping.set_text(f'Outside view: XY {counts[0]}/{len(p)}, XZ {counts[1]}/{len(p)} | '
                               f'undefined density (grey): {undefined} | no local gas-stability classification')
        return [*self.scatters, self.heading, self.accuracy_text, self.measurements, self.clipping]


def research_figure(data, report):
    import matplotlib.pyplot as plt
    t = data['time_s']/YEAR_S
    snapshots, a = report['snapshots'], report['accuracy']
    def values(key):
        return np.array([s[key] if s[key] is not None else np.nan for s in snapshots])
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), constrained_layout=True)
    ax = axes.flat
    ax[0].plot(t, values('energy'))
    ax[0].set_ylabel('Total energy [J]')
    ax[1].plot(t, a['energy_error'], label='signed current error')
    ax[1].plot(t, a['running_peak_energy_error'], label='running peak |error|')
    for sign in (-1, 1):
        ax[1].axhline(sign*a['energy_budget'], color='red', linestyle='--')
    ax[1].set_ylabel('ΔE / |E₀|; dashed = budget')
    for key, label in [('momentum_residual', '|ΔP|/(MV)'), ('angular_momentum_residual', '|ΔL|/(MRV)'), ('com_residual', 'COM residual/R')]:
        ax[2].plot(t, a[key], label=label)
    ax[2].set_ylabel('Normalized conservation residual')
    for key in ('r50_m', 'r90_m', 'extent_m'):
        ax[3].plot(t, values(key)/AU_M, label=key.replace('_m', ''))
    ax[3].set_ylabel('COM-centred radii [AU]')
    ax[4].plot(t, values('mean_density_half_kg_m3'), label='actual mass within r50 / volume')
    ax[4].set_ylabel('Enclosed-volume density [kg/m³]')
    for key in ('virial_ratio_U', 'virial_ratio_W'):
        ax[5].plot(t, values(key), label=key)
    ax[5].set_ylabel('2 K_internal / |U| or |W|')
    for key in ('kinetic', 'kinetic_internal', 'kinetic_bulk'):
        ax[6].plot(t, values(key), label=key)
    ax[6].set_ylabel('Kinetic energy components [J]')
    for key in ('aperture_mass_fraction', 'inflow_mass_fraction', 'contraction_half'):
        ax[7].plot(t, values(key), label=key)
    ax[7].set_ylabel('Mass fractions / geometric contraction')
    ax[8].plot(t, values('diameter_half_over_jeans'))
    ax[8].axhline(1, color='grey', linestyle='--')
    ax[8].set_ylabel('Hypothetical (2 r50)/λ_J; no pressure')
    for axis in axes.flat:
        axis.set_xlabel('Elapsed time [yr]')
        axis.grid(alpha=.2)
        if axis.get_legend_handles_labels()[0]:
            axis.legend(fontsize=6)
    fig.suptitle(f'Collisionless diagnostics | Energy: {a["energy_status"]}; peak={a["peak_energy_error"]:.3%}\n'
                 'Criteria apply only to saved interval; no stellar-formation or full-duration validation', fontsize=11)
    return fig


def render(directory, output, *, view='com', animation=False, fps=5, dpi=120):
    if isinstance(fps, bool) or not isinstance(fps, int) or fps < 1:
        raise ValueError('fps must be a positive integer')
    if isinstance(dpi, bool) or not isinstance(dpi, int) or dpi < 1:
        raise ValueError('dpi must be a positive integer')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data, report = load_dataset(directory)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    display = ProjectionView(data, report, view=view)
    peak = int(np.argmax(np.abs(report['accuracy']['energy_error'])))
    for name, index in [('initial', 0), ('peak_error', peak), ('final', len(data['time_s'])-1)]:
        display.update(index)
        display.figure.savefig(output/f'{name}.png', dpi=dpi)
    fig = research_figure(data, report)
    fig.savefig(output/'research_diagnostics.png', dpi=dpi+10)
    plt.close(fig)
    if animation:
        from matplotlib.animation import FuncAnimation, PillowWriter
        indices = np.unique(np.linspace(0, len(data['time_s'])-1, min(21, len(data['time_s'])), dtype=int))
        movie = FuncAnimation(display.figure, display.update, frames=indices, blit=False)
        movie.save(output/'diagnostic_preview.gif', writer=PillowWriter(fps=fps), dpi=round(dpi*2/3))
    metadata = {'view': view, 'fps': fps, 'dpi': dpi,
                'requested_output_path': str(output.resolve()),
                'half_width_au': display.half_width,
                'density_scale': 'fixed logarithmic over ALL saved frames',
                'density_vmin_kg_m3': display.norm.vmin if display.norm else None,
                'density_vmax_kg_m3': display.norm.vmax if display.norm else None,
                'frame_selection': 'all particles; initial, peak energy error, final',
                'animation_frames': len(indices) if animation else 0,
                'visualization_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'source_diagnostics_sha256': hashlib.sha256((Path(directory)/'diagnostics.json').read_bytes()).hexdigest()}
    (output/'visualization.json').write_text(json.dumps(metadata, indent=2)+'\n')
    plt.close(display.figure)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--view', choices=('com', 'fixed'), default='com')
    parser.add_argument('--animation', action='store_true', help='Save a short <=21-frame diagnostic GIF')
    args = parser.parse_args()
    render(args.data, args.output, view=args.view, animation=args.animation)
