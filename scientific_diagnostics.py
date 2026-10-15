"""Phase 2 measurements for a collisionless cloud; SI numeric inputs only.

No inferred gas temperature, pressure support, star-formation stages, or changes
to particle states. Undefined estimates are represented explicitly, never floored.
"""

import json
from pathlib import Path

import numpy as np

from nbody import compute_acceleration, diagnostics as conserved_quantities


def _numeric(value):
    if hasattr(value, 'unit'):
        raise TypeError('Convert unit-bearing values explicitly to SI first')
    value = np.asarray(value, dtype=float)
    if not np.isfinite(value).all():
        raise ValueError('Inputs must be finite')
    return value


def _cloud(pos, mass, vel=None):
    pos, mass = _numeric(pos), _numeric(mass)
    if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) == 0:
        raise ValueError('Positions must have nonempty shape (N,3)')
    if mass.shape != (len(pos),) or np.any(mass <= 0):
        raise ValueError('Masses must be positive with shape (N,)')
    if vel is None:
        return pos, mass
    vel = _numeric(vel)
    if vel.shape != pos.shape:
        raise ValueError('Velocities must match positions')
    return pos, mass, vel


def _positive(value, name, zero=False):
    value = _numeric(value)
    if value.ndim != 0 or (value < 0 if zero else value <= 0):
        raise ValueError(f'{name} must be a finite {"nonnegative" if zero else "positive"} scalar')
    return float(value)


def centre_of_mass(pos, mass, vel=None):
    """Mass-weighted position [m], or (position, velocity [m/s])."""
    if vel is None:
        pos, mass = _cloud(pos, mass)
        return np.sum(mass[:, None]*pos, axis=0)/mass.sum()
    pos, mass, vel = _cloud(pos, mass, vel)
    return (np.sum(mass[:, None]*pos, axis=0)/mass.sum(),
            np.sum(mass[:, None]*vel, axis=0)/mass.sum())


def enclosed_mass_radii(pos, mass, fractions=(.1, .5, .9)):
    """Smallest COM-centred particle radius enclosing at least q of total mass.

    Discrete step quantiles, no interpolation. All particles on a tied shell are
    enclosed; actual enclosed mass can exceed the requested fraction.
    """
    pos, mass = _cloud(pos, mass)
    fractions = _numeric(fractions)
    if fractions.ndim != 1 or len(fractions) == 0 or np.any((fractions <= 0) | (fractions > 1)):
        raise ValueError('Fractions must be a nonempty vector in (0,1]')
    radius = np.linalg.norm(pos-centre_of_mass(pos, mass), axis=1)
    order = np.argsort(radius, kind='stable')
    cumulative = np.cumsum(mass[order])
    return radius[order][np.searchsorted(cumulative, fractions*cumulative[-1], side='left')]


def enclosed_mass(pos, mass, radius):
    """Mass [kg] within an inclusive COM-centred aperture of radius [m]."""
    pos, mass = _cloud(pos, mass)
    radius = _positive(radius, 'radius', zero=True)
    r = np.linalg.norm(pos-centre_of_mass(pos, mass), axis=1)
    return float(mass[r <= radius].sum())


def cloud_geometry(pos, mass):
    """RMS, discrete r10/r50/r90, extent, and actual enclosed-r50 mean density.

    rho_half=M(r<=r50)/(4*pi*r50^3/3). The extent is not a gas boundary.
    """
    pos, mass = _cloud(pos, mass)
    r = np.linalg.norm(pos-centre_of_mass(pos, mass), axis=1)
    r10, r50, r90 = enclosed_mass_radii(pos, mass)
    mh = float(mass[r <= r50].sum())
    return {'rms_radius_m': float(np.sqrt(np.sum(mass*r*r)/mass.sum())),
            'r10_m': float(r10), 'r50_m': float(r50), 'r90_m': float(r90),
            'extent_m': float(r.max()), 'mass_within_r50_kg': mh,
            'mean_density_half_kg_m3': mh/(4*np.pi*r50**3/3) if r50 > 0 else None,
            'concentration_r50_over_r90': float(r50/r90) if r90 > 0 else None}


def local_particle_density(pos, mass, k=16):
    """Neighbour mass-density proxy [kg/m^3], neighbour radius [m].

    rho_i=(k-1)/k * M_k / [(4*pi/3) d_k^3], excluding self. Equal-distance
    neighbours on the kth shell share fractional weight to total k neighbours.
    For equal-mass homogeneous Poisson sampling the (k-1)/k factor removes the
    inverse-volume bias. Unequal masses require locally representative marks;
    finite-boundary and density-gradient biases are NOT corrected. This is not SPH.
    Zero d_k -> NaN estimate, explicitly counted by snapshot(); no softening floor.
    """
    pos, mass = _cloud(pos, mass)
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 2 <= k < len(pos):
        raise ValueError('k must be an integer satisfying 2 <= k < N')
    rho, support = np.full(len(pos), np.nan), np.empty(len(pos))
    for i in range(len(pos)):
        distances = np.linalg.norm(pos-pos[i], axis=1)
        distances[i] = np.inf
        h = float(np.partition(distances, k-1)[k-1])
        support[i] = h
        if h == 0:
            continue
        interior, shell = distances < h, distances == h
        neighbour_mass = mass[interior].sum() + (k-int(interior.sum()))/int(shell.sum())*mass[shell].sum()
        rho[i] = (k-1)/k * neighbour_mass / (4*np.pi*h**3/3)
    return rho, support


def jeans_scales(density, sound_speed, *, G):
    """Hypothetical uniform-gas lambda_J [m], M_J [kg], t_ff [s].

    M_J uses a sphere of radius lambda_J/2. No pressure force exists in this model.
    """
    density = _positive(density, 'density')
    sound_speed = _positive(sound_speed, 'sound_speed')
    G = _positive(G, 'G')
    length = float(np.sqrt(np.pi*sound_speed**2/(G*density)))
    return {'jeans_length_m': length, 'jeans_mass_kg': float(np.pi/6*density*length**3),
            'free_fall_s': float(np.sqrt(3*np.pi/(32*G*density)))}


def virial_measurements(pos, vel, mass, softening, *, G):
    """Bulk/internal K and both 2K_internal/|U| and 2K_internal/|W|.

    W=sum m(r-Rcm).a=-G sum_i<j m_i*m_j*r_ij^2/(r_ij^2+eps^2)^1.5.
    W!=U for nonzero softening. Ratios do not certify equilibrium.
    """
    pos, mass, vel = _cloud(pos, mass, vel)
    d = conserved_quantities(pos, vel, mass, softening, G=G)
    internal = float(.5*np.sum(mass[:, None]*(vel-d['com_velocity'])**2))
    bulk = float(.5*mass.sum()*np.dot(d['com_velocity'], d['com_velocity']))
    acc = compute_acceleration(pos, mass, softening, G=G)
    W = float(np.sum(mass[:, None]*(pos-d['com'])*acc))
    return {**d, 'kinetic_internal': internal, 'kinetic_bulk': bulk,
            'virial_W': W, 'virial_ratio_U': 2*internal/abs(d['potential']) if d['potential'] != 0 else None,
            'virial_ratio_W': 2*internal/abs(W) if W != 0 else None}


def radial_flow(pos, vel, mass):
    """COM-frame radial velocities; negative is inflow. r=0 is undefined.

    Inflow fraction uses total mass; mean radial velocity uses only defined mass.
    """
    pos, mass, vel = _cloud(pos, mass, vel)
    com, vcom = centre_of_mass(pos, mass, vel)
    dr, dv = pos-com, vel-vcom
    r = np.linalg.norm(dr, axis=1)
    valid = r > 0
    vr = np.full(len(pos), np.nan)
    vr[valid] = np.sum(dr[valid]*dv[valid], axis=1)/r[valid]
    return {'radial_velocity_m_s': vr,
            'mean_radial_velocity_m_s': float(np.sum(mass[valid]*vr[valid])/mass[valid].sum()) if valid.any() else None,
            'inflow_mass_fraction': float(mass[vr < 0].sum()/mass.sum()),
            'radial_velocity_defined_mass_fraction': float(mass[valid].sum()/mass.sum())}


def contraction(radius, reference_radius):
    """1-r/r_initial: geometric contraction; negative for expansion; not clipped."""
    radius = _positive(radius, 'radius', zero=True)
    reference_radius = _positive(reference_radius, 'reference_radius')
    return 1-radius/reference_radius


def snapshot(pos, vel, mass, *, G, softening, sound_speed, k,
             reference_half_radius, aperture_radius):
    """Measurements and qualification flags; no physical stage classification."""
    reference_half_radius = _positive(reference_half_radius, 'reference_half_radius', zero=True)
    aperture_radius = _positive(aperture_radius, 'aperture_radius', zero=True)
    sound_speed = _positive(sound_speed, 'sound_speed')
    G = _positive(G, 'G')
    geometry = cloud_geometry(pos, mass)
    d = virial_measurements(pos, vel, mass, softening, G=G)
    rho, support = local_particle_density(pos, mass, k)
    flow = radial_flow(pos, vel, mass)
    mean_density = geometry['mean_density_half_kg_m3']
    jeans = jeans_scales(mean_density, sound_speed, G=G) if mean_density is not None else {
        'jeans_length_m': None, 'jeans_mass_kg': None, 'free_fall_s': None}
    warnings = ['Collisionless particle density proxy; no hydrodynamic gas density or pressure forces.',
                'Jeans scales are hypothetical uniform-gas indicators, not a cause or test of this collapse.',
                'Neighbour density has uncorrected boundary/gradient bias; radii about COM can span fragments.']
    if mean_density is None:
        warnings.append('Zero half-mass volume: global density and Jeans scales undefined.')
    if reference_half_radius == 0:
        warnings.append('Zero reference half-mass radius: geometric contraction undefined.')
    invalid = int(np.count_nonzero(~np.isfinite(rho)))
    if invalid:
        warnings.append(f'{invalid} local density estimates undefined (zero neighbour volume).')
    summary = {**geometry, **d, **jeans,
               'jeans_characteristic_length_m': 2*geometry['r50_m'],
               'diameter_half_over_jeans': 2*geometry['r50_m']/jeans['jeans_length_m'] if mean_density is not None else None,
               'contraction_half': contraction(geometry['r50_m'], reference_half_radius) if reference_half_radius > 0 else None,
               'aperture_radius_m': aperture_radius,
               'aperture_mass_kg': enclosed_mass(pos, mass, aperture_radius),
               'aperture_mass_fraction': enclosed_mass(pos, mass, aperture_radius)/np.sum(mass),
               'mean_radial_velocity_m_s': flow['mean_radial_velocity_m_s'],
               'inflow_mass_fraction': flow['inflow_mass_fraction'],
               'radial_velocity_defined_mass_fraction': flow['radial_velocity_defined_mass_fraction'],
               'undefined_density_count': invalid, 'warnings': warnings}
    return summary, {'local_density_kg_m3': rho, 'neighbour_radius_m': support,
                     'radial_velocity_m_s': flow['radial_velocity_m_s']}


def accuracy_history(snapshots, times, *, mass_total, reference_radius, G, energy_budget=.01):
    """Numerical error metrics over exactly the supplied sample times; never extrapolate."""
    times = _numeric(times)
    if times.ndim != 1 or len(times) != len(snapshots) or len(times) == 0 or np.any(np.diff(times) <= 0):
        raise ValueError('Matching, nonempty, strictly increasing sample times required')
    M = _positive(mass_total, 'mass_total')
    R = _positive(reference_radius, 'reference_radius')
    G = _positive(G, 'G')
    budget = _positive(energy_budget, 'energy_budget')
    first = snapshots[0]
    if first['energy'] == 0:
        raise ValueError('Relative energy error undefined for zero initial energy')
    V = np.sqrt(G*M/R)
    energies = _numeric([d['energy'] for d in snapshots])
    delta = (energies-first['energy'])/abs(first['energy'])
    p = _numeric([d['momentum'] for d in snapshots])
    l = _numeric([d['angular_momentum'] for d in snapshots])
    cm = _numeric([d['com'] for d in snapshots])
    residual = cm-cm[0]-(times-times[0])[:, None]*np.asarray(first['com_velocity'])
    return {'energy_error': delta, 'running_peak_energy_error': np.maximum.accumulate(abs(delta)),
            'peak_energy_error': float(np.max(abs(delta))), 'energy_budget': budget,
            'energy_status': 'FAIL' if np.max(abs(delta)) > budget else 'PASS (sampled interval only)',
            'momentum_residual': np.linalg.norm(p-p[0], axis=1)/(M*V),
            'angular_momentum_residual': np.linalg.norm(l-l[0], axis=1)/(M*R*V),
            'com_residual': np.linalg.norm(residual, axis=1)/R,
            'sample_start_s': float(times[0]), 'sample_end_s': float(times[-1]),
            'scope': 'Energy budget on supplied samples only; no full-duration or physical-model validation.'}


def json_ready(value):
    """Strict JSON: undefined floating estimates become null, not nonstandard NaN."""
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def save_diagnostics(path, payload):
    """Exclusive creation protects previous research data."""
    with Path(path).open('x') as stream:
        json.dump(json_ready(payload), stream, indent=2, allow_nan=False)
        stream.write('\n')
