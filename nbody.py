"""Direct softened Newtonian gravity and fixed-step KDK, without rendering.

Boundaries accept finite floating-point arrays in SI: position [m], velocity
[m/s], mass [kg], acceleration [m/s^2], softening [m], dt [s], G [m^3/kg/s^2].
Consistent dimensionless inputs are also allowed with the corresponding G.
Unit-bearing quantities must be explicitly converted by the caller. No global
physical parameters, random state, trajectory corrections, or timestep selection.
"""

import numpy as np


def _array(value, name):
    if hasattr(value, "unit"):
        raise TypeError(f"{name}: explicitly convert quantities to consistent numeric units")
    result = np.asarray(value, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite values")
    return result


def _scalar(value, name, nonnegative=False):
    result = _array(value, name)
    if result.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    result = float(result)
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _state(pos, mass):
    pos, mass = _array(pos, "pos"), _array(mass, "mass")
    if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) == 0:
        raise ValueError("pos must have nonempty shape (N, 3)")
    if mass.shape != (len(pos),) or (mass <= 0).any():
        raise ValueError("mass must have shape (N,) and be strictly positive")
    return pos, mass


def _vectors(value, shape, name):
    value = _array(value, name)
    if value.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    return value


def compute_acceleration(pos, mass, softening, *, G):
    """Plummer-softened gravity [m/s^2]; preserves legacy summation order.

    a_i = G sum_(j != i) m_j (r_j-r_i) / (|r_j-r_i|^2+eps^2)^(3/2).
    Zero G is a no-gravity reference. Coincident distinct particles at eps=0
    raise FloatingPointError rather than concealing the singularity.
    """
    pos, mass = _state(pos, mass)
    softening = _scalar(softening, "softening", nonnegative=True)
    G = _scalar(G, "G", nonnegative=True)
    acc = np.zeros_like(pos)
    if G == 0:
        return acc
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        for i in range(len(pos)):
            dr = pos - pos[i]
            r2 = np.sum(dr**2, axis=1) + softening**2
            r2[i] = np.inf
            inv_r3 = r2**(-1.5)
            acc[i] = G * np.sum(mass[:, None] * dr * inv_r3[:, None], axis=0)
    return acc


def leapfrog_step(pos, vel, acc, dt, mass, softening, *, G):
    """One KDK step. Caller supplies acceleration consistent with the input state.

    Use constant dt for the usual second-order symplectic map. Negative dt is
    permitted for reversibility tests. Returns new arrays without mutating inputs.
    """
    pos, mass = _state(pos, mass)
    vel = _vectors(vel, pos.shape, "vel")
    acc = _vectors(acc, pos.shape, "acc")
    dt = _scalar(dt, "dt")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        vel_half = vel + 0.5 * dt * acc
        pos_new = pos + dt * vel_half
        acc_new = compute_acceleration(pos_new, mass, softening, G=G)
        vel_new = vel_half + 0.5 * dt * acc_new
    return pos_new, vel_new, acc_new


def potential_energy(pos, mass, softening, *, G):
    """U = -G sum_(i<j) m_i m_j / sqrt(r_ij^2+eps^2) [J]. No self terms."""
    pos, mass = _state(pos, mass)
    softening = _scalar(softening, "softening", nonnegative=True)
    G = _scalar(G, "G", nonnegative=True)
    if G == 0:
        return 0.0
    potential = 0.0
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        for i in range(len(pos) - 1):
            dr = pos[i + 1:] - pos[i]
            distance = np.sqrt(np.sum(dr**2, axis=1) + softening**2)
            potential -= G * mass[i] * np.sum(mass[i + 1:] / distance)
    return float(potential)


def diagnostics(pos, vel, mass, softening, *, G):
    """Inertial-frame K, U, E [J], P [kg m/s], L about origin [kg m^2/s],
    centre of mass [m], and COM velocity [m/s]. No values are corrected.
    Angular momentum is sum(r_i cross m_i v_i); bulk kinetic energy is retained.
    """
    pos, mass = _state(pos, mass)
    vel = _vectors(vel, pos.shape, "vel")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        kinetic = float(0.5 * np.sum(mass[:, None] * vel**2))
        potential = potential_energy(pos, mass, softening, G=G)
        momentum = np.sum(mass[:, None] * vel, axis=0)
        angular = np.sum(np.cross(pos, mass[:, None] * vel), axis=0)
        com = np.sum(mass[:, None] * pos, axis=0) / np.sum(mass)
    return {"kinetic": kinetic, "potential": potential, "energy": kinetic + potential,
            "momentum": momentum, "angular_momentum": angular, "com": com,
            "com_velocity": momentum / np.sum(mass)}
