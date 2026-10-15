"""Explicit astronomical-input to SI conversions; Astropy constants are pinned by provenance."""

import numpy as np
from astropy import constants as const
from astropy import units as u


def legacy_si_constants():
    """Match the exact conversion expressions used by the preserved script."""
    return (float(const.G.to(u.m**3 / (u.kg * u.s**2)).value),
            float((1.0 * u.M_sun).to(u.kg).value),
            float((1.0 * u.pc).to(u.m).value))


def free_fall_time(mass_kg, radius_m, gravitational_constant_si):
    density = mass_kg / (4/3 * np.pi * radius_m**3)
    return float(np.sqrt(3 * np.pi / (32 * gravitational_constant_si * density)))
