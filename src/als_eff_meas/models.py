"""Modelling side: facet reflectance, roughness, and PCGrate results.

These are what a measured efficiency curve gets compared against. Nothing here
touches beamline data.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from .constants import HC_EV_NM

__all__ = [
    "nevot_croce",
    "read_pcgrate_xml",
    "refl_angle",
    "reflectance",
    "reflectance_interp",
]


def _fresnel(incidence_angle, n_complex, polarization):
    """Fresnel reflectance from the angle *from normal* and complex index.

    ``n_complex`` follows the ``n' + ik`` convention (``k > 0`` absorbing).
    ``als_functions.reflectance`` used ``n' - ik`` instead; since ``R = |r|^2``
    and conjugating the index conjugates ``r``, the reflectance is identical.
    Only the sign convention differed, not the number.
    """
    sqrt_term = np.sqrt(1 - (np.sin(incidence_angle) / n_complex) ** 2)
    r_s = (np.cos(incidence_angle) - n_complex * sqrt_term) / (
        np.cos(incidence_angle) + n_complex * sqrt_term
    )
    r_p = (sqrt_term - n_complex * np.cos(incidence_angle)) / (
        sqrt_term + n_complex * np.cos(incidence_angle)
    )
    R_s = (np.conjugate(r_s) * r_s).real
    R_p = (np.conjugate(r_p) * r_p).real
    if polarization == "s":
        return R_s
    if polarization == "p":
        return R_p
    return 0.5 * (R_s + R_p)


def reflectance(graze_angle_deg, refr_index_data, polarization="NP"):
    """Fresnel reflectance of a flat surface at a **grazing** angle.

    ``graze_angle_deg`` is measured from the *surface*, as everywhere else in
    this package; it is converted to the angle from normal internally.

    ``refr_index_data`` is the four-column CXRO-derived table from
    :func:`als_eff_meas.io.cxro_to_n_k`: ``energy, wavelength, n, k``. Returns
    ``(energy, R)`` on that table's own energy grid.

    .. warning::

       ``wavelength_order_efficiency.reflectance`` took the angle **from the
       normal** — it used ``incidence_angle_deg * pi/180`` with no ``90 -``
       — while ``als_functions.reflectance`` took it from the surface. The two
       functions have the same name and the same signature and disagree by 90
       degrees. Use :func:`reflectance_interp` for the from-normal behaviour.
    """
    incidence_angle = (90.0 - graze_angle_deg) * np.pi / 180.0
    n_complex = refr_index_data[:, 2] + 1j * refr_index_data[:, 3]
    return refr_index_data[:, 0], _fresnel(incidence_angle, n_complex, polarization)


def reflectance_interp(incidence_angle_deg, refr_index_data, wav_nm, polarization="NP"):
    """Fresnel reflectance from the **normal**, interpolated onto ``wav_nm``.

    The ``wavelength_order_efficiency`` convention. Interpolating the optical
    constants onto the measured wavelength grid is what makes this the useful
    variant for overlaying a model on data.
    """
    incidence_angle = incidence_angle_deg * np.pi / 180.0
    data_wav = refr_index_data[:, 1]
    n = interp1d(data_wav, refr_index_data[:, 2])(wav_nm)
    k = interp1d(data_wav, refr_index_data[:, 3])(wav_nm)
    return _fresnel(incidence_angle, n + 1j * k, polarization)


def refl_angle(gamma_deg, delta_deg, alpha_deg):
    """Incidence angle on a blazed groove facet, in degrees.

    ``zeta = arcsin(sin(gamma) cos(delta - alpha))`` — the projection of the
    cone half-angle onto the facet, where ``delta`` is the facet orientation and
    ``alpha`` the azimuthal incidence angle. This is the angle that goes into
    :func:`reflectance` for a blazed grating: the facet, not the substrate, is
    what the beam sees.
    """
    gamma = gamma_deg * np.pi / 180.0
    delta = delta_deg * np.pi / 180.0
    alpha = alpha_deg * np.pi / 180.0
    zeta = np.arcsin(np.sin(gamma) * np.cos(delta - alpha))
    return round((180.0 / np.pi) * zeta, 3)


def nevot_croce(energy_range, graze_angle_deg, refr_index_data, roughness_nm):
    """Nevot-Croce roughness damping factor.

    ``exp(-4 sigma^2 k0^2 sin(theta) Re sqrt(n^2 - cos^2 theta))``, which
    suppresses specular reflectance as surface roughness becomes comparable to
    the wavelength. ``roughness_nm`` is the RMS roughness sigma.
    """
    theta = graze_angle_deg * np.pi / 180.0
    k0 = 2 * np.pi / (HC_EV_NM / np.asarray(energy_range))
    data_energy = refr_index_data[:, 0]
    n_complex = refr_index_data[:, 2] + 1j * refr_index_data[:, 3]
    sqrt_term = np.sqrt(n_complex**2 - np.cos(theta) ** 2).real
    sampled = np.array(
        [sqrt_term[np.argmin(abs(e - data_energy))] for e in np.asarray(energy_range)]
    )
    return np.exp(-4.0 * (roughness_nm**2) * (k0**2) * np.sin(theta) * sampled)


def read_pcgrate_xml(path: Path, polarization: str = "NP"):
    """Read modelled per-order efficiency from a PCGrate XML results file.

    Returns ``(orders, efficiency)``, both ``dict`` keyed by wavelength index:
    ``orders[i]`` is the list of order numbers computed at scan point ``i`` and
    ``efficiency[i]`` the matching efficiencies.

    .. note::

       The original ``wavelength_order_efficiency.efficiency`` had
       ``eff_str == 'Efficiency_TM'`` — a comparison where an assignment was
       meant — so requesting TM raised ``UnboundLocalError`` and TM was never
       modelled. Fixed here; it cannot change any past result, because the
       broken path could not run.
    """
    attribute = {
        "TE": "Efficiency_TE",
        "TM": "Efficiency_TM",
        "NP": "Efficiency",
    }.get(polarization, "Efficiency")

    root = ET.parse(path).getroot()
    orders: dict = defaultdict(list)
    efficiency: dict = defaultdict(list)
    for i, scan_point in enumerate(root[1]):
        for child in scan_point[0]:
            orders[i].append(int(child.get("Order_number")))
            efficiency[i].append(float(child.get(attribute)))
    return dict(orders), dict(efficiency)
