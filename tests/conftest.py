"""Shared fixtures.

Two kinds of test live here. The *baseline* tests need the private beamline tree
and skip without it. The rest run on **synthetic scans** written into ``tmp_path``
by the :func:`write_scan` fixture — those exercise the real reading, centroiding
and geometry code with no data at all, which is what makes CI worth running.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from als_eff_meas import fit_geometry, load_campaign, measure_efficiency

BASELINE_CAMPAIGN = "ogre_laminar_I_oct21"


@pytest.fixture(scope="session")
def campaign():
    c = load_campaign(BASELINE_CAMPAIGN)
    if not c.data.is_dir():
        pytest.skip(f"beamline data not available at {c.data}; set ALS_DATA_ROOT")
    return c


@pytest.fixture(scope="session")
def geometry(campaign):
    with warnings.catch_warnings():
        # ODR probes |x - x0| > R while converging, which makes arc_func take
        # the square root of a negative number. Harmless and expected.
        warnings.simplefilter("ignore", RuntimeWarning)
        return fit_geometry(campaign)


@pytest.fixture(scope="session")
def table(campaign, geometry):
    return measure_efficiency(campaign, geometry)


@pytest.fixture
def write_scan(tmp_path):
    """Write a scan file in the beamline's format and return its path.

    One timestamp header line, then four whitespace-separated columns:
    position, flux, I0, energy. Only the first two are ever read, but the
    others are written so the fixture exercises the same parsing as real data.
    """

    def _write(name, position, flux, energy_eV=500.0, directory=None):
        position = np.asarray(position, dtype=float)
        flux = np.asarray(flux, dtype=float)
        directory = directory or tmp_path
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / name
        lines = ["Sunday, November 05, 2017 11:00:39 AM"]
        lines += [
            f"{p:.6f}\t{f:.6E}\t{4.122031e-1:.6E}\t{energy_eV:.3f}"
            for p, f in zip(position, flux, strict=True)
        ]
        path.write_text("\n".join(lines) + "\n")
        return path

    return _write


@pytest.fixture
def gaussian_scan(write_scan):
    """A scan file holding one Gaussian peak of known centre.

    Symmetric, so the flux-weighted centroid should land on ``centre`` exactly.
    """

    def _write(name, centre=0.0, amplitude=1.0, sigma=0.5, span=6.0, n=241, **kwargs):
        x = np.linspace(centre - span, centre + span, n)
        flux = amplitude * np.exp(-((x - centre) ** 2) / (2 * sigma**2))
        return write_scan(name, x, flux, **kwargs)

    return _write


def parse_results_file(path):
    """Read a ``name value+/-error`` results file into ``{name: (nominal, err)}``.

    This is the format ``np.savetxt(..., fmt='%.18s')`` produces for an object
    array of ufloats, and therefore the format of ``arc_params.txt`` and
    ``arc_angles.txt``.
    """
    out = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        *name, value = line.rsplit(" ", 1)
        nominal, _, err = value.partition("+/-")
        out[" ".join(name)] = (float(nominal), float(err))
    return out
