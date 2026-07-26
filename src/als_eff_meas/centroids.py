"""Locating peaks in a scan.

The default centroid is model-independent: a flux-weighted mean over a window
around each peak. That choice matters — an order sitting on a sloping background
is not a Gaussian, and fitting one to it biases the position.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from uncertainties import unumpy as unp

from .io import read_scan
from .stats import nominal

__all__ = [
    "theta_centroids_wmean",
    "trap_odr_centroid",
    "wmean_centroid",
    "x_centroids_wmean",
]


def wmean_centroid(upos, uflux):
    """Flux-weighted mean position, ``sum(x*I) / sum(I)``.

    Everything in the window contributes, background included, so the window
    width chosen by the caller is part of the measurement.
    """
    return np.sum(np.array(upos) * np.array(uflux)) / np.sum(np.array(uflux))


def _window(index, min_index, max_index, step, ux, uflux):
    """The +/-``step`` window around ``index``, clipped at the scan ends.

    Ported from the repeated block in ``als_functions.x_centroids_wmean``. Note
    the clipping is asymmetric: near the low end the window becomes
    ``[min_index, index+step)``, so it is wider on the high side. Preserved as
    written.
    """
    if (index - min_index) < step:
        return ux[min_index : index + step], uflux[min_index : index + step]
    if (max_index - index) < step:
        return ux[index - step : max_index], uflux[index - step : max_index]
    return ux[index - step : index + step], uflux[index - step : index + step]


def x_centroids_wmean(x_file: Path, meas_x, peak_step: int):
    """Centroid each peak in a linear scan, starting from rough guesses.

    Two passes, both from the original:

    1. Around each guess in ``meas_x``, take a +/-``peak_step`` window and find
       the sample of maximum flux. The guess only has to be close enough to put
       the true peak inside the window.
    2. Re-centre the window on that maximum and take the flux-weighted mean.

    The second pass is what makes the result insensitive to the guess: two
    different guesses that snap to the same maximum give an identical centroid.

    Returns an array of ufloats, one per entry in ``meas_x``.
    """
    ux, uflux = read_scan(x_file, "linear")
    x = unp.nominal_values(ux)
    min_index = np.argmin(abs(x - min(x)))
    max_index = np.argmin(abs(x - max(x)))
    step = peak_step

    max_pos = []
    for guess in meas_x:
        index = np.argmin(abs(guess - x))
        ux_range, uflux_range = _window(index, min_index, max_index, step, ux, uflux)
        flux_range = unp.nominal_values(uflux_range)
        max_pos.append(ux_range[np.argmin(max(flux_range) - flux_range)])

    ucentroids = []
    for pos in max_pos:
        index = np.argmin(abs(nominal(pos) - x))
        ux_range, uflux_range = _window(index, min_index, max_index, step, ux, uflux)
        ucentroids.append(wmean_centroid(ux_range, uflux_range))
    return np.array(ucentroids)


def theta_centroids_wmean(file_list):
    """Flux-weighted centroid of each goniometer scan, over the whole scan.

    No windowing here — a theta scan through a single order is one broad peak,
    so the whole scan is the window.
    """
    ucentroids = []
    for path in file_list:
        utheta, uflux = read_scan(path, "angle")
        ucentroids.append(wmean_centroid(utheta, uflux))
    return ucentroids


def _trap_func(B, x):
    from astropy.modeling import models

    return models.Trapezoid1D(amplitude=B[0], x_0=B[1], width=B[2], slope=B[3])(x) + B[4]


def trap_odr_centroid(
    theta,
    flux,
    theta_err,
    flux_err,
    *,
    centroid_guess,
    width_guess,
    slope_guess=9.1,
):
    """Model-dependent alternative: fit a trapezoid by orthogonal distance regression.

    Requires the optional ``astropy`` extra.

    .. note::

       The original ``als_functions.trap_odr_centroid`` derived
       ``centroid_guess`` and ``width_guess`` by calling ``find_theta_center``
       and ``find_theta_width`` — **neither of which is defined anywhere in the
       codebase**. The function raised ``NameError`` on any call, so it never
       contributed to a published result. The guesses are required arguments
       here rather than silently invented.

    A trapezoid is the right shape for a flat-topped beam clipped by a finite
    detector, but :func:`wmean_centroid` is what the pipeline actually uses.
    """
    from scipy import odr as sp_odr

    model = sp_odr.Model(_trap_func)
    data = sp_odr.RealData(theta, flux, sx=theta_err, sy=flux_err)
    out = sp_odr.ODR(
        data,
        model,
        beta0=[max(flux), centroid_guess, width_guess, slope_guess, 0.0],
    ).run()
    centroid = np.round(out.beta[1], 4)
    centroid_err = np.round(out.sd_beta[1], 4)
    theta_fit = np.linspace(theta[0], theta[-1], 1000)
    return centroid, centroid_err, theta_fit, _trap_func(out.beta, theta_fit)
