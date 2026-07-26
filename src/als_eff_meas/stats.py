"""Uncertainty bookkeeping shared by every stage of the reduction."""

from __future__ import annotations

import numpy as np
from uncertainties import ufloat


def list_tot_err(ufloat_list):
    """Combine a list of correlated measurements into one value.

    The total error is the scatter between the nominal values (treated as a
    systematic) added in quadrature with the mean of the individual measurement
    errors::

        sigma_tot = sqrt(std(nominals)**2 + mean(ufloats).std_dev**2)

    Ported verbatim from ``als_functions.list_tot_err``.
    """
    nominal_list = [item.nominal_value for item in ufloat_list]
    list_nominal = np.mean(nominal_list)
    sys_err = np.std(nominal_list)
    meas_err = np.mean(ufloat_list).std_dev
    tot_err = np.sqrt(sys_err**2 + meas_err**2)
    return ufloat(list_nominal, tot_err)


def nominal(value):
    """Nominal value of a ufloat, or the value unchanged if it is a plain number.

    Used wherever a ufloat is only being *compared* or used to pick an array
    index. ``uncertainties`` deprecated ``__lt__``/``__gt__``/``__abs__`` on
    ``AffineScalarFunc``; those operators worked on nominal values, so doing it
    explicitly is identical arithmetic that will not break on a future release.
    """
    return value.nominal_value if hasattr(value, "nominal_value") else value
