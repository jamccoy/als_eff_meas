"""Turning arc scans into absolute diffraction efficiency.

Data files are laid down in strict pairs: even index is the direct beam, odd
index is the dispersed arc at the same energy. Efficiency is the ratio of the
two, both corrected for the dark level and for a per-measurement noise floor::

    efficiency = (order_flux - order_noise - dark) / (direct_flux - dir_noise - dark)
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp

from .centroids import wmean_centroid, x_centroids_wmean
from .constants import HC_EV_NM, WAVELENGTH_REL_ERR
from .geometry import predict_orders_new
from .io import gather_file, read_scan
from .stats import list_tot_err, nominal

__all__ = ["EfficiencyTable", "eff_data", "measure_eff_three_pt"]


@dataclass
class EfficiencyTable:
    """Efficiency indexed by diffraction order.

    ``energy[n]`` and ``efficiency[n]`` are parallel lists for order ``n``.
    ``groove_spacing_nm`` is re-derived from the measured dispersion as a
    cross-check on the value assumed during the reduction.
    """

    energy: dict = field(default_factory=lambda: defaultdict(list))
    efficiency: dict = field(default_factory=lambda: defaultdict(list))
    groove_spacing_nm: object = None

    @property
    def orders(self) -> list[int]:
        return sorted(self.efficiency)

    def exclude(self, pairs) -> EfficiencyTable:
        """Drop measurements by ``(order, energy)``, returning a new table.

        This replaces the roughly fifty hardcoded ``np.delete`` calls in the
        original scripts, which removed points by array index. An index only
        means anything if every earlier deletion is also applied, in order,
        which is why those blocks could not safely be reordered or partially
        commented out.
        """
        drop = {(int(n), float(e)) for n, e in pairs}
        out = EfficiencyTable(groove_spacing_nm=self.groove_spacing_nm)
        for n in self.efficiency:
            for e, eff in zip(self.energy[n], self.efficiency[n], strict=True):
                if (int(n), float(e)) not in drop:
                    out.energy[n].append(e)
                    out.efficiency[n].append(eff)
        return out

    def to_csv(self, path: Path) -> None:
        """Write ``order, energy_eV, efficiency, efficiency_err``, sorted."""
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["order", "energy_eV", "efficiency", "efficiency_err"])
            for n in self.orders:
                rows = sorted(
                    zip(self.energy[n], self.efficiency[n], strict=True), key=lambda r: r[0]
                )
                for e, eff in rows:
                    w.writerow([n, e, eff.nominal_value, eff.std_dev])


def _three_point_sum(uflux, index):
    """Sum the peak sample and its two neighbours.

    Three points rather than one because the stage step is comparable to the
    order width — a single sample lands at a different place on the profile for
    each order, biasing the ratio.
    """
    return uflux[index - 1] + uflux[index] + uflux[index + 1]


def measure_eff_three_pt(
    data_index,
    energy,
    dark,
    scan_dir: Path,
    d_nm,
    throw,
    zero_x,
    center_x,
    radius,
    pos_factor,
    neg_factor,
    peak_step,
    offset,
    order_noise,
    dir_noise,
):
    """Absolute efficiency at one energy, from the direct/arc scan pair.

    Returns ``(efficiency, propagating_orders, dispersion)``, with the
    efficiency an array of ufloats parallel to the order list.
    """
    wav_nm = (HC_EV_NM / (energy[data_index])) * ufloat(1, WAVELENGTH_REL_ERR)
    direct_scan = gather_file(2 * data_index, scan_dir)
    arc_scan = gather_file(2 * data_index + 1, scan_dir)

    udir_x, udir_flux = read_scan(direct_scan, "linear")
    dir_centroid = wmean_centroid(udir_x, udir_flux)
    dir_x_nom = unp.nominal_values(udir_x)
    dir_index = np.argmin(abs(nominal(dir_centroid) - dir_x_nom))
    meas_dir_flux = _three_point_sum(udir_flux, dir_index)

    uarc_x, uarc_flux = read_scan(arc_scan, "linear")
    arc_x_nom = unp.nominal_values(uarc_x)
    order_pos, prop_orders = predict_orders_new(
        wav_nm, d_nm, throw, zero_x, center_x, radius, uarc_x, pos_factor, neg_factor, offset
    )
    measured_centroids = x_centroids_wmean(arc_scan, unp.nominal_values(order_pos), peak_step)

    arc_lo = nominal(center_x) - nominal(radius)
    arc_hi = nominal(center_x) + nominal(radius)

    meas_dif_flux = []
    for centroid in measured_centroids:
        c = nominal(centroid)
        index = np.argmin(abs(c - arc_x_nom))
        # An order whose centroid falls off the arc did not propagate; it is
        # recorded at the dark level, i.e. zero efficiency, not dropped --
        # dropping it would silently shorten the order list for that energy.
        if arc_lo <= c <= arc_hi:
            meas_dif_flux.append(_three_point_sum(uarc_flux, index))
        else:
            meas_dif_flux.append(dark)

    disp = list_tot_err(np.diff(measured_centroids))
    efficiency = (np.array(meas_dif_flux) - np.array(order_noise) - dark) / (
        meas_dir_flux - dir_noise - dark
    )
    return efficiency, prop_orders, disp


def eff_data(
    energy_list,
    dark_list,
    data_dir_list,
    d_nm,
    throw,
    zero_x,
    x_0,
    R,
    pos_factor_list,
    neg_factor_list,
    peak_step,
    order_noise_list,
    dir_noise_list,
    offset=0,
) -> EfficiencyTable:
    """Run the efficiency measurement over every energy in every gain group.

    Measurements are grouped by detector gain, since each group has its own dark
    scan and its own directory. Results are re-indexed by diffraction order,
    because that is how efficiency is plotted and compared to a model.
    """
    energy_tot, eff_tot, prop_orders_tot, disp_tot = [], [], [], []
    for i in range(len(energy_list)):
        for j in range(len(energy_list[i])):
            energy_tot.append(energy_list[i][j])
            eff, prop_orders, disp = measure_eff_three_pt(
                j,
                energy_list[i],
                dark_list[i],
                data_dir_list[i],
                d_nm,
                throw,
                zero_x,
                x_0,
                R,
                pos_factor_list[i],
                neg_factor_list[i],
                peak_step,
                offset,
                order_noise_list[i][j],
                dir_noise_list[i][j],
            )
            eff_tot.append(eff)
            prop_orders_tot.append(prop_orders)
            disp_tot.append(disp)

    table = EfficiencyTable()
    for k in range(len(energy_tot)):
        for l in range(len(prop_orders_tot[k])):
            n = int(prop_orders_tot[k][l])
            table.efficiency[n].append(eff_tot[k][l])
            table.energy[n].append(float(energy_tot[k]))

    table.groove_spacing_nm = list_tot_err(
        (HC_EV_NM * throw) / (np.array(energy_tot) * np.array(disp_tot))
    )
    return table
