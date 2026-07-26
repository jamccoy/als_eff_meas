"""The two reduction stages, as functions you call rather than modules you import.

The original scripts performed the geometry fit *as a side effect of being
imported*, and did the same for the efficiency measurement. Anything that wanted
a single value from either had to re-run the whole thing.
"""

from __future__ import annotations

from dataclasses import dataclass

from .centroids import wmean_centroid
from .config import Campaign, EnergyGroup
from .efficiency import EfficiencyTable, eff_data
from .geometry import (
    ArcFit,
    ConeAngles,
    arc_angles,
    check_arc_angles,
    combine_xy_data,
    fit_arc_to_circle,
)
from .io import dark_val, gather_file, gather_files, read_scan
from .stats import list_tot_err

__all__ = ["Geometry", "fit_geometry", "measure_efficiency"]


@dataclass(frozen=True)
class Geometry:
    """Result of the arc-fitting stage.

    ``checks`` holds each angle re-derived from the others. Comparing it against
    ``angles`` is the standing consistency test on the fit -- with the caveat,
    documented in :func:`als_eff_meas.geometry.check_arc_angles` and in
    ``docs/GEOMETRY.md``, that the ``eta`` check is not a meaningful identity.
    """

    arc: ArcFit
    angles: ConeAngles
    checks: ConeAngles
    throw: object
    zero_x: object
    zero_y: object
    direct_x: object

    @property
    def summary(self) -> dict:
        return {
            "radius": self.arc.radius,
            "x center": self.arc.x_center,
            "y center": self.arc.y_center,
            **self.angles.degrees,
        }

    @property
    def check_summary(self) -> dict:
        """The re-derived angles, in degrees. ``roll`` has no independent check."""
        deg = self.checks.degrees
        return {name: deg[name] for name in ("eta", "gamma", "alpha", "yaw")}


def _direct_beam_x(campaign: Campaign):
    """Centroid of the direct beam, averaged over the direct-beam scans."""
    directory = campaign.data / campaign.direct_beam_dir
    centroids = []
    for name in campaign.direct_beam_files:
        upos, uflux = read_scan(directory / name, "linear")
        centroids.append(wmean_centroid(upos, uflux))
    return list_tot_err(centroids)


def fit_geometry(campaign: Campaign) -> Geometry:
    """Map the diffraction arc and recover the mount angles.

    Centroid the orders at each energy, fit a circle to the resulting point
    cloud, then invert the arc parameters against the direct-beam and
    zero-order positions.
    """
    campaign.check_data_present()
    scans_dir = campaign.yaw_scans_dir

    x_files = [gather_file(s.x_scan, scans_dir) for s in campaign.yaw_scans]
    y_files = [gather_files(s.y_scans, scans_dir) for s in campaign.yaw_scans]
    meas_x = [s.meas_x for s in campaign.yaw_scans]
    orders = [s.orders for s in campaign.yaw_scans]
    energy = [s.energy_eV for s in campaign.yaw_scans]

    x_udata, y_udata, x_uzero, y_uzero = combine_xy_data(
        x_files, y_files, meas_x, orders, energy, campaign.throw_mm, campaign.peak_step
    )
    arc = fit_arc_to_circle(x_udata, y_udata, *campaign.arc_guess)

    zero_x, zero_y = list_tot_err(x_uzero), list_tot_err(y_uzero)
    direct_x = _direct_beam_x(campaign)

    angles = arc_angles(
        arc.radius,
        campaign.throw_mm,
        arc.x_center - direct_x,
        arc.x_center - zero_x,
        arc.y_center - zero_y,
    )
    return Geometry(
        arc=arc,
        angles=angles,
        checks=check_arc_angles(angles),
        throw=campaign.throw_mm,
        zero_x=zero_x,
        zero_y=zero_y,
        direct_x=direct_x,
    )


def _noise(group: EnergyGroup, attr: str, index: int):
    """Per-measurement noise floor, defaulting to zero.

    A scalar 0.0 broadcasts against however many orders propagate at that
    energy, which is what the original achieved with hand-written lists of
    zeros of the right length.
    """
    values = getattr(group, attr)
    if values is None:
        return 0.0
    return values[index]


def measure_efficiency(campaign: Campaign, geometry: Geometry) -> EfficiencyTable:
    """Measure absolute efficiency for every energy in every gain group.

    Any ``[efficiency] exclude`` entries in the campaign are applied last, so
    the unfiltered result is always what the raw data says.
    """
    campaign.check_data_present()
    groups = campaign.energy_groups

    table = eff_data(
        energy_list=[g.energies for g in groups],
        dark_list=[dark_val(g.dark_dig, campaign.dark_dir(g)) for g in groups],
        data_dir_list=[campaign.group_dir(g) for g in groups],
        d_nm=campaign.groove_spacing_nm,
        throw=geometry.throw,
        zero_x=geometry.zero_x,
        x_0=geometry.arc.x_center,
        R=geometry.arc.radius,
        pos_factor_list=[g.pos_factor for g in groups],
        neg_factor_list=[g.neg_factor for g in groups],
        peak_step=campaign.peak_step,
        order_noise_list=[
            [_noise(g, "order_noise", j) for j in range(len(g.energies))] for g in groups
        ],
        dir_noise_list=[
            [_noise(g, "dir_noise", j) for j in range(len(g.energies))] for g in groups
        ],
        offset=campaign.offset,
    )
    if campaign.exclusions:
        table = table.exclude(campaign.exclusions)
    return table
