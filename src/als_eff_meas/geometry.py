"""Arc fitting and the off-plane cone geometry.

In the extreme off-plane (conical) mount the diffracted orders do not land on a
line — they land on an arc, the intersection of the diffraction cone with the
detector plane. Recovering the mount geometry therefore happens in two steps:
fit the arc, then invert the arc parameters for the angles.

See ``docs/GEOMETRY.md`` for the relations and sign conventions.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import odr as sp_odr
from uncertainties import ufloat
from uncertainties import unumpy as unp

from .centroids import theta_centroids_wmean, wmean_centroid, x_centroids_wmean
from .constants import HC_EV_NM, REFERENCE_THROW_MM
from .io import read_scan
from .stats import list_tot_err, nominal

__all__ = [
    "ArcFit",
    "ConeAngles",
    "arc_angles",
    "arc_func",
    "arc_given_yaw",
    "check_arc_angles",
    "combine_xy_data",
    "find_zero_order_scans",
    "fit_arc_to_circle",
    "measure_avg_throw",
    "measure_disp",
    "measure_graze_angle",
    "predict_orders_new",
    "radius_curv_throw",
    "write_arc_angles",
    "write_arc_params",
]

_DEG = 180.0 / np.pi


def _scalar(value):
    """Unwrap the 0-d object array ``unp.arcsin`` returns for scalar input.

    ``uncertainties.unumpy`` functions always return arrays, so calling one on a
    single ufloat gives a 0-d array wrapping it. The value is identical; this
    just keeps :class:`ConeAngles` fields usable as plain ufloats.
    """
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return value.item()
    return value


@dataclass(frozen=True)
class ArcFit:
    """Best-fit circle through the order positions. All values in mm."""

    radius: object
    x_center: object
    y_center: object


@dataclass(frozen=True)
class ConeAngles:
    """Mount geometry recovered from the arc, in **radians**.

    ``gamma`` is the cone half-angle, ``alpha`` the azimuthal incidence angle,
    ``roll`` the grating roll about the beam, ``eta`` the graze angle, and
    ``yaw`` the yaw. Use :attr:`degrees` for the values as reported.
    """

    gamma: object
    alpha: object
    roll: object
    eta: object
    yaw: object

    @property
    def degrees(self) -> dict:
        return {
            "gamma": _DEG * self.gamma,
            "alpha": _DEG * self.alpha,
            "roll": _DEG * self.roll,
            "eta": _DEG * self.eta,
            "yaw": _DEG * self.yaw,
        }


def arc_func(B, x):
    """Upper semicircle of radius ``B[0]`` centred at ``(B[1], B[2])``."""
    return np.sqrt(B[0] ** 2 - (x - B[1]) ** 2) + B[2]


def arc_given_yaw(yaw_deg, eta_deg, throw):
    """Forward geometry: predict the arc radius from yaw and graze angle."""
    yaw = yaw_deg * np.pi / 180.0
    eta = eta_deg * np.pi / 180.0
    alpha = np.arctan(np.sin(yaw) / np.tan(eta))
    gamma = np.arcsin(np.sin(eta) / np.cos(alpha))
    return alpha, gamma, np.round(throw * np.sin(gamma), 2)


def find_zero_order_scans(y_file_list, meas_x, energy, orders):
    """Pick the theta scans that sit on zero order, for the throw measurement.

    Note the test is ``orders[i][0] == 0`` — the *first* order listed for an
    energy, not any of them. A campaign whose order list starts at -1 yields
    nothing here. Preserved as written; the throw is in practice supplied by the
    campaign rather than measured (see :func:`radius_curv_throw`).
    """
    theta_zero_files = []
    for i in range(len(energy)):
        if orders[i][0] == 0:
            zero_index = np.argmin(abs(np.array(meas_x[i])))
            theta_zero_files.append(y_file_list[i][zero_index])
    return theta_zero_files


def combine_xy_data(x_file_list, y_file_list, meas_x, orders, energy, throw, peak_step):
    """Build the (x, y) point cloud of order positions across all energies.

    ``x`` is the stage position of the order; ``y`` is ``L sin(theta)``, the
    goniometer angle projected onto the detector plane.

    The angles are passed through ``arcsin((235/L) sin(theta))`` first. With the
    usual ``L = 235 +/- 2`` mm the nominal value is unchanged, but the ratio is
    ``1.000 +/- 0.009``, so this step **inflates the angular uncertainty** by the
    fractional throw error. That is its only effect and it is deliberate: the
    throw error is correlated between x and y and this is where it enters.

    Returns ``(x_udata, y_udata, x_uzero, y_uzero)`` — the flattened point cloud
    and the zero-order positions per energy.
    """
    tot_x, tot_y, x_uzero, y_uzero = [], [], [], []
    for i in range(len(energy)):
        x_ucentroids = x_centroids_wmean(x_file_list[i], meas_x[i], peak_step)
        tot_x.append(x_ucentroids.tolist())

        theta = (np.pi / 180.0) * np.array(theta_centroids_wmean(y_file_list[i]))
        theta_fixed = unp.arcsin((REFERENCE_THROW_MM / throw) * unp.sin(theta))
        y_ucentroids = throw * unp.sin(theta_fixed)
        tot_y.append(y_ucentroids.tolist())

        for j in range(len(x_ucentroids)):
            if orders[i][j] == 0:
                zero_index = np.argmin(abs(unp.nominal_values(x_ucentroids)))
                x_uzero.append(x_ucentroids[zero_index])
                y_uzero.append(y_ucentroids[zero_index])

    x_udata = [item for sublist in tot_x for item in sublist]
    y_udata = [item for sublist in tot_y for item in sublist]
    return x_udata, y_udata, x_uzero, y_uzero


def fit_arc_to_circle(x_udata, y_udata, R_guess, x_0_guess, y_0_guess) -> ArcFit:
    """Fit the order positions to a circle by orthogonal distance regression.

    ODR rather than least squares because both coordinates carry error: x from
    the stage, y from the goniometer *and* the throw.
    """
    x_data, y_data = unp.nominal_values(x_udata), unp.nominal_values(y_udata)
    x_err, y_err = unp.std_devs(x_udata), unp.std_devs(y_udata)
    data = sp_odr.RealData(x_data, y_data, sx=x_err, sy=y_err)
    out = sp_odr.ODR(data, sp_odr.Model(arc_func), beta0=[R_guess, x_0_guess, y_0_guess]).run()
    return ArcFit(
        radius=ufloat(out.beta[0], out.sd_beta[0]),
        x_center=ufloat(out.beta[1], out.sd_beta[1]),
        y_center=ufloat(out.beta[2], out.sd_beta[2]),
    )


def arc_angles(R, throw, delta_x_dir, delta_x_zero, delta_y_zero) -> ConeAngles:
    """Invert the arc parameters for the mount angles. Returns radians.

    ``delta_x_dir`` is the direct beam's offset from the arc centre,
    ``delta_x_zero`` and ``delta_y_zero`` the zero order's::

        gamma = arcsin(R / L)
        alpha = arcsin(delta_x_dir / R)
        roll  = arcsin(delta_x_zero / R) - alpha
        eta   = arcsin(|delta_y_zero| / L * cos(roll))
        yaw   = arcsin(delta_x_dir / (L cos(roll) cos(eta)))

    .. warning::

       The ``eta`` expression groups as ``(|dy| / L) * cos(roll)``, **not**
       ``|dy| / (L cos(roll))``, which is what the symmetry with ``yaw`` would
       lead you to expect. This is how the original is written and it is what
       produced the published angles, so it is preserved exactly. The two differ
       by ``cos^2(roll)`` inside the arcsin — for the October 2021 roll of
       -0.6 deg that is 1 part in 10^4, well inside the quoted error, but it
       would not be for a large roll.
    """
    gamma = unp.arcsin(R / throw)
    alpha = unp.arcsin(delta_x_dir / R)
    roll = unp.arcsin(delta_x_zero / R) - alpha
    abs_dy = delta_y_zero if nominal(delta_y_zero) >= 0 else -delta_y_zero
    eta = unp.arcsin(abs_dy / throw * unp.cos(roll))
    yaw = unp.arcsin(delta_x_dir / (throw * unp.cos(roll) * unp.cos(eta)))
    return ConeAngles(
        gamma=_scalar(gamma),
        alpha=_scalar(alpha),
        roll=_scalar(roll),
        eta=_scalar(eta),
        yaw=_scalar(yaw),
    )


def check_arc_angles(angles: ConeAngles) -> ConeAngles:
    """Re-derive each angle from the others as a consistency check.

    .. warning::

       ``eta_check`` is computed as ``arctan(arcsin(yaw) / tan(alpha))``. Taking
       the arcsine of an angle in radians is dimensionally meaningless; the
       identity being checked is almost certainly ``arctan(sin(yaw)/tan(alpha))``,
       matching the ``alpha_check`` line directly below it. Preserved verbatim,
       because this function's output was written to ``arc_angles_check.txt``
       alongside the published numbers. **Treat ``eta`` disagreement here as
       uninformative**; the other three checks are sound.
    """
    eta_ = unp.arctan(unp.arcsin(angles.yaw) / unp.tan(angles.alpha))
    gamma_ = unp.arcsin(unp.sin(angles.eta) / unp.cos(angles.alpha))
    alpha_ = unp.arctan(unp.sin(angles.yaw) / unp.tan(angles.eta))
    yaw_ = unp.arcsin(unp.tan(angles.alpha) * unp.tan(angles.eta))
    return ConeAngles(
        gamma=_scalar(gamma_),
        alpha=_scalar(alpha_),
        roll=angles.roll,
        eta=_scalar(eta_),
        yaw=_scalar(yaw_),
    )


def predict_orders_new(
    wavelength_nm,
    groove_spacing_nm,
    throw,
    zero_order_x,
    center_x,
    radius,
    arc_x,
    pos_factor,
    neg_factor,
    offset,
):
    """Predict where orders -20..20 land, keeping only those on the measured arc.

    Positions come from the small-angle dispersion ``x_n = x_0 + n lambda L / d``
    rather than the full conical grating equation. ``pos_factor`` and
    ``neg_factor`` scale the two sides independently and ``offset`` shifts the
    positive orders, both to absorb roll.
    """
    arc_x_nom = unp.nominal_values(arc_x)
    lo, hi = min(arc_x_nom), max(arc_x_nom)
    center, R = nominal(center_x), nominal(radius)

    order_positions, propagating_orders = [], []
    for n in range(-20, 21):
        arg = (n * wavelength_nm * throw) / groove_spacing_nm + zero_order_x - center_x
        if abs(nominal(arg)) < R:
            if n <= 0:
                order_x = (
                    zero_order_x + neg_factor * (n * wavelength_nm * throw) / groove_spacing_nm
                )
            else:
                order_x = (
                    zero_order_x
                    + pos_factor * (n * wavelength_nm * throw) / groove_spacing_nm
                    - offset
                )
            x = nominal(order_x)
            if (x < hi) and (x > lo) and (x < center + R) and (x > center - R):
                order_positions.append(order_x)
                propagating_orders.append(n)
    return np.array(order_positions), np.array(propagating_orders)


def measure_graze_angle(file_list):
    """Graze angle from paired direct / zero-order theta scans.

    ``eta = (theta_zero - theta_direct) / 2`` per pair, then combined. Each entry
    of ``file_list`` is a ``(direct, zero_order)`` pair of paths.
    """
    meas_graze_deg, direct_angle_deg = [], []
    for direct_path, zero_path in file_list:
        upos, uflux = read_scan(direct_path, "angle")
        direct_theta = wmean_centroid(upos, uflux)
        direct_angle_deg.append(direct_theta)
        upos, uflux = read_scan(zero_path, "angle")
        zero_theta = wmean_centroid(upos, uflux)
        meas_graze_deg.append(0.5 * (zero_theta - direct_theta))
    return list_tot_err(meas_graze_deg), list_tot_err(direct_angle_deg)


def measure_disp(x_file, meas_x, energy, peak_step):
    """Average and specific dispersion from the spacing between order centroids."""
    ucentroids = x_centroids_wmean(x_file, meas_x, peak_step)
    avg_disp = list_tot_err(np.diff(ucentroids))
    return avg_disp, avg_disp / (HC_EV_NM / float(energy))


def radius_curv_throw(theta_scan: Path, curv_param, detector_width):
    """Measure the throw from the angular width of the flat-topped direct beam.

    Walks the discrete circumradius along the scan to find where the profile
    turns over at the beam edges, converts the angular width to a throw given
    the detector width.

    .. warning::

       This never worked in practice. Both ``als.py:17`` and the campaign files
       carry the comment "this doesn't actually work...", and every published
       run instead assumed ``L = 235 +/- 2`` mm. When the curvature threshold
       finds no edges the original silently returned ``ufloat(235., 0.5)`` — a
       *different* uncertainty from the assumed value, arrived at by failure.
       This port raises instead, so a failed measurement cannot be mistaken for
       a successful one.
    """
    utheta, uflux = read_scan(theta_scan, "angle")
    theta = unp.nominal_values(utheta)
    flux = unp.nominal_values(uflux)
    theta_err = unp.std_devs(utheta)

    R = []
    for i in range(1, len(theta) - 1):
        A = np.array([theta[i - 1], flux[i - 1]])
        B = np.array([theta[i], flux[i]])
        C = np.array([theta[i + 1], flux[i + 1]])
        a = np.linalg.norm(B - A)
        b = np.linalg.norm(C - B)
        c = np.linalg.norm(C - A)
        # Three collinear samples -- the normal case along a flat top or a flat
        # background -- give a degenerate triangle and an infinite circumradius,
        # i.e. zero curvature. That is the right answer, not an error.
        with np.errstate(divide="ignore", invalid="ignore"):
            R.append((a * b * c) / np.sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)))
    with np.errstate(divide="ignore"):
        curv = 1.0 / np.array([np.inf, *R, np.inf])
    curv = np.nan_to_num(curv, nan=0.0, posinf=0.0)

    index = np.where(np.signbit(float(curv_param) - curv))[0]
    if not np.any(index):
        raise ValueError(
            f"throw could not be measured from {Path(theta_scan).name}: no samples exceed "
            f"curv_param={curv_param}. Lower it, omit this file, or supply the throw "
            f"in the campaign."
        )
    start_index = index[0]
    metric = curv[index] * curv[index + 1] * curv[index + 2] * curv[index + 3]
    end_index = index[np.sort(metric.argsort()[:2])][1]
    theta1 = ufloat(theta[start_index], theta_err[start_index])
    theta2 = ufloat(theta[end_index], theta_err[end_index])
    delta_theta = (np.pi / 180.0) * np.diff([theta1, theta2])[0]
    return detector_width / delta_theta


def measure_avg_throw(file_list, curv_param, detector_width):
    """Average :func:`radius_curv_throw` over several scans."""
    return list_tot_err([radius_curv_throw(f, curv_param, detector_width) for f in file_list])


def write_arc_params(path: Path, fit: ArcFit) -> None:
    """Write ``arc_params.txt`` in the original's format, for continuity."""
    data = np.array(
        [["radius", "x center", "y center"], [fit.radius, fit.x_center, fit.y_center]]
    ).T
    np.savetxt(path, data, fmt="%.18s")


def write_arc_angles(path: Path, angles: ConeAngles) -> None:
    """Write ``arc_angles.txt`` in the original's format, for continuity."""
    deg = angles.degrees
    names = ["gamma", "alpha", "roll", "eta", "yaw"]
    data = np.array([names, [deg[n] for n in names]]).T
    np.savetxt(path, data, fmt="%.18s")
