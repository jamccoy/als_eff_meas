"""Reading scans, finding files, and the instrument routines built on them.

Everything here runs on synthetic scans, so it needs no beamline data.
"""

from __future__ import annotations

import numpy as np
import pytest
from uncertainties import ufloat
from uncertainties import unumpy as unp

from als_eff_meas.constants import FLUX_REL_ERR, POSITION_ERR_DEG, POSITION_ERR_MM
from als_eff_meas.geometry import (
    find_zero_order_scans,
    measure_avg_throw,
    measure_disp,
    measure_graze_angle,
    radius_curv_throw,
)
from als_eff_meas.io import dark_val, gather_file, gather_files, read_scan

# --- read_scan ---------------------------------------------------------------


def test_read_scan_reads_the_first_two_columns(write_scan):
    path = write_scan("s_000001.dat", [0.0, 1.0, 2.0], [10.0, 20.0, 30.0])
    upos, uflux = read_scan(path, "linear")
    assert unp.nominal_values(upos) == pytest.approx([0.0, 1.0, 2.0])
    assert unp.nominal_values(uflux) == pytest.approx([10.0, 20.0, 30.0])


def test_read_scan_assigns_the_assumed_errors(write_scan):
    path = write_scan("s_000001.dat", [0.0, 1.0], [10.0, 20.0])
    _, uflux = read_scan(path, "linear")
    assert unp.std_devs(uflux) == pytest.approx([FLUX_REL_ERR * 10.0, FLUX_REL_ERR * 20.0])

    upos, _ = read_scan(path, "linear")
    assert unp.std_devs(upos) == pytest.approx([POSITION_ERR_MM] * 2)
    upos, _ = read_scan(path, "angle")
    assert unp.std_devs(upos) == pytest.approx([POSITION_ERR_DEG] * 2)


def test_read_scan_keeps_flux_errors_positive_for_negative_flux(write_scan):
    """Dark-corrected flux can go negative; its uncertainty must not."""
    path = write_scan("s_000001.dat", [0.0, 1.0], [-5.0, 20.0])
    _, uflux = read_scan(path, "linear")
    assert (unp.std_devs(uflux) >= 0).all()


def test_read_scan_rejects_an_unknown_data_type(write_scan):
    path = write_scan("s_000001.dat", [0.0, 1.0], [1.0, 2.0])
    with pytest.raises(ValueError, match=r"linear.*angle"):
        read_scan(path, "radial")


# --- gather_file -------------------------------------------------------------


def test_gather_file_matches_on_substring(write_scan):
    for n in (10, 11, 20):
        write_scan(f"yaw_scans_0000{n}.dat", [0.0], [1.0])
    assert gather_file(11, write_scan("t.dat", [0.0], [1.0]).parent).name == "yaw_scans_000011.dat"


def test_gather_file_is_deterministic_when_several_names_match(write_scan):
    """Substring matching is ambiguous by design; sorting at least makes it stable.

    ``3`` matches both of these. The documented behaviour is the first in sorted
    order, and a campaign relying on that is a campaign with a bug.
    """
    write_scan("scan_000031.dat", [0.0], [1.0])
    directory = write_scan("scan_000003.dat", [0.0], [1.0]).parent
    assert gather_file(3, directory).name == "scan_000003.dat"


def test_gather_file_ignores_other_extensions(write_scan):
    write_scan("scan_000030.txt", [0.0], [1.0])
    directory = write_scan("scan_000030.dat", [0.0], [1.0]).parent
    assert gather_file(30, directory).suffix == ".dat"


def test_gather_file_raises_a_useful_error_when_nothing_matches(write_scan):
    directory = write_scan("scan_000030.dat", [0.0], [1.0]).parent
    with pytest.raises(FileNotFoundError, match="99"):
        gather_file(99, directory)


def test_gather_files_preserves_the_requested_order(write_scan):
    for n in (41, 42, 43):
        write_scan(f"scan_0000{n}.dat", [0.0], [1.0])
    directory = write_scan("t.dat", [0.0], [1.0]).parent
    got = gather_files([43, 41, 42], directory)
    assert [p.name for p in got] == ["scan_000043.dat", "scan_000041.dat", "scan_000042.dat"]


# --- dark_val ----------------------------------------------------------------


def test_dark_val_is_the_median_not_the_mean(write_scan):
    """A spike in a dark scan must not move the dark level."""
    directory = write_scan("dark_000001.dat", [0, 1, 2, 3, 4], [1.0, 1.0, 1.0, 1.0, 500.0]).parent
    assert dark_val(1, directory).nominal_value == pytest.approx(1.0)


def test_dark_val_averages_the_middle_two_when_even(write_scan):
    directory = write_scan("dark_000001.dat", [0, 1, 2, 3], [1.0, 2.0, 4.0, 8.0]).parent
    assert dark_val(1, directory).nominal_value == pytest.approx(3.0)


def test_dark_val_carries_an_uncertainty(write_scan):
    directory = write_scan("dark_000001.dat", [0, 1, 2], [10.0, 20.0, 30.0]).parent
    assert dark_val(1, directory).std_dev > 0


# --- routines built on scans -------------------------------------------------


def test_measure_graze_angle_halves_the_zero_order_offset(gaussian_scan):
    """eta = (theta_zero - theta_direct) / 2, per direct/zero-order pair."""
    pairs = [
        (gaussian_scan("d1.dat", centre=0.0), gaussian_scan("z1.dat", centre=3.2)),
        (gaussian_scan("d2.dat", centre=0.0), gaussian_scan("z2.dat", centre=3.4)),
    ]
    eta, direct = measure_graze_angle(pairs)
    assert eta.nominal_value == pytest.approx(1.65, abs=1e-6)
    assert direct.nominal_value == pytest.approx(0.0, abs=1e-6)
    # The two pairs disagree by 0.1 deg in eta, so the scatter must show up.
    assert eta.std_dev > 0.04


def test_measure_disp_recovers_a_known_order_spacing(write_scan):
    """Dispersion is the mean gap between adjacent order centroids."""
    x = np.linspace(-2.0, 10.0, 601)
    flux = sum(np.exp(-((x - c) ** 2) / (2 * 0.15**2)) for c in (0.0, 2.5, 5.0, 7.5))
    path = write_scan("arc.dat", x, flux)
    avg_disp, spec_disp = measure_disp(path, [0.0, 2.5, 5.0, 7.5], 500, peak_step=10)
    assert avg_disp.nominal_value == pytest.approx(2.5, abs=1e-3)
    assert spec_disp.nominal_value == pytest.approx(avg_disp.nominal_value / (1239.84193 / 500))


def _flat_top(width_deg=1.0, span=3.0, n=601):
    """A flat-topped beam profile with sharp shoulders, as a theta scan."""
    theta = np.linspace(-span, span, n)
    return theta, np.where(abs(theta) <= width_deg / 2, 1.0, 1e-4)


def test_radius_curv_throw_scales_with_detector_width(write_scan):
    """throw = detector_width / delta_theta, so it is linear in the width."""
    theta, flux = _flat_top()
    path = write_scan("throw.dat", theta, flux)
    one = radius_curv_throw(path, curv_param=0.5, detector_width=ufloat(10.0, 0.01))
    two = radius_curv_throw(path, curv_param=0.5, detector_width=ufloat(20.0, 0.02))
    assert one.nominal_value > 0
    assert two.nominal_value == pytest.approx(2 * one.nominal_value)


def test_radius_curv_throw_raises_rather_than_inventing_a_throw(write_scan):
    """The original silently returned 235 +/- 0.5 mm when the fit found no edges.

    That is a *different* uncertainty from the assumed value, reached by
    failure, and indistinguishable downstream from a real measurement.
    """
    theta = np.linspace(-3, 3, 601)
    path = write_scan("smooth.dat", theta, np.exp(-(theta**2)))
    with pytest.raises(ValueError, match="curv_param"):
        radius_curv_throw(path, curv_param=1e9, detector_width=ufloat(10.0, 0.01))


def test_measure_avg_throw_combines_several_scans(write_scan):
    theta, flux = _flat_top()
    paths = [write_scan(f"t{i}.dat", theta, flux) for i in range(3)]
    got = measure_avg_throw(paths, curv_param=0.5, detector_width=ufloat(10.0, 0.01))
    single = radius_curv_throw(paths[0], curv_param=0.5, detector_width=ufloat(10.0, 0.01))
    assert got.nominal_value == pytest.approx(single.nominal_value)


# --- find_zero_order_scans ---------------------------------------------------


def test_find_zero_order_scans_picks_the_centroid_nearest_zero():
    files = [["a", "b", "c"]]
    got = find_zero_order_scans(
        files, meas_x=[[-3.6, 0.11, 3.85]], energy=[500], orders=[[0, 1, 2]]
    )
    assert got == ["b"]


def test_find_zero_order_scans_tests_only_the_first_order_listed():
    """``orders[i][0] == 0`` -- a list starting at -1 yields nothing, even
    though zero order is present. Preserved from the original."""
    files = [["a", "b", "c"]]
    got = find_zero_order_scans(
        files, meas_x=[[-3.6, 0.11, 3.85]], energy=[500], orders=[[-1, 0, 1]]
    )
    assert got == []
