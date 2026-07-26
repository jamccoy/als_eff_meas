"""Data-independent tests. These run on a clean clone with no beamline data."""

from __future__ import annotations

import numpy as np
import pytest
from uncertainties import ufloat

from als_eff_meas import list_campaigns, load_campaign
from als_eff_meas.centroids import wmean_centroid
from als_eff_meas.efficiency import EfficiencyTable
from als_eff_meas.geometry import (
    ArcFit,
    arc_angles,
    arc_func,
    arc_given_yaw,
    check_arc_angles,
    fit_arc_to_circle,
    write_arc_angles,
    write_arc_params,
)
from als_eff_meas.io import cxro_to_n_k, read_cxro_table
from als_eff_meas.models import (
    nevot_croce,
    read_pcgrate_xml,
    refl_angle,
    reflectance,
    reflectance_interp,
)
from als_eff_meas.stats import list_tot_err, nominal
from conftest import parse_results_file


def test_campaign_loads_without_touching_disk(tmp_path, monkeypatch):
    """Loading a campaign must not read beamline data.

    The original config modules ran ``dark_val()`` and ``gather_file()`` at
    import time, so merely *naming* a dataset opened files.
    """
    monkeypatch.setenv("ALS_DATA_ROOT", str(tmp_path / "does-not-exist"))
    c = load_campaign("ogre_laminar_I_oct21")
    assert c.name == "ogre_laminar_I_oct21"
    assert c.throw_mm.nominal_value == 235.0
    assert c.throw_mm.std_dev == 2.0
    assert c.groove_spacing_nm.nominal_value == 158.0
    assert [s.energy_eV for s in c.yaw_scans] == [500, 600]
    assert len(c.exclusions) == 15
    with pytest.raises(FileNotFoundError, match="ALS_DATA_ROOT"):
        c.check_data_present()


def test_energy_group_arange_matches_numpy():
    c = load_campaign("ogre_laminar_I_oct21")
    assert c.energy_groups[0].energies == np.arange(440, 560, 20).tolist()
    assert c.energy_groups[1].energies == np.arange(560, 1220, 20).tolist()


def test_campaigns_are_discoverable():
    assert "ogre_laminar_I_oct21" in list_campaigns()


def test_unknown_campaign_names_the_alternatives():
    with pytest.raises(FileNotFoundError, match="ogre_laminar_I_oct21"):
        load_campaign("no_such_campaign")


def test_list_tot_err_combines_scatter_and_measurement_error():
    """sqrt(std(nominals)^2 + mean(errors)^2), the original's definition."""
    values = [ufloat(1.0, 0.1), ufloat(2.0, 0.1), ufloat(3.0, 0.1)]
    got = list_tot_err(values)
    assert got.nominal_value == pytest.approx(2.0)
    # The measurement term is the error on the *mean* of three independent
    # values, i.e. 0.1/sqrt(3), not 0.1.
    expected = np.sqrt(np.std([1.0, 2.0, 3.0]) ** 2 + (0.1 / np.sqrt(3)) ** 2)
    assert got.std_dev == pytest.approx(expected)


def test_list_tot_err_of_identical_values_is_just_the_measurement_error():
    got = list_tot_err([ufloat(5.0, 0.2)] * 4)
    assert got.nominal_value == pytest.approx(5.0)
    assert got.std_dev == pytest.approx(0.2)


def test_nominal_passes_plain_numbers_through():
    assert nominal(3.5) == 3.5
    assert nominal(ufloat(3.5, 0.1)) == 3.5


def test_wmean_centroid_of_a_symmetric_peak_is_its_centre():
    x = np.linspace(-5, 5, 101)
    flux = np.exp(-(x**2) / 2)
    got = wmean_centroid(
        np.array([ufloat(v, 0.01) for v in x]),
        np.array([ufloat(v, 0.001 * v) for v in flux]),
    )
    assert got.nominal_value == pytest.approx(0.0, abs=1e-9)


def test_wmean_centroid_follows_a_shifted_peak():
    x = np.linspace(-5, 15, 201)
    flux = np.exp(-((x - 4.0) ** 2) / 2)
    got = wmean_centroid(
        np.array([ufloat(v, 0.01) for v in x]),
        np.array([ufloat(v, 0.001 * v) for v in flux]),
    )
    assert got.nominal_value == pytest.approx(4.0, abs=1e-6)


def test_fit_arc_recovers_a_known_circle():
    """A circle in, the same circle out."""
    R, x0, y0 = 7.5, 3.0, 5.0
    x = np.linspace(x0 - 5, x0 + 5, 40)
    y = arc_func([R, x0, y0], x)
    ux = [ufloat(v, 0.01) for v in x]
    uy = [ufloat(v, 0.01) for v in y]
    fit = fit_arc_to_circle(ux, uy, 50.0, 0.0, 47.0)
    assert fit.radius.nominal_value == pytest.approx(R, abs=1e-6)
    assert fit.x_center.nominal_value == pytest.approx(x0, abs=1e-6)
    assert fit.y_center.nominal_value == pytest.approx(y0, abs=1e-6)


def test_arc_angles_gamma_is_arcsin_R_over_throw():
    angles = arc_angles(
        ufloat(7.62, 0.08),
        ufloat(235.0, 2.0),
        ufloat(1.0, 0.01),
        ufloat(1.0, 0.01),
        ufloat(-6.8, 0.1),
    )
    assert angles.gamma.nominal_value == pytest.approx(np.arcsin(7.62 / 235.0))


def test_arc_angles_uses_the_absolute_zero_order_offset():
    """The sign of delta_y_zero must not change eta."""
    args = (ufloat(7.62, 0.08), ufloat(235.0, 2.0), ufloat(1.0, 0.01), ufloat(1.0, 0.01))
    positive = arc_angles(*args, ufloat(6.8, 0.1))
    negative = arc_angles(*args, ufloat(-6.8, 0.1))
    assert positive.eta.nominal_value == pytest.approx(negative.eta.nominal_value)
    assert positive.eta.std_dev == pytest.approx(negative.eta.std_dev)


def test_arc_given_yaw_inverts_arc_angles():
    """Forward and inverse geometry agree on the cone half-angle."""
    yaw_deg, eta_deg, throw = 0.85, 1.66, 235.0
    _alpha, gamma, radius = arc_given_yaw(yaw_deg, eta_deg, throw)
    assert radius == pytest.approx(throw * np.sin(gamma), abs=0.005)
    assert np.degrees(gamma) == pytest.approx(np.degrees(np.arcsin(radius / throw)), abs=0.01)


def test_reflectance_is_high_at_grazing_and_low_near_normal():
    """Total external reflection is the whole reason for a grazing-incidence mount."""
    energy = np.linspace(400, 1200, 50)
    table = np.column_stack(
        [energy, 1239.84193 / energy, 1 - 1e-3 * np.ones_like(energy), 1e-4 * np.ones_like(energy)]
    )
    _, grazing = reflectance(1.0, table)
    _, steep = reflectance(45.0, table)
    assert grazing.mean() > 0.9
    assert steep.mean() < 0.01


def test_reflectance_polarizations_bracket_the_unpolarized_case():
    energy = np.linspace(400, 1200, 20)
    table = np.column_stack(
        [energy, 1239.84193 / energy, 1 - 1e-3 * np.ones_like(energy), 1e-4 * np.ones_like(energy)]
    )
    _, s = reflectance(10.0, table, "s")
    _, p = reflectance(10.0, table, "p")
    _, np_ = reflectance(10.0, table, "NP")
    assert np.allclose(np_, 0.5 * (s + p))
    assert (s >= p).all()


def test_refl_angle_reduces_to_gamma_when_facet_faces_the_beam():
    assert refl_angle(1.86, 27.0, 27.0) == pytest.approx(1.86, abs=1e-3)


def test_nevot_croce_is_unity_without_roughness_and_damps_with_it():
    energy = np.linspace(400, 1200, 20)
    table = np.column_stack(
        [energy, 1239.84193 / energy, 1 - 1e-3 * np.ones_like(energy), 1e-4 * np.ones_like(energy)]
    )
    assert np.allclose(nevot_croce(energy, 1.5, table, 0.0), 1.0)
    damped = nevot_croce(energy, 1.5, table, 2.0)
    assert (damped < 1.0).all()
    assert damped[-1] < damped[0], "damping must worsen towards shorter wavelength"


def test_efficiency_table_csv_round_trip(tmp_path):
    t = EfficiencyTable()
    t.energy[1] = [500.0, 460.0]
    t.efficiency[1] = [ufloat(0.25, 0.01), ufloat(0.26, 0.02)]
    t.energy[-1] = [500.0]
    t.efficiency[-1] = [ufloat(0.1, 0.005)]
    out = tmp_path / "eff.csv"
    t.to_csv(out)
    rows = [line.split(",") for line in out.read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["-1", "1", "1"]
    assert [float(r[1]) for r in rows] == [500.0, 460.0, 500.0]


def _consistent_arc(gamma_deg=1.857, alpha_deg=27.2, throw=235.0):
    """Arc inputs for a physically self-consistent, zero-roll configuration.

    The check identities only hold if the inputs describe a real geometry, so
    the offsets are derived from gamma and alpha rather than chosen freely::

        R  = L sin(gamma)                 dx = R sin(alpha)
        dy = -L sin(eta),  where  sin(eta) = sin(gamma) cos(alpha)

    Zero roll means the zero-order and direct offsets coincide.
    """
    gamma, alpha = np.radians(gamma_deg), np.radians(alpha_deg)
    R = throw * np.sin(gamma)
    dx = R * np.sin(alpha)
    dy = -throw * np.sin(gamma) * np.cos(alpha)  # = -L sin(eta)
    return (
        ufloat(R, 0.08),
        ufloat(throw, 2.0),
        ufloat(dx, 0.01),
        ufloat(dx, 0.01),
        ufloat(dy, 0.1),
    )


def test_check_arc_angles_recovers_the_angles_it_was_given():
    """gamma, alpha and yaw are genuine identities and round-trip exactly.

    eta is not -- ``check_arc_angles`` takes the arcsine of an angle in radians
    -- so it is excluded here and asserted against the 2021 baseline instead,
    in test_geometry_baseline.py.
    """
    angles = arc_angles(*_consistent_arc())
    checks = check_arc_angles(angles)
    for name in ("gamma", "alpha", "yaw"):
        assert getattr(checks, name).nominal_value == pytest.approx(
            getattr(angles, name).nominal_value, rel=1e-9
        ), f"{name} does not round-trip"


def test_check_arc_angles_round_trip_needs_a_consistent_geometry():
    """Guard on the test above: inconsistent offsets must *not* round-trip.

    If they did, the identities would be vacuous and the check worthless.
    """
    angles = arc_angles(
        ufloat(7.62, 0.08),
        ufloat(235.0, 2.0),
        ufloat(1.0, 0.01),
        ufloat(1.0, 0.01),
        ufloat(-6.8, 0.1),
    )
    checks = check_arc_angles(angles)
    assert checks.gamma.nominal_value != pytest.approx(angles.gamma.nominal_value, rel=1e-3)


def test_check_arc_angles_leaves_roll_alone():
    """Roll has no independent identity to check it against."""
    angles = arc_angles(*_consistent_arc())
    assert check_arc_angles(angles).roll is angles.roll


def test_write_and_reparse_arc_params(tmp_path):
    """The writers must produce the format the baseline parser reads."""
    fit = ArcFit(ufloat(7.617118, 0.077862), ufloat(3.49954, 0.03751), ufloat(5.04234, 0.09994))
    path = tmp_path / "arc_params.txt"
    write_arc_params(path, fit)
    got = parse_results_file(path)
    assert got == {"radius": (7.62, 0.08), "x center": (3.50, 0.04), "y center": (5.04, 0.10)}


def test_write_and_reparse_arc_angles(tmp_path):
    angles = arc_angles(
        ufloat(7.62, 0.08),
        ufloat(235.0, 2.0),
        ufloat(1.0, 0.01),
        ufloat(1.0, 0.01),
        ufloat(-6.8, 0.1),
    )
    path = tmp_path / "arc_angles.txt"
    write_arc_angles(path, angles)
    got = parse_results_file(path)
    assert set(got) == {"gamma", "alpha", "roll", "eta", "yaw"}
    assert got["gamma"][0] == pytest.approx(np.degrees(angles.gamma.nominal_value), abs=0.01)


def test_cxro_to_n_k_converts_delta_beta(tmp_path):
    """n = 1 - delta, k = beta, wavelength = hc / E."""
    energy = np.array([500.0, 600.0, 700.0])
    delta, beta = np.array([1e-3, 8e-4, 6e-4]), np.array([1e-4, 9e-5, 8e-5])
    src = tmp_path / "Au_CXRO.dat"
    src.write_text(
        "header\nheader\nheader\n"
        + "\n".join(f"{e} {d} {b}" for e, d, b in zip(energy, delta, beta, strict=True))
        + "\n"
    )
    table = cxro_to_n_k(src, "Au")
    assert table[:, 0] == pytest.approx(energy)
    assert table[:, 1] == pytest.approx(1239.84193 / energy)
    assert table[:, 2] == pytest.approx(1 - delta)
    assert table[:, 3] == pytest.approx(beta)
    assert not list(tmp_path.glob("Au_optical_constants*")), "must not write unless asked"


def test_cxro_to_n_k_writes_only_when_given_an_output_directory(tmp_path):
    src = tmp_path / "Au_CXRO.dat"
    src.write_text("h\nh\nh\n500.0 1e-3 1e-4\n600.0 8e-4 9e-5\n")
    out = tmp_path / "out"
    out.mkdir()
    cxro_to_n_k(src, "Au", out_dir=out)
    assert (out / "Au_optical_constants_n_k.txt").is_file()


def test_read_cxro_table_skips_the_three_header_lines(tmp_path):
    src = tmp_path / "t.dat"
    src.write_text("a\nb\nc\n500.0 1e-3 1e-4\n600.0 8e-4 9e-5\n")
    assert read_cxro_table(src).shape == (2, 3)


PCGRATE_XML = """<?xml version="1.0"?>
<Results>
  <Problem name="test"/>
  <ScanPoints>
    <Point>
      <Orders>
        <Order Order_number="-1" Efficiency="0.10" Efficiency_TE="0.11" Efficiency_TM="0.12"/>
        <Order Order_number="0" Efficiency="0.50" Efficiency_TE="0.51" Efficiency_TM="0.52"/>
      </Orders>
    </Point>
    <Point>
      <Orders>
        <Order Order_number="-1" Efficiency="0.20" Efficiency_TE="0.21" Efficiency_TM="0.22"/>
        <Order Order_number="0" Efficiency="0.60" Efficiency_TE="0.61" Efficiency_TM="0.62"/>
      </Orders>
    </Point>
  </ScanPoints>
</Results>
"""


def test_read_pcgrate_xml_indexes_orders_by_scan_point(tmp_path):
    path = tmp_path / "r.xml"
    path.write_text(PCGRATE_XML)
    orders, efficiency = read_pcgrate_xml(path)
    assert orders == {0: [-1, 0], 1: [-1, 0]}
    assert efficiency == {0: [0.10, 0.50], 1: [0.20, 0.60]}


@pytest.mark.parametrize("polarization,expected", [("TE", 0.11), ("TM", 0.12), ("NP", 0.10)])
def test_read_pcgrate_xml_selects_the_polarization(tmp_path, polarization, expected):
    """TM was unreachable in the original -- a comparison where an assignment
    was meant meant it raised UnboundLocalError."""
    path = tmp_path / "r.xml"
    path.write_text(PCGRATE_XML)
    _, efficiency = read_pcgrate_xml(path, polarization)
    assert efficiency[0][0] == pytest.approx(expected)


def test_reflectance_interp_agrees_with_reflectance_on_the_same_angle():
    """The two take the angle from opposite references; 90 - x must reconcile them."""
    energy = np.linspace(400, 1200, 40)
    table = np.column_stack(
        [energy, 1239.84193 / energy, 1 - 1e-3 * np.ones_like(energy), 1e-4 * np.ones_like(energy)]
    )
    graze = 1.5
    _, from_surface = reflectance(graze, table)
    from_normal = reflectance_interp(90.0 - graze, table, 1239.84193 / energy)
    assert from_normal == pytest.approx(from_surface, rel=1e-9)


def test_exclude_returns_a_new_table_and_leaves_the_original_alone():
    t = EfficiencyTable()
    t.energy[1] = [440.0, 460.0, 480.0]
    t.efficiency[1] = [ufloat(0.1, 0.01), ufloat(0.2, 0.01), ufloat(0.3, 0.01)]
    filtered = t.exclude([(1, 460)])
    assert [round(e) for e in filtered.energy[1]] == [440, 480]
    assert [round(e) for e in t.energy[1]] == [440, 460, 480]
