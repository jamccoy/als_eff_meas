"""Pin the geometry stage to the saved output of the October 2021 run.

``October2021/Data/OGRE_laminar_I/results/`` holds what the original scripts
produced, written by ``np.savetxt(..., fmt='%.18s')`` — so the file itself is
the record, to the precision it was printed at. Twelve numbers in all: three
arc parameters, five mount angles, and the four re-derived consistency values.
If a change to this package moves any of them, it changed the science.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import parse_results_file

# radius, x center, y center from arc_params.txt; the five angles from
# arc_angles.txt. Reproduced here so the test states its own expectations.
ARC_PARAMS = {
    "radius": (7.62, 0.08),
    "x center": (3.50, 0.04),
    "y center": (5.04, 0.10),
}
ARC_ANGLES = {
    "gamma": (1.857, 0.025),
    "alpha": (27.2, 0.5),
    "roll": (-0.60, 0.22),
    "eta": (1.660, 0.029),
    "yaw": (0.850, 0.013),
}

#: From arc_angles_check.txt -- each angle re-derived from the others.
ARC_ANGLES_CHECK = {
    "eta": (1.652, 0.026),
    "gamma": (1.867, 0.033),
    "alpha": (27.1, 0.5),
    "yaw": (0.854, 0.023),
}


def _assert_matches(got, expected, label):
    nominal, err = expected
    digits = len(str(nominal).split(".")[1])
    assert np.round(got.nominal_value, digits) == pytest.approx(nominal), (
        f"{label}: {got.nominal_value} does not round to the baseline {nominal}"
    )
    err_digits = len(str(err).split(".")[1])
    assert np.round(got.std_dev, err_digits) == pytest.approx(err), (
        f"{label} uncertainty: {got.std_dev} does not round to the baseline {err}"
    )


@pytest.mark.parametrize("name,expected", {**ARC_PARAMS, **ARC_ANGLES}.items())
def test_matches_2021_baseline(geometry, name, expected):
    """Each value reproduces the saved run at the precision it was written."""
    _assert_matches(geometry.summary[name], expected, name)


@pytest.mark.parametrize("name,expected", ARC_ANGLES_CHECK.items())
def test_consistency_check_matches_2021_baseline(geometry, name, expected):
    """The re-derived angles reproduce arc_angles_check.txt too.

    This pins ``check_arc_angles``, which the pipeline now runs on every fit.
    Note it asserts reproduction, not correctness: ``eta`` here comes from an
    expression that takes the arcsine of an angle in radians and is not a
    meaningful identity. See docs/GEOMETRY.md.
    """
    _assert_matches(geometry.check_summary[name], expected, f"{name} check")


def test_consistency_check_agrees_with_the_fit_within_errors(geometry):
    """gamma, alpha and yaw are real identities: fit and check must overlap.

    A disagreement here means the arc fit and the angle inversion have drifted
    apart -- the failure this check exists to catch. ``eta`` is excluded for the
    reason above.
    """
    for name in ("gamma", "alpha", "yaw"):
        fitted, checked = geometry.summary[name], geometry.check_summary[name]
        separation = abs(fitted.nominal_value - checked.nominal_value)
        combined = np.hypot(fitted.std_dev, checked.std_dev)
        assert separation <= combined, (
            f"{name}: fit {fitted} and check {checked} differ by more than "
            f"their combined uncertainty"
        )


def test_baseline_files_still_say_what_this_test_claims(campaign):
    """Guard against the on-disk baseline drifting away from the values above.

    Without this, someone could regenerate the results directory and the test
    would keep passing against stale hardcoded numbers.
    """
    results = campaign.data / "results"
    on_disk = {
        **parse_results_file(results / "arc_params.txt"),
        **parse_results_file(results / "arc_angles.txt"),
    }
    for name, expected in {**ARC_PARAMS, **ARC_ANGLES}.items():
        assert on_disk[name] == expected, f"{name} changed on disk: {on_disk[name]} != {expected}"

    checks = parse_results_file(results / "arc_angles_check.txt")
    for name, expected in ARC_ANGLES_CHECK.items():
        assert checks[name] == expected, (
            f"{name} check changed on disk: {checks[name]} != {expected}"
        )


def test_geometry_is_independent_of_working_directory(campaign, geometry, tmp_path, monkeypatch):
    """The original wrote results relative to the CWD captured at import time.

    ``als_functions.py:17`` did ``home_dir = os.getcwd()`` and every function
    chdir'd around it, so importing from the wrong place sent output elsewhere.
    Nothing here should care where it runs.
    """
    import warnings

    from als_eff_meas import fit_geometry

    monkeypatch.chdir(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        again = fit_geometry(campaign)

    for name, value in geometry.summary.items():
        assert again.summary[name].nominal_value == value.nominal_value
        assert again.summary[name].std_dev == value.std_dev
    assert not list(tmp_path.iterdir()), "reduction wrote files into the working directory"
