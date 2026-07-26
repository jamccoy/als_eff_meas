"""Pin the efficiency stage to ``results/measured/efficiency_order_*.txt``.

Those files are the per-order efficiency the original scripts wrote, at eight
decimal places. Every point this package produces must land on them.
"""

from __future__ import annotations

import numpy as np
import pytest

#: The saved files print 8 decimals, so agreement below this is exact.
PRINT_PRECISION = 1e-7

#: From results/measured_groove_spacing.txt. Re-derived from the measured
#: dispersion, so it is an independent check on the assumed 158.0 +/- 0.5 nm.
GROOVE_SPACING = (157.8, 2.7)


def _saved_orders(campaign):
    base = campaign.data / "results" / "measured"
    if not base.is_dir():
        pytest.skip(f"no saved efficiency baseline at {base}")
    out = {}
    for path in base.glob("efficiency_order_*.txt"):
        n = int(path.stem.rsplit("_", 1)[1])
        out[n] = {
            round(float(a)): float(b)
            for a, b in (line.split() for line in path.read_text().splitlines() if line.strip())
        }
    return out


def test_every_saved_point_is_reproduced(campaign, table):
    """All 249 measured points, across nine orders, to print precision."""
    saved = _saved_orders(campaign)
    assert saved, "baseline directory is empty"

    compared = 0
    worst, worst_at = 0.0, None
    for n, points in saved.items():
        assert n in table.efficiency, f"order {n} is in the baseline but not in the result"
        got = dict(zip((round(e) for e in table.energy[n]), table.efficiency[n], strict=True))
        for energy, expected in points.items():
            assert energy in got, (
                f"order {n} at {energy} eV is in the baseline but was not measured"
            )
            diff = abs(got[energy].nominal_value - expected)
            if diff > worst:
                worst, worst_at = diff, (n, energy)
            compared += 1
    assert worst < PRINT_PRECISION, f"largest deviation {worst:.2e} at order/energy {worst_at}"
    assert compared == 249, f"expected 249 baseline points, compared {compared}"


def test_no_extra_points_survive_the_exclusions(campaign, table):
    """After exclusions the result must not contain points the baseline dropped.

    This is what pins the ``[efficiency] exclude`` list. Order 7 is exempt: it
    was measured but never plotted, so it has no baseline file.
    """
    saved = _saved_orders(campaign)
    for n, points in saved.items():
        got = {round(e) for e in table.energy[n]}
        assert got == set(points), (
            f"order {n}: extra {sorted(got - set(points))}, missing {sorted(set(points) - got)}"
        )


def test_order_7_measured_but_not_in_baseline(campaign, table):
    """Order 7 propagates at the top two energies and is correctly absent from the baseline."""
    saved = _saved_orders(campaign)
    assert 7 not in saved
    assert sorted(round(e) for e in table.energy[7]) == [1180, 1200]
    assert 7 not in campaign.plot_orders


def test_groove_spacing_is_rederived_consistently(table):
    """The dispersion-derived spacing reproduces measured_groove_spacing.txt.

    It agrees with the 158.0 +/- 0.5 nm assumed during the reduction, which is
    the check that the orders were identified correctly.
    """
    nominal, err = GROOVE_SPACING
    assert np.round(table.groove_spacing_nm.nominal_value, 1) == pytest.approx(nominal)
    assert np.round(table.groove_spacing_nm.std_dev, 1) == pytest.approx(err)


def test_exclusions_are_applied_by_value_not_index(table):
    """``exclude`` targets (order, energy) pairs, so order does not matter.

    The original ``np.delete(..., 0)`` blocks only worked if every earlier
    deletion had already been applied; reordering or commenting one out
    silently removed the wrong point.
    """
    shuffled = [(-1, 460), (5, 900), (-2, 880)]
    a = table.exclude(shuffled)
    b = table.exclude(list(reversed(shuffled)))
    for n in a.efficiency:
        assert [round(e) for e in a.energy[n]] == [round(e) for e in b.energy[n]]
    assert 460 not in {round(e) for e in a.energy[-1]}
    assert 900 not in {round(e) for e in a.energy[5]}
    assert 880 not in {round(e) for e in a.energy[-2]}
