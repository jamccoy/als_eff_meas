"""Smoke tests for the plotting helpers.

These assert that the figures build and that ``save`` writes where it is told —
not that anything looks right. The point of `plots` being a separate module is
that nothing in the reduction depends on it, so this is all the coverage it
needs.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest
from uncertainties import ufloat

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from als_eff_meas.efficiency import EfficiencyTable
from als_eff_meas.geometry import ArcFit, arc_func
from als_eff_meas.plots import (
    mark_centroid,
    plot_arc_fit,
    plot_efficiency,
    plot_scan,
    save,
)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def scan():
    x = np.linspace(-3, 3, 51)
    flux = np.exp(-(x**2))
    return (
        np.array([ufloat(v, 0.01) for v in x]),
        np.array([ufloat(v, 0.001 * v) for v in flux]),
    )


def test_plot_scan_draws_one_errorbar_series(scan):
    ax = plot_scan(*scan)
    assert len(ax.containers) == 1


def test_mark_centroid_adds_a_line_at_the_centroid(scan):
    ax = plot_scan(*scan)
    mark_centroid(ufloat(0.4, 0.02), scan[1], ax=ax)
    assert any(line.get_xdata()[0] == pytest.approx(0.4) for line in ax.lines)


def test_plot_arc_fit_draws_the_arc_and_its_centre():
    x = np.linspace(-2, 8, 30)
    fit = ArcFit(ufloat(7.62, 0.08), ufloat(3.5, 0.04), ufloat(5.04, 0.1))
    y = arc_func([7.62, 3.5, 5.04], x)
    ax = plot_arc_fit(
        [ufloat(v, 0.01) for v in x],
        [ufloat(v, 0.01) for v in y],
        fit,
    )
    assert ax.get_xlabel() and ax.get_ylabel()
    assert any(line.get_label() == "best fit arc" for line in ax.lines)


def test_plot_efficiency_draws_one_series_per_order():
    t = EfficiencyTable()
    for n in (-1, 0, 1):
        t.energy[n] = [500.0, 460.0, 480.0]
        t.efficiency[n] = [ufloat(0.2, 0.01)] * 3
    ax = plot_efficiency(t)
    assert len(ax.containers) == 3
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == ["order -1", "order 0", "order 1"]


def test_plot_efficiency_honours_an_explicit_order_list():
    t = EfficiencyTable()
    for n in (-1, 0, 1):
        t.energy[n] = [500.0]
        t.efficiency[n] = [ufloat(0.2, 0.01)]
    ax = plot_efficiency(t, orders=[0, 99])  # 99 is absent and must be skipped
    assert len(ax.containers) == 1


def test_save_creates_missing_directories(tmp_path, scan):
    plot_scan(*scan)
    out = tmp_path / "figures" / "nested" / "scan.png"
    save(out)
    assert out.is_file() and out.stat().st_size > 0


def test_save_closes_the_figure(tmp_path, scan):
    plot_scan(*scan)
    save(tmp_path / "a.png")
    assert not plt.get_fignums()
