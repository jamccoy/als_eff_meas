"""Diagnostic plotting, kept out of the compute path.

In the original, ``extract_flux_profile`` saved a PNG on *every* read — which is
why the data tree holds roughly twelve thousand of them — and the centroid
routines wrote figures as a side effect of returning numbers. Nothing here is
called unless you call it.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp

from .geometry import arc_func

__all__ = ["mark_centroid", "plot_arc_fit", "plot_efficiency", "plot_scan", "save"]


def plot_scan(upos, uflux, ax=None, label="experimental data"):
    """Errorbar plot of one scan."""
    ax = ax or plt.gca()
    ax.errorbar(
        unp.nominal_values(upos),
        unp.nominal_values(uflux),
        xerr=unp.std_devs(upos),
        yerr=unp.std_devs(uflux),
        fmt=".",
        ecolor="b",
        label=label,
    )
    return ax


def mark_centroid(ucentroid, uflux, color="r", label="centroid position", ax=None):
    """Vertical line at a centroid, with its uncertainty as a horizontal bar."""
    ax = ax or plt.gca()
    ax.axvline(x=ucentroid.nominal_value, color=color, linestyle="--", label=label)
    ax.errorbar(
        ucentroid.nominal_value,
        0.5 * max(unp.nominal_values(uflux)),
        xerr=ucentroid.std_dev,
        yerr=0,
        ecolor=color,
    )
    return ax


def plot_arc_fit(x_udata, y_udata, fit, ax=None):
    """The order point cloud with the best-fit arc and its radius uncertainty."""
    ax = ax or plt.gca()
    plot_scan(x_udata, y_udata, ax=ax, label="order positions")
    R, x0, y0 = fit.radius, fit.x_center, fit.y_center
    xs = np.linspace(
        x0.nominal_value - 0.95 * R.nominal_value,
        x0.nominal_value + 0.95 * R.nominal_value,
        10001,
    )
    beta = np.array([R.nominal_value, x0.nominal_value, y0.nominal_value])
    ax.plot(xs, arc_func(beta, xs), "k-", label="best fit arc")
    ax.fill_between(
        xs,
        arc_func(np.array([R.nominal_value - R.std_dev, beta[1], beta[2]]), xs),
        arc_func(np.array([R.nominal_value + R.std_dev, beta[1], beta[2]]), xs),
        facecolor="black",
        alpha=0.2,
    )
    ax.errorbar(
        x0.nominal_value,
        y0.nominal_value,
        xerr=x0.std_dev,
        yerr=y0.std_dev,
        fmt="o",
        ecolor="k",
        label="best fit arc centre",
    )
    ax.set_xlabel("$x$ (mm)")
    ax.set_ylabel("$y$ (mm)")
    ax.set_aspect("equal", "datalim")
    ax.legend()
    return ax


def plot_efficiency(table, orders=None, ax=None):
    """Efficiency against energy, one series per diffraction order."""
    ax = ax or plt.gca()
    for n in orders if orders is not None else table.orders:
        if n not in table.efficiency:
            continue
        e = np.array(table.energy[n], dtype=float)
        eff = np.array(table.efficiency[n])
        order = np.argsort(e)
        ax.errorbar(
            e[order],
            unp.nominal_values(eff)[order],
            yerr=unp.std_devs(eff)[order],
            fmt=".-",
            label=f"order {n}",
        )
    ax.set_xlabel("energy (eV)")
    ax.set_ylabel("absolute efficiency")
    ax.legend(fontsize="small")
    return ax


def save(path: Path, fig=None, close: bool = True) -> None:
    """Save and close, creating the directory if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    (fig or plt.gcf()).savefig(path, bbox_inches="tight")
    if close:
        plt.close(fig or plt.gcf())
