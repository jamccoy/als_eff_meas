# als_eff_meas

Reduction of Advanced Light Source beamline scans to grating diffraction efficiency.

A grazing-incidence grating in the extreme off-plane (conical) mount disperses orders
onto an **arc**, not a line. Measuring efficiency therefore means finding that arc
before you can find the orders on it. This package does both:

```
photodiode scans (.dat)
        |
        v
  centroids  ------>  arc fit (radius, centre)
        |                     |
        |                     v
        |             cone geometry: gamma, alpha, roll, eta, yaw
        v                     |
  order flux / direct flux <--+
        |
        v
  efficiency vs. order vs. energy
```

## The measurement

The beamline monochromator sets photon energy; a photodiode on an x-translation stage
sweeps across the diffraction arc, and a goniometer sweeps in the orthogonal direction,
at a throw distance L ~ 235 mm. Data files are laid down in strict **direct/arc pairs**
— even-numbered file is the direct beam, odd-numbered the dispersed arc — and

```
efficiency = (order_flux - order_noise - dark) / (direct_flux - dir_noise - dark)
```

with the dark taken as the median of a shutter-closed scan and the flux summed over the
peak sample and its two neighbours. Uncertainties are carried through every step with
the `uncertainties` package.

The mount angles are recovered from the fitted arc radius, the throw, and the offsets of
the direct and zero-order spots from the arc centre. See [docs/GEOMETRY.md](docs/GEOMETRY.md)
for the relations, the sign conventions, and the two places where the original code does
something that looks wrong and is preserved anyway.

## Install

Requires Python 3.12 (pinned in `.python-version`).

```bash
python3 -m venv .venv && .venv/bin/pip install -e '.[dev]'
```

## Point it at the data

The raw beamline data is **private and not in this repository**. Tell the package where
it lives:

```bash
export ALS_DATA_ROOT=~/Documents/ALS
```

The expected layout and file format are in [docs/DATA.md](docs/DATA.md). Without
`ALS_DATA_ROOT` the package still imports and the data-independent tests still pass.

## Use

```python
from als_eff_meas import load_campaign, fit_geometry, measure_efficiency

campaign = load_campaign("ogre_laminar_I_oct21")
geometry = fit_geometry(campaign)          # arc fit + cone angles
table    = measure_efficiency(campaign, geometry)
table.to_csv("efficiency.csv")
```

```python
>>> f"{geometry.summary['radius']:.3f} mm"
'7.617+/-0.078 mm'
>>> f"{geometry.angles.degrees['yaw']:.3f} deg"
'0.850+/-0.013 deg'
>>> table.orders
[-2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
>>> f"{table.groove_spacing_nm:.1f} nm"          # re-derived, vs. 158.0+/-0.5 assumed
'157.8+/-2.7 nm'
```

Each campaign is a TOML file in `src/als_eff_meas/campaigns/` holding energies, scan
file digits, rough order positions, dark directories, the assumed throw, and the list of
bad points to exclude. Adding a dataset means adding a file, not editing code.

So far only `ogre_laminar_I_oct21` is ported — it is the one with a verified baseline.
`list_campaigns()` shows what is available, and
[docs/VERIFICATION.md](docs/VERIFICATION.md) describes how to port another without
breaking it.

## Verify

```bash
.venv/bin/python -m pytest
```

The reduction is pinned to the saved output of the last real run (October 2021, OGRE
laminar I): **twelve arc, angle and consistency values, and all 249 per-order efficiency
points**, reproduce to the precision they were originally written at. That baseline is
what makes a refactor safe — see [docs/VERIFICATION.md](docs/VERIFICATION.md), which
also records two discrepancies the comparison turned up in the original exclusion logic.

> **A green CI badge does not mean the science was verified.** CI has no beamline data,
> so it runs the 58 data-independent tests — the physics, the file parsing, the
> configuration — and *skips* the 20 that compare against the baseline. To actually
> verify the reduction you need the data, with `ALS_DATA_ROOT` set, locally.

## Data availability

The raw ALS scans are not public and are not distributed with this package. What is here
is the reduction: given data in the layout described in [docs/DATA.md](docs/DATA.md),
it produces absolute efficiency per order per energy. Nothing about the method is
specific to the gratings measured — a campaign file is a directory layout, a set of
energies and an assumed throw.

Requests for the underlying measurements should go to the author.

## Provenance

This package descends from analysis scripts written between 2016 and 2021, kept in a
private archive. Those scripts produced the published results; this is a restructuring of
them, not a reimplementation, and the numerics were carried across unchanged. Where the
original does something that looks wrong and was preserved anyway, the docstring says so
and says why — see `arc_angles`, `check_arc_angles` and `gather_file`.

## License

BSD-3-Clause.
