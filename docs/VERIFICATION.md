# Verification

The point of this repository is that the reduction was restructured without changing the
science. This is the evidence.

```bash
export ALS_DATA_ROOT=~/Documents/ALS
.venv/bin/python -m pytest
```

34 tests. The data-dependent ones skip cleanly when `ALS_DATA_ROOT` is missing, so a
clean clone still passes on a machine without the private data — but then it is proving
nothing about the science. **Run these with the data before trusting a refactor.**

## The baseline

October 2021, OGRE laminar grating I — the most recent campaign, and the last one the
original reduction was run against. Its saved outputs in
`$ALS_DATA_ROOT/October2021/Data/OGRE_laminar_I/results/` are what those scripts
produced, and they are treated here as ground truth.

They have to be, because **the original scripts cannot be re-run**: they hardcode data
paths that no longer exist, and part of the chain is Python 2 importing from Python 3.
The saved outputs are the only surviving record of what the code did, so the comparison
is against them rather than against a fresh run.

### Stage 1 — geometry (`arc_params.txt`, `arc_angles.txt`)

Eight numbers, each reproduced to the precision it was written at:

| | baseline | this package |
|---|---|---|
| radius | 7.62 +/- 0.08 mm | 7.6171 +/- 0.0779 |
| x centre | 3.50 +/- 0.04 mm | 3.4995 +/- 0.0375 |
| y centre | 5.04 +/- 0.10 mm | 5.0423 +/- 0.0999 |
| gamma | 1.857 +/- 0.025 deg | 1.8575 +/- 0.0247 |
| alpha | 27.2 +/- 0.5 deg | 27.2152 +/- 0.4874 |
| roll | -0.60 +/- 0.22 deg | -0.5986 +/- 0.2240 |
| eta | 1.660 +/- 0.029 deg | 1.6605 +/- 0.0287 |
| yaw | 0.850 +/- 0.013 deg | 0.8498 +/- 0.0132 |

`tests/test_geometry_baseline.py` also re-reads the files on disk and checks they still
say what the test claims, so the baseline cannot drift out from under a passing test.

### Stage 2 — efficiency (`results/measured/efficiency_order_*.txt`)

**249 points across nine orders, every one matching to better than 1e-7** — the saved
files print eight decimals, so that is exact agreement. The largest deviation is
5.0e-9, at order 4, 740 eV.

The independently re-derived groove spacing lands on **157.8 +/- 2.7 nm**, matching
`measured_groove_spacing.txt` exactly and agreeing with the 158.0 +/- 0.5 nm assumed
during the reduction. That agreement is the check that the orders were identified
correctly: if the arc fit or the order prediction were wrong, the dispersion would not
reproduce the spacing.

## Two things the comparison uncovered

Both are recorded in `src/als_eff_meas/campaigns/ogre_laminar_I_oct21.toml`.

### The exclusion comments do not match the exclusion code

`als_eff.py` removed bad points with fifteen blocks like:

```python
#delete n=-2 at 900 eV
order_no = -2; energy_index = 9
energy_dict[order_no] = np.delete(energy_dict[order_no], 0).tolist()
```

Every block deletes index **0**, not the energy named in its own comment. That happens
to be correct for orders -1, 3, 4, 5 and 6, because those deletions were listed in
ascending energy and always hit the lowest remaining point. **For order -2 it is not.**
The comments say 900, 920 and 940 eV; the code removed 880, 900 and 920.

The exclusion list in the campaign file is the set that actually ran, recovered by
diffing a full unfiltered run against `results/measured/`. It is expressed as
`(order, energy)` pairs, which have no ordering dependence — the original blocks were
only correct as a sequence, so commenting one out silently removed the wrong point from
every block after it.

### `offset = -0.05` changes exactly one point

`als_eff.py:19` shifts the predicted position of the positive orders to absorb roll. It
is almost always irrelevant, because the centroid search re-finds the true peak either
way. With `offset = 0` all 249 points still match except one: **order 6 at 1080 eV**, a
weak order near the edge of the arc, where 0.05 mm is enough to select a different local
maximum. That single point also moves the re-derived groove spacing from 157.8 to
158.0 nm.

It is a good illustration of where this reduction is fragile: not in the arithmetic, but
in which sample gets called the peak.

## Regenerating the efficiency baseline

Should you need a fresh CSV to diff against:

```python
from als_eff_meas import load_campaign, fit_geometry, measure_efficiency

campaign = load_campaign("ogre_laminar_I_oct21")
table = measure_efficiency(campaign, fit_geometry(campaign))
table.to_csv("baselines/ogre_laminar_I_oct21.csv")
```

`baselines/` is gitignored — it derives from private data. The eight geometry numbers
are small enough to live in the test file itself; a full efficiency table is not.

## Adding a campaign

Only OGRE laminar I is ported so far. To add another, translate that sample's original
configuration module into a TOML file and check it against that sample's own `results/`
directory the same way. Two cautions from doing it once:

- Those config modules have blocks inside `'''` quotes. Only what was *outside* them was
  active at the time of the saved run, so porting a commented-out energy group will not
  reproduce the baseline.
- Derive the exclusions by diffing against `results/measured/`, not by reading the
  `#delete` comments. See above.

Both are reasons to treat a new campaign as unverified until its own baseline
reproduces, rather than assuming the translation was faithful.
