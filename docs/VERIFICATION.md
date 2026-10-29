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

## The second and third campaigns: TASTE, December 2018

`taste_dec18_round3a` (440-800 eV) and `taste_dec18_round3b` (80-420 eV) are the two
halves of the campaign published as McCoy, McEntaffer & Miles, ApJ 891, 114 (2020). They
are separate campaigns because each has its own arc fit: the low-energy half was taken
with the beamline's order sorter inserted, which moved the beam on the grating, and the
paper says so.

**296 points, every one matching to better than 5e-9** -- 208 in round 3a across orders
-4 to 8, and 88 in round 3b across -2 to 5. Round 3a's eight geometry numbers reproduce
at the precision they were saved:

| | baseline | this package |
|---|---|---|
| radius | 7.04 +/- 0.12 mm | 7.0359 +/- 0.1242 |
| x centre | 2.73 +/- 0.05 mm | 2.7331 +/- 0.0478 |
| y centre | 6.20 +/- 0.14 mm | 6.1971 +/- 0.1425 |
| gamma | 1.731 +/- 0.032 deg | 1.7304 +/- 0.0323 |
| alpha | -23.4 +/- 0.6 deg | 23.3925 +/- 0.6098 |
| roll | -1.17 +/- 0.04 deg | 1.1744 +/- 0.0417 |
| eta | 1.56 +/- 0.04 deg | 1.5571 +/- 0.0363 |
| yaw | -0.687 +/- 0.012 deg | 0.6873 +/- 0.0125 |

The sign difference on alpha, roll and yaw is the legacy writer's, not a disagreement:
`arc_angles` here returns the magnitudes the original printed with a leading minus. Both
halves re-derive the groove spacing from their own dispersion -- 401 +/- 8 nm in round 3a
and 382 +/- 14 nm in round 3b -- matching each round's
`measured_groove_spacing.txt` exactly, against the 400.0 +/- 0.5 nm assumed.

### Three things this port had to recover, and one it could not

**The script on disk is not the record of the published run.** For OGRE the surviving
driver matched the saved output. Here it does not, three times over: `als_eff_3a.py`
carries thirty `#delete` blocks where the saved output reflects twenty-two,
`als_eff_3b.py` carries six that were never applied at all, and
`december_2018_taste3b.py` sets `R_correct = 1.05 # 0.85` where 1.05 cannot even
complete the reduction. All three were edited after the run that produced
`results/`, which is itself byte-identical to `results_3a_rev23`, the last of that
round's twenty-one saved revisions and the state the paper was submitted from a week
later. The campaign files hold what reproduces the saved output, and say where they
depart from the code.

**`R_correct` had to become a knob.** `radius_correction` under `[efficiency]` scales
the fitted radius between the arc fit and the order prediction, exactly where
`R = R_correct*R` sat in the driver. It is 1.0 for OGRE, which is why it did not exist
before. It is not cosmetic: it decides how many orders are counted as propagating, and
the original's hand-written per-energy noise lists have one entry per order, so a wrong
value does not give wrong numbers -- it raises a broadcast error. That makes those lists
an accidental, and strict, check on it.

**The throw is 233.0 mm, which is nobody's published number.** Table 1 gives
L = 232.0 +/- 1.4 mm, measured; `als.py` assumed 235.0 +/- 2.0 for every campaign. Both
fail the order-count check -- 232.0 admits a tenth order at 520-560 eV, 235.0 drops a
fifteenth at 740-800 eV -- and 233.0 reproduces all eight geometry numbers and all 296
efficiencies. It sits inside Table 1's uncertainty. Round 3a's efficiencies are
insensitive to it across 232.1-233.9; round 3b's are not, and pick out 233.0.

**Round 3b's direct beam is lost.** Its configuration module defines no direct-beam
block, and the `dir_x` its saved alpha and yaw imply, -0.183 +/- 0.005, matches no scan
surviving in that round. The campaign uses the round's own first direct-beam pair, at 80
and 100 eV, which gives alpha 19.24 against the saved 19.4 and yaw 0.729 against 0.730 --
close, and still not it. This costs nothing: `direct_x` feeds only alpha and yaw, and
neither enters the efficiency reduction. It is recorded because the alpha this campaign
reports should not be quoted as the mount. The paper's Table 1 is the mount.

## Adding a campaign

Three campaigns are ported. To add another, translate that sample's original
configuration module into a TOML file and check it against that sample's own `results/`
directory the same way. Cautions from doing it three times:

- Those config modules have blocks inside `'''` quotes. Only what was *outside* them was
  active at the time of the saved run, so porting a commented-out energy group will not
  reproduce the baseline.
- Derive the exclusions by diffing against `results/measured/`, not by reading the
  `#delete` comments. See above.

- The driver script may have been edited after the run that produced `results/`. Where
  the two disagree, the saved output is the record and the script is a later draft.
- Where a sample saved numbered revisions, check which one `results/` equals before
  treating it as the published state.

All of these are reasons to treat a new campaign as unverified until its own baseline
reproduces, rather than assuming the translation was faithful.
