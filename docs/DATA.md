# The data

**Raw beamline data is private and is not in this repository.** It stays in its own
tree, and this repo lives outside that tree so that a stray `git add` cannot reach it.
The `.gitignore` blocks `*.dat`, `*.abs`, `*.script`, `data/` and friends as a second
line of defence.

Point the package at the data with an environment variable:

```bash
export ALS_DATA_ROOT=~/Documents/ALS
```

Every campaign path in `src/als_eff_meas/campaigns/*.toml` is relative to it, so no
absolute path appears in any source file. The scripts this package replaces hardcoded
paths belonging to four different machines, one of which pointed at a directory that
was not where the data lived.

## Layout

```
$ALS_DATA_ROOT/
└── <Campaign>/                       e.g. October2021, June2017, December2018
    ├── Data/
    │   └── <sample>/                 e.g. OGRE_laminar_I, TASTE, silicon_master
    │       ├── yaw_scans/            arc mapping: one x scan + one theta scan per order
    │       ├── zero_order_scans/     direct-beam reference scans
    │       ├── darks/                shutter-closed scans, one directory per gain
    │       │   ├── dark_seven/       10^7 gain
    │       │   ├── dark_eight/       10^8 gain
    │       │   └── dark_nine/
    │       ├── energy_scans/
    │       │   └── data_<range>eV-<gain>_gain/
    │       └── results/              outputs; the verification baseline lives here
    └── Figures/
        └── <sample>/
```

Directories are never created by the reduction — they must already exist.

## File format

Scans are `.dat` (older campaigns use `.abs`): whitespace-delimited, one header line
holding a timestamp, then four columns.

```
Sunday, November 05, 2017 11:00:39 AM
0.000000	1.501584E-2	4.122031E-1	499.968
```

| Column | Meaning |
|---|---|
| 0 | stage position (mm) for an x scan, goniometer angle (deg) for a theta scan |
| 1 | photodiode flux |
| 2 | I0, the beam-current monitor |
| 3 | actual monochromator energy (eV) |

Only columns 0 and 1 are read — by this package and by the original alike. **The
efficiency is therefore not normalised by I0**; it relies on the direct-beam scan in
each pair being taken close enough in time that the ring current has not moved. Column 3
is likewise unused, so the *nominal* energy from the campaign file is what labels each
measurement, not the monochromator readback.

## File pairing

Within an energy-scan directory, files are ordered **direct, arc, direct, arc**: the
scan for energy index `i` is at digits `2i` (direct beam) and `2i+1` (dispersed arc).

Files are found by **substring match on the filename**, not by index — `gather_file(30,
…)` returns the first sorted `.dat` whose name contains `"30"`. This is fragile by
construction, and the original code's own comment called it "finicky". It is preserved
because changing it changes which file each measurement reads. When adding a campaign,
make sure the digits you list are unambiguous within their directory.

## Reference data

`reference/` is for public CXRO optical-constant tables (`energy, delta, beta`, three
header lines), which are not measurements and may be committed. Nothing else belongs
there.

## Artefacts in the data tree

The original wrote into the data directories as a side effect of computing: a
diagnostic PNG for every scan read (roughly twelve thousand of them across the tree) and
a `test.txt` dropped into the data directory on every centroid call. If you are working
with a tree that predates this package, expect to find both. This package writes nothing
unless asked; see `als_eff_meas.plots`.
