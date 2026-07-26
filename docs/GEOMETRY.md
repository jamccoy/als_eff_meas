# Geometry and conventions

This document is **normative**. Read it before changing anything in `geometry.py` or
filing a bug about a sign.

## The mount

A grating in the extreme off-plane (conical) mount, at grazing incidence. The beam
arrives nearly along the grooves rather than across them, so the diffracted orders lie
on a **cone** about the groove direction. Where that cone meets the detector plane, the
orders lie on an **arc** — not a line, which is why the reduction has to fit a circle
before it can identify anything.

```
        y  (goniometer, projected: y = L sin theta)
        ^
        |        . -- .            <- orders on the arc
        |     .          .
        |   .              .
        | .        x        .      <- (x_c, y_c), the fitted arc centre
        +--------------------> x  (translation stage, mm)
```

## Symbols

| Symbol | Name | Meaning |
|---|---|---|
| `L` | throw | grating to detector distance, ~235 mm |
| `R` | arc radius | radius of the fitted circle, mm |
| `(x_c, y_c)` | arc centre | centre of the fitted circle, mm |
| `gamma` | cone half-angle | half-opening angle of the diffraction cone |
| `eta` | graze angle | angle between the beam and the grating surface |
| `alpha` | azimuthal incidence | in-cone azimuth of the incident beam |
| `beta_n` | azimuthal diffraction | in-cone azimuth of order *n* |
| `yaw` | yaw | rotation of the grating about the surface normal |
| `roll` | roll | rotation of the grating about the beam |
| `d` | groove spacing | nm; 158.0 +/- 0.5 for the OGRE laminar grating |

**Angles are radians inside `ConeAngles`; degrees everywhere they are reported.** Use
`ConeAngles.degrees`. Graze angles are measured **from the surface**, not from the
normal — with the single exception of `models.reflectance_interp`, which is documented
where it lives.

## Relations

Forward (`geometry.arc_given_yaw`), predicting the arc from the mount:

```
alpha = arctan( sin(yaw) / tan(eta) )
gamma = arcsin( sin(eta) / cos(alpha) )
R     = L sin(gamma)
```

Inverse (`geometry.arc_angles`), recovering the mount from the fitted arc. `dx_dir` is
the direct beam's offset from the arc centre, `dx_zero` and `dy_zero` the zero order's:

```
gamma = arcsin( R / L )
alpha = arcsin( dx_dir / R )
roll  = arcsin( dx_zero / R ) - alpha
eta   = arcsin( |dy_zero| / L * cos(roll) )
yaw   = arcsin( dx_dir / (L cos(roll) cos(eta)) )
```

> **The `eta` line groups as `(|dy_zero| / L) * cos(roll)`, not `|dy_zero| / (L
> cos(roll))`.** The symmetry with the `yaw` line below it suggests the latter was
> intended. The former is what the original computed and what produced the published
> angles, so it is preserved exactly. The two differ by `cos^2(roll)` inside the
> arcsine: negligible at the October 2021 roll of -0.6 degrees (1 part in 10^4, far
> inside the quoted error), but not negligible for a large roll. If you ever measure a
> grating with several degrees of roll, revisit this before trusting `eta`.

## Grating equation

In the conical mount:

```
sin(alpha) + sin(beta_n) = n lambda / (d sin(gamma))
```

`predict_orders` and `predict_orders_alt` use this directly. The routine actually in
service, `predict_orders_new`, instead uses the small-angle dispersion along the arc:

```
x_n = x_0 + n lambda L / d
```

with independent scale factors for the positive and negative orders and an `offset` on
the positive side, both absorbing roll. This is an approximation, and it is adequate
only because the predicted position is a *starting guess*: `x_centroids_wmean` re-finds
the true peak within `peak_step` samples. A prediction that lands on the wrong peak is
the one failure mode that matters, and it is why `peak_step` is a campaign parameter.

## Consistency check

`check_arc_angles` re-derives each angle from the others:

```
gamma_check = arcsin( sin(eta) / cos(alpha) )
alpha_check = arctan( sin(yaw) / tan(eta) )
yaw_check   = arcsin( tan(alpha) tan(eta) )
eta_check   = arctan( arcsin(yaw) / tan(alpha) )     <- see below
```

> **`eta_check` takes the arcsine of an angle in radians**, which is meaningless. By
> symmetry with `alpha_check` directly above it, `sin(yaw)` was almost certainly meant.
> It is preserved verbatim because its output went into `arc_angles_check.txt` alongside
> the published numbers. **Treat a disagreement in `eta_check` as uninformative.** The
> other three checks are sound; for October 2021 they agree with the fitted values to
> well within their errors:
>
> | | fitted | check |
> |---|---|---|
> | gamma | 1.857 +/- 0.025 | 1.867 +/- 0.033 |
> | alpha | 27.2 +/- 0.5 | 27.1 +/- 0.5 |
> | yaw | 0.850 +/- 0.013 | 0.854 +/- 0.023 |

## Uncertainty

Every quantity is a `ufloat` and errors propagate automatically. The inputs are assumed,
not measured (`constants.py`): 0.01 mm on stage position, 0.01 deg on goniometer angle,
0.1% on flux, 0.02% on wavelength.

Repeated measurements are combined by `stats.list_tot_err`, which adds the *scatter
between* the nominal values to the mean of the individual errors in quadrature. When the
repeats disagree by more than their individual errors, the scatter dominates — the
combination is deliberately conservative.

The throw enters twice, and this is worth understanding. `combine_xy_data` passes the
goniometer angles through `arcsin((235/L) sin(theta))` before projecting them. With `L =
235 +/- 2` mm the nominal value is unchanged, but the ratio is `1.000 +/- 0.009`, so the
step exists purely to **inflate the angular uncertainty** by the fractional throw error
before the arc fit sees it. That is intentional: the throw error is common to `x` and
`y` and this is where it is injected.

## What is assumed rather than measured

- **The throw.** `radius_curv_throw` was written to measure it from the angular width of
  the flat-topped direct beam. It never worked; the campaign files carry the comment
  "this doesn't actually work..." and every published run assumed 235 +/- 2 mm. The
  original silently fell back to `235 +/- 0.5` on failure — a *different* uncertainty,
  arrived at by failing. This package raises instead.
- **The groove spacing.** 158.0 +/- 0.5 nm, carried through the original reduction as a
  bare constant; the computation that had produced it was already commented out by the
  time of the last run, so its provenance does not survive in the code that used it. The
  reduction re-derives it from the measured dispersion as a cross-check and gets
  157.8 +/- 2.7 nm, which agrees.
