"""Physical constants and the instrument's assumed measurement errors.

These values are lifted from the original 2021 analysis scripts, where they
appeared as bare literals inside the functions that used them.
"""

#: Planck constant times c, in eV nm. Converts photon energy to wavelength.
HC_EV_NM = 1239.84193

#: Assumed stage position error, mm. ``extract_flux_profile`` used 0.01 for both
#: the linear and the angular branch -- the two branches were identical.
POSITION_ERR_MM = 0.01

#: Assumed goniometer position error, degrees.
POSITION_ERR_DEG = 0.01

#: Assumed flux error as a fraction of flux (0.1%).
FLUX_REL_ERR = 0.001

#: Assumed fractional error on the monochromator energy, applied to the
#: wavelength in the efficiency measurement as ``ufloat(1, 0.0002)``.
WAVELENGTH_REL_ERR = 0.0002

#: Throw used to normalise goniometer angles in ``combine_xy_data``. This is a
#: reference value baked into the original code as the literal ``235.``; it is
#: *not* the same quantity as the campaign's assumed throw, even when the two
#: happen to be numerically equal. See ``geometry.combine_xy_data``.
REFERENCE_THROW_MM = 235.0
