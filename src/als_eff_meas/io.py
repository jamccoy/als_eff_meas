"""Reading beamline scans and CXRO optical-constant tables.

Every function here takes an explicit :class:`~pathlib.Path`. The original code
did ``os.chdir(data_dir) ... os.chdir(home_dir)`` around each read, with
``home_dir`` captured at import time, so results landed wherever the process
happened to be when the module was first imported. Nothing in this package
depends on the current working directory.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from uncertainties import unumpy as unp

from .constants import FLUX_REL_ERR, POSITION_ERR_DEG, POSITION_ERR_MM

__all__ = [
    "cxro_to_n_k",
    "dark_val",
    "gather_file",
    "gather_files",
    "read_cxro_table",
    "read_scan",
]


def gather_file(file_dig, data_dir: Path, extension: str = ".dat") -> Path:
    """Return the first scan file whose *name contains* ``str(file_dig)``.

    .. warning::

       This is substring matching on the whole filename, not an index. Asking
       for ``30`` matches ``..._000030.dat`` but would equally match
       ``scan_300.dat`` or a file whose *date* contains "30". The original code
       carried the comment "finicky" about exactly this. The behaviour is
       preserved deliberately -- changing it changes which file each
       measurement reads -- so campaign digits must stay unambiguous.

    Files are sorted before matching, so the choice among several matches is at
    least deterministic.
    """
    data_dir = Path(data_dir)
    limited = sorted(f.name for f in data_dir.iterdir() if extension in f.name)
    matches = [f for f in limited if str(file_dig) in f]
    if not matches:
        raise FileNotFoundError(
            f"no file matching {file_dig!r} with extension {extension!r} in {data_dir}"
        )
    return data_dir / matches[0]


def gather_files(file_digs, data_dir: Path, extension: str = ".dat") -> list[Path]:
    """Apply :func:`gather_file` to each digit in turn, preserving order.

    Note the original did *not* sort in this variant while ``gather_file`` did.
    Sorting here too makes the two consistent; it changes the result only when a
    digit is ambiguous, which is already a configuration error.
    """
    return [gather_file(d, data_dir, extension) for d in file_digs]


def read_scan(path: Path, data_type: str = "linear"):
    """Read one photodiode scan and attach the assumed measurement errors.

    The file is whitespace-delimited with a single header line (a timestamp)::

        Sunday, November 05, 2017 11:00:39 AM
        0.000000   1.501584E-2   4.122031E-1   499.968

    Column 0 is stage position (mm) or goniometer angle (deg); column 1 is
    flux. Columns 2 and 3 (I0 and the actual monochromator energy) are read by
    neither this code nor the original.

    Returns ``(upos, uflux)`` as ``uncertainties`` arrays. Unlike the original
    ``extract_flux_profile``, this writes no diagnostic PNG -- see
    :mod:`als_eff_meas.plots`.
    """
    if data_type not in ("linear", "angle"):
        raise ValueError(f"data_type must be 'linear' or 'angle', got {data_type!r}")
    data = np.loadtxt(path, skiprows=1)
    pos, flux = data[:, 0], data[:, 1]
    pos_err = (POSITION_ERR_MM if data_type == "linear" else POSITION_ERR_DEG) * np.ones(len(pos))
    flux_err = FLUX_REL_ERR * flux
    return unp.uarray(pos, pos_err), unp.uarray(flux, abs(flux_err))


def dark_val(start_dig, darks_dir: Path):
    """Median flux of a shutter-closed scan, used as the dark level.

    Median rather than mean, so a stray spike in a dark scan does not shift the
    subtraction.

    Implemented by ordering on nominal values and then combining the middle
    element(s) as ufloats. ``np.median`` on an object array does exactly this,
    but reaches the ordering through the ``AffineScalarFunc`` comparison
    operators that ``uncertainties`` has deprecated.
    """
    uflux = read_scan(gather_file(start_dig, darks_dir), "linear")[1]
    order = np.argsort(unp.nominal_values(uflux))
    n = len(order)
    if n % 2:
        return uflux[order[n // 2]]
    return 0.5 * (uflux[order[n // 2 - 1]] + uflux[order[n // 2]])


def read_cxro_table(path: Path):
    """Read a CXRO delta/beta table: ``energy(eV), delta, beta``, 3 header lines."""
    return np.genfromtxt(fname=str(path), skip_header=3)


def cxro_to_n_k(path: Path, material: str, out_dir: Path | None = None) -> np.ndarray:
    """Convert a CXRO delta/beta table to ``energy, wavelength, n, k``.

    ``n = 1 - delta`` and ``k = beta``. Returns the array; writes
    ``<material>_optical_constants_n_k.txt`` only if ``out_dir`` is given (the
    original always wrote, into whatever directory was current).
    """
    from .constants import HC_EV_NM

    raw = read_cxro_table(path)
    energy, delta, beta = raw[:, 0], raw[:, 1], raw[:, 2]
    table = np.array([energy, HC_EV_NM / energy, 1 - delta, beta]).transpose()
    if out_dir is not None:
        out = Path(out_dir) / f"{material}_optical_constants_n_k.txt"
        np.savetxt(
            out,
            table,
            header="complex index of refraction from CXRO\nenergy(eV)\twavelength(nm)\tn\tk",
        )
    return table
