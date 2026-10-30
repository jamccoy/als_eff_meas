"""Campaigns as data.

The original stored each experiment's parameters in a Python module that was
consumed with ``from ogre_laminar_I_Oct21 import *`` — and importing it ran
``dark_val()`` and ``gather_file()``, so *reading the configuration touched the
disk*. Selecting a dataset meant commenting out five import lines and
uncommenting a sixth.

Here a campaign is a TOML file. Loading one reads no beamline data; the paths
are resolved but not opened. Adding a dataset means adding a file.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from endstation.paths import data_root as _data_root
from endstation.paths import require_dir
from uncertainties import ufloat

__all__ = ["Campaign", "EnergyGroup", "YawScan", "data_root", "list_campaigns", "load_campaign"]

_CAMPAIGN_DIR = Path(__file__).parent / "campaigns"


#: Where the beamline tree lives when ``$ALS_DATA_ROOT`` says nothing.
DEFAULT_DATA_ROOT = "~/Documents/ALS"


def data_root() -> Path:
    """Root of the private beamline data tree.

    From ``$ALS_DATA_ROOT``, defaulting to ``~/Documents/ALS``. The data is not
    in this repository and never will be.

    Deliberately does not check that the directory exists -- that is
    :meth:`Campaign.check_data_present`'s job. Separating them is what lets a
    caller ask "where would it be?" in order to name it in a skip message,
    without first having to handle an exception.
    """
    return _data_root("ALS_DATA_ROOT", DEFAULT_DATA_ROOT)


def _uf(pair) -> object:
    """``[nominal, std_dev]`` -> ufloat. A bare number means an exact value."""
    if isinstance(pair, (int, float)):
        return ufloat(float(pair), 0.0)
    nominal, std = pair
    return ufloat(float(nominal), float(std))


def _energies(spec) -> list:
    """Accept an explicit list or ``{arange = [start, stop, step]}``."""
    if isinstance(spec, dict):
        start, stop, step = spec["arange"]
        return np.arange(start, stop, step).tolist()
    return list(spec)


@dataclass(frozen=True)
class YawScan:
    """One energy's worth of arc-mapping scans.

    ``x_scan`` maps the orders along the arc; one ``y_scan`` per order maps it
    in goniometer angle. ``meas_x`` are rough starting positions — they only
    have to be within ``peak_step`` samples of the true peak.
    """

    energy_eV: float
    x_scan: object
    y_scans: list
    meas_x: list
    orders: list


@dataclass(frozen=True)
class EnergyGroup:
    """A set of energies sharing one detector gain, and so one dark scan."""

    name: str
    directory: str
    energies: list
    dark_dir: str
    dark_dig: object
    pos_factor: float = 1.0
    neg_factor: float = 1.0


@dataclass(frozen=True)
class Campaign:
    """Everything needed to reduce one sample from one beamtime."""

    name: str
    description: str
    data_dir: str
    figure_dir: str
    throw_mm: object
    groove_spacing_nm: object
    peak_step: int
    arc_guess: tuple
    yaw_scans: list[YawScan] = field(default_factory=list)
    energy_groups: list[EnergyGroup] = field(default_factory=list)
    direct_beam_dir: str = ""
    direct_beam_files: list = field(default_factory=list)
    detector_width_mm: object = None
    offset: float = 0.0
    #: Scales the fitted arc radius between the arc fit and the order
    #: prediction, where ``R = R_correct*R`` sat in the original drivers. It
    #: therefore changes how many orders are counted as propagating without
    #: touching the recovered angles. 1.0 leaves the fit alone.
    radius_correction: float = 1.0
    #: Name of the per-group noise-floor file, relative to :attr:`data`. The
    #: floors were read off the raw scans, so they are beamline-derived and
    #: live in the private tree rather than in this repository; a campaign
    #: that names one cannot be reduced without it. Empty means no floors,
    #: which is a true zero and not a missing file.
    noise_file: str = ""
    exclusions: list = field(default_factory=list)
    plot_orders: list = field(default_factory=list)

    @property
    def data(self) -> Path:
        return data_root() / self.data_dir

    @property
    def yaw_scans_dir(self) -> Path:
        return self.data / "yaw_scans"

    @property
    def results_dir(self) -> Path:
        return self.data / "results"

    @property
    def noise_path(self) -> Path:
        return self.data / self.noise_file

    def group_dir(self, group: EnergyGroup) -> Path:
        return self.data / group.directory

    def dark_dir(self, group: EnergyGroup) -> Path:
        return self.data / group.dark_dir

    def check_data_present(self) -> None:
        """Raise with a useful message if the private data tree is missing.

        Naming the variable is the whole point: ``FileNotFoundError: [Errno 2]``
        on a path the user never typed is unactionable, and every project here
        had written its own version of this message.
        """
        require_dir(self.data, env_var="ALS_DATA_ROOT", hint="docs/DATA.md")


def load_campaign(name_or_path) -> Campaign:
    """Load a campaign by name from ``campaigns/``, or from an explicit path."""
    path = Path(name_or_path)
    if not path.suffix:
        path = _CAMPAIGN_DIR / f"{name_or_path}.toml"
    if not path.is_file():
        raise FileNotFoundError(
            f"no campaign {name_or_path!r}; available: {', '.join(list_campaigns())}"
        )
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    inst = raw["instrument"]
    arc = raw["arc_fit"]
    return Campaign(
        name=raw["name"],
        description=raw.get("description", ""),
        data_dir=raw["paths"]["data"],
        figure_dir=raw["paths"].get("figures", ""),
        throw_mm=_uf(inst["throw_mm"]),
        groove_spacing_nm=_uf(inst["groove_spacing_nm"]),
        peak_step=int(inst["peak_step"]),
        detector_width_mm=_uf(inst["detector_width_mm"]) if "detector_width_mm" in inst else None,
        arc_guess=(arc["radius_guess"], arc["x_center_guess"], arc["y_center_guess"]),
        yaw_scans=[
            YawScan(
                energy_eV=s["energy_eV"],
                x_scan=s["x_scan"],
                y_scans=list(s["y_scans"]),
                meas_x=list(s["meas_x"]),
                orders=list(s["orders"]),
            )
            for s in raw.get("yaw_scan", [])
        ],
        energy_groups=[
            EnergyGroup(
                name=g.get("name", g["dir"]),
                directory=g["dir"],
                energies=_energies(g["energies"]),
                dark_dir=g["dark_dir"],
                dark_dig=g.get("dark_dig", 1),
                pos_factor=float(g.get("pos_factor", 1.0)),
                neg_factor=float(g.get("neg_factor", 1.0)),
            )
            for g in raw.get("energy_group", [])
        ],
        direct_beam_dir=raw.get("direct_beam", {}).get("dir", ""),
        direct_beam_files=list(raw.get("direct_beam", {}).get("files", [])),
        offset=float(raw.get("efficiency", {}).get("offset", 0.0)),
        radius_correction=float(raw.get("efficiency", {}).get("radius_correction", 1.0)),
        noise_file=str(raw.get("efficiency", {}).get("noise_file", "")),
        exclusions=[tuple(p) for p in raw.get("efficiency", {}).get("exclude", [])],
        plot_orders=list(raw.get("efficiency", {}).get("plot_orders", [])),
    )


def list_campaigns() -> list[str]:
    """Names of the campaigns shipped with the package."""
    return sorted(p.stem for p in _CAMPAIGN_DIR.glob("*.toml"))
