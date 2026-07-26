"""Reduction of ALS beamline scans to grating diffraction efficiency.

Two stages, in order::

    campaign = load_campaign("ogre_laminar_I_oct21")
    geometry = fit_geometry(campaign)            # where is the arc?
    table    = measure_efficiency(campaign, geometry)   # what is on it?

The raw data is private and lives outside this repository; point at it with
``$ALS_DATA_ROOT``.
"""

from .config import Campaign, data_root, list_campaigns, load_campaign
from .efficiency import EfficiencyTable
from .geometry import ArcFit, ConeAngles
from .pipeline import Geometry, fit_geometry, measure_efficiency

__version__ = "0.1.0"

__all__ = [
    "ArcFit",
    "Campaign",
    "ConeAngles",
    "EfficiencyTable",
    "Geometry",
    "data_root",
    "fit_geometry",
    "list_campaigns",
    "load_campaign",
    "measure_efficiency",
]
