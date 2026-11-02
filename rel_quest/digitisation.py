"""The uncertainty of a value read off a paper's figure, from the figure's calibration.

A calibration (outputs/verification/ai_check/calibrations/<key>.json; format in
outputs/verification/ai_check/CALIBRATION.md) maps a figure's pixels to data:
two ticks per axis, optionally one set per panel. A value read off a figure is
good to about one marker radius, `marker_px` pixels at the calibration's dpi.

- On a log D axis (or a log D value plotted linearly) that radius is a fixed
  number of dex: sigma(log D) = marker_px * dex-per-pixel.
- On a linear D axis it is a fixed delta-D, so sigma(log D) = dD / (D ln 10)
  grows as D shrinks. Once dD reaches D the figure cannot tell D from zero and
  says only D < D + dD: the value becomes an upper bound (censor '<').

These numbers fill a cell's `sigma_log_D` with `sigma_source = 'measured'`, or
its `censor`/bound with `value_origin = 'censored'` (see the record schema).
"""

from __future__ import annotations

import math
import re
from typing import Dict, Optional

# Marker radius assumed when a calibration does not record one: about 0.9 mm,
# a typical journal symbol, i.e. 7 px at 200 dpi (measured on s41365-017-0229-4
# Fig. 3, whose markers are 14 px across).
DEFAULT_MARKER_PX = 7.0
DEFAULT_DPI = 200


LOCATOR = re.compile(r"\b(Fig(?:ure)?\.?\s*S?\d+[a-z]?|Table\s*S?I?\s*\d+)", re.I)


def figure_key(doi: Optional[str], locator: Optional[str]) -> str:
    """The key of a figure or table: calibrations/<key>.json, review/figures/<key>.png."""
    paper = re.sub(r"[^A-Za-z0-9._-]+", "_", doi) if doi else "no_doi"
    return f"{paper}__{re.sub(r'[^A-Za-z0-9]+', '_', locator or 'unknown').strip('_')}"


def cited_locator(locator, ai_evidence) -> Optional[str]:
    """The figure or table an entry cites: its locator, else the first one the AI's evidence names."""
    if locator is not None and str(locator).strip() and str(locator) != "nan":
        return str(locator).strip()
    m = LOCATOR.search(str(ai_evidence or ""))
    if not m:
        return None
    return re.sub(r"^Fig(?:ure)?\.?\s*", "Figure ", re.sub(r"\s+", " ", m.group(1)), flags=re.I)


def panel_for(cal: Dict, element: str, entry_id: str) -> Optional[Dict]:
    """The calibrated panel that plots this element of this entry (the whole figure when it has no panels)."""
    for panel in cal.get("panels") or [cal]:
        if ("elements" not in panel or element in panel["elements"]) and \
                ("entries" not in panel or entry_id in panel["entries"]):
            return panel
    return None


def marker_px(cal: Dict, panel: Dict) -> float:
    """Marker radius in pixels: the panel's, else the figure's, else the default scaled to the dpi."""
    if panel.get("marker_px") is not None:
        return float(panel["marker_px"])
    if cal.get("marker_px") is not None:
        return float(cal["marker_px"])
    return DEFAULT_MARKER_PX * cal.get("dpi", DEFAULT_DPI) / DEFAULT_DPI


def _units_per_px(axis: Dict) -> float:
    """Data units per pixel along an axis: dex for a log axis, the plotted quantity for a linear one."""
    (p1, v1), (p2, v2) = axis["ticks"]
    if axis.get("scale") == "log":
        v1, v2 = math.log10(v1), math.log10(v2)
    return abs(v2 - v1) / abs(p2 - p1)


def digitisation_sigma(y_axis: Dict, log_d: float, radius_px: float,
                       reading_px: Optional[float] = None, sigma_floor: float = 0.0) -> Dict:
    """{'sigma_log_D', 'censor', 'log_bound'} for a log D value read off this y axis.

    A readable value has a sigma and no censor; one below the axis resolution
    (linear D axis, dD >= D, dD = one marker radius) has censor '<', the bound
    log10(D + dD), and no sigma. `reading_px` (default: the marker radius) is how
    precisely a marker's position is known — finer on a vector figure read from
    the PDF's own coordinates — and `sigma_floor` the smallest sigma reported.
    """
    units = _units_per_px(y_axis)
    bound_step, read_step = radius_px * units, (radius_px if reading_px is None else reading_px) * units
    if y_axis["quantity"] == "logD" or y_axis.get("scale") == "log":
        return {"sigma_log_D": max(read_step, sigma_floor), "censor": None, "log_bound": None}
    d = 10 ** log_d
    if bound_step >= d:
        return {"sigma_log_D": None, "censor": "<", "log_bound": math.log10(d + bound_step)}
    return {"sigma_log_D": max(read_step / (d * math.log(10)), sigma_floor), "censor": None, "log_bound": None}


def figure_uncertainty(cal: Dict, element: str, entry_id: str, log_d: float) -> Optional[Dict]:
    """digitisation_sigma for one value of one entry on a calibrated figure; None when no panel plots it."""
    if cal.get("overlay") is False:
        return None
    panel = panel_for(cal, element, entry_id)
    if panel is None or "y" not in panel:
        return None
    reading = panel.get("reading_px", cal.get("reading_px"))
    floor = panel.get("sigma_floor", cal.get("sigma_floor", 0.0))
    return digitisation_sigma(panel["y"], log_d, marker_px(cal, panel), reading, floor)
