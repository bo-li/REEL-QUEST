# Reading values off a figure

A value read off a figure is only as good as the axis it was measured
against. Never estimate a point by eye from the whole page: calibrate the
axes in pixels, locate each marker's centre in pixels, and let
`scripts/figure.py read` convert. The calibration file is kept with the
contribution, so the reviewer can redraw your readings on the page and see
whether they sit on the markers.

`F` below is `python .agents/skills/reel-quest-contribute/scripts/figure.py`.

## Steps

1. **Render the page** at 200 dpi (300 for a small figure):
   `F render papers/<paper>.pdf <page> --dpi 200`. All pixels below are
   pixels of this render, (0, 0) at the top left, y growing downwards.
2. **Identify the series.** From the caption and legend: which curve or
   symbol is which extractant and condition, what the axes plot (D or log D,
   log or linear scale; x the element, [HNO3], [ligand], temperature ...).
   Curves from other papers (cited in the caption) are not this paper's.
3. **Calibrate both axes.** Crop around each axis
   (`F crop <png> x0 y0 x1 y1 --zoom 4 --grid 20`) and find the centre pixel
   of two labelled ticks per axis, far apart. On a log axis use decades
   (0.1 and 100), not minor ticks.
4. **Locate every marker.** Crop around each cluster with a fine grid and
   read the centre of each marker (the middle of its symbol, not an error
   bar's end). Measure one isolated marker's half-width: that is
   `marker_px`, how precisely a value can be read.
5. **Write the calibration**, `contributions/figures/<doi-slug>/<locator>.json`
   (locator with spaces as `_`, e.g. `Figure_2.json`):

   ```json
   {"pdf": "papers/<doi-slug>.pdf", "page": 4, "dpi": 200,
    "locator": "Figure 2", "crop": [130, 480, 800, 980],
    "y": {"quantity": "D", "scale": "log", "ticks": [[512.0, 1000], [940.0, 0.1]]},
    "x": {"quantity": "element", "ticks": [[205.0, 57], [770.0, 71]]},
    "marker_px": 7,
    "series": {"A0001": {"La": [205.4, 601.2], "Ce": [245.9, 640.0]}},
    "note": "open squares: 0.1 M L2 in n-dodecane from 3 M HNO3; circles are ref. 14's data"}
   ```

   - `y.quantity`: `D` when the axis is labelled D, `logD` when it is log D;
     `scale`: `log` for a logarithmic axis, else `linear`.
   - `x.quantity`: `element` (ticks in atomic numbers: La 57 ... Lu 71) or
     the condition field the axis plots (`Acid_Concentration_M`, ...).
   - `series`: one entry per condition id of the contribution, each element
     with its marker's `[x, y]`.
   - Several panels with their own axes: `"panels": [{"y": ..., "x": ...,
     "series": ...}, ...]` instead of top-level `y`, `x`, `series`.
   - A figure drawn as vectors can be read more precisely than its markers
     suggest; leave that to the reviewer and keep `marker_px`.
6. **Read it**:
   `F read contributions/figures/<doi-slug>/Figure_2.json --overlay outputs/contributions/<doi-slug>_Figure_2.png`.
   It prints the cells, in the contribution's grammar, with `sigma_log_D`
   from the axis and `marker_px` (a value below a linear axis's resolution
   comes out as a `<` bound), and warns on stderr when a point filed as one
   element sits at another element's position, or gives the x value each
   point sits at on a condition axis: check it against the condition.
7. **Look at the overlay.** Every red cross must sit on its marker's centre.
   Fix the calibration and read again until they do. Then copy the cells
   into the contribution and name the calibration file in `why`.

An axis broken in the middle, a figure with no tick labels, or markers
hidden behind each other cannot be read reliably: read what you can, give
the rest as `Unsure` in `note`, and say why.
