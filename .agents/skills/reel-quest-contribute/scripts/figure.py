#!/usr/bin/env python3
"""Read distribution ratios off a paper's figure, reproducibly.

    python figure.py text <pdf> [<page> ...]
        The PDF's text layer, page by page ("===== page N ====="), to find the
        experimental section, compound names, tables and captions. It can be
        garbled (tables, scans): read values from the rendered page instead.
    python figure.py render <pdf> <page> [<page> ...] [--dpi 200] [--out DIR]
        PNG of each page (1-based), printed one path per line.
        Default DIR: outputs/contributions/pages/ in the repository.
    python figure.py crop <png> <x0> <y0> <x1> <y1> [--zoom 3] [--grid 50] [--out PATH]
        An enlarged region with a labelled pixel grid in the page's own
        coordinates, to locate axis ticks and marker centres.
    python figure.py read <calibration.json> [--overlay PATH]
        The cells each series' points give, as JSON on stdout (the
        contribution's cells grammar), with sigma_log_D from the figure's
        resolution (rel_quest.digitisation). Warnings go to stderr. --overlay
        draws the readings on the page for a visual check.

A calibration (references/figures.md) maps pixels of one render to data:

    {"pdf": "papers/x.pdf", "page": 3, "dpi": 200, "locator": "Figure 2",
     "crop": [x0, y0, x1, y1],
     "y": {"quantity": "D" | "logD", "scale": "log" | "linear", "ticks": [[py, value], [py, value]]},
     "x": {"quantity": "element" | "<a condition field>", "scale": "linear" | "log",
           "ticks": [[px, value], [px, value]]},
     "marker_px": 7,
     "series": {"A0001": {"La": [px, py], "Ce": [px, py]}}}

A figure with several panels gives "panels": [{"y": ..., "x": ..., "series": ...}]
instead of top-level y, x and series. On an element axis, x values are
atomic numbers (La = 57 ... Lu = 71, Y = 39, Sc = 21).
"""

import argparse
import json
import math
import sys
from pathlib import Path


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "rel_quest" / "__init__.py").is_file():
            return parent
    sys.exit("figure.py must run inside a REEL-QUEST checkout")


ROOT = repository_root()
sys.path.insert(0, str(ROOT))

from rel_quest.digitisation import digitisation_sigma, marker_px  # noqa: E402

Z = {"Sc": 21, "Y": 39, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
     "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71}
BY_Z = {z: e for e, z in Z.items()}


def resolve(path) -> Path:
    path = Path(path)
    if path.is_absolute() or path.exists():
        return path
    return ROOT / path


def render_page(pdf, page: int, dpi: int):
    """The page (1-based) as a PIL image at `dpi`."""
    import pypdfium2 as pdfium
    doc = pdfium.PdfDocument(str(resolve(pdf)))
    if not 1 <= page <= len(doc):
        raise SystemExit(f"page {page} is outside 1-{len(doc)} of {pdf}")
    return doc[page - 1].render(scale=dpi / 72).to_pil().convert("RGB")


def cmd_text(args) -> int:
    from pypdf import PdfReader
    reader = PdfReader(str(resolve(args.pdf)))
    for page in args.pages or range(1, len(reader.pages) + 1):
        if not 1 <= page <= len(reader.pages):
            raise SystemExit(f"page {page} is outside 1-{len(reader.pages)} of {args.pdf}")
        print(f"===== page {page} =====")
        print(reader.pages[page - 1].extract_text() or "")
    return 0


def cmd_render(args) -> int:
    out = Path(args.out) if args.out else ROOT / "outputs" / "contributions" / "pages"
    out.mkdir(parents=True, exist_ok=True)
    for page in args.pages:
        image = render_page(args.pdf, page, args.dpi)
        path = out / f"{Path(args.pdf).stem}_p{page}.png"
        image.save(path)
        print(path)
    return 0


def cmd_crop(args) -> int:
    from PIL import Image, ImageDraw
    image = Image.open(args.png).convert("RGB")
    x0, y0, x1, y1 = args.box
    zoom = args.zoom
    crop = image.crop((x0, y0, x1, y1)).resize(((x1 - x0) * zoom, (y1 - y0) * zoom), Image.NEAREST)
    if args.grid:
        draw = ImageDraw.Draw(crop)
        for x in range((x0 // args.grid + 1) * args.grid, x1, args.grid):
            draw.line([((x - x0) * zoom, 0), ((x - x0) * zoom, crop.height)], fill=(0, 160, 255), width=1)
            draw.text(((x - x0) * zoom + 2, 2), str(x), fill=(0, 100, 220))
        for y in range((y0 // args.grid + 1) * args.grid, y1, args.grid):
            draw.line([(0, (y - y0) * zoom), (crop.width, (y - y0) * zoom)], fill=(0, 160, 255), width=1)
            draw.text((2, (y - y0) * zoom + 2), str(y), fill=(0, 100, 220))
    out = Path(args.out) if args.out else Path(args.png).with_name(
        f"{Path(args.png).stem}_crop_{x0}_{y0}_{x1}_{y1}.png")
    crop.save(out)
    print(out)
    return 0


def axis_value(axis: dict, pixel: float) -> float:
    """The data value at `pixel` along a two-tick axis (log10-linear on a log scale)."""
    (p1, v1), (p2, v2) = axis["ticks"]
    if p1 == p2:
        raise SystemExit("an axis's two ticks must be at different pixels")
    log = axis.get("scale") == "log"
    if log:
        v1, v2 = math.log10(v1), math.log10(v2)
    value = v1 + (pixel - p1) * (v2 - v1) / (p2 - p1)
    return 10 ** value if log else value


def read_point(y_axis: dict, py: float, radius: float, reading_px=None, sigma_floor=0.0) -> dict:
    """One cell: log_D with its sigma, or an upper bound when the figure cannot resolve it."""
    value = axis_value(y_axis, py)
    if y_axis["quantity"] == "logD":
        log_d = value
    elif value > 0:
        log_d = math.log10(value)
    else:
        log_d = None
    if log_d is None:          # at or below a linear D axis's zero: only a bound is known
        step = radius * abs(y_axis["ticks"][1][1] - y_axis["ticks"][0][1]) / \
            abs(y_axis["ticks"][1][0] - y_axis["ticks"][0][0])
        return {"censor": "<", "log_bound": round(math.log10(max(value, 0) + step), 4)}
    sigma = digitisation_sigma(y_axis, log_d, radius, reading_px, sigma_floor)
    if sigma["censor"]:
        return {"censor": sigma["censor"], "log_bound": round(sigma["log_bound"], 4)}
    return {"log_D": round(log_d, 4), "sigma_log_D": round(sigma["sigma_log_D"], 4), "value_origin": "digitised"}


def panels(cal: dict):
    return cal.get("panels") or [cal]


def cmd_read(args) -> int:
    cal = json.loads(Path(args.calibration).read_text(encoding="utf-8"))
    out, notes = {}, []
    for index, panel in enumerate(panels(cal)):
        if "y" not in panel:
            raise SystemExit(f"panel {index}: no y axis")
        radius = marker_px(cal, panel)
        reading = panel.get("reading_px", cal.get("reading_px"))
        floor = panel.get("sigma_floor", cal.get("sigma_floor", 0.0))
        x_axis = panel.get("x")
        for series, points in panel.get("series", {}).items():
            cells = out.setdefault(series, {})
            for element, (px, py) in points.items():
                if element not in Z:
                    raise SystemExit(f"{series}: {element!r} is not an element symbol")
                if element in cells:
                    raise SystemExit(f"{series}: {element} is read twice")
                cells[element] = read_point(panel["y"], py, radius, reading, floor)
                if x_axis:
                    x = axis_value(x_axis, px)
                    if x_axis.get("quantity") == "element":
                        nearest = BY_Z.get(round(x))
                        if nearest != element:
                            notes.append(f"{series}: the point filed as {element} sits at Z = {x:.2f} "
                                         f"({nearest or 'no element'}) on the element axis")
                    else:
                        notes.append(f"{series} {element}: x = {x:.4g} ({x_axis.get('quantity', 'x')})")
    for note in notes:
        print(note, file=sys.stderr)
    if args.overlay:
        overlay(cal, out, Path(args.overlay))
        print(f"overlay: {args.overlay}", file=sys.stderr)
    print(json.dumps(out, indent=1))
    return 0


def overlay(cal: dict, cells: dict, path: Path) -> None:
    """The page with each reading marked (red cross, element label) and the ticks (blue)."""
    from PIL import ImageDraw
    image = render_page(cal["pdf"], cal["page"], cal.get("dpi", 200))
    draw = ImageDraw.Draw(image)
    for panel in panels(cal):
        for axis, vertical in (("y", True), ("x", False)):
            for pixel, value in (panel.get(axis) or {}).get("ticks", []):
                if vertical:
                    draw.line([(0, pixel), (image.width, pixel)], fill=(0, 120, 255), width=1)
                else:
                    draw.line([(pixel, 0), (pixel, image.height)], fill=(0, 120, 255), width=1)
        for series, points in panel.get("series", {}).items():
            for element, (px, py) in points.items():
                r = 6
                draw.line([(px - r, py), (px + r, py)], fill=(220, 0, 0), width=2)
                draw.line([(px, py - r), (px, py + r)], fill=(220, 0, 0), width=2)
                draw.text((px + r + 1, py - r - 10), f"{series} {element}", fill=(220, 0, 0))
    if cal.get("crop"):
        image = image.crop(tuple(cal["crop"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("text", help="the PDF's text, page by page")
    p.add_argument("pdf")
    p.add_argument("pages", type=int, nargs="*")
    p.set_defaults(func=cmd_text)
    p = sub.add_parser("render", help="render PDF pages to PNG")
    p.add_argument("pdf")
    p.add_argument("pages", type=int, nargs="+")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--out")
    p.set_defaults(func=cmd_render)
    p = sub.add_parser("crop", help="an enlarged region with a pixel grid")
    p.add_argument("png")
    p.add_argument("box", type=int, nargs=4, metavar="X0 Y0 X1 Y1")
    p.add_argument("--zoom", type=int, default=3)
    p.add_argument("--grid", type=int, default=50, help="grid spacing in page pixels (0: none)")
    p.add_argument("--out")
    p.set_defaults(func=cmd_crop)
    p = sub.add_parser("read", help="cells from a calibration")
    p.add_argument("calibration")
    p.add_argument("--overlay")
    p.set_defaults(func=cmd_read)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
