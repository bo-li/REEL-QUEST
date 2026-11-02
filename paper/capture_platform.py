"""Capture the inputs of the platform figures (Fig. 2, Fig. S3) from the running software.

Writes paper/figures/platform/:
  assistant.png        the dashboard's AI assistant answering QUESTION (cropped after its table)
  browser_table.png, browser_detail.png
                       the ligand browser: the table the sidebar search filters, and LIGAND's
                       detail view (structure, identifiers, trend fit and its statistics)
  sf_selectors.png, sf_plot.png
                       the separation-factor viewer (Fig. S3): its two element selectors set to
                       PAIR, and the log SF plot they drive
  mcp_exchange.json    one real MCP tool call (get_best_separators Eu/La) and its result
  capture.json         when, from which commits, which model answered, and the full answer

Nothing here is drawn or typed by hand: the dashboard is started from this checkout,
the assistant answers with the model configured in .streamlit/secrets.toml (free mode,
OPENROUTER_FREE_API_KEY), and the MCP server is spawned over stdio from --mcp-repo.
paper/figures.py composes the figure from these files.

--app-root runs the dashboard from another tree: the public export, whose ids and
controls carry no source labels, is what a reader will run. --only browser (or sf_viewer)
recaptures that part alone and keeps the other captures and their record in capture.json:
    python scripts/export_public_repo.py --out <dir>
    python paper/capture_platform.py --only browser --app-root <dir> --app-python <python with streamlit>

Run with an interpreter that has Playwright (pip install playwright; playwright install
chromium), e.g.
    python paper/capture_platform.py --app-python python --mcp-python /path/to/python-with-mcp
--app-python runs the dashboard (default: python on PATH); --mcp-python needs the mcp SDK
(the `mcp` extra). The assistant's answer is model-generated: check its numbers against
agent_tools.get_best_separators after every capture.
"""

import argparse
import json
import os
import re
import socket
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paper" / "figures" / "platform"
DATABASE = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
QUESTION = "What are the best ligands for separating Eu from La?"
# One of the dashboard's checked free presets, pinned so the answer does not depend on which
# model the "openrouter/free" router happens to pick (its answers vary in quality).
# Free endpoints are retired without notice (nemotron-3-super's was, 2026-10-20).
FREE_MODEL = "nvidia/nemotron-3-ultra-550b-a55b:free"
MCP_CALL = ("get_best_separators", {"element1": "Eu", "element2": "La", "limit": 3})
LIGAND = "TODGA"  # the browser's example: the most-measured diglycolamide
# The viewer's example pair: the industrially hard adjacent pair, measured together under the
# most conditions after Eu/La, which Fig. 5c already shows.
PAIR = ("Nd", "Pr")
ELEMENT_NAMES = {"Nd": "Neodymium", "Pr": "Praseodymium", "Eu": "Europium", "La": "Lanthanum"}
HIDE_CHROME = "header[data-testid='stHeader'], [data-testid='stToolbar'] {visibility: hidden !important;}"
SCALE = 3  # device pixels per CSS pixel
# Viewport widths (CSS px). At 1550 px with the sidebar collapsed, the browser table shows
# every column and fills the figure's full width; at 1100 px the chat keeps each table row
# on one line, which keeps the half-width panel short.
DASHBOARD_WIDTH = 1550
ASSISTANT_WIDTH = 1100
GRID_ROW = 35  # px per row of st.dataframe (header included)

MCP_CLIENT = """
import asyncio, json, sys
from mcp import Client
from mcp.client.stdio import StdioServerParameters

async def main(name, arguments, env):
    params = StdioServerParameters(command=sys.executable, args=["-m", "rel_quest.mcp_server"], env=env)
    async with Client(params, mode="legacy", read_timeout_seconds=120) as client:
        tools = await client.list_tools()
        result = await client.call_tool(name, arguments)
        print(json.dumps({"server": client.server_info.name, "tools": [t.name for t in tools.tools],
                          "call": name, "arguments": arguments, "is_error": bool(result.is_error),
                          "result": json.loads(result.content[0].text)}))

name, arguments, env = json.loads(sys.argv[1])
asyncio.run(main(name, arguments, env))
"""


def git_commit(repo):
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_dashboard(python, root=ROOT):
    port = free_port()
    proc = subprocess.Popen([python, "-m", "streamlit", "run", "app.py", "--server.port", str(port),
                             "--server.headless", "true", "--browser.gatherUsageStats", "false"],
                            cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    url = f"http://localhost:{port}"
    for _ in range(120):
        try:
            urllib.request.urlopen(url, timeout=2)
            return proc, url
        except OSError:
            time.sleep(1)
    proc.terminate()
    raise RuntimeError("the dashboard did not start")


def clip(page, top, bottom, pad=8):
    """Viewport clip spanning the main column between two y positions (CSS px)."""
    main = page.locator("[data-testid='stMainBlockContainer']").bounding_box()
    return {"x": main["x"], "y": max(top - pad, 0), "width": main["width"], "height": bottom - top + pad}


def open_dashboard(browser, url, width, height=2400):
    page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=SCALE)
    page.goto(url, wait_until="networkidle")
    page.wait_for_selector("text=Database Overview", timeout=180_000)
    page.add_style_tag(content=HIDE_CHROME)
    page.wait_for_timeout(3000)
    return page


def collapse_sidebar(page):
    page.locator("[data-testid='stSidebarCollapseButton'] button").first.click(force=True)
    page.wait_for_timeout(2000)


def capture_assistant(browser, url, free_model=FREE_MODEL):
    page = open_dashboard(browser, url, ASSISTANT_WIDTH)
    page.get_by_text("AI Assistant").first.click()
    page.wait_for_selector("text=REEL-QUEST AI Assistant", timeout=60_000)
    page.wait_for_timeout(2000)
    collapse_sidebar(page)  # at this width it overlays the page
    if free_model:  # choose the model in the dashboard's own Settings dialog
        page.get_by_role("button", name="Settings").first.click()
        dialog = page.get_by_role("dialog")
        dialog.locator("[data-testid='stSelectbox']").first.click()
        page.get_by_role("option", name=free_model, exact=True).click()
        dialog.get_by_role("button", name="Save").click()
        page.wait_for_selector(f"text={free_model}", timeout=60_000)
        page.wait_for_timeout(2000)
    model = page.get_by_text("Free models ·").first.inner_text() if page.get_by_text("Free models ·").count() else None
    if model is None:
        raise RuntimeError("free mode is not offered: set OPENROUTER_FREE_API_KEY in .streamlit/secrets.toml")
    box = page.get_by_placeholder("Ask about ligands")
    box.fill(QUESTION)
    box.press("Enter")
    page.wait_for_timeout(5000)
    page.wait_for_selector("text=Thinking", state="detached", timeout=300_000)
    page.wait_for_timeout(4000)
    messages = page.locator("[data-testid='stChatMessage']")
    first, last = messages.first, messages.last
    first.evaluate("e => e.scrollIntoView({block: 'start'})")
    page.wait_for_timeout(1500)
    top, end = first.bounding_box(), last.bounding_box()
    # Crop after the answer's (last) table when it has one, to keep the panel short;
    # capture.json keeps the whole answer.
    tables = last.locator("table")
    has_table = tables.count() > 0
    if has_table:
        end = tables.last.bounding_box()
    page.screenshot(path=OUT / "assistant.png", clip=clip(page, top["y"], end["y"] + end["height"] + 6))
    answer = last.inner_text()
    page.close()
    # The agent reports a failed call (retired model, rate limit, bad key) as a ⚠️ message;
    # never let one become the figure.
    if "⚠️" in answer or not has_table:
        raise RuntimeError(f"the assistant did not answer with a table:\n{answer[:500]}")
    return {"model": model, "question": QUESTION, "answer": answer}


def search(page, name):
    box = page.get_by_label("🔍 Search ligand name")
    box.fill(name)
    box.press("Enter")
    page.wait_for_timeout(6000)


def browser_table(page):
    head = page.get_by_text("Ligand Browser (").first
    head.evaluate("e => e.scrollIntoView({block: 'start'})")
    page.wait_for_timeout(3000)
    top = head.bounding_box()["y"]  # the first table below the heading (the detail view adds others)
    return head, min((b for b in (e.bounding_box() for e in page.locator("[data-testid='stDataFrame']").all())
                      if b and b["y"] > top), key=lambda b: b["y"])


def capture_browser(browser, url, ligand=LIGAND):
    """The ligand browser: the table the sidebar search filters, with the ligand's row
    selected, and the detail view that selection opens."""
    page = open_dashboard(browser, url, DASHBOARD_WIDTH)
    search(page, ligand)
    # the search box keeps the sidebar busy: click its collapse button directly
    page.locator("[data-testid='stSidebarCollapseButton'] button").first.evaluate("e => e.click()")
    page.wait_for_timeout(2000)
    head, _ = browser_table(page)
    n = int(re.search(r"\((\d+) ligands", head.inner_text()).group(1))
    title = page.get_by_text(f"Selectivity Trend: {ligand} (")
    for row in range(n):  # the table keeps the database's order: try each row, never a fixed index
        _, grid = browser_table(page)
        page.mouse.click(grid["x"] + 18, grid["y"] + GRID_ROW * (row + 1) + GRID_ROW / 2)
        page.wait_for_timeout(6000)
        if title.count():
            break
    else:
        raise RuntimeError(f"no row of the browser table opens {ligand}")
    head, grid = browser_table(page)
    page.screenshot(path=OUT / "browser_table.png",
                    clip=clip(page, head.bounding_box()["y"], grid["y"] + grid["height"]))
    top = page.get_by_text("Structure", exact=True).first
    top.evaluate("e => e.scrollIntoView({block: 'start'})")
    page.wait_for_timeout(5000)
    # The two columns end at different heights: the left one above "Chemical Identifiers"
    # (structure, SMILES, family, synonyms), the right one after the trend statistics and the
    # data-quality line. No straight cut does both, so the left column below its end is blanked.
    from PIL import Image
    box_top = lambda e: e.bounding_box()["y"]  # noqa: E731
    main = page.locator("[data-testid='stMainBlockContainer']")  # the sidebar has a "Selectivity Trend" too
    quality = main.get_by_text("Data Quality:").first.bounding_box()
    identifiers = main.get_by_text("Chemical Identifiers", exact=True).first.bounding_box()
    right = main.get_by_text("Selectivity Trend", exact=True).first.bounding_box()["x"]
    if not 0 < identifiers["y"] - box_top(top) < quality["y"] - box_top(top) or right <= 0:
        raise RuntimeError("the detail view's layout changed: check the browser_detail crop")
    box = clip(page, top.bounding_box()["y"], quality["y"] + quality["height"] + 6)
    path = OUT / "browser_detail.png"
    page.screenshot(path=path, clip=box)
    shot = Image.open(path).convert("RGB")
    cut = [round(SCALE * v) for v in (0, identifiers["y"] - 10 - box["y"], right - 12 - box["x"], box["height"])]
    shot.paste((255, 255, 255), cut)
    shot.save(path)
    page.close()
    return {"ligand": ligand}


def capture_sf_viewer(browser, url, pair=PAIR):
    """Two crops of the separation-factor viewer with `pair` selected: its element selectors,
    and the plot they drive. The options between them (x-axis choice, solvent filter) are left
    out to keep the figure short."""
    page = open_dashboard(browser, url, DASHBOARD_WIDTH)
    page.locator("[data-testid='stSidebarCollapseButton'] button").first.evaluate("e => e.click()")
    page.wait_for_timeout(2000)
    for side, element in zip(("ln1", "ln2"), pair):
        # each selector's buttons carry the key f"{side}_{element}"; a click reruns the app
        page.locator(f".st-key-{side}_{element} button[data-testid^='stBaseButton']").first.evaluate("e => e.click()")
        page.wait_for_selector(f"text={side.capitalize()}: {ELEMENT_NAMES[element]} ({element})", timeout=60_000)
        page.wait_for_timeout(4000)
    if not page.get_by_text(f"logSF[{pair[0]}/{pair[1]}]").count():
        raise RuntimeError(f"the viewer did not switch to {pair[0]}/{pair[1]}")
    head = page.get_by_text("Select First Element").first
    head.evaluate("e => e.scrollIntoView({block: 'start'})")
    page.wait_for_timeout(4000)
    top = head.bounding_box()["y"]
    caption = page.get_by_text(f"Ln2: {ELEMENT_NAMES[pair[1]]}").first.bounding_box()
    page.screenshot(path=OUT / "sf_selectors.png", clip=clip(page, top, caption["y"] + caption["height"] + 4))
    chart = next(b for b in (e.bounding_box() for e in page.locator("[data-testid='stPlotlyChart']").all())
                 if b and b["y"] > caption["y"])
    page.screenshot(path=OUT / "sf_plot.png", clip=clip(page, chart["y"], chart["y"] + chart["height"]))
    page.close()
    return {"pair": "/".join(pair)}


def capture_mcp(python, repo):
    name, arguments = MCP_CALL
    env = dict(os.environ, PYTHONPATH=str(repo), RELQUEST_DB=str(DATABASE))
    run = subprocess.run([python, "-c", MCP_CLIENT, json.dumps([name, arguments, env])],
                         cwd=repo, capture_output=True, text=True, timeout=300)
    if run.returncode:
        raise RuntimeError(f"MCP client failed:\n{run.stderr[-2000:]}")
    exchange = json.loads(run.stdout.strip().splitlines()[-1])
    (OUT / "mcp_exchange.json").write_text(json.dumps(exchange, indent=2, ensure_ascii=False) + "\n",
                                           encoding="utf-8")
    return exchange


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--app-python", default="python", help="interpreter that runs the dashboard")
    parser.add_argument("--mcp-python", default=sys.executable, help="interpreter with the mcp SDK")
    parser.add_argument("--mcp-repo", default=str(ROOT), help="checkout holding rel_quest/mcp_server.py")
    parser.add_argument("--free-model", default=FREE_MODEL,
                        help="free-mode model for the assistant ('' keeps the dashboard default)")
    parser.add_argument("--app-root", default=str(ROOT),
                        help="tree the dashboard runs from (e.g. the public export)")
    parser.add_argument("--only", choices=["browser", "sf_viewer"], help="recapture this part alone")
    args = parser.parse_args()
    from playwright.sync_api import sync_playwright

    OUT.mkdir(parents=True, exist_ok=True)
    app_root = Path(args.app_root).resolve()
    record = OUT / "capture.json"
    meta = json.loads(record.read_text(encoding="utf-8")) if args.only and record.exists() else {}
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    if not args.only:
        mcp_repo = Path(args.mcp_repo).resolve()
        exchange = capture_mcp(args.mcp_python, mcp_repo)
    proc, url = start_dashboard(args.app_python, app_root)
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            if not args.only:
                assistant = capture_assistant(browser, url, args.free_model)
            parts = {name: fn(browser, url) for name, fn in
                     (("browser", capture_browser), ("sf_viewer", capture_sf_viewer))
                     if args.only in (None, name)}
            browser.close()
    finally:
        proc.terminate()
    provenance = app_root / "RELEASE_PROVENANCE.json"
    release = json.loads(provenance.read_text(encoding="utf-8")) if provenance.exists() else {}
    for name, details in parts.items():
        meta[name] = {"captured": now, "app_root_commit": git_commit(app_root),
                      "exported_from": (release.get("dev_commit") or "")[:7] or None,
                      "exported_tree_clean": release.get("dev_tree_clean"), **details}
    if not args.only:
        meta.update({"captured": now, "commit": git_commit(ROOT), "mcp_repo_commit": git_commit(mcp_repo),
                     "database": DATABASE.name, "assistant": assistant,
                     "mcp_call": {"tool": exchange["call"], "arguments": exchange["arguments"]}})
    record.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    files = {"browser": "browser_table.png, browser_detail.png", "sf_viewer": "sf_selectors.png, sf_plot.png"}
    print(f"wrote {OUT.relative_to(ROOT)}: " + ("" if args.only else "assistant.png, mcp_exchange.json, ")
          + ", ".join(files[name] for name in parts) + ", capture.json")
    if not args.only:
        print(f"assistant: {assistant['model']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
