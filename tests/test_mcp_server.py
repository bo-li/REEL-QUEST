"""End-to-end test of the MCP server (rel_quest/mcp_server.py).

The server is launched as a real subprocess and spoken to over stdio with the
MCP client SDK, in the classic `initialize` handshake that Claude Desktop and
Claude Code use. Every answer is checked against the raw JSONL, read here
without rel_quest, so a bug shared by the server and LigandDatabase cannot
confirm itself.

Ways the server can fail, written before it existed:

 1. Anything printed to stdout corrupts the protocol, since stdout *is* the
    channel (LigandDatabase used to print "Loaded N ligands" there). The
    handshake fails.
 2. The tool list drifts from the vocabulary: a trend or mechanism enum that
    misses a value makes those ligands unreachable (the CLI rejected Peaked
    and U_Shaped this way).
 3. Plot tools leak through. They return Plotly figures, which no MCP client
    can render.
 4. An unknown ligand or an invalid enum value comes back as a normal result
    instead of an error, and the client model reads it as an answer.
 5. Condition parameters come back empty: agent_tools read `acid_name`,
    `solvent`, ... while records store `Acid_Name`, `Solvent_Name`, ...
 6. A condition reaches the client without the DOI of the paper that
    measured it. The data are CC BY 4.0, and attribution has to survive.
 7. A separation factor divides medians taken from different conditions or
    papers, instead of pairing both elements within one condition.
 8. A bound (D > x, D < x) or a mixed-extractant condition is pooled into a
    median as if it were a value.
 9. Y and Sc queries return nothing (element_coverage is La-Lu only).
10. Pointing the server at a data file outside the checkout (a non-editable
    install plus the Zenodo JSONL) fails because something still reads
    data/ or vocab/ relative to the package. KNOWN GAP: utils.load_vocabulary,
    plotting and schema read vocab/ beside the package, so this is an
    expected failure until that is decided (docs/MCP_SERVER.md).

Over HTTP (`--transport http`, the hosted mode for Hugging Face Spaces), also
written before the code:

H1. The HTTP transport answers differently from stdio. Every check above
    reruns over HTTP (TestMCPServerOverHTTP).
H2. --port (or $PORT) is ignored, and the server lands on a port the proxy
    does not route to.
H3. Behind the Spaces proxy the Host header is <space>.hf.space. A server
    bound to 0.0.0.0 must accept it. Bound to 127.0.0.1 it must refuse it:
    the SDK's DNS-rebinding guard, which is right for a laptop.
H4. A free Space sleeps and restarts. A stateful server then answers the
    client's old Mcp-Session-Id with 404, so it must be stateless.
H5. Replies come back as SSE streams, which a buffering proxy can hold back.
    They must be plain JSON.
H6. Opening the Space URL in a browser shows a 404. GET / must explain the
    server and give the connector URL, with https:// when the proxy
    terminates TLS (X-Forwarded-Proto). That URL must be the one to paste.
H7. Only one protocol era works. Both the classic initialize handshake and
    the 2026-07-28 discovery ("auto") must connect.

On the Hugging Face Space (space/app.py: a Gradio demo, with this server
mounted at /mcp of the same app), also written before the code:

S1. The mounted server answers differently from the standalone one. Every
    check above reruns against /mcp of the Gradio app (TestMCPServerInSpaceApp).
S2. The session manager's lifespan never starts inside Gradio, so every /mcp
    request fails ("Task group is not initialized").
S3. Gradio's own MCP mode switches on (GRADIO_MCP_SERVER=true in the Space's
    settings) and adds a second, different tool set at /gradio_api/mcp.
    Gradio 6.30 cannot build its server against mcp 2.x and quietly turns the
    mode off; the app also turns it off, for a Gradio that can.
S4. Gradio shadows /mcp (route order, or SSR's Node front proxy) and answers
    it with its own 404 page.
S5. Bound to 0.0.0.0, as on a Space, the DNS-rebinding guard still refuses
    the Space's host name with 421.
S6. The Gradio page is replaced by the MCP landing page, or shows a connector
    URL other than the Space's https://<SPACE_HOST>/mcp.
S7. A demo tab disagrees with the MCP tool it is built on, or answers an
    unknown ligand with an empty table instead of an error.
S8. The bundle misses a file the app reads (vocab/, the database, a module)
    and the Space crashes on import, or it carries files beyond those. The
    app is launched from the bundle scripts/build_space.py assembles, not
    from the checkout.
S9. The Space card or requirements are refused or misread: a
    short_description over 60 characters, an sdk_version other than the
    Gradio tested here, a Python ZeroGPU does not run, or a pin on gradio,
    spaces or huggingface_hub (the platform installs those).
S11. The REST API (rel_quest/api.py, mounted at /api) is missing from the
    Space, or answers differently there from the standalone app on the same
    data (TestRESTAPIInSpaceApp).
S12. Mounted at /api it shadows Gradio (the page, /config, /gradio_api) or
    /mcp, or Gradio shadows it.
S13. Its documentation does not load behind Gradio.
S14. Gradio's middleware drops the API's CORS header, so a browser tool on
    another site cannot call it.
S15. The demo does not tell a visitor where the REST API is: the "Use with an
    AI assistant" tab must give https://<SPACE_HOST>/api/docs.
S10. ZeroGPU will not start a Space without a @spaces.GPU function. Only the
    live Space shows this (the spaces package is a no-op anywhere else); see
    docs/MCP_SERVER.md.

Artifacts: outputs/mcp_e2e_transcript.json (stdio),
outputs/mcp_e2e_transcript_http.json and outputs/mcp_e2e_transcript_space.json,
every request and response.
"""

import asyncio
import json
import math
import os
import shutil
import statistics
import sys
import socket
import subprocess
import tempfile
import time
import unittest
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB_FILE = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
TRANSCRIPT = ROOT / "outputs" / "mcp_e2e_transcript.json"

try:
    from mcp import Client
    from mcp.client.stdio import StdioServerParameters
except ImportError:  # pragma: no cover - optional extra
    Client = None


def raw_records(path=DB_FILE):
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def raw_by_name(records, name):
    return next(r for r in records if r["ligand"]["Name"] == name)


def censored(cell):
    return cell.get("censor") in (">", "<")


def payload(result):
    """The JSON a tool returned, from structured content or the text block."""
    if result.structured_content is not None:
        data = result.structured_content
        # A non-object return is wrapped as {"result": ...}
        return data["result"] if set(data) == {"result"} else data
    return json.loads(result.content[0].text)


def server_params(env=None, cwd=ROOT, pythonpath=None):
    full_env = dict(os.environ)
    full_env.pop("RELQUEST_DB", None)
    full_env.update(env or {})
    if pythonpath:
        full_env["PYTHONPATH"] = str(pythonpath)
    return StdioServerParameters(
        command=sys.executable,
        args=["-m", "rel_quest.mcp_server"],
        cwd=str(cwd),
        env=full_env,
    )


class Session:
    """Run one client session and record the transcript."""

    def __init__(self, params, log, mode="legacy"):
        self.params = params
        self.log = log
        self.mode = mode

    async def run(self, steps):
        async with Client(self.params, mode=self.mode, read_timeout_seconds=120) as client:
            self.log.append({"event": "initialize",
                             "server": client.server_info.name if client.server_info else None,
                             "instructions": client.instructions})
            return await steps(client, self)

    async def call(self, client, name, arguments=None):
        result = await client.call_tool(name, arguments or {})
        entry = {"call": name, "arguments": arguments or {}, "is_error": bool(result.is_error)}
        if result.is_error:
            entry["error"] = result.content[0].text if result.content else None
        else:
            entry["result"] = payload(result)
        self.log.append(entry)
        return result


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
class TestMCPServerEndToEnd(unittest.TestCase):

    transcript = TRANSCRIPT

    @classmethod
    def target(cls):
        """What the client connects to: a stdio subprocess here, a URL over HTTP."""
        return server_params()

    @classmethod
    def setUpClass(cls):
        cls.records = raw_records()
        cls.log = []
        cls.results = asyncio.run(Session(cls.target(), cls.log).run(cls._exercise))

    @classmethod
    def tearDownClass(cls):
        cls.transcript.parent.mkdir(exist_ok=True)
        cls.transcript.write_text(json.dumps(cls.log, indent=1, ensure_ascii=False, default=str))

    @staticmethod
    async def _exercise(client, session):
        out = {"instructions": client.instructions}
        tools = await client.list_tools()
        out["tools"] = {t.name: t for t in tools.tools}
        session.log.append({"call": "tools/list", "result": sorted(out["tools"])})

        resources = await client.list_resources()
        out["resources"] = [str(r.uri) for r in resources.resources]
        vocab = await client.read_resource("reel-quest://vocab")
        out["vocab"] = json.loads(vocab.contents[0].text)
        session.log.append({"call": "resources/read", "uri": "reel-quest://vocab", "result": out["vocab"]})

        call = session.call
        out["stats"] = await call(client, "get_database_statistics")
        out["details"] = await call(client, "get_ligand_details", {"name": "TODGA"})
        out["details_mix"] = await call(client, "get_ligand_details", {"name": "T2EHDGA"})
        out["conditions"] = await call(client, "get_ligand_conditions", {"name": "TODGA", "limit": 500})
        out["conditions_eu"] = await call(client, "get_ligand_conditions",
                                          {"name": "TODGA", "element": "Eu", "limit": 500})
        out["by_name"] = await call(client, "search_ligands_by_name", {"name": "BTBP"})
        out["sf"] = await call(client, "get_best_separators",
                               {"element1": "Eu", "element2": "La", "limit": 500})
        out["y"] = await call(client, "get_ligands_for_element", {"element": "Y", "limit": 500})
        out["peaked"] = await call(client, "search_by_trend", {"trend": "Peaked", "limit": 500})
        out["solvation"] = await call(client, "search_by_mechanism",
                                      {"mechanism": "Solvation", "limit": 500})
        out["dga"] = await call(client, "search_by_subtype", {"subtype": "DGA", "limit": 500})
        out["compare"] = await call(client, "compare_ligands", {"names": ["TODGA", "No-Such-Ligand"]})
        out["unknown"] = await call(client, "get_ligand_details", {"name": "No-Such-Ligand-XYZ"})
        out["bad_trend"] = await call(client, "search_by_trend", {"trend": "Upward"})
        return out

    # --- 1. handshake / 3. tool surface / 2. vocabulary --------------------

    def test_handshake_and_instructions(self):
        instructions = self.results["instructions"] or ""
        for needle in ("CC BY 4.0", "DOI", "TDdDGA", "bound"):
            self.assertIn(needle, instructions)

    def test_tool_surface(self):
        expected = {
            "get_database_statistics", "search_ligands_by_name", "get_ligand_details",
            "get_ligand_conditions", "search_by_trend", "search_by_mechanism",
            "search_by_subtype", "get_best_separators", "compare_ligands",
            "get_ligands_for_element",
        }
        self.assertEqual(set(self.results["tools"]), expected)
        for name, tool in self.results["tools"].items():
            self.assertFalse(name.startswith("plot_"))
            self.assertTrue(tool.description, f"{name} has no description")
            self.assertTrue(tool.annotations and tool.annotations.read_only_hint,
                            f"{name} is not marked read-only")

    def test_enums_cover_every_stored_value(self):
        def enum(tool, param):
            schema = self.results["tools"][tool].input_schema
            prop = schema["properties"][param]
            if "enum" in prop:
                return set(prop["enum"])
            ref = prop["$ref"].split("/")[-1]
            return set(schema["$defs"][ref]["enum"])

        stored = lambda key, field: {r[key][field] for r in self.records if r[key].get(field)}
        self.assertEqual(enum("search_by_trend", "trend"),
                         {"Ascending", "Reverse", "Flat", "Mixed", "Peaked", "U_Shaped",
                          "Insufficient_Data"})
        self.assertLessEqual(stored("selectivity", "trend_classification"),
                             enum("search_by_trend", "trend"))
        self.assertLessEqual(stored("ligand", "Mechanism"), enum("search_by_mechanism", "mechanism"))
        self.assertLessEqual(stored("ligand", "Subtype"), enum("search_by_subtype", "subtype"))

        vocab = self.results["vocab"]
        self.assertIn("reel-quest://vocab", self.results["resources"])
        self.assertEqual(sum(vocab["subtypes"].values()), len(self.records))
        self.assertIn("Y", vocab["elements"])

    # --- 4. errors ---------------------------------------------------------

    def test_unknown_name_and_bad_enum_are_errors(self):
        self.assertTrue(self.results["unknown"].is_error)
        self.assertTrue(self.results["bad_trend"].is_error)
        compare = payload(self.results["compare"])
        self.assertEqual(compare["not_found"], ["No-Such-Ligand"])

    # --- counts ------------------------------------------------------------

    def test_statistics_match_the_file(self):
        stats = payload(self.results["stats"])
        self.assertEqual(stats["total_ligands"], len(self.records))
        peaked = payload(self.results["peaked"])
        self.assertEqual(peaked["total_found"],
                         sum(r["selectivity"]["trend_classification"] == "Peaked" for r in self.records))
        solvation = payload(self.results["solvation"])
        self.assertEqual(solvation["total_found"],
                         sum(r["ligand"]["Mechanism"] == "Solvation" for r in self.records))
        dga = payload(self.results["dga"])
        self.assertEqual(dga["total_found"], sum(r["ligand"]["Subtype"] == "DGA" for r in self.records))
        by_name = payload(self.results["by_name"])
        self.assertTrue(by_name and all("btbp" in hit["name"].lower() for hit in by_name))

    # --- 5. condition parameters / 6. DOIs ---------------------------------

    def test_details_carry_real_condition_parameters(self):
        details = payload(self.results["details"])
        raw = raw_by_name(self.records, "TODGA")
        self.assertEqual(details["name"], "TODGA")
        first, raw_first = details["sample_conditions"][0], raw["experimental_data"][0]
        self.assertEqual(first["acid"], raw_first["Acid_Name"])
        self.assertEqual(first["acid_conc_M"], raw_first["Acid_Concentration_M"])
        self.assertEqual(first["solvent"], raw_first["Solvent_Name"])
        self.assertEqual(first["DOI"], raw_first["reference"]["DOI"])
        for cond in details["sample_conditions"]:
            self.assertIsNotNone(cond["acid"])
            self.assertIsNotNone(cond["DOI"])

    def test_conditions_are_complete_and_cited(self):
        raw = raw_by_name(self.records, "TODGA")
        conds = payload(self.results["conditions"])
        self.assertEqual(conds["total_conditions"], len(raw["experimental_data"]))
        self.assertEqual(len(conds["conditions"]), len(raw["experimental_data"]))
        by_id = {c["condition_id"]: c for c in raw["experimental_data"]}
        for cond in conds["conditions"]:
            source = by_id[cond["condition_id"]]
            self.assertEqual(cond["reference"]["DOI"], source["reference"]["DOI"])
            self.assertTrue(cond["reference"]["DOI"])
            self.assertEqual(cond["is_mixture"], bool(source.get("Co_Extractants")))
            self.assertEqual(set(cond["D"]), set(source["distribution_ratios"]))
            for elem, cell in source["distribution_ratios"].items():
                self.assertEqual(cond["D"][elem]["censor"], cell.get("censor") if censored(cell) else None)
                self.assertAlmostEqual(cond["D"][elem]["log_D"], cell["log_value"], places=6)

        eu = payload(self.results["conditions_eu"])
        self.assertEqual(eu["total_conditions"],
                         sum("Eu" in c["distribution_ratios"] for c in raw["experimental_data"]))

    # --- 7. paired separation factors / 8. bounds and mixtures --------------

    def test_separation_factors_are_paired_within_conditions(self):
        sf = payload(self.results["sf"])
        expected = {}
        for record in self.records:
            pairs = []
            for cond in record["experimental_data"]:
                dr = cond["distribution_ratios"]
                if cond.get("Co_Extractants") or "Eu" not in dr or "La" not in dr:
                    continue
                if censored(dr["Eu"]) or censored(dr["La"]):
                    continue
                if dr["Eu"]["log_value"] is None or dr["La"]["log_value"] is None:
                    continue  # 19 empty cells, e.g. BAPhen (Pramanik 2024)
                pairs.append(dr["Eu"]["log_value"] - dr["La"]["log_value"])
            if pairs:
                expected[record["ligand"]["Name"]] = (statistics.median(pairs), len(pairs))

        got = {row["name"]: row for row in sf["top_separators"]}
        self.assertEqual(set(got), set(expected))
        for name, (median, n) in expected.items():
            self.assertAlmostEqual(got[name]["median_log_SF"], median, delta=0.001, msg=name)
            self.assertEqual(got[name]["n_paired_conditions"], n, msg=name)
        ranks = [abs(row["median_log_SF"]) for row in sf["top_separators"]]
        self.assertEqual(ranks, sorted(ranks, reverse=True))

    def test_medians_leave_out_bounds_and_mixtures(self):
        raw = raw_by_name(self.records, "T2EHDGA")
        self.assertTrue(any(c.get("Co_Extractants") for c in raw["experimental_data"]))
        details = payload(self.results["details_mix"])
        for elem, summary in details["median_logD_by_element"].items():
            values = [c["distribution_ratios"][elem]["log_value"]
                      for c in raw["experimental_data"]
                      if not c.get("Co_Extractants") and elem in c["distribution_ratios"]
                      and not censored(c["distribution_ratios"][elem])]
            self.assertAlmostEqual(summary["median_logD"], statistics.median(values),
                                   delta=0.001, msg=elem)
            self.assertEqual(summary["n_conditions"], len(values), msg=elem)

    # --- 9. Y and Sc --------------------------------------------------------

    def test_yttrium_is_queryable(self):
        y = payload(self.results["y"])
        expected = sum(any("Y" in c["distribution_ratios"] for c in r["experimental_data"])
                       for r in self.records)
        self.assertGreater(expected, 0)
        self.assertEqual(y["total_ligands_with_data"], expected)


# --- HTTP ---------------------------------------------------------------------

def outside_address():
    """This machine's non-loopback IPv4 address, or None. No packet is sent:
    connecting a UDP socket only picks the outgoing interface."""
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
            sock.connect(("192.0.2.1", 9))  # TEST-NET-1, never routed
            address = sock.getsockname()[0]
    except OSError:
        return None
    return None if address.startswith("127.") else address


def free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def http(method, url, body=None, headers=None, timeout=60):
    """(status, headers, text) for one request; HTTP errors are answers too."""
    data = json.dumps(body).encode() if body is not None else None
    request = urllib.request.Request(url, data=data, method=method, headers=headers or {})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.headers, response.read().decode()
    except urllib.error.HTTPError as error:
        return error.code, error.headers, error.read().decode()


JSON_RPC_HEADERS = {"Content-Type": "application/json",
                    "Accept": "application/json, text/event-stream"}
INITIALIZE = {"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
    "protocolVersion": "2025-06-18", "capabilities": {},
    "clientInfo": {"name": "raw-http-check", "version": "0"}}}


class HTTPServer:
    """`python -m rel_quest.mcp_server --transport http` as a subprocess."""

    cwd = ROOT

    def __init__(self, host="127.0.0.1"):
        self.port = free_port()
        self.base = f"http://127.0.0.1:{self.port}"
        # stderr to a file: uvicorn's access log would fill an undrained pipe.
        self.log = tempfile.TemporaryFile(mode="w+")
        self.proc = subprocess.Popen(
            self.command(host), cwd=str(self.cwd), env=self.environment(host),
            stdout=self.log, stderr=self.log,
        )
        deadline = time.time() + 180  # RDKit + the database load slowly on /mnt/c
        while time.time() < deadline:
            if self.proc.poll() is not None:
                break
            try:
                if http("GET", self.base + "/", timeout=5)[0] == 200:
                    return
            except OSError:
                pass
            time.sleep(1)
        output = self.output()
        self.stop()
        raise RuntimeError(f"HTTP server did not come up:\n{output[-3000:]}")

    def command(self, host):
        return [sys.executable, "-m", "rel_quest.mcp_server", "--transport", "http",
                "--host", host, "--port", str(self.port)]

    def environment(self, host):
        env = dict(os.environ)
        env.pop("RELQUEST_DB", None)
        return env

    def output(self):
        self.log.seek(0)
        return self.log.read()

    def stop(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(10)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
        self.log.close()


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
class TestMCPServerOverHTTP(TestMCPServerEndToEnd):
    """H1: every check of the stdio class, inherited, over streamable HTTP."""

    transcript = ROOT / "outputs" / "mcp_e2e_transcript_http.json"

    @classmethod
    def target(cls):
        return cls.server.base + "/mcp"

    @classmethod
    def setUpClass(cls):
        cls.server = HTTPServer()  # H2: on the port it was given
        try:
            super().setUpClass()
        except BaseException:
            cls.server.stop()
            raise

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.server.stop()

    def test_loopback_binding_refuses_a_foreign_host(self):
        """H3, the local half: bound to 127.0.0.1 the DNS-rebinding guard holds."""
        status, _, _ = http("POST", self.target(), INITIALIZE,
                            {**JSON_RPC_HEADERS, "Host": "reel-quest.hf.space"})
        self.assertEqual(status, 421)

    def test_stateless_a_stale_session_id_still_works(self):
        """H4: the request a client sends after the Space slept and restarted."""
        status, _, text = http("POST", self.target(),
                               {"jsonrpc": "2.0", "id": 7, "method": "tools/list"},
                               {**JSON_RPC_HEADERS, "Mcp-Session-Id": "stale-from-before-restart",
                                "MCP-Protocol-Version": "2025-06-18"})
        self.assertEqual(status, 200, text)
        names = {tool["name"] for tool in json.loads(text)["result"]["tools"]}
        self.assertIn("get_ligand_conditions", names)

    def test_replies_are_plain_json(self):
        """H5"""
        status, headers, text = http("POST", self.target(), INITIALIZE, JSON_RPC_HEADERS)
        self.assertEqual(status, 200, text)
        self.assertTrue(headers["Content-Type"].startswith("application/json"),
                        headers["Content-Type"])
        self.assertEqual(json.loads(text)["id"], 1)

    def test_landing_page_gives_the_connector_url(self):
        """H6, direct: the absolute /mcp URL of this server."""
        status, headers, text = http("GET", self.server.base + "/")
        self.assertEqual(status, 200)
        self.assertTrue(headers["Content-Type"].startswith("text/html"))
        self.assertIn(self.target(), text)
        self.assertIn("CC BY 4.0", text)
        self.assertIn(str(len(self.records)), text)

    def test_modern_discovery_connects(self):
        """H7: the classic handshake is what the inherited checks use."""
        log = []

        async def steps(client, session):
            return {t.name for t in (await client.list_tools()).tools}

        tools = asyncio.run(Session(self.target(), log, mode="auto").run(steps))
        self.assertIn("get_best_separators", tools)


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
class TestMCPServerBehindAProxy(unittest.TestCase):
    """H3 and H6 as the Spaces proxy sees them: bound to 0.0.0.0, addressed
    by the Space's host name, TLS terminated in front.

    Requests go to the non-loopback address, as the proxy's do. uvicorn
    trusts X-Forwarded-* from 127.0.0.1 only, so over loopback it would fix
    the scheme itself and hide a landing page that ignores the header."""

    @classmethod
    def setUpClass(cls):
        address = outside_address()
        if address is None:
            raise unittest.SkipTest("no non-loopback interface to stand in for the proxy")
        cls.server = HTTPServer(host="0.0.0.0")
        cls.server.base = f"http://{address}:{cls.server.port}"

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()

    PROXY = {"Host": "reel-quest.hf.space", "X-Forwarded-Proto": "https"}

    def test_accepts_the_space_host_name(self):
        status, _, text = http("POST", self.server.base + "/mcp", INITIALIZE,
                               {**JSON_RPC_HEADERS, **self.PROXY})
        self.assertEqual(status, 200, text)
        self.assertEqual(json.loads(text)["result"]["serverInfo"]["name"], "reel-quest")

    def test_landing_page_advertises_the_https_url(self):
        status, _, text = http("GET", self.server.base + "/", headers=self.PROXY)
        self.assertEqual(status, 200)
        self.assertIn("https://reel-quest.hf.space/mcp", text)
        self.assertNotIn("http://reel-quest.hf.space", text)


# --- Hugging Face Space ---------------------------------------------------------

# A Python with gradio, spaces and space/requirements.txt, like the Space's
# (docs/MCP_SERVER.md, "Testing"). The bundle checks run without it.
SPACE_PYTHON = os.environ.get("RELQUEST_SPACE_PYTHON")
if SPACE_PYTHON and os.sep in SPACE_PYTHON:
    # The app runs from inside its bundle, so a path relative to where pytest
    # started (the docs' build/space-venv/bin/python) must be made absolute,
    # without resolving the venv's symlink to its base interpreter.
    SPACE_PYTHON = str(Path(SPACE_PYTHON).absolute())
elif SPACE_PYTHON:
    SPACE_PYTHON = shutil.which(SPACE_PYTHON) or SPACE_PYTHON
SPACE = ROOT / "space"
PLATFORM_PACKAGES = {"gradio", "spaces", "huggingface_hub", "huggingface-hub"}


def build_space(dest):
    """The bundle `hf upload` sends, assembled by the script that sends it."""
    proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_space.py"),
                           "--out", str(dest)], capture_output=True, text=True, timeout=300)
    if proc.returncode:
        raise RuntimeError(f"build_space.py failed:\n{proc.stdout}\n{proc.stderr}")
    return dest


def front_matter(path):
    """The YAML header of a Space card, as flat key: value pairs."""
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "---", "the card must open with front matter"
    card = {}
    for line in lines[1:lines.index("---", 1)]:
        key, _, value = line.partition(":")
        if value.strip():
            card[key.strip()] = value.strip().strip('"')
    return card


def gradio_call(base, api_name, *args):
    """(event, data) of one call to a Gradio endpoint: POST starts the job,
    GET streams its events; the last is "complete" or "error"."""
    status, _, text = http("POST", f"{base}/gradio_api/call/{api_name}", {"data": list(args)},
                           {"Content-Type": "application/json"})
    if status != 200:
        raise AssertionError(f"{api_name}: {status} {text[:500]}")
    event_id = json.loads(text)["event_id"]
    _, _, stream = http("GET", f"{base}/gradio_api/call/{api_name}/{event_id}", timeout=120)
    event = data = None
    for line in stream.splitlines():
        if line.startswith("event:"):
            event = line.partition(":")[2].strip()
        elif line.startswith("data:"):
            data = line.partition(":")[2].strip()
    return event, json.loads(data) if data else None


class SpaceApp(HTTPServer):
    """space/app.py as the Space runs it: `python app.py` inside the bundle,
    with Gradio's address and port, and the Space's host name, from the
    environment."""

    SPACE_HOST = "bolivu-reel-quest.hf.space"

    def __init__(self, bundle, host="0.0.0.0"):
        self.cwd = bundle
        super().__init__(host)

    def command(self, host):
        return [SPACE_PYTHON, "app.py"]

    def environment(self, host):
        env = super().environment(host)
        env.pop("PYTHONPATH", None)  # the bundle's own rel_quest, not the checkout's
        env.update({
            "GRADIO_SERVER_NAME": host,
            "GRADIO_SERVER_PORT": str(self.port),
            "SPACE_HOST": self.SPACE_HOST,
            "GRADIO_MCP_SERVER": "true",  # S3: as if set in the Space's settings
            "GRADIO_ANALYTICS_ENABLED": "False",
        })
        return env


class TestSpaceBundle(unittest.TestCase):
    """S8 and S9: what scripts/build_space.py hands to `hf upload`."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.bundle = build_space(Path(cls.tmp.name) / "space")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_bundle_holds_what_the_app_reads_and_nothing_else(self):
        top = {p.name for p in self.bundle.iterdir()}
        self.assertEqual(top, {"app.py", "README.md", "requirements.txt", "rel_quest", "vocab",
                               "data", "LICENSE", "LICENSE-DATA", "CITATION.cff"})
        self.assertEqual({p.name for p in (self.bundle / "data").iterdir()}, {DB_FILE.name})
        self.assertEqual((self.bundle / "data" / DB_FILE.name).read_bytes(), DB_FILE.read_bytes())
        for folder, pattern in (("rel_quest", "*.py"), ("vocab", "*.json")):
            self.assertEqual({p.name for p in (self.bundle / folder).iterdir()},
                             {p.name for p in (ROOT / folder).glob(pattern)}, folder)
        self.assertFalse([p for p in self.bundle.rglob("*") if p.suffix == ".pyc"
                          or p.name == "__pycache__"])

    def test_card_and_requirements_are_what_the_platform_accepts(self):
        card = front_matter(self.bundle / "README.md")
        self.assertEqual(card["sdk"], "gradio")
        self.assertEqual(card["app_file"], "app.py")
        self.assertLessEqual(len(card["short_description"]), 60)
        self.assertIn(card["python_version"], ("3.10", "3.12", "3.10.13", "3.12.12"))
        self.assertEqual(card["license"], "cc-by-4.0")
        if SPACE_PYTHON:
            tested = subprocess.run([SPACE_PYTHON, "-c", "import gradio; print(gradio.__version__)"],
                                    capture_output=True, text=True, check=True).stdout.strip()
            self.assertEqual(card["sdk_version"], tested)

        requirements = [line.split("#")[0].strip() for line in
                        (self.bundle / "requirements.txt").read_text().splitlines()]
        names = {line.split(">")[0].split("=")[0].split("<")[0].strip().lower()
                 for line in requirements if line}
        self.assertFalse(names & PLATFORM_PACKAGES)
        self.assertIn("mcp", names)


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
@unittest.skipUnless(SPACE_PYTHON, "set RELQUEST_SPACE_PYTHON to a Python with gradio "
                                   "(docs/MCP_SERVER.md, Testing)")
class TestMCPServerInSpaceApp(TestMCPServerEndToEnd):
    """S1: every check of the stdio class, against /mcp of the Space's Gradio
    app, launched from the bundle as the Space launches it."""

    transcript = ROOT / "outputs" / "mcp_e2e_transcript_space.json"

    @classmethod
    def target(cls):
        return cls.server.base + "/mcp"

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        try:
            cls.server = SpaceApp(build_space(Path(cls.tmp.name) / "space"))
        except BaseException:
            cls.tmp.cleanup()
            raise
        try:
            super().setUpClass()
        except BaseException:
            cls.server.stop()
            cls.tmp.cleanup()
            raise

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.server.stop()
        cls.tmp.cleanup()

    PROXY = {"Host": SpaceApp.SPACE_HOST, "X-Forwarded-Proto": "https"}

    def test_space_host_name_is_accepted(self):
        """S5"""
        status, _, text = http("POST", self.target(), INITIALIZE, {**JSON_RPC_HEADERS, **self.PROXY})
        self.assertEqual(status, 200, text)
        self.assertEqual(json.loads(text)["result"]["serverInfo"]["name"], "reel-quest")

    def test_stateless_and_plain_json_inside_gradio(self):
        """H4 and H5 through Gradio's middleware."""
        status, headers, text = http("POST", self.target(),
                                     {"jsonrpc": "2.0", "id": 7, "method": "tools/list"},
                                     {**JSON_RPC_HEADERS, "Mcp-Session-Id": "stale-from-before-restart",
                                      "MCP-Protocol-Version": "2025-06-18"})
        self.assertEqual(status, 200, text)
        self.assertTrue(headers["Content-Type"].startswith("application/json"), headers["Content-Type"])
        self.assertEqual(len(json.loads(text)["result"]["tools"]), 10)

    def test_mcp_path_belongs_to_the_mcp_server(self):
        """S4: a GET gets the MCP transport's event stream, as from the
        standalone server, not Gradio's page. Only the headers are read: the
        stream stays open."""
        request = urllib.request.Request(self.target(), headers={"Accept": "text/event-stream"})
        with urllib.request.urlopen(request, timeout=30) as response:
            status, content_type = response.status, response.headers.get("Content-Type", "")
        self.assertEqual(status, 200)
        self.assertTrue(content_type.startswith("text/event-stream"), content_type)

    def test_gradio_mcp_mode_stays_off(self):
        """S3: GRADIO_MCP_SERVER=true is set for this app; no second server appears."""
        status, _, text = http("POST", self.server.base + "/gradio_api/mcp", INITIALIZE, JSON_RPC_HEADERS)
        self.assertEqual(status, 404, text[:300])

    def test_page_is_the_demo_with_the_space_connector_url(self):
        """S6"""
        status, headers, page = http("GET", self.server.base + "/")
        self.assertEqual(status, 200)
        self.assertTrue(headers["Content-Type"].startswith("text/html"))
        self.assertNotIn("REEL-QUEST MCP server</title>", page)
        status, _, config = http("GET", self.server.base + "/config")
        self.assertEqual(status, 200)
        self.assertIn(f"https://{SpaceApp.SPACE_HOST}/mcp", config)

    # --- S7: the demo tabs ---------------------------------------------------

    def test_separation_tab_matches_the_tool(self):
        event, data = gradio_call(self.server.base, "best_separators", "Eu", "La",
                                  len(self.records))
        self.assertEqual(event, "complete", data)
        _, table, plot = data
        tool = payload(self.results["sf"])["top_separators"]
        rows = [dict(zip(table["headers"], row)) for row in table["data"]]
        self.assertEqual([r["Ligand"] for r in rows], [t["name"] for t in tool])
        self.assertEqual([r["median log SF"] for r in rows], [t["median_log_SF"] for t in tool])
        self.assertEqual([r["paired conditions"] for r in rows], [t["n_paired_conditions"] for t in tool])
        self.assertTrue(plot)

    def test_ligand_tab_matches_the_tool(self):
        event, data = gradio_call(self.server.base, "ligand", "TODGA", "All")
        self.assertEqual(event, "complete", data)
        summary, table, plot = data
        self.assertIn("TODGA", summary)
        raw = raw_by_name(self.records, "TODGA")["experimental_data"]
        for doi in {c["reference"]["DOI"] for c in raw}:  # 6: every paper, as a link
            self.assertIn(f"(https://doi.org/{doi})", summary)
        rows = [dict(zip(table["headers"], row)) for row in table["data"]]
        self.assertEqual([r["Condition"] for r in rows], [c["condition_id"] for c in raw])
        for row, cond in zip(rows, raw):
            self.assertIn(cond["reference"]["DOI"], row["DOI"])
            for elem, cell in cond["distribution_ratios"].items():
                if cell["log_value"] is None:
                    continue
                shown = row[elem]
                self.assertTrue(shown.startswith(cell["censor"]) if censored(cell)
                                else shown[0] in "-0123456789", (elem, shown))
                # the stored log D, to the table's two decimals
                self.assertEqual(shown.lstrip("<> "), f"{cell['log_value']:.2f}", elem)
        self.assertTrue(plot)

        event, data = gradio_call(self.server.base, "ligand", "TODGA", "Y")
        self.assertEqual(event, "complete", data)
        self.assertEqual(len(data[1]["data"]),
                         sum("Y" in c["distribution_ratios"] for c in raw))

    def test_unknown_ligand_is_an_error(self):
        event, _ = gradio_call(self.server.base, "ligand", "No-Such-Ligand-XYZ", "All")
        self.assertEqual(event, "error")

    def test_element_tab_matches_the_tool(self):
        event, data = gradio_call(self.server.base, "ligands_for_element", "Y",
                                  len(self.records))
        self.assertEqual(event, "complete", data)
        summary, table = data
        tool = payload(self.results["y"])
        self.assertEqual([row[0] for row in table["data"]], [r["name"] for r in tool["results"]])
        self.assertIn(str(tool["total_ligands_with_data"]), summary)


@unittest.skipUnless(SPACE_PYTHON, "set RELQUEST_SPACE_PYTHON to a Python with gradio "
                                   "(docs/MCP_SERVER.md, Testing)")
class TestRESTAPIInSpaceApp(unittest.TestCase):
    """S11-S15: the REST API inside the Space's Gradio app, launched from the bundle."""

    PATHS = ["/api/v1/health", "/api/v1/stats", "/api/v1/vocab",
             "/api/v1/ligands?subtype=DGA&elements=Eu,La&limit=200",
             "/api/v1/ligands?elements=Y", "/api/v1/ligands/TODGA",
             "/api/v1/ligands/no-such-ligand"]

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        try:
            cls.server = SpaceApp(build_space(Path(cls.tmp.name) / "space"))
        except BaseException:
            cls.tmp.cleanup()
            raise

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        cls.tmp.cleanup()

    def get(self, path, headers=None):
        status, response_headers, text = http("GET", self.server.base + path, headers=headers or {})
        return status, response_headers, text

    def test_answers_as_the_standalone_app(self):                    # S11
        from fastapi.testclient import TestClient
        from rel_quest.api import create_app
        standalone = TestClient(create_app())
        for path in self.PATHS:
            with self.subTest(path):
                status, _, text = self.get(path)
                expected = standalone.get(path)
                self.assertEqual(status, expected.status_code, text[:300])
                self.assertEqual(json.loads(text), expected.json())

    def test_gradio_and_mcp_keep_their_paths(self):                  # S12
        status, headers, _ = self.get("/")
        self.assertEqual(status, 200)
        self.assertTrue(headers["Content-Type"].startswith("text/html"))
        self.assertEqual(self.get("/config")[0], 200)
        self.assertEqual(self.get("/gradio_api/info")[0], 200)
        event, data = gradio_call(self.server.base, "ligands_for_element", "Y", 5)
        self.assertEqual(event, "complete", data)
        status, _, text = http("POST", self.server.base + "/mcp", INITIALIZE, JSON_RPC_HEADERS)
        self.assertEqual(status, 200, text[:300])

    def test_docs_load_behind_gradio(self):                          # S13
        status, _, page = self.get("/api/docs")
        self.assertEqual(status, 200)
        self.assertIn("/api/v1/openapi.json", page)
        status, _, spec = self.get("/api/v1/openapi.json")
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(spec)["info"]["title"], "REEL-QUEST API")

    def test_cross_origin_reads_are_allowed(self):                   # S14
        status, headers, _ = self.get("/api/v1/health", {"Origin": "https://example.org"})
        self.assertEqual(status, 200)
        self.assertEqual(headers.get("Access-Control-Allow-Origin"), "*")

    def test_demo_points_to_the_api(self):                           # S15
        _, _, config = self.get("/config")
        self.assertIn(f"https://{SpaceApp.SPACE_HOST}/api/docs", config)


# The deployed Space, e.g. RELQUEST_SPACE_URL=https://bolivu-reel-quest.hf.space,
# with HF_TOKEN=$(hf auth token) while the Space is private.
LIVE_SPACE = os.environ.get("RELQUEST_SPACE_URL", "").rstrip("/")


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
@unittest.skipUnless(LIVE_SPACE, "set RELQUEST_SPACE_URL to the deployed Space (docs/MCP_SERVER.md)")
class TestMCPServerOnTheLiveSpace(TestMCPServerEndToEnd):
    """S1, S4 and S6 through Hugging Face's proxy, and the only check of S10:
    a Space that ZeroGPU refused to start answers nothing at all."""

    transcript = ROOT / "outputs" / "mcp_e2e_transcript_live_space.json"

    @classmethod
    def headers(cls):
        token = os.environ.get("HF_TOKEN")
        return {"Authorization": f"Bearer {token}"} if token else {}

    @classmethod
    def target(cls):
        from mcp.client.streamable_http import streamable_http_client
        from mcp.shared._httpx_utils import create_mcp_http_client

        return streamable_http_client(LIVE_SPACE + "/mcp",
                                      http_client=create_mcp_http_client(headers=cls.headers()))

    def test_mcp_path_reaches_python_not_a_front_proxy(self):
        """S4: with SSR on, Node fronts the app and /mcp would get its page."""
        request = urllib.request.Request(LIVE_SPACE + "/mcp",
                                         headers={"Accept": "text/event-stream", **self.headers()})
        with urllib.request.urlopen(request, timeout=60) as response:
            status, content_type = response.status, response.headers.get("Content-Type", "")
        self.assertEqual(status, 200)
        self.assertTrue(content_type.startswith("text/event-stream"), content_type)

    def test_page_shows_its_own_connector_url(self):
        """S6"""
        status, _, config = http("GET", LIVE_SPACE + "/config", headers=self.headers())
        self.assertEqual(status, 200)
        host = urllib.parse.urlsplit(LIVE_SPACE).netloc
        self.assertIn(f"https://{host}/mcp", config)


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
class TestMCPServerStdout(unittest.TestCase):
    """1. The SDK client skips stray stdout lines, so the session tests above
    pass even when the server prints to stdout. Other clients may not, so check
    the raw stream: every line must be a JSON-RPC message."""

    def test_stdout_carries_only_json_rpc(self):
        import subprocess

        params = server_params()
        initialize = {"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "raw-stdout-check", "version": "0"}}}
        proc = subprocess.run(
            [params.command, *params.args], cwd=params.cwd, env=params.env,
            input=json.dumps(initialize) + "\n", capture_output=True, text=True, timeout=180,
        )
        lines = [line for line in proc.stdout.splitlines() if line.strip()]
        self.assertTrue(lines, f"no response; stderr:\n{proc.stderr[-2000:]}")
        for line in lines:
            message = json.loads(line)  # a print() here raises
            self.assertEqual(message.get("jsonrpc"), "2.0")
        self.assertEqual(json.loads(lines[0]).get("id"), 1)


@unittest.skipIf(Client is None, "mcp SDK not installed (pip install -e '.[mcp]')")
class TestMCPServerOutsideTheCheckout(unittest.TestCase):
    """10. The package copied away from data/ and vocab/, data from RELQUEST_DB:
    what a non-editable install plus the Zenodo JSONL looks like."""

    # Known gap, not a regression: vocab/ (elements.json, solvents.json) is
    # still read beside the package. Remove this marker when the package can
    # be pointed at a data directory; pytest then reports it as passing.
    @unittest.expectedFailure
    def test_runs_on_a_data_file_named_by_env(self):
        records = raw_records()
        keep = [r for r in records if r["ligand"]["Name"] in ("TODGA", "T2EHDGA", "TBP")]
        self.assertEqual(len(keep), 3)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            shutil.copytree(ROOT / "rel_quest", tmp / "site" / "rel_quest",
                            ignore=shutil.ignore_patterns("__pycache__"))
            data = tmp / "zenodo.jsonl"
            data.write_text("".join(json.dumps(r) + "\n" for r in keep))
            params = server_params(env={"RELQUEST_DB": str(data)}, cwd=tmp,
                                   pythonpath=tmp / "site")
            log = []

            calls = [
                ("get_database_statistics", {}),
                ("search_ligands_by_name", {"name": "DGA"}),
                ("get_ligand_details", {"name": "TODGA"}),
                ("get_ligand_conditions", {"name": "TBP", "element": "Eu"}),
                ("search_by_trend", {"trend": "Ascending"}),
                ("search_by_mechanism", {"mechanism": "Solvation"}),
                ("search_by_subtype", {"subtype": "DGA"}),
                ("get_best_separators", {"element1": "Eu", "element2": "La"}),
                ("compare_ligands", {"names": ["TODGA", "T2EHDGA"]}),
                ("get_ligands_for_element", {"element": "Eu"}),
            ]

            async def steps(client, session):
                tools = await client.list_tools()
                results = [await session.call(client, name, args) for name, args in calls]
                return {t.name for t in tools.tools}, results

            tools, results = asyncio.run(Session(params, log).run(steps))
        self.assertEqual(tools, {name for name, _ in calls})
        for (name, _), result in zip(calls, results):
            self.assertFalse(result.is_error, f"{name}: {result.content}")
        self.assertEqual(payload(results[0])["total_ligands"], 3)
        self.assertEqual(payload(results[6])["total_found"], 2)  # TODGA, T2EHDGA


if __name__ == "__main__":
    unittest.main()
