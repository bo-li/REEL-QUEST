"""End-to-end test of the Claude Code skill (skills/reel-quest, bundled by
scripts/build_skill.py).

The MCP server answers ten fixed questions with the processing enforced. The
skill hands a code-running agent the database and rel_quest's canonical
helpers, so it can ask new questions without reinventing the processing.
Every check runs on a bundle built into a temporary directory, from another
working directory, with PYTHONPATH scrubbed: what a user gets after
installing it, not the checkout.

Ways the skill can fail, written before it existed:

K1. The bundle misses something the code reads (vocab/, the database, a
    module), or carries more than the skill needs, and works only inside
    the checkout.
K2. Claude Code refuses or never triggers the skill: a name that is not
    lowercase-hyphen within 64 characters, a description over 1024
    characters or without the words a chemist would use, a body too long to
    load cheaply.
K3. A recipe in SKILL.md rots: a helper is renamed or changes shape, and the
    instructions an agent copies no longer run, or print a number that is
    not the database's. Every python block runs here, from the bundle, and
    its numbers are checked against the raw JSONL read without rel_quest.
K4. A recipe drops the processing it exists to teach: pools a bound or a
    mixed-extractant condition, or pairs elements across conditions.
K5. query.py, the quick-lookup command, answers differently from the
    agent_tools function of the same name (which the MCP tests check
    against the raw file), prints anything but JSON on success, or exits 0
    on an unknown ligand or an invalid vocabulary value.
K6. A real Claude Code session never loads the skill, or loads it and gets
    the number wrong: for a question the fixed tools answer, for one only
    custom code answers, and for one that does not name the skill. Opt-in
    (RELQUEST_SKILL_LIVE=1): it runs `claude -p` three times.

Artifacts: outputs/skill_recipes_transcript.json (every recipe's output) and,
with K6, outputs/skill_live_transcript.json (the three sessions' events).
"""

import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB_FILE = ROOT / "data" / "sx_ligands_merged_v14_flat.jsonl"
SOURCE = ROOT / "skills" / "reel-quest"
RECIPES_TRANSCRIPT = ROOT / "outputs" / "skill_recipes_transcript.json"
LIVE_TRANSCRIPT = ROOT / "outputs" / "skill_live_transcript.json"
TRENDS = {"Ascending", "Reverse", "Flat", "Mixed", "Peaked", "U_Shaped", "Insufficient_Data"}


def raw_records():
    with open(DB_FILE, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def censored(cell):
    return cell.get("censor") in (">", "<")


def measured(cell):
    return cell.get("log_value") is not None and not censored(cell)


def build_skill(dest):
    proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_skill.py"), "--out", str(dest)],
                          capture_output=True, text=True, timeout=300)
    if proc.returncode:
        raise RuntimeError(f"build_skill.py failed:\n{proc.stdout}\n{proc.stderr}")
    return dest


def clean_env():
    env = dict(os.environ)
    for key in ("PYTHONPATH", "RELQUEST_DB"):
        env.pop(key, None)
    return env


def front_matter(text):
    lines = text.splitlines()
    assert lines[0] == "---", "SKILL.md must open with front matter"
    end = lines.index("---", 1)
    fields = {}
    for line in lines[1:end]:
        key, _, value = line.partition(":")
        fields[key.strip()] = value.strip()
    return fields, "\n".join(lines[end + 1:])


def python_blocks(body):
    """(label, code) of every ```python block; the label is its `# Recipe:` line."""
    blocks = []
    for code in re.findall(r"```python\n(.*?)```", body, flags=re.S):
        label = re.search(r"^# Recipe: (.+)$", code, flags=re.M)
        blocks.append((label.group(1).strip() if label else None, code))
    return blocks


# --- ground truth, from the raw file ------------------------------------------

def paired_log_sf(record, e1, e2, keep=lambda cond: True):
    """log D(e1) - log D(e2) within each non-mixture condition measuring both, no bounds."""
    out = []
    for cond in record["experimental_data"]:
        dr = cond["distribution_ratios"]
        if cond.get("Co_Extractants") or e1 not in dr or e2 not in dr or not keep(cond):
            continue
        if measured(dr[e1]) and measured(dr[e2]):
            out.append((dr[e1]["log_value"] - dr[e2]["log_value"], cond))
    return out


def strong_nitric(cond):
    return cond.get("Acid_Name") == "HNO3" and (cond.get("Acid_Concentration_M") or 0) >= 3


class SkillBundle:
    """A bundle in a temporary directory, and a separate working directory."""

    def __init__(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = build_skill(Path(self.tmp.name) / "reel-quest")
        self.cwd = Path(self.tmp.name) / "elsewhere"
        self.cwd.mkdir()

    def python(self, code, timeout=300):
        code = code.replace("SKILL_DIR", str(self.dir))  # what SKILL.md tells the agent to do
        return subprocess.run([sys.executable, "-"], input=code, cwd=self.cwd, env=clean_env(),
                              capture_output=True, text=True, timeout=timeout)

    def query(self, *args):
        return subprocess.run([sys.executable, str(self.dir / "scripts" / "query.py"), *args],
                              cwd=self.cwd, env=clean_env(), capture_output=True, text=True,
                              timeout=300)

    def cleanup(self):
        self.tmp.cleanup()


class TestSkillBundle(unittest.TestCase):
    """K1, K2"""

    @classmethod
    def setUpClass(cls):
        cls.bundle = SkillBundle()

    @classmethod
    def tearDownClass(cls):
        cls.bundle.cleanup()

    def test_bundle_holds_what_the_skill_reads_and_nothing_else(self):
        root = self.bundle.dir
        self.assertEqual({p.name for p in root.iterdir()},
                         {"SKILL.md", "scripts", "reference", "requirements.txt", "rel_quest",
                          "vocab", "data", "LICENSE", "LICENSE-DATA", "CITATION.cff"})
        self.assertEqual({p.name for p in (root / "data").iterdir()}, {DB_FILE.name})
        self.assertEqual((root / "data" / DB_FILE.name).read_bytes(), DB_FILE.read_bytes())
        for folder, pattern in (("rel_quest", "*.py"), ("vocab", "*.json")):
            self.assertEqual({p.name for p in (root / folder).iterdir()},
                             {p.name for p in (ROOT / folder).glob(pattern)}, folder)
        self.assertEqual({p.name for p in (root / "scripts").iterdir()}, {"query.py"})
        self.assertFalse([p for p in root.rglob("*") if p.suffix == ".pyc" or p.name == "__pycache__"])

    def test_front_matter_and_size_are_what_claude_code_loads(self):
        fields, body = front_matter((self.bundle.dir / "SKILL.md").read_text(encoding="utf-8"))
        self.assertEqual(set(fields), {"name", "description"})
        self.assertRegex(fields["name"], r"^[a-z0-9]+(-[a-z0-9]+)*$")
        self.assertLessEqual(len(fields["name"]), 64)
        self.assertEqual(fields["name"], self.bundle.dir.name)
        description = fields["description"]
        self.assertLessEqual(len(description), 1024)
        self.assertNotRegex(description, r"[<>]")
        for word in ("rare-earth", "lanthanide", "extractant", "separation factor",
                     "distribution ratio", "TODGA"):
            self.assertIn(word.lower(), description.lower(), word)
        self.assertLess(len(body.splitlines()), 500)


class TestSkillRecipes(unittest.TestCase):
    """K3 and K4: every python block in SKILL.md, run as an agent would run it."""

    @classmethod
    def setUpClass(cls):
        cls.records = raw_records()
        cls.by_name = {r["ligand"]["Name"]: r for r in cls.records}
        cls.bundle = SkillBundle()
        _, body = front_matter((cls.bundle.dir / "SKILL.md").read_text(encoding="utf-8"))
        cls.blocks = python_blocks(body)
        cls.outputs = {}
        for label, code in cls.blocks:
            proc = cls.bundle.python(code)
            cls.outputs[label] = proc
        RECIPES_TRANSCRIPT.parent.mkdir(exist_ok=True)
        RECIPES_TRANSCRIPT.write_text(json.dumps(
            [{"recipe": label, "code": code, "exit": cls.outputs[label].returncode,
              "stdout": cls.outputs[label].stdout, "stderr": cls.outputs[label].stderr[-2000:]}
             for label, code in cls.blocks], indent=1, ensure_ascii=False))

    @classmethod
    def tearDownClass(cls):
        cls.bundle.cleanup()

    CHECKED = {"Look up a ligand", "Paired separation factors under chosen conditions",
               "Cite the papers behind a result", "One element across ligands",
               "Recompute a selectivity trend on a subset"}

    def out(self, label):
        proc = self.outputs[label]
        self.assertEqual(proc.returncode, 0, f"{label}:\n{proc.stderr[-2000:]}")
        return proc.stdout

    def test_every_recipe_is_labelled_runs_and_is_checked(self):
        labels = [label for label, _ in self.blocks]
        self.assertNotIn(None, labels, "a python block without a '# Recipe:' line")
        self.assertEqual(set(labels), self.CHECKED)
        for label in labels:
            self.out(label)
            self.assertIn("SKILL_DIR", dict(self.blocks)[label])  # runs from the skill, not a checkout

    def test_look_up_a_ligand(self):
        raw = self.by_name["TODGA"]
        out = self.out("Look up a ligand")
        self.assertIn(raw["ligand"]["Subtype"], out)
        self.assertIn(raw["selectivity"]["trend_classification"], out)
        self.assertRegex(out, rf"\b{len(raw['experimental_data'])} conditions")

    def test_paired_separation_factors(self):
        pairs = paired_log_sf(self.by_name["TODGA"], "Eu", "La", strong_nitric)
        self.assertGreater(len(pairs), 1)
        median = statistics.median(v for v, _ in pairs)
        out = self.out("Paired separation factors under chosen conditions")
        self.assertIn(f"median log SF(Eu/La) = {median:.3f} over {len(pairs)} conditions", out)

    def test_cite_the_papers_behind_a_result(self):
        pairs = paired_log_sf(self.by_name["TODGA"], "Eu", "La", strong_nitric)
        dois = {cond["reference"]["DOI"] for _, cond in pairs}
        out = self.out("Cite the papers behind a result")
        shown = set(re.findall(r"10\.\d{4,9}/\S+", out))
        self.assertEqual(shown, dois)

    def test_one_element_across_ligands(self):
        expected = {}
        for record in self.records:
            values = [c["distribution_ratios"]["Y"]["log_value"] for c in record["experimental_data"]
                      if not c.get("Co_Extractants") and "Y" in c["distribution_ratios"]
                      and measured(c["distribution_ratios"]["Y"])]
            if values:
                expected[record["ligand"]["Name"]] = statistics.median(values)
        out = self.out("One element across ligands")
        for name, median in expected.items():
            self.assertRegex(out, rf"{re.escape(name)}\s+{median:.2f}", name)

    def test_recompute_a_selectivity_trend(self):
        out = self.out("Recompute a selectivity trend on a subset")
        full = re.search(r"all conditions: (\w+)", out).group(1)
        subset = re.search(r"HNO3 >= 3 M: (\w+)", out).group(1)
        self.assertEqual(full, self.by_name["TODGA"]["selectivity"]["trend_classification"])
        self.assertIn(subset, TRENDS)


class TestSkillQuery(unittest.TestCase):
    """K5: query.py against agent_tools, both run on the same file."""

    CALLS = [
        ("get_database_statistics", {}),
        ("search_ligands_by_name", {"name": "BTBP"}),
        ("get_ligand_details", {"name": "TODGA"}),
        ("get_ligand_conditions", {"name": "TBP", "element": "Eu"}),
        ("search_by_trend", {"trend": "Peaked"}),
        ("search_by_mechanism", {"mechanism": "Solvation"}),
        ("search_by_subtype", {"subtype": "DGA"}),
        ("get_best_separators", {"element1": "Eu", "element2": "La"}),
        ("compare_ligands", {"names": ["TODGA", "T2EHDGA"]}),
        ("get_ligands_for_element", {"element": "Y"}),
    ]

    @classmethod
    def setUpClass(cls):
        cls.bundle = SkillBundle()

    @classmethod
    def tearDownClass(cls):
        cls.bundle.cleanup()

    def test_each_tool_answers_as_agent_tools_does(self):
        from rel_quest import agent_tools
        from rel_quest.core import LigandDatabase

        db = LigandDatabase()
        direct = {
            "get_database_statistics": lambda: agent_tools.get_database_statistics(db),
            "search_ligands_by_name": lambda: agent_tools.search_ligands_by_name(db, "BTBP"),
            "get_ligand_details": lambda: agent_tools.get_ligand_details(db, "TODGA"),
            "get_ligand_conditions": lambda: agent_tools.get_ligand_conditions(db, "TBP", "Eu"),
            "search_by_trend": lambda: agent_tools.search_by_trend(db, "Peaked"),
            "search_by_mechanism": lambda: agent_tools.search_by_mechanism(db, "Solvation"),
            "search_by_subtype": lambda: agent_tools.search_by_extractant_type(db, "DGA"),
            "get_best_separators": lambda: agent_tools.get_best_separators(db, "Eu", "La"),
            "compare_ligands": lambda: agent_tools.compare_ligands(db, ["TODGA", "T2EHDGA"]),
            "get_ligands_for_element": lambda: agent_tools.get_ligands_for_element(db, "Y"),
        }
        for tool, args in self.CALLS:
            proc = self.bundle.query(tool, json.dumps(args))
            self.assertEqual(proc.returncode, 0, f"{tool}: {proc.stderr[-1000:]}")
            got = json.loads(proc.stdout)  # stdout is JSON and nothing else
            self.assertEqual(got, json.loads(json.dumps(direct[tool](), default=str)), tool)

    def test_vocabulary_lists_the_valid_values(self):
        proc = self.bundle.query("vocab")
        self.assertEqual(proc.returncode, 0, proc.stderr[-1000:])
        vocab = json.loads(proc.stdout)
        self.assertEqual(set(vocab["trends"]), TRENDS)
        self.assertIn("Y", vocab["elements"])
        self.assertIn("DGA", vocab["subtypes"])

    def test_errors_exit_nonzero_with_the_reason_on_stderr(self):
        unknown = self.bundle.query("get_ligand_details", '{"name": "No-Such-Ligand-XYZ"}')
        self.assertEqual(unknown.returncode, 1)
        self.assertIn("No-Such-Ligand-XYZ", unknown.stderr)
        self.assertEqual(unknown.stdout, "")

        bad = self.bundle.query("search_by_trend", '{"trend": "Upward"}')
        self.assertEqual(bad.returncode, 2)
        self.assertIn("Peaked", bad.stderr)  # names the valid values

        tool = self.bundle.query("plot_selectivity_trend", '{"ligand_name": "TODGA"}')
        self.assertEqual(tool.returncode, 2)


@unittest.skipUnless(os.environ.get("RELQUEST_SKILL_LIVE") == "1" and shutil.which("claude"),
                     "set RELQUEST_SKILL_LIVE=1 (runs `claude -p` three times)")
class TestSkillInClaudeCode(unittest.TestCase):
    """K6: the skill as a project skill of an empty project, asked three questions."""

    @classmethod
    def setUpClass(cls):
        cls.records = raw_records()
        cls.tmp = tempfile.TemporaryDirectory()
        project = Path(cls.tmp.name) / "project"
        build_skill(project / ".claude" / "skills" / "reel-quest")
        cls.sessions = {}
        for key, prompt in cls.QUESTIONS.items():
            proc = subprocess.run(
                ["claude", "-p", prompt, "--output-format", "stream-json", "--verbose",
                 "--allowedTools", "Skill", "Read", "Glob", "Grep", "Bash(python:*)", "Bash(python3:*)"],
                cwd=project, capture_output=True, text=True, timeout=900)
            events = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
            cls.sessions[key] = {"prompt": prompt, "exit": proc.returncode, "events": events,
                                 "stderr": proc.stderr[-2000:]}
        LIVE_TRANSCRIPT.parent.mkdir(exist_ok=True)
        LIVE_TRANSCRIPT.write_text(json.dumps(cls.sessions, indent=1, ensure_ascii=False, default=str))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    QUESTIONS = {
        "fixed": "Using the reel-quest skill: which ligand has the largest median Eu/La "
                 "separation factor, and what is its median log SF?",
        "custom": "Using the reel-quest skill: for TODGA in nitric acid at 3 M or more, what is "
                  "the median log separation factor of Eu over La, and over how many paired "
                  "conditions? Give the number to three decimals.",
        "unnamed": "How many extractants have distribution ratios measured for yttrium in "
                   "rare-earth solvent extraction? Name them.",
    }

    def answer(self, key):
        session = self.sessions[key]
        self.assertEqual(session["exit"], 0, session["stderr"])
        results = [e for e in session["events"] if e.get("type") == "result"]
        return results[-1].get("result", "")

    def used_skill(self, key):
        for event in self.sessions[key]["events"]:
            for block in (event.get("message") or {}).get("content") or []:
                if (isinstance(block, dict) and block.get("type") == "tool_use"
                        and block.get("name") == "Skill"
                        and "reel-quest" in json.dumps(block.get("input"))):
                    return True
        return False

    def test_fixed_question(self):
        self.assertTrue(self.used_skill("fixed"))
        answer = self.answer("fixed")
        self.assertIn("DBPhthA", answer)
        self.assertIn("4.3", answer)

    def test_custom_analysis(self):
        by_name = {r["ligand"]["Name"]: r for r in self.records}
        pairs = paired_log_sf(by_name["TODGA"], "Eu", "La", strong_nitric)
        median = statistics.median(v for v, _ in pairs)
        self.assertTrue(self.used_skill("custom"))
        answer = self.answer("custom")
        self.assertIn(f"{median:.3f}", answer)
        self.assertRegex(answer, rf"\b{len(pairs)}\b")

    def test_skill_triggers_unnamed(self):
        self.assertTrue(self.used_skill("unnamed"))
        expected = [r["ligand"]["Name"] for r in self.records
                    if any("Y" in c["distribution_ratios"] and measured(c["distribution_ratios"]["Y"])
                           for c in r["experimental_data"])]
        answer = self.answer("unnamed")
        for name in expected:
            self.assertIn(name, answer)


if __name__ == "__main__":
    unittest.main()
