"""Checks on deploy.sh, the documented one-command launcher.

docs/DEPLOYMENT_CHECKLIST.md tells users to run `./deploy.sh`, and for a
long time that command could not work at all: the file carried CRLF line
endings, so the kernel looked for an interpreter literally named
"/bin/bash\r" and reported "bad interpreter"; git recorded mode 100644, so
it was not executable even once the endings were repaired; and its
preflight checks still required the removed `rex_db` package and the v6
database (PROJECT_REVIEW_2026-07-31.md finding 3).

Nothing in the suite exercised it, because running it means building a
virtualenv and installing rdkit. These tests instead check the properties
whose absence made it unrunnable — cheaply, so they run on every commit.
"""

import re
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEPLOY = ROOT / "deploy.sh"


class TestDeployScript(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not DEPLOY.exists():  # pragma: no cover
            raise unittest.SkipTest("deploy.sh not present")
        cls.raw = DEPLOY.read_bytes()
        cls.text = cls.raw.decode()

    def test_has_no_crlf_line_endings(self):
        """A CRLF shebang makes the script unrunnable, not merely untidy."""
        self.assertNotIn(b"\r\n", self.raw,
                         "deploy.sh contains CRLF; the shebang will fail with "
                         "'/bin/bash^M: bad interpreter'")

    def test_line_endings_are_pinned_by_gitattributes(self):
        """Pinning is what stops a Windows checkout reintroducing CRLF."""
        gitattributes = ROOT / ".gitattributes"
        self.assertTrue(gitattributes.exists(), ".gitattributes is missing")
        self.assertRegex(gitattributes.read_text(), r"\*\.sh\s+text\s+eol=lf")

    def test_is_executable_in_the_git_index(self):
        """chmod +x locally is not enough; git must record mode 100755."""
        result = subprocess.run(
            ["git", "ls-files", "-s", "deploy.sh"],
            cwd=ROOT, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout.strip():  # pragma: no cover
            self.skipTest("deploy.sh is not tracked by git here")
        mode = result.stdout.split()[0]
        self.assertEqual(mode, "100755",
                         f"git records mode {mode}; ./deploy.sh will not run")

    def test_parses_as_valid_bash(self):
        result = subprocess.run(["bash", "-n", str(DEPLOY)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0,
                         f"bash -n failed: {result.stderr}")

    def test_help_exits_cleanly_without_side_effects(self):
        result = subprocess.run(["bash", str(DEPLOY), "--help"],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--port", result.stdout)

    def test_rejects_an_invalid_port(self):
        result = subprocess.run(["bash", str(DEPLOY), "--port", "notanumber"],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 1)
        self.assertIn("--port must be a number", result.stdout + result.stderr)

    def test_rejects_an_unknown_option(self):
        result = subprocess.run(["bash", str(DEPLOY), "--bogus"],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 1)

    def test_references_no_retired_paths(self):
        """The preflight must not require things that were deleted."""
        retired = ["rex_db", "sx_ligands_merged_v6", "test_ux_improvements"]
        for name in retired:
            with self.subTest(name=name):
                self.assertNotIn(name, self.text,
                                 f"deploy.sh still references retired '{name}'")

    def test_every_required_file_it_checks_actually_exists(self):
        """The REQUIRED_FILES list is the preflight; it must be truthful."""
        block = re.search(r"REQUIRED_FILES=\((.*?)\)", self.text, re.S)
        self.assertIsNotNone(block, "REQUIRED_FILES array not found")
        paths = re.findall(r'"([^"]+)"', block.group(1))
        self.assertGreater(len(paths), 3, "suspiciously short preflight list")
        for path in paths:
            with self.subTest(path=path):
                self.assertTrue((ROOT / path).exists(),
                                f"deploy.sh requires '{path}', which does not exist")

    def test_does_not_branch_on_a_dead_status_check(self):
        """`set -e` aborts first, so `if [ $? -eq 0 ]` can never be false.

        The old script had two such blocks, both purporting to report an
        installation or database-load failure that could not reach them.
        """
        self.assertIn("set -e", self.text)
        self.assertNotRegex(
            self.text, r"if\s+\[\s+\$\?\s+-eq\s+0\s+\]",
            "unreachable $? branch under `set -e`; use `if ! cmd; then`")


class TestDeploymentDocs(unittest.TestCase):
    """The checklist is the entry point users read before running anything."""

    def test_checklist_does_not_teach_retired_paths(self):
        checklist = ROOT / "docs" / "DEPLOYMENT_CHECKLIST.md"
        if not checklist.exists():  # pragma: no cover
            self.skipTest("checklist not present")
        text = checklist.read_text()
        for name in ("rex_db", "sx_ligands_merged_v6", "test_ux_improvements"):
            with self.subTest(name=name):
                self.assertNotIn(name, text,
                                 f"checklist still instructs users to copy '{name}'")


if __name__ == "__main__":
    unittest.main()
