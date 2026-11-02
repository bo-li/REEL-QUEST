"""Packaging metadata tests.

Dependencies are declared twice on purpose: pyproject.toml makes the project
installable, while requirements.txt stays a flat list because the Dockerfile
copies it on its own layer (`COPY requirements.txt` before the source) so the
dependency layer stays cached. These tests stop the two drifting apart.
"""

import re
import unittest
from pathlib import Path

try:  # stdlib from 3.11; tomli is present on 3.10 as a pytest dependency
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - depends on interpreter
    try:
        import tomli as tomllib
    except ModuleNotFoundError:  # pragma: no cover
        tomllib = None

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"
REQUIREMENTS = ROOT / "requirements.txt"


def parse_requirement(line):
    """Split 'package>=1.2.3' into ('package', '>=1.2.3')."""
    match = re.match(r"^([A-Za-z0-9._-]+)\s*(.*)$", line.strip())
    return (match.group(1).lower(), match.group(2).replace(" ", "")) if match else None


def read_requirements():
    entries = {}
    for line in REQUIREMENTS.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parsed = parse_requirement(line)
        if parsed:
            entries[parsed[0]] = parsed[1]
    return entries


class TestPackagingMetadata(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if tomllib is None:  # pragma: no cover
            raise unittest.SkipTest("no TOML parser available")
        cls.pyproject = tomllib.loads(PYPROJECT.read_text())

    def test_pyproject_exists_and_declares_the_package(self):
        project = self.pyproject["project"]
        self.assertEqual(project["name"], "reel-quest")
        self.assertIn("rel_quest", self.pyproject["tool"]["setuptools"]["packages"])

    def test_version_is_sourced_from_the_package(self):
        """Version must not be duplicated in pyproject."""
        import rel_quest

        self.assertIn("version", self.pyproject["project"]["dynamic"])
        attr = self.pyproject["tool"]["setuptools"]["dynamic"]["version"]["attr"]
        self.assertEqual(attr, "rel_quest.__version__")
        self.assertRegex(rel_quest.__version__, r"^\d+\.\d+")

    def test_dependencies_match_requirements_txt(self):
        """Every runtime/dev dependency is declared identically in both files."""
        declared = {}
        for spec in self.pyproject["project"]["dependencies"]:
            name, constraint = parse_requirement(spec)
            declared[name] = constraint
        for extra in self.pyproject["project"].get("optional-dependencies", {}).values():
            for spec in extra:
                name, constraint = parse_requirement(spec)
                declared[name] = constraint

        required = read_requirements()

        self.assertEqual(
            set(declared) - set(required), set(),
            "declared in pyproject.toml but missing from requirements.txt",
        )
        self.assertEqual(
            set(required) - set(declared), set(),
            "listed in requirements.txt but missing from pyproject.toml",
        )
        for name, constraint in sorted(required.items()):
            with self.subTest(package=name):
                self.assertEqual(
                    declared[name], constraint,
                    f"{name}: pyproject says '{declared[name]}', "
                    f"requirements.txt says '{constraint}'",
                )

    def test_pytest_config_lives_in_pyproject(self):
        """pytest.ini was folded in; two config sources would be ambiguous."""
        self.assertFalse(
            (ROOT / "pytest.ini").exists(),
            "pytest.ini is back — pytest config belongs in pyproject.toml",
        )
        options = self.pyproject["tool"]["pytest"]["ini_options"]
        self.assertEqual(options["testpaths"], ["tests"])
        self.assertIn("manual", options["norecursedirs"])

    def test_lock_pins_every_requirement(self):
        """requirements.lock (Dockerfile, locked CI job) must pin every
        package requirements.txt names, at a version its bound admits."""
        lock = {}
        for line in (ROOT / "requirements.lock").read_text().splitlines():
            match = re.match(r"^([A-Za-z0-9._-]+)==(\S+)", line.strip())
            if match:
                lock[match.group(1).lower().replace("_", "-")] = match.group(2)
        for name, constraint in sorted(read_requirements().items()):
            key = name.replace("_", "-")
            with self.subTest(package=name):
                self.assertIn(key, lock, f"{name} is not pinned in requirements.lock")
                bound = re.match(r">=([0-9][0-9.]*)", constraint)
                if bound:
                    pinned = tuple(int(x) for x in re.findall(r"\d+", lock[key])[:3])
                    lower = tuple(int(x) for x in bound.group(1).split(".")[:3])
                    self.assertGreaterEqual(pinned, lower)

    def test_requires_python_covers_the_ci_matrix(self):
        self.assertEqual(self.pyproject["project"]["requires-python"], ">=3.10")


if __name__ == "__main__":
    unittest.main()
