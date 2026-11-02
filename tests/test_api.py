"""Tests for the REST API (rel_quest/api.py).

Counts are asserted against the live LigandDatabase rather than hardcoded,
so the tests keep passing when records are added — what they pin is that the
API agrees with the Python API, not a snapshot of the data.

Two regressions are pinned deliberately:
  - the trend filter accepts all 7 vocabulary values (Peaked and U_Shaped
    are rejected by the CLI's --trend; the API must not repeat that);
  - the element filter reads condition-level measurements, so Y and Sc are
    queryable even though selectivity.element_coverage omits them
    (PROJECT_REVIEW_2026-07-31.md, finding 2).

Mounted under a prefix (the Hugging Face Space serves it at /api of its Gradio
app, space/app.py), written before the code:
  M1. The mounted app answers differently from the standalone one.
  M2. Its documentation does not load: the Swagger page asks for an OpenAPI
      URL without the mount prefix, or the root redirect leaves the mount.
  M3. It loads a second copy of the database instead of the one it is given.
"""

import unittest

try:
    from fastapi.testclient import TestClient
    HAVE_FASTAPI = True
except ImportError:  # pragma: no cover
    HAVE_FASTAPI = False

TRENDS = ["Ascending", "Reverse", "Flat", "Mixed", "Insufficient_Data",
          "Peaked", "U_Shaped"]


@unittest.skipUnless(HAVE_FASTAPI, "fastapi not installed")
class TestRestApi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from rel_quest import LigandDatabase
        from rel_quest.api import create_app

        cls.db = LigandDatabase()
        cls.client = TestClient(create_app())

    # -- infrastructure -----------------------------------------------------

    def test_health_reports_live_counts(self):
        response = self.client.get("/api/v1/health")
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body["status"], "ok")
        self.assertEqual(body["ligands"], len(self.db.get_all_ligands()))
        self.assertEqual(
            body["measurements"],
            sum(len(l.experimental_conditions) for l in self.db.get_all_ligands()),
        )

    def test_openapi_spec_is_served(self):
        response = self.client.get("/api/v1/openapi.json")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["info"]["title"], "REEL-QUEST API")

    def test_root_redirects_to_docs(self):
        response = self.client.get("/", follow_redirects=False)
        self.assertIn(response.status_code, (302, 307))
        self.assertEqual(response.headers["location"], "/api/docs")

    def test_stats_matches_python_api(self):
        response = self.client.get("/api/v1/stats")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["total_ligands"],
                         self.db.get_statistics()["total_ligands"])

    # -- listing and filters ------------------------------------------------

    def test_list_pagination_is_consistent(self):
        total = len(self.db.get_all_ligands())
        first = self.client.get("/api/v1/ligands",
                                params={"limit": 10, "offset": 0}).json()
        second = self.client.get("/api/v1/ligands",
                                 params={"limit": 10, "offset": 10}).json()
        self.assertEqual(first["total"], total)
        self.assertEqual(first["count"], 10)
        names_first = {r["name"] for r in first["results"]}
        names_second = {r["name"] for r in second["results"]}
        self.assertEqual(names_first & names_second, set(),
                         "pages at different offsets must not overlap")

    def test_subtype_filter_matches_python_api(self):
        for subtype in ("DGA", "Malonamide", "Phosphine_oxide"):
            with self.subTest(subtype=subtype):
                body = self.client.get("/api/v1/ligands",
                                       params={"subtype": subtype}).json()
                self.assertEqual(body["total"],
                                 len(self.db.filter_by_subtype(subtype)))

    def test_all_seven_trend_values_are_accepted(self):
        """The CLI rejects Peaked/U_Shaped; the API must accept all 7."""
        for trend in TRENDS:
            with self.subTest(trend=trend):
                response = self.client.get("/api/v1/ligands",
                                           params={"trend": trend})
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.json()["total"],
                                 len(self.db.filter_by_trend(trend)))

    def test_element_filter_reads_conditions_not_coverage(self):
        """Y and Sc must be queryable (element_coverage omits them)."""
        for element in ("Y", "Sc"):
            with self.subTest(element=element):
                body = self.client.get("/api/v1/ligands",
                                       params={"elements": element}).json()
                self.assertGreaterEqual(body["total"], 1)
                self.assertIn("TODGA", {r["name"] for r in body["results"]})

    def test_element_filter_requires_all_elements(self):
        both = self.client.get("/api/v1/ligands",
                               params={"elements": "Eu,La"}).json()["total"]
        eu_only = self.client.get("/api/v1/ligands",
                                  params={"elements": "Eu"}).json()["total"]
        self.assertGreater(both, 0)
        self.assertLessEqual(both, eu_only)

    def test_logd_range_filter(self):
        body = self.client.get("/api/v1/ligands",
                               params={"logd_min": 2.0}).json()
        self.assertGreater(body["total"], 0)
        # Every returned ligand must actually have a measurement >= 2.0.
        for result in body["results"]:
            ligand = self.db.query_by_name(result["name"])[0]
            values = [
                data.get("log_value")
                for condition in ligand.experimental_conditions
                for data in (condition.get("distribution_ratios") or {}).values()
                if isinstance(data, dict) and data.get("log_value") is not None
            ]
            self.assertTrue(any(v >= 2.0 for v in values),
                            f"{result['name']} has no log D >= 2.0")

    def test_invalid_filter_value_returns_422_with_allowed_list(self):
        response = self.client.get("/api/v1/ligands",
                                   params={"subtype": "NotASubtype"})
        self.assertEqual(response.status_code, 422)
        detail = response.json()["detail"]
        self.assertEqual(detail["param"], "subtype")
        self.assertIn("DGA", detail["allowed"])

    def test_vocab_counts_sum_to_total(self):
        body = self.client.get("/api/v1/vocab").json()
        total = len(self.db.get_all_ligands())
        self.assertEqual(sum(body["trends"].values()), total)
        self.assertEqual(sum(body["mechanisms"].values()), total)

    # -- detail lookup ------------------------------------------------------

    def test_detail_by_exact_name(self):
        body = self.client.get("/api/v1/ligands/TODGA").json()
        self.assertEqual(body["name"], "TODGA")
        self.assertEqual(body["matched_by"], "name")
        self.assertEqual(body["subtype"], "DGA")
        self.assertGreater(len(body["experimental_data"]), 0)

    def test_detail_resolves_alternative_names(self):
        """Old names from renames/merges must still resolve (CLAUDE.md promise)."""
        with_alternates = next(
            l for l in self.db.get_all_ligands()
            if l.alternative_names
            and not self.db.query_by_name(l.alternative_names[0])
        )
        alternate = with_alternates.alternative_names[0]
        response = self.client.get(f"/api/v1/ligands/{alternate}")
        self.assertEqual(response.status_code, 200, f"{alternate!r} did not resolve")
        body = response.json()
        self.assertEqual(body["name"], with_alternates.name)
        self.assertIn(body["matched_by"], ("alternative_name", "fuzzy_name"))

    def test_detail_unknown_name_is_404_with_suggestions(self):
        response = self.client.get("/api/v1/ligands/definitely-not-a-ligand")
        self.assertEqual(response.status_code, 404)
        self.assertIn("suggestions", response.json()["detail"])

    def test_detail_elements_measured_in_periodic_order(self):
        body = self.client.get("/api/v1/ligands/TODGA").json()
        measured = body["elements_measured"]
        self.assertIn("Y", measured)
        self.assertIn("Sc", measured)
        self.assertLess(measured.index("Sc"), measured.index("Y"),
                        "Sc (Z=21) must sort before Y (Z=39)")
        self.assertLess(measured.index("Y"), measured.index("La"),
                        "Y (Z=39) must sort before La (Z=57)")


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(HAVE_FASTAPI, "fastapi not installed")
class TestRestApiMounted(unittest.TestCase):
    """M1-M3: create_app(mount="") inside another app at /api, as the Space mounts it."""

    PATHS = ["/api/v1/health", "/api/v1/stats", "/api/v1/vocab",
             "/api/v1/ligands?subtype=DGA&elements=Eu,La&limit=200",
             "/api/v1/ligands?trend=Peaked", "/api/v1/ligands/TODGA",
             "/api/v1/ligands/no-such-ligand", "/api/v1/ligands?trend=Wobbly"]

    @classmethod
    def setUpClass(cls):
        from starlette.applications import Starlette
        from starlette.routing import Mount
        from rel_quest import LigandDatabase
        from rel_quest.api import create_app

        cls.db = LigandDatabase()
        cls.inner = create_app(db=cls.db, mount="")
        cls.mounted = TestClient(Starlette(routes=[Mount("/api", app=cls.inner)]))
        cls.standalone = TestClient(create_app(db=cls.db))

    def test_same_answers_as_the_standalone_app(self):              # M1
        for path in self.PATHS:
            with self.subTest(path):
                a, b = self.mounted.get(path), self.standalone.get(path)
                self.assertEqual(a.status_code, b.status_code)
                self.assertEqual(a.json(), b.json())

    def test_docs_and_redirect_stay_under_the_mount(self):          # M2
        docs = self.mounted.get("/api/docs")
        self.assertEqual(docs.status_code, 200)
        self.assertIn("/api/v1/openapi.json", docs.text)
        spec = self.mounted.get("/api/v1/openapi.json")
        self.assertEqual(spec.status_code, 200)
        self.assertIn("/v1/ligands", " ".join(spec.json()["paths"]))
        self.assertEqual(self.mounted.get("/api/redoc").status_code, 200)
        redirect = self.mounted.get("/api/", follow_redirects=False)
        self.assertIn(redirect.status_code, (302, 307))
        self.assertEqual(redirect.headers["location"], "/api/docs")

    def test_uses_the_database_it_is_given(self):                   # M3
        self.assertIs(self.inner.state.db, self.db)

