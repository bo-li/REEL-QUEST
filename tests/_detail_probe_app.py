"""Streamlit probe app: renders ONLY the ligand detail view, for one row.

The detail view is reachable in the real dashboard only by selecting a row in
the browser table, which AppTest cannot drive (st.dataframe selection state is
not settable). This probe calls the render function directly with a stub
selection event so the ~440-line detail path can be exercised in CI.

Driven by test_dashboard_smoke.py; the row index comes from
st.session_state["probe_row"]. Not collected by pytest (filename is not test_*).
"""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import app  # noqa: E402  (path set up above)


class _Selection:
    """Stands in for st.dataframe's selection object."""

    def __init__(self, rows):
        self.rows = rows
        self.columns = []


class _Event:
    def __init__(self, rows):
        self.selection = _Selection(rows)


db = app.load_database()
ligands = db.get_all_ligands()

row = int(st.session_state.get("probe_row", 0))

if row < 0:
    # Negative row means "no row selected" — exercises the placeholder branch.
    app._render_ligand_detail(ligands, _Event([]))
else:
    row = min(row, len(ligands) - 1)
    app._render_ligand_detail(ligands, _Event([row]))
