#!/usr/bin/env python3
"""Focused invariants for the curated database build."""

from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DATABASE_DIR = REPO_ROOT / "code" / "database"
if str(DATABASE_DIR) not in sys.path:
    sys.path.insert(0, str(DATABASE_DIR))

import build_enz_cpd_seq_database as build  # noqa: E402


def test_bare_and_partially_annotated_substrates_are_retained():
    compounds = build.parse_compound_str(
        "cellulose; unnamed polymer | CID:NA",
        role="acceptor",
        index="1",
    )
    by_name = {compound.canonical_name: compound for compound in compounds}

    assert set(by_name) == {"cellulose", "unnamed polymer"}
    assert all(compound.is_acceptor for compound in by_name.values())
    assert by_name["cellulose"].pubchem_id == build.NA


def test_duplicate_activity_indices_report_all_source_rows():
    records = [
        {"INDEX": "10001"},
        {"INDEX": "10002"},
        {"INDEX": "10001"},
        {"INDEX": ""},
    ]

    assert build.duplicate_activity_indices(records) == {
        "10001": [1, 3],
        build.NA: [4],
    }
