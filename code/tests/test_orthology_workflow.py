#!/usr/bin/env python3
"""Small contract tests for the family-neutral FuncPred-OG workflow."""

from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
ORTHOLOGY_DIR = REPO_ROOT / "code" / "orthology"
if str(ORTHOLOGY_DIR) not in sys.path:
    sys.path.insert(0, str(ORTHOLOGY_DIR))

from build_orthology_database import build_index, load_orthogroups  # noqa: E402
from predict_with_funcpred_og import transfer_predictions  # noqa: E402


def test_orthofinder_table_is_indexed_and_characterized_members_are_marked(tmp_path):
    table = tmp_path / "Orthogroups.tsv"
    table.write_text(
        "Orthogroup\tspecies_a\tspecies_b\n"
        "OG0001\tseq_a1, seq_a2\tseq_b1\n"
        "OG0002\tseq_a3\tseq_b2\n",
        encoding="utf-8",
    )

    groups = load_orthogroups(table)
    index = build_index(groups, {"seq_a1", "seq_b2", "absent"})

    assert index["member_to_orthogroup"]["seq_b1"] == "OG0001"
    assert index["orthogroups"]["OG0001"]["characterized_members"] == ["seq_a1"]
    assert index["orthogroups"]["OG0002"]["characterized_members"] == ["seq_b2"]


def test_funcpred_og_reports_support_and_abstention():
    orthology = {
        "member_to_orthogroup": {"ref_a": "OG1", "ref_b": "OG1", "ref_c": "OG2"},
        "orthogroups": {
            "OG1": {"characterized_members": ["ref_a", "ref_b"]},
            "OG2": {"characterized_members": []},
        },
    }
    labels = {
        "ref_a": {("substrate", "class_x"): 1},
        "ref_b": {("substrate", "class_x"): 0},
    }
    hits = {
        "query_1": {"target": "ref_a", "percent_identity": "72", "bits": "100"},
        "query_2": {"target": "ref_c", "percent_identity": "51", "bits": "90"},
    }

    rows = transfer_predictions(hits, orthology, labels, min_support=0.5)

    assert rows[0]["score"] == "0.5"
    assert rows[0]["predicted"] == 1
    assert rows[0]["n_evidence"] == 2
    assert rows[1]["status"] == "no_characterized_labels_in_orthogroup"
