#!/usr/bin/env python3
"""Integration checks for the family-neutral starter tables and builder."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DATABASE_DIR = REPO_ROOT / "code" / "database"
if str(DATABASE_DIR) not in sys.path:
    sys.path.insert(0, str(DATABASE_DIR))

import build_funczyme_database as build  # noqa: E402


def read_header(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return set(next(csv.reader(handle, delimiter="\t")))


def test_starter_tables_expose_fields_consumed_by_the_builder():
    activity_fields = read_header(REPO_ROOT / "data/templates/activity_data.tsv")
    compound_fields = read_header(REPO_ROOT / "data/templates/compound_data.tsv")

    assert {
        "INDEX",
        "SPECIES",
        "ENZYME_COMMON_NAME",
        "SUBSTRATE",
        "MANUAL_SEQUENCE",
        "FAMILY_MEMBERSHIP_REFERENCE",
        "FAMILY_MEMBERSHIP_CALL",
    } <= activity_fields
    assert {
        "ROW_ID",
        "COMPOUND_NAME",
        "CID",
        "INCHIKEY",
        "CHEBI",
        "SMILES",
    } <= compound_fields


def test_empty_starter_tables_build_an_empty_database(tmp_path: Path):
    output_db = tmp_path / "funczymedb.json"
    output_fasta = tmp_path / "sequences.fa"
    command = [
        sys.executable,
        str(DATABASE_DIR / "build_funczyme_database.py"),
        "--compound_file",
        str(REPO_ROOT / "data/templates/compound_data.tsv"),
        "--activity_file",
        str(REPO_ROOT / "data/templates/activity_data.tsv"),
        "--output_db",
        str(output_db),
        "--output_fasta",
        str(output_fasta),
        "--person",
        "test",
        "--skip_kingdom_check",
    ]

    subprocess.run(command, check=True, capture_output=True, text=True)

    database = json.loads(output_db.read_text(encoding="utf-8"))
    assert database["enzymes"] == {}
    assert database["compounds"] == {}
    assert database["reactions"] == {}
    assert database["log"]["build_database"]["person"] == "test"
    assert output_fasta.read_text(encoding="utf-8") == ""


def test_curation_and_membership_provenance_are_retained():
    compound = build.Compound("example substrate")
    compound_registry = {compound.canonical_name: compound}
    name_lookup = {compound.canonical_name: compound}
    enzymes = {}
    nonmatches = {role: set() for role in ("SUBSTRATE", "DONOR", "ACCEPTOR", "PRODUCT")}
    record = {
        "INDEX": "source-1",
        "SPECIES": "Genus species",
        "KINGDOM": "example clade",
        "ENZYME_COMMON_NAME": "example enzyme",
        "ENZYME_FULL_NAME": "example enzyme full name",
        "SUBSTRATE": "example substrate",
        "GENERAL_ENZYME_FAMILY": "example family",
        "FAMILY_MEMBERSHIP_REFERENCE": "example reference",
        "FAMILY_MEMBERSHIP_CALL": "accepted",
        "SOURCE_DATASET": "example dataset",
        "ACQUISITION_METHOD": "manual review",
        "CURATION_STATUS": "reviewed",
        "EVIDENCE_TYPE": "biochemical assay",
        "OTHER_COMMENTS": "example note",
        "CURATED_BY": "AB",
    }

    build.ingest_enzyme_record(
        record,
        enzymes,
        compound_registry,
        name_lookup,
        nonmatches,
    )

    enzyme = next(iter(enzymes.values()))
    reaction = next(iter(enzyme.reactions))
    assert enzyme.additional_metadata["GENERAL_ENZYME_FAMILY"] == "example family"
    assert enzyme.additional_metadata["FAMILY_MEMBERSHIP_REFERENCE"] == "example reference"
    assert enzyme.additional_metadata["FAMILY_MEMBERSHIP_CALL"] == "accepted"
    assert reaction.reaction_metadata["SOURCE_DATASET"] == "example dataset"
    assert reaction.reaction_metadata["CURATION_STATUS"] == "reviewed"
    assert reaction.reaction_metadata["CURATED_BY"] == "AB"
