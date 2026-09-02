#!/usr/bin/env python3
"""Build a family-neutral orthogroup index from an OrthoFinder table."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def load_orthogroups(path: Path) -> dict[str, list[str]]:
    """Return ``orthogroup -> members`` from OrthoFinder's Orthogroups.tsv."""
    groups: dict[str, list[str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "Orthogroup":
            raise ValueError("The first Orthogroups.tsv column must be 'Orthogroup'.")
        for row_number, row in enumerate(reader, start=2):
            if not row or not row[0].strip():
                continue
            group = row[0].strip()
            if group in groups:
                raise ValueError(f"Duplicate orthogroup {group!r} on row {row_number}.")
            members: list[str] = []
            for cell in row[1:]:
                members.extend(item.strip() for item in cell.split(",") if item.strip())
            groups[group] = sorted(set(members))
    return groups


def load_characterized_ids(database_path: Path) -> set[str]:
    with database_path.open(encoding="utf-8") as handle:
        database = json.load(handle)
    return set(database.get("enzymes", {}))


def build_index(groups: dict[str, list[str]], characterized: set[str]) -> dict:
    member_to_group: dict[str, str] = {}
    records: dict[str, dict] = {}
    for group, members in sorted(groups.items()):
        for member in members:
            previous = member_to_group.setdefault(member, group)
            if previous != group:
                raise ValueError(
                    f"Sequence {member!r} occurs in both {previous!r} and {group!r}."
                )
        records[group] = {
            "members": members,
            "characterized_members": sorted(characterized.intersection(members)),
        }
    return {"orthogroups": records, "member_to_orthogroup": member_to_group}


def annotate_database(database_path: Path, output_path: Path, mapping: dict[str, str]) -> None:
    with database_path.open(encoding="utf-8") as handle:
        database = json.load(handle)
    for enzyme_id, record in database.get("enzymes", {}).items():
        record.setdefault("orthology_data", {})["orthogroup"] = mapping.get(enzyme_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(database, handle, indent=2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Index OrthoFinder orthogroups and mark characterized members."
    )
    parser.add_argument("--orthogroups", type=Path, required=True,
                        help="OrthoFinder Orthogroups.tsv file.")
    parser.add_argument("--database", type=Path, required=True,
                        help="FuncZymeDB JSON containing characterized enzymes.")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output orthology JSON.")
    parser.add_argument("--annotated-database", type=Path,
                        help="Optional copy of the database annotated with orthogroups.")
    parser.add_argument("--person", required=True,
                        help="Name or initials of the person running the build.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    groups = load_orthogroups(args.orthogroups)
    characterized = load_characterized_ids(args.database)
    result = build_index(groups, characterized)
    result["log"] = {
        "build_orthology_database": {
            "time": datetime.now(timezone.utc).isoformat(),
            "person": args.person,
            "orthogroups": str(args.orthogroups),
            "database": str(args.database),
            "n_orthogroups": len(groups),
            "n_members": len(result["member_to_orthogroup"]),
            "n_characterized_members": sum(
                len(group["characterized_members"])
                for group in result["orthogroups"].values()
            ),
        }
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    if args.annotated_database:
        annotate_database(
            args.database,
            args.annotated_database,
            result["member_to_orthogroup"],
        )


if __name__ == "__main__":
    main()
