#!/usr/bin/env python3
"""Transfer functional labels from characterized members of an orthogroup."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


SEARCH_COLUMNS = ("query", "target", "percent_identity", "alignment_length", "evalue", "bits")


def run_mmseqs(query: Path, reference: Path, output: Path, threads: int) -> None:
    with tempfile.TemporaryDirectory(prefix="funczyme-mmseqs-") as temp_dir:
        command = [
            "mmseqs", "easy-search", str(query), str(reference), str(output), temp_dir,
            "--threads", str(threads), "--format-output", "query,target,pident,alnlen,evalue,bits",
        ]
        subprocess.run(command, check=True)


def best_hits(path: Path) -> dict[str, dict[str, str]]:
    hits: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, fieldnames=SEARCH_COLUMNS, delimiter="\t")
        for row in reader:
            query = row["query"]
            score = float(row["bits"])
            if query not in hits or score > float(hits[query]["bits"]):
                hits[query] = row
    return hits


def load_labels(path: Path) -> dict[str, dict[tuple[str, str], int]]:
    labels: dict[str, dict[tuple[str, str], int]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "task", "label", "value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Labels must contain columns: {sorted(required)}")
        for row in reader:
            value = int(row["value"])
            if value not in (0, 1):
                raise ValueError("Label values must be 0 or 1.")
            labels[row["sequence_id"]][(row["task"], row["label"])] = value
    return labels


def transfer_predictions(
    hits: dict[str, dict[str, str]], orthology: dict, labels: dict, min_support: float,
    query_ids: list[str] | None = None,
) -> list[dict]:
    mapping = orthology["member_to_orthogroup"]
    groups = orthology["orthogroups"]
    rows: list[dict] = []
    for query in sorted(query_ids if query_ids is not None else hits):
        hit = hits.get(query)
        if hit is None:
            rows.append({
                "sequence_id": query, "orthogroup": "", "task": "", "label": "",
                "score": "", "predicted": "", "n_evidence": 0,
                "best_reference": "", "percent_identity": "", "status": "no_reference_hit",
            })
            continue
        group_id = mapping.get(hit["target"])
        members = groups.get(group_id, {}).get("characterized_members", [])
        observations: dict[tuple[str, str], list[int]] = defaultdict(list)
        for member in members:
            for key, value in labels.get(member, {}).items():
                observations[key].append(value)
        if not observations:
            rows.append({
                "sequence_id": query, "orthogroup": group_id or "", "task": "",
                "label": "", "score": "", "predicted": "", "n_evidence": 0,
                "best_reference": hit["target"], "percent_identity": hit["percent_identity"],
                "status": "no_characterized_labels_in_orthogroup",
            })
            continue
        for (task, label), values in sorted(observations.items()):
            score = sum(values) / len(values)
            rows.append({
                "sequence_id": query, "orthogroup": group_id, "task": task,
                "label": label, "score": f"{score:.8g}",
                "predicted": int(score >= min_support), "n_evidence": len(values),
                "best_reference": hit["target"], "percent_identity": hit["percent_identity"],
                "status": "ok",
            })
    return rows


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                if identifier:
                    identifiers.append(identifier)
    return identifiers


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run FuncPred-OG by orthogroup placement and label transfer."
    )
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True,
                        help="Family reference FASTA using IDs indexed in the orthology JSON.")
    parser.add_argument("--orthology", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--search-results", type=Path,
                        help="Reuse an MMseqs result TSV instead of running a search.")
    parser.add_argument("--min-support", type=float, default=0.5)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 <= args.min_support <= 1:
        raise SystemExit("--min-support must be between 0 and 1.")
    search_path = args.search_results
    if search_path is None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        search_path = args.output.with_suffix(".mmseqs.tsv")
        run_mmseqs(args.query_fasta, args.reference_fasta, search_path, args.threads)
    with args.orthology.open(encoding="utf-8") as handle:
        orthology = json.load(handle)
    rows = transfer_predictions(
        best_hits(search_path), orthology, load_labels(args.labels), args.min_support,
        query_ids=fasta_ids(args.query_fasta),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sequence_id", "orthogroup", "task", "label", "score", "predicted",
        "n_evidence", "best_reference", "percent_identity", "status",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
