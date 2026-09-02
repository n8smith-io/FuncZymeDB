"""Shared I/O and validation for the family-neutral FuncPred-AI workflow."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np


def load_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(parts)
                identifier = line[1:].split()[0]
                if not identifier or identifier in records:
                    raise ValueError(f"Missing or duplicate FASTA identifier: {identifier!r}")
                parts = []
            elif identifier is None:
                raise ValueError("Sequence data precedes the first FASTA header.")
            else:
                parts.append(line)
    if identifier is not None:
        records[identifier] = "".join(parts)
    return records


def save_embeddings(path: Path, sequence_ids: list[str], matrix: np.ndarray, model: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        sequence_ids=np.asarray(sequence_ids, dtype=str),
        embeddings=np.asarray(matrix, dtype=np.float32),
        model=np.asarray(model),
    )


def load_embeddings(path: Path) -> tuple[list[str], np.ndarray, str]:
    with np.load(path, allow_pickle=False) as payload:
        ids = payload["sequence_ids"].astype(str).tolist()
        matrix = np.asarray(payload["embeddings"], dtype=np.float32)
        model = str(payload["model"].item())
    if matrix.ndim != 2 or matrix.shape[0] != len(ids):
        raise ValueError("Embedding archive must contain an n_sequences by n_features matrix.")
    if len(ids) != len(set(ids)):
        raise ValueError("Embedding sequence identifiers must be unique.")
    return ids, matrix, model


def load_labels(path: Path) -> tuple[dict[tuple[str, str], dict[str, int]], dict[str, str]]:
    tasks: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    groups: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "task", "label", "value", "evaluation_group"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Labels must contain columns: {sorted(required)}")
        for row_number, row in enumerate(reader, start=2):
            sequence_id = row["sequence_id"].strip()
            key = (row["task"].strip(), row["label"].strip())
            if not sequence_id or not all(key):
                raise ValueError(f"Blank sequence_id, task, or label on row {row_number}.")
            value = int(row["value"])
            if value not in (0, 1):
                raise ValueError(f"Label value on row {row_number} is not 0 or 1.")
            if sequence_id in tasks[key]:
                raise ValueError(f"Duplicate label observation on row {row_number}.")
            tasks[key][sequence_id] = value
            group = row["evaluation_group"].strip() or sequence_id
            previous = groups.setdefault(sequence_id, group)
            if previous != group:
                raise ValueError(f"Conflicting evaluation groups for {sequence_id!r}.")
    return dict(tasks), groups


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
