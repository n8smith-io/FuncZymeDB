#!/usr/bin/env python3
"""Apply a trained FuncPred-AI model bundle to a compatible embedding archive."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import joblib

from workflow_utils import load_embeddings


def predict(bundle: dict, ids: list[str], matrix) -> list[dict]:
    if matrix.shape[1] != bundle["feature_count"]:
        raise ValueError("Embedding feature count does not match the model bundle.")
    threshold = float(bundle.get("threshold", 0.5))
    rows: list[dict] = []
    for (task, label), model in sorted(bundle["models"].items()):
        scores = model.predict_proba(matrix)[:, 1]
        rows.extend(
            {"sequence_id": sequence_id, "task": task, "label": label,
             "score": f"{score:.10g}", "threshold": threshold,
             "predicted": int(score >= threshold)}
            for sequence_id, score in zip(ids, scores)
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run FuncPred-AI predictions.")
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ids, matrix, embedding_model = load_embeddings(args.embeddings)
    bundle = joblib.load(args.models)
    if embedding_model != bundle["embedding_model"]:
        raise SystemExit(
            f"Embedding model {embedding_model!r} does not match bundle "
            f"{bundle['embedding_model']!r}."
        )
    rows = predict(bundle, ids, matrix)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sequence_id", "task", "label", "score", "threshold", "predicted"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
