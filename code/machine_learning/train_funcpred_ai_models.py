#!/usr/bin/env python3
"""Train one reproducible binary classifier per family-defined substrate label."""

from __future__ import annotations

import argparse
import csv
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from workflow_utils import load_embeddings, load_labels, sha256


def fit_models(
    sequence_ids: list[str], matrix: np.ndarray, tasks: dict, groups: dict[str, str],
    folds: int, seed: int,
) -> tuple[dict, list[dict], list[dict]]:
    index = {sequence_id: row for row, sequence_id in enumerate(sequence_ids)}
    models: dict[tuple[str, str], LogisticRegression] = {}
    metrics: list[dict] = []
    predictions: list[dict] = []
    for key, observed in sorted(tasks.items()):
        ids = [sequence_id for sequence_id in sequence_ids if sequence_id in observed]
        y = np.asarray([observed[sequence_id] for sequence_id in ids], dtype=int)
        x = matrix[[index[sequence_id] for sequence_id in ids]]
        group_values = np.asarray([groups[sequence_id] for sequence_id in ids])
        if len(np.unique(y)) < 2:
            raise ValueError(f"Task {key!r} has only one observed class.")
        unique_groups = np.unique(group_values)
        n_splits = min(folds, len(unique_groups))
        if n_splits < 2:
            raise ValueError(f"Task {key!r} requires at least two evaluation groups.")
        out_of_fold = np.full(len(ids), np.nan)
        for train_rows, test_rows in GroupKFold(n_splits=n_splits).split(x, y, group_values):
            if len(np.unique(y[train_rows])) < 2:
                raise ValueError(
                    f"Task {key!r} has a cross-validation fold with one training class; "
                    "merge groups or collect more observations."
                )
            fold_model = LogisticRegression(
                class_weight="balanced", max_iter=2000, random_state=seed
            ).fit(x[train_rows], y[train_rows])
            out_of_fold[test_rows] = fold_model.predict_proba(x[test_rows])[:, 1]
        model = LogisticRegression(
            class_weight="balanced", max_iter=2000, random_state=seed
        ).fit(x, y)
        models[key] = model
        metrics.append({
            "task": key[0], "label": key[1], "n_observations": len(ids),
            "n_positive": int(y.sum()), "n_groups": len(unique_groups),
            "n_splits": n_splits,
            "average_precision": float(average_precision_score(y, out_of_fold)),
            "roc_auc": float(roc_auc_score(y, out_of_fold)),
        })
        predictions.extend(
            {"sequence_id": sequence_id, "task": key[0], "label": key[1],
             "value": int(value), "out_of_fold_score": f"{score:.10g}",
             "evaluation_group": groups[sequence_id]}
            for sequence_id, value, score in zip(ids, y, out_of_fold)
        )
    return models, metrics, predictions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train group-aware FuncPred-AI models.")
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.folds < 2:
        raise SystemExit("--folds must be at least 2.")
    sequence_ids, matrix, embedding_model = load_embeddings(args.embeddings)
    tasks, groups = load_labels(args.labels)
    models, metrics, predictions = fit_models(
        sequence_ids, matrix, tasks, groups, args.folds, args.seed
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bundle = {
        "schema_version": 1,
        "model_type": "funcpred_ai_binary_logistic",
        "embedding_model": embedding_model,
        "feature_count": matrix.shape[1],
        "threshold": 0.5,
        "models": models,
    }
    joblib.dump(bundle, args.output_dir / "funcpred_ai_models.joblib")
    with (args.output_dir / "evaluation.json").open("w", encoding="utf-8") as handle:
        json.dump({"metrics": metrics}, handle, indent=2)
    with (args.output_dir / "out_of_fold_predictions.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        fields = ["sequence_id", "task", "label", "value", "out_of_fold_score", "evaluation_group"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(predictions)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "scikit_learn": sklearn.__version__,
        "embedding_model": embedding_model, "folds": args.folds, "seed": args.seed,
        "embedding_sha256": sha256(args.embeddings), "labels_sha256": sha256(args.labels),
    }
    with (args.output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
