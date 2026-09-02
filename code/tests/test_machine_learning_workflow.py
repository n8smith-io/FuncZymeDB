#!/usr/bin/env python3
"""Contract tests for FuncPred-AI training and prediction."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
ML_DIR = REPO_ROOT / "code" / "machine_learning"
if str(ML_DIR) not in sys.path:
    sys.path.insert(0, str(ML_DIR))

from predict_with_funcpred_ai import predict  # noqa: E402
from train_funcpred_ai_models import fit_models  # noqa: E402
from workflow_utils import load_embeddings, save_embeddings  # noqa: E402


def test_embedding_archive_round_trip(tmp_path):
    path = tmp_path / "embeddings.npz"
    expected = np.arange(12, dtype=np.float32).reshape(3, 4)
    save_embeddings(path, ["a", "b", "c"], expected, "test_model")

    ids, observed, model = load_embeddings(path)

    assert ids == ["a", "b", "c"]
    np.testing.assert_array_equal(observed, expected)
    assert model == "test_model"


def test_group_aware_training_produces_deployable_predictions():
    ids = [f"seq_{index}" for index in range(12)]
    values = np.asarray([index % 2 for index in range(12)], dtype=int)
    matrix = np.column_stack((values, 1 - values, np.arange(12), np.ones(12))).astype(float)
    tasks = {("activity", "label_a"): dict(zip(ids, values))}
    groups = {sequence_id: f"group_{index}" for index, sequence_id in enumerate(ids)}

    models, metrics, out_of_fold = fit_models(ids, matrix, tasks, groups, folds=3, seed=7)
    bundle = {
        "models": models,
        "feature_count": matrix.shape[1],
        "embedding_model": "test_model",
        "threshold": 0.5,
    }
    rows = predict(bundle, ids, matrix)

    assert len(metrics) == 1
    assert metrics[0]["n_splits"] == 3
    assert len(out_of_fold) == len(ids)
    assert len(rows) == len(ids)
    assert {row["task"] for row in rows} == {"activity"}
