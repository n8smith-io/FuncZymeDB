# Methods — FuncPred-AI

## Purpose

FuncPred-AI trains one binary logistic-regression model for each substrate label
defined by an enzyme-family instance. Protein features are mean-pooled ESMC
embeddings. Evaluation uses held-out biological groups, not random sequence
rows, to reduce information leakage among closely related proteins.

This template deliberately does not prescribe a substrate ontology. Copy
`data/templates/substrate_labels.tsv` to `data/curated/` and define appropriate
substrate tasks and labels for the family. Every sequence/label pair must have
an explicit binary value; an untested activity must not silently be treated as
inactive.

## Generate embeddings

```bash
python code/machine_learning/generate_protein_embeddings.py \
  --fasta results/database/sequences.fa \
  --model esmc_600m \
  --device cuda \
  --output results/machine_learning/characterized_embeddings.npz
```

The same ESMC model must be used during training and prediction. Model weights
download on first use and are cached by the ESM package.

## Train and evaluate

```bash
python code/machine_learning/train_funcpred_ai_models.py \
  --embeddings results/machine_learning/characterized_embeddings.npz \
  --labels data/curated/substrate_labels.tsv \
  --output-dir results/machine_learning/trained_models \
  --folds 5 \
  --seed 42
```

The training stage writes a deployable model bundle, out-of-fold predictions,
per-label average precision and ROC AUC, and a manifest containing versions,
parameters, and input hashes. `evaluation_group` should normally be an
orthogroup or a stricter sequence-similarity cluster. Select the grouping rule
before model comparison and document it for the instance.

## Predict

Generate embeddings for query proteins with the same embedding script and
model, then run:

```bash
python code/machine_learning/predict_with_funcpred_ai.py \
  --embeddings results/prediction/query_embeddings.npz \
  --models results/machine_learning/trained_models/funcpred_ai_models.joblib \
  --output results/prediction/funcpred_ai.tsv
```

The default 0.5 threshold is a transparent starting point, not a universally
validated biological cutoff. An instance should calibrate or select thresholds
using only training folds and report uncertainty and applicability limits.
