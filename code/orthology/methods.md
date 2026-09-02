# Methods — FuncPred-OG

## Purpose

FuncPred-OG provides an interpretable orthology-based baseline and prediction
route. OrthoFinder defines orthogroups over a family-wide sequence set. New
queries are placed by their best MMseqs2 hit to that reference, and labels are
transferred only from characterized members of the assigned orthogroup.

## Build the orthology index

Run OrthoFinder with characterized and uncharacterized family members together,
using stable FASTA identifiers that match the FuncZymeDB enzyme identifiers.
Then run from the repository root:

```bash
python code/orthology/build_orthology_database.py \
  --orthogroups data/orthology/Orthogroups.tsv \
  --database results/database/funczymedb.json \
  --output results/orthology/orthology.json \
  --annotated-database results/database/funczymedb.orthology.json \
  --person AB
```

The OrthoFinder command, version, input FASTA files, and parameters are part of
the instance's scientific methods and must be recorded alongside its results.

## Predict

```bash
python code/orthology/predict_with_funcpred_og.py \
  --query-fasta data/prediction/query.fa \
  --reference-fasta data/orthology/family_reference.fa \
  --orthology results/orthology/orthology.json \
  --labels data/curated/substrate_labels.tsv \
  --output results/prediction/funcpred_og.tsv
```

The output reports the assigned orthogroup, best reference hit, percent
identity, number of characterized observations, label support fraction, and
prediction. Queries with no characterized labels in their assigned group
remain explicit abstentions.

Use leakage-resistant withheld groups when evaluating this transfer method.
Never count a query's own label as transfer evidence.
