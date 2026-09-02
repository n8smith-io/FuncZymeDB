# FuncZymeDB template

FuncZymeDB is a starter repository for building a curated enzyme, reaction,
compound, and sequence database for one protein family and for developing
reproducible functional-prediction workflows around it.

This repository is intentionally family-neutral. It contains runnable tools
for database construction, orthology-based prediction (FuncPred-OG), and
protein-language-model prediction (FuncPred-AI), together with empty input
tables and workflow contracts. It contains no biological records, trained
models, family-membership rules, publication figures, or analysis results.

## Start a new instance

Create a new repository with GitHub's **Use this template** button. In the new
repository:

1. Replace the placeholder values in `environment/instance.yaml`.
2. Add curated records to copies of the tables in `data/templates/`.
3. Define and document a reproducible protein-family membership rule.
4. Define family-appropriate substrate tasks in `substrate_labels.tsv`.
5. Run or adapt the supplied orthology and model-training stages as needed.
6. Add dependencies to `environment/environment.yml` and tests to `code/tests/`.
7. Update this README with the instance's scope, setup, run order, data
   provenance, maintainers, and citation information.

Do not commit generated results, model caches, downloaded databases, secrets,
or machine-specific paths. See `data/README.md` and `.gitignore` for the default
data policy.

## Core workflow contracts

A FuncZymeDB instance should make the following boundaries explicit:

1. **Curation** — publication-specific activity records, compound vocabulary,
   sequence identifiers, evidence, and curator provenance.
2. **Family membership** — the reference model or set, software version,
   thresholds, accepted sequences, and rejected sequences.
3. **Substrate labels** — family-appropriate substrate classes and the rules
   used to derive them from curated records.
4. **Evaluation** — leakage-resistant grouping, held-out data, metrics,
   uncertainty, baselines, and model-selection criteria.
5. **Prediction** — versioned inputs and models, applicability limits, and
   auditable outputs for uncharacterized sequences.

The supplied stages document their inputs, outputs, assumptions, provenance,
and exact run commands in neighboring `methods.md` files. An instance may run
only the database stage; FuncPred-OG and FuncPred-AI are optional downstream
workflows, not requirements for using the database.

## Repository layout

```text
code/            workflow code and tests
data/templates/  tracked, empty curation-table templates
data/            curated inputs and ignored generated intermediates
docs/            tracked design and reproducibility documentation
environment/     instance configuration and environment specification
metadata/        capsule and release metadata guidance
results/         generated figures, reports, and tables
scratch/         exploratory and large intermediate work
```

The top-level scientific directories follow the Code Ocean capsule model. Keep
source and tests in `code/`, inputs in `data/`, reproducible outputs in
`results/`, disposable work in `scratch/`, and dependency specifications in
`environment/`.

## Included workflows

- `build_funczyme_database.py` builds the JSON database and sequence FASTA
  from the curated tables, with optional NCBI retrieval, manual sequences,
  lineage validation, sequence caching, and duplicate merging.
- `update_funczyme_database.py` previews and applies small approved curation
  updates without rebuilding the complete database.
- `export_compound_classifications.py` exports compound classifications from a
  built database into a TSV.
- `extend_species_lineages.py` adds missing species to a lineage table for
  subsequent manual curation.
- `build_orthology_database.py` indexes OrthoFinder groups and annotates the
  characterized database; `predict_with_funcpred_og.py` performs auditable
  orthogroup-based label transfer.
- `generate_protein_embeddings.py`, `train_funcpred_ai_models.py`, and
  `predict_with_funcpred_ai.py` provide an end-to-end, family-neutral ESMC and
  group-aware logistic-regression workflow.

See the `methods.md` file in each stage directory for schemas, commands, and
validation expectations. Publication-specific visualization is intentionally
excluded: figures belong in an instance or publication repository where their
scientific questions and result contracts are defined.

## Environment

Create the supplied environment to run the full workflow:

```bash
conda env create -f environment/environment.yml
conda activate funczymedb
cp .env.example .env  # only when local overrides are needed
```

ESMC weights download on first use. GPU acceleration is optional; add the CUDA
package appropriate to the target system when needed. Keep the environment file
and the tested environment in sync.

## License and citation

Code and documentation are provided under the MIT License. `CITATION.cff`
currently points to the associated preprint and should be updated when the work
is formally published. Each family instance should also identify its own
authors and release. Curated datasets may require additional source-specific
attribution; record that provenance in the instance documentation.
