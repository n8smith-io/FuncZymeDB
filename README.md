# FuncZymeDB template

FuncZymeDB is a starter repository for building a curated enzyme, reaction,
compound, and sequence database for one protein family and for developing
reproducible functional-prediction workflows around it.

This repository is intentionally family-neutral. It contains the reusable
database builder and maintenance utilities, project conventions, empty
curation tables, and workflow contracts, but no biological records, trained
models, family-membership rules, publication figures, or analysis results.

## Start a new instance

Create a new repository with GitHub's **Use this template** button. In the new
repository:

1. Replace the placeholder values in `environment/instance.yaml`.
2. Add curated records to copies of the tables in `data/templates/`.
3. Define and document a reproducible protein-family membership rule.
4. Adapt or add only the workflow stages needed for the instance under `code/`.
5. Add dependencies to `environment/environment.yml` and tests to `code/tests/`.
6. Update this README with the instance's scope, setup, run order, data
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
3. **Functional labels** — family-appropriate labels and the rules used to
   derive them from curated records.
4. **Evaluation** — leakage-resistant grouping, held-out data, metrics,
   uncertainty, baselines, and model-selection criteria.
5. **Prediction** — versioned inputs and models, applicability limits, and
   auditable outputs for uncharacterized sequences.

Stage-specific code belongs in clearly named subdirectories of `code/`. Every
stage should document its inputs, outputs, assumptions, provenance, and exact
run command in a neighboring `methods.md` or README.

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

## Included database tools

- `build_enz_cpd_seq_database.py` builds the JSON database and sequence FASTA
  from the curated tables, with optional NCBI retrieval, manual sequences,
  lineage validation, sequence caching, and duplicate merging.
- `update_enz_cpd_db_from_tsv.py` previews and applies small approved curation
  updates without rebuilding the complete database.
- `add_npclassifier_to_compound_tsv.py` exports compound classifications from a
  built database into a TSV.
- `merge_species_lineages.py` adds missing species to a lineage table for
  subsequent manual curation.

See `code/database/methods.md` for schemas and commands. Family membership,
orthology, modeling, prediction, and visualization should be added only after
their family-specific scientific contracts are defined.

## Environment

The starter environment is deliberately small:

```bash
conda env create -f environment/environment.yml
conda activate funczymedb
cp .env.example .env  # only when local overrides are needed
```

Add and pin scientific or command-line dependencies when an instance actually
uses them. Keep `environment/environment.yml` and the tested environment in
sync.

## License and citation

Code and documentation are provided under the MIT License. Update
`CITATION.cff` with the instance authors and preferred citation before a
release. Curated datasets may require additional source-specific attribution;
record that provenance in the instance documentation.
