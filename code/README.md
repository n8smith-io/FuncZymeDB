# Code layout

The template includes three family-neutral, independently runnable stages:

- `database/` builds and maintains FuncZymeDB;
- `orthology/` builds the orthogroup index and runs FuncPred-OG;
- `machine_learning/` embeds proteins, trains models, and runs FuncPred-AI.

`tests/` stays under `code/` because tests are executable project code. Prefer
descriptive snake-case script names over acronyms or numbered filenames unless
a directory is an explicitly ordered pipeline.

An instance may add stages such as `family_membership`, `structure`, or
`phylogeny`. Publication-specific figures should remain outside the reusable
template.

Each implemented stage should include adjacent documentation that records:

- scientific purpose and assumptions;
- tracked inputs and generated outputs;
- software and reference-data versions;
- parameters, thresholds, random seeds, and responsible person;
- an exact command runnable from the repository root;
- tests or validation checks.

Scripts should read from `data/`, write generated artifacts to `results/` or
ignored data subdirectories, and avoid machine-specific absolute paths.
