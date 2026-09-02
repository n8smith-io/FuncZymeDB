# Code layout

Add only the stages required by this FuncZymeDB instance. Prefer descriptive
subdirectory and script names over numbered filenames unless a directory is an
explicitly ordered pipeline.

Suggested stage names are `database`, `family_membership`, `orthology`,
`structure`, `phylogeny`, `modeling`, `prediction`, and `visualization`.
These are contracts, not mandatory modules.

Each implemented stage should include adjacent documentation that records:

- scientific purpose and assumptions;
- tracked inputs and generated outputs;
- software and reference-data versions;
- parameters, thresholds, random seeds, and responsible person;
- an exact command runnable from the repository root;
- tests or validation checks.

Scripts should read from `data/`, write generated artifacts to `results/` or
ignored data subdirectories, and avoid machine-specific absolute paths.
