# Methods — database construction

## Purpose

The database stage converts curated activity and compound tables into a JSON
database of enzymes, publication-specific reactions, compounds, identifiers,
and provenance. It can also build the corresponding sequence FASTA.

The stage is family-neutral. `SUBSTRATE` is the universal reaction input;
`DONOR`, `ACCEPTOR`, and `PRODUCT` are optional, non-exclusive roles used only
when they fit the enzyme family.

## Inputs

Copy the empty tables from `data/templates/` to `data/curated/`. Activity
`INDEX` values are required, unique, stable source-record identifiers. Compound
names may be bare names or semicolon-separated aliases, with the canonical name
first. Chemical structures and identifiers are optional.

If lineage validation is used, provide `data/curated/species_lineages.tsv` with
`BINOMIAL_NAME` and `CLADE1` through `CLADE8` columns. If manually curated
sequences are used, FASTA headers must identify activity-table enzymes as
described by `build_funczyme_database.py --help`.

## Build

Run from the repository root:

```bash
mkdir -p results/database
python code/database/build_funczyme_database.py \
  --compound_file data/curated/compound_data.tsv \
  --activity_file data/curated/activity_data.tsv \
  --species_lineages data/curated/species_lineages.tsv \
  --output_db results/database/funczymedb.json \
  --output_fasta results/database/sequences.fa \
  --person AB
```

In a Code Ocean reproducible run, the working directory is `/code`; use the
same interface with capsule-relative paths:

```bash
mkdir -p ../results/database
python database/build_funczyme_database.py \
  --compound_file ../data/curated/compound_data.tsv \
  --activity_file ../data/curated/activity_data.tsv \
  --species_lineages ../data/curated/species_lineages.tsv \
  --output_db ../results/database/funczymedb.json \
  --output_fasta ../results/database/sequences.fa \
  --person AB
```

Add `--fetch_seqs --ncbi_email <email>` for NCBI sequence retrieval. Use
`--seq_cache` to avoid refetching unchanged sequences, `--manual_fasta` for
curated sequences, and `--merge_seqs` with a tracked merge log when duplicate
sequence names require an irreproducible human choice.

## Maintenance utilities

Use `update_funczyme_database.py` for small, approval-driven curation fixes.
It previews compound-centered change packages before modifying a copy of the
database. Run each utility with `--help` for its complete interface.

`export_compound_classifications.py` exports cached compound classifications
from a built database. `extend_species_lineages.py` appends missing species with
blank lineage fields for manual curation.

## Reproducibility assumptions

- Preserve curated tables and manual merge decisions as tracked inputs.
- Record the responsible person and retain the database `log` block.
- Pin the environment and external reference-data versions.
- Treat fetched databases and web services as changing external dependencies;
  freeze release inputs and reuse caches when reproducing a release.
- Write generated database snapshots and FASTA files under `results/`.
