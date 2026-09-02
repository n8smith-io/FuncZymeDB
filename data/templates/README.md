# Curation table templates

Copy these empty tab-separated tables into `data/curated/` and populate them
for the instance. Do not edit the templates in place.

## Activity records

Each row in `activity_data.tsv` represents one publication-specific source
record. `INDEX` must be unique and stable. `SUBSTRATE` is the universal input;
`ACCEPTOR`, `DONOR`, and `PRODUCT` are optional, non-exclusive roles. Multiple
values are separated with semicolons. Record accession, evidence, acquisition,
review, and curator provenance wherever available.

`FAMILY_MEMBERSHIP_REFERENCE` and `FAMILY_MEMBERSHIP_CALL` record the evidence
used to include a sequence in this instance. They do not prescribe a particular
database or boundary method.

## Compound records

`compound_data.tsv` standardizes canonical names and aliases. Chemical
identifiers and structures are optional because not every substrate is a
discrete small molecule. State how structures and aliases were verified.

## Substrate labels

`substrate_labels.tsv` is the explicit contract between curation and both
FuncPred methods. Each row records a binary substrate observation for one
sequence, task, and family-defined substrate label. Tasks may represent
different substrate roles, ontology levels, or prediction endpoints. Both
positive (`1`) and experimentally meaningful negative (`0`) observations must
be explicit; do not silently convert missing or untested activities to zero.
`evaluation_group` defines the unit kept together during cross-validation,
normally an orthogroup or a stricter sequence-similarity cluster.
