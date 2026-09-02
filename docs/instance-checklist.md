# Instance checklist

## Scope

- [ ] Define the protein family and intended prediction tasks.
- [ ] Name maintainers and the supported contact channel.
- [ ] Replace every `REPLACE_ME` value in `config/instance.yaml`.

## Curation and provenance

- [ ] Copy and populate the empty curation tables.
- [ ] Define stable record and sequence identifiers.
- [ ] Record source, acquisition method, review status, and curator.
- [ ] Document licensing and redistribution constraints.
- [ ] Freeze release inputs with checksums or archived identifiers.

## Reproducibility

- [ ] Document the family-membership rule and rejected sequences.
- [ ] Pin software, reference data, parameters, and random seeds.
- [ ] Separate curated inputs from generated outputs.
- [ ] Add validation tests and a root-relative run guide.
- [ ] Record provenance in machine-readable outputs.

## Evaluation and release

- [ ] Define leakage-resistant groups before model selection.
- [ ] Include simple baselines and uncertainty estimates.
- [ ] State applicability limits and unsupported uses.
- [ ] Update the README, citation metadata, changelog, and release tag.
- [ ] Confirm no secrets, local paths, caches, or large artifacts are tracked.
