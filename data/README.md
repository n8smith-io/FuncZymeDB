# Data policy

`data/templates/` contains empty tracked schemas. Copy them to `data/curated/`
when starting an instance, then explicitly add only small, irreplaceable
curated inputs whose licensing permits redistribution.

Generated intermediates, logs, downloaded resources, model weights, and large
arrays are ignored by default. Record how to regenerate them, including source
versions and content hashes. Use an external data archive or object store for
large release artifacts; do not bypass repository limits by committing them to
ordinary Git history.

Before committing any data, verify:

- redistribution is permitted;
- personal or confidential information is absent;
- source and curator provenance are recorded;
- identifiers are stable;
- generated and curated files are clearly separated.
