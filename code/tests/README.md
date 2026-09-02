# Tests

Database, orthology, and machine-learning tests live together here because they
are executable project code. Keep fixtures small. Tests should
cover schema validation, stable identifiers, provenance records, deterministic
splits, and failure at invalid stage boundaries. Run tests with:

```bash
python -m pytest
```
