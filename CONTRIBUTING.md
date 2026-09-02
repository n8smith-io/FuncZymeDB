# Contributing

Open an issue before a substantial change so its scope, scientific assumptions,
and data implications can be agreed upon. Keep pull requests focused and use a
short, imperative commit subject.

For code changes:

1. Create a branch from `main`.
2. Update adjacent methods or run documentation.
3. Add or update tests under `code/tests/` and run them in the environment
   specified by `environment/environment.yml`.
4. Confirm generated data, credentials, and machine-specific paths are absent.
5. Describe the inputs, outputs, validation, and compatibility impact in the
   pull request.

Curated-data changes must identify their source, acquisition method, reviewer,
and any licensing constraints. Never rewrite an existing stable record
identifier to give it a new meaning.
