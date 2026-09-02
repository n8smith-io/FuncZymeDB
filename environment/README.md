# Environment

`environment.yml` defines the portable local environment and
`instance.yaml` records the instance-level scientific configuration. Keep both
in sync with the code and tests.

Code Ocean represents the capsule image with files in this directory, commonly
including a generated `Dockerfile` and optional `postInstall` script. Add those
files when exporting or publishing a capsule; do not install dependencies at
runtime from analysis scripts.
