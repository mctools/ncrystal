Standalone test scripts
-----------------------

Standalone tests in this directory are pure python scripts which can be invoked
directly without being in a specific environment. They are mostly intended for
quick tests of source-code consistency.

Examples:

- That all files contain appropriate license boiler-plate.
- That directories do not contain any unexpected files.
- That files pass static linter checks.
- That the NCrystal version encoded in various places is consistent.
- That various pyproject.toml files have consistent meta-data.
- That source files does not contain white-space errors (e.g. non-unix line
  endings, non-utf8 chars, tabs, etc.).
- ...
