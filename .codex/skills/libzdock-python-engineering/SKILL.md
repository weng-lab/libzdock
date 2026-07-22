---
name: libzdock-python-engineering
description: Engineering rules for libzdock's uv-managed Python package, including packaging, layout, typing, tests, tooling, documentation, Make integration, and CI. Use for changes under python/ or any workflow that installs, builds, validates, or publishes the zdock Python distribution.
license: BSD-2-Clause
metadata:
  author: Arjan van der Velde, Weng Lab
  copyright: Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
  version: 1.0.0
---

# libzdock Python Engineering

## Normative Language

The keywords **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, and **MAY** are normative.

- **MUST** indicates an absolute requirement.
- **MUST NOT** indicates an absolute prohibition.
- **SHOULD** indicates a strong recommendation that may be departed from only for a documented reason.
- **SHOULD NOT** indicates a discouraged action that may be taken only for a documented reason.
- **MAY** indicates an optional action.

Interpret these keywords according to RFC 2119 and RFC 8174 when they appear in uppercase.

## Project And Package Layout

- `python/pyproject.toml` MUST be the sole package and tool configuration file.
- `python/uv.lock` MUST be the sole Python lockfile and MUST remain committed.
- Dependency management and project commands MUST run through uv. Agents MUST NOT restore `setup.py`, ad hoc pip workflows, Poetry, or a second lockfile.
- Runtime code MUST remain under `python/src/zdock/` and tests under `python/tests/`.
- Tests MUST import the installed `zdock` package. They MUST NOT modify `sys.path` to import source files directly.
- Hatchling MUST remain the build backend unless a concrete packaging requirement justifies a change.
- The distribution version SHOULD track the repository release version.
- Supported Python versions MUST be declared in `requires-python`; Black MUST target the lowest supported version.

## API And Implementation

- The public import MUST remain `from zdock import ...`.
- Parser failures and unsupported format operations SHOULD raise `ZDOCKError` or a more specific documented exception.
- Public functions, classes, methods, reusable test helpers, and tests MUST have complete type annotations and informative docstrings.
- Return types MUST be precise. Functions MUST NOT silently return `None` when their contract implies a result.
- Filesystem code SHOULD use `pathlib.Path`.
- Subprocess calls MUST use explicit argument lists, `shell=False`, captured output, and bounded timeouts.
- Mutable implementation details SHOULD NOT be exposed when an immutable public view is sufficient.
- `Any`, local-import cycle workarounds, and type or lint suppressions SHOULD be avoided. Findings SHOULD be fixed structurally.

## Tooling

- Development dependencies MUST include current compatible versions of isort, Black, mypy, and pylint in uv's `dev` dependency group.
- isort, Black, mypy, and pylint settings MUST live in `python/pyproject.toml`.
- mypy MUST run in strict mode.
- Quality checks MUST run against both `src` and `tests`.
- Python-related pre-commit hooks MUST use the uv-managed environment.
- Generated `.venv`, `.uv-cache`, `__pycache__`, build, distribution, and tool-cache directories MUST NOT be treated as source.

## Tests And Documentation

- Behavior changes MUST include tests for the happy path and the relevant failure boundary.
- Every test MUST state its invariant in its docstring.
- Parser tests SHOULD cover ordinary ZDOCK, M-ZDOCK, malformed records, mixed formats, and limiting behavior as applicable.
- CLI integration tests MUST invoke the native binaries built by GNU Make and use a timeout.
- Markdown prose MUST NOT be artificially hard-wrapped to satisfy source-code line limits.
- `python/README.md` MUST document uv setup, test, quality, build, and import workflows.

## Validation

Run Python checks from the repository root:

```console
uv --directory python lock --check
uv --directory python run isort --check-only src tests
uv --directory python run black --check src tests
uv --directory python run mypy src tests
PYLINTHOME=python/.pylint.d uv --directory python run pylint src tests
uv --directory python run python -m unittest discover -s tests -v
uv --directory python build
```

The equivalent GNU Make entry points are `make python-test` and `make python-check`. Agents changing Make or CI integration MUST verify those entry points as well.

Before handoff, scan for unjustified suppressions:

```console
rg -n "type:\s*ignore|pylint:\s*disable|#\s*noqa|#\s*pyright|#\s*mypy|pragma:\s*no cover" python/src python/tests python/pyproject.toml
```
