<!--
Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
SPDX-License-Identifier: BSD-2-Clause
-->

# Repository Instructions

## Normative Language

The keywords **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, and **MAY** are normative.

- **MUST** indicates an absolute requirement.
- **MUST NOT** indicates an absolute prohibition.
- **SHOULD** indicates a strong recommendation that may be departed from only for a documented reason.
- **SHOULD NOT** indicates a discouraged action that may be taken only for a documented reason.
- **MAY** indicates an optional action.

Interpret these keywords according to RFC 2119 and RFC 8174 when they appear in uppercase. Lowercase uses retain their ordinary meaning.

## Required Skills

Before substantive work in this repository, agents MUST load and follow the skills relevant to the task:

`$libzdock-project-principles`

`$libzdock-python-engineering`

`$git-commit-hygiene`

Native code, build, CI, dependency, licensing, documentation, and cross-language changes MUST use `$libzdock-project-principles`.

Changes under `python/` or to Python-related Make, CI, or pre-commit integration MUST additionally use `$libzdock-python-engineering`.

Staging, committing, tagging, or pushing MUST use `$git-commit-hygiene`.

## Repository Shape

- `GNUmakefile` is the authoritative portable build for native tools and tests.
- `Makefile` exists only to direct BSD make users to GNU make.
- `include/` contains public native headers.
- `src/` contains native library and command-line implementation.
- `src/libpdb++/`, `include/pdb++.h`, and `doc/man/pdb++.3` contain upstream-derived code and notices.
- `test/` contains native tests and shared docking fixtures.
- `python/` is an independent uv-managed Python project with an installed-package test layout.
- `.codex/skills/` contains the durable engineering rules for this repository.

Agents MUST preserve unrelated user changes and the fixed-format data under `test/data/`.
