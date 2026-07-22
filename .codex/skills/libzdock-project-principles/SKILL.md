---
name: libzdock-project-principles
description: Project-specific engineering rules for libzdock native code, GNU Make portability, dependencies, licensing, CI, tests, and repository structure. Use for any C++ change, build or workflow change, dependency update, source branding change, native test work, or cross-language repository maintenance.
---

# libzdock Project Principles

## Normative Language

The keywords **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, and **MAY** are normative.

- **MUST** indicates an absolute requirement.
- **MUST NOT** indicates an absolute prohibition.
- **SHOULD** indicates a strong recommendation that may be departed from only for a documented reason.
- **SHOULD NOT** indicates a discouraged action that may be taken only for a documented reason.
- **MAY** indicates an optional action.

Interpret these keywords according to RFC 2119 and RFC 8174 when they appear in uppercase.

## Project Boundary

libzdock provides native PDB and ZDOCK/M-ZDOCK parsing, coordinate transformations, command-line tools, and a separate Python parser package.

Changes SHOULD remain small and explicit. New abstraction layers MUST have a concrete maintenance, safety, or reuse benefit.

## Build And Portability

- The project MUST use GNU Make. It MUST NOT introduce CMake unless the user explicitly reverses this decision.
- `GNUmakefile` MUST remain the authoritative build. `Makefile` MUST remain a clear BSD-make redirect to `gmake`.
- The build MUST support 64-bit Linux, macOS, and FreeBSD without assuming GNU userland outside GNU Make.
- Linux documentation and CI MAY use `make`. macOS and BSD documentation and CI SHOULD use `gmake` when the platform make is not GNU Make.
- Compiler, archiver, flags, and linker settings MUST remain overridable from the command line.
- Native code MUST build warning-free with the configured GCC, Clang, and Apple Clang warning set.
- Catch2 and Eigen MUST remain pinned submodules unless the user explicitly chooses another dependency mechanism.
- Dependency updates MUST include a clean rebuild rather than relying on stale objects.

## Licensing And Branding

- Project-owned files MUST use the short copyright and SPDX marker from `copyright.txt`.
- The project license MUST remain BSD 2-Clause unless the user explicitly changes it.
- `LICENSE` and `python/LICENSE` MUST carry the complete BSD 2-Clause text and current copyright range.
- Upstream `libpdb++` notices MUST remain intact. Agents MUST NOT relabel upstream-derived files as solely BSD-2-Clause.
- Fixed-format fixtures and vendored submodules MUST remain excluded from automated license and whitespace rewriting.
- `.pre-commit-config.yaml` SHOULD enforce branding, syntax, whitespace, spelling, and project quality rules.

## Native Correctness

- Object lifetime and copy behavior MUST use C++ value semantics rather than byte-clearing non-trivial objects.
- `zdock::PDB` copying MUST preserve the selected record graph and shared record identity without re-running user filters.
- Fixed-size record formatting MUST be bounded. Tests SHOULD force truncation and maximum-width paths, not only ordinary inputs.
- Intentional switch fallthrough MUST be explicit to both readers and compilers.
- M-ZDOCK transformation order and reverse-rotation compatibility MUST be treated as observable behavior.
- Coordinate and transformation tests SHOULD use independently calculated expected values rather than restating production expressions.

## Testing Standards

- Every test case MUST include a nearby comment or docstring that states the invariant being tested.
- Tests MUST assert externally meaningful behavior. They MUST NOT merely describe the current implementation.
- Regression tests SHOULD fail when the motivating defect is reintroduced.
- Boundary-sensitive code SHOULD use adversarial values, malformed input, and exact edge conditions.
- Identity, ownership, copy, and model tests MUST distinguish equality of values from sharing of objects.
- Changes to native code or shared behavior MUST run a clean build, native tests, and Python integration tests.
- Safety-sensitive native changes SHOULD also run AddressSanitizer and UndefinedBehaviorSanitizer when supported locally.

## Validation

Run the normal repository workflow with the platform's GNU Make command:

```console
make clean
make -j4 all
make cpp-test
make python-test
make python-check
```

On systems where GNU Make is installed as `gmake`, substitute `gmake` consistently.

Before handoff, agents MUST run `git diff --check` and the relevant pre-commit hooks. Generated build output and caches MUST NOT be committed.
