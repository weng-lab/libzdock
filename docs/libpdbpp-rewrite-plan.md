<!--
Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
SPDX-License-Identifier: BSD-2-Clause
-->

# libpdb++ Rewrite Plan

Target: replace the C-style string manipulation inside `libpdb++` with a modern,
class-based implementation, without changing observable behavior until each
change is explicitly pinned and validated.

Reference specification: [wwPDB File Format v3.3+](https://www.wwpdb.org/documentation/file-format)
(record-by-record column layouts). Upstream libpdb++ implements the Feb 1992
Brookhaven format plus PDBRUN 5/6 USER records and the zlab `ZDATOM` extension.

## Status

Work proceeds on branch `pdbcpprewrite`; this document is the living plan and
work log for that effort.

| Phase | State |
| ----- | ----- |
| 0 — Document and pin current behavior | Characterization suite complete: `test/pdb_records.cpp` (21 test cases, 330 assertions) plus extended `test/formatting.cpp`; full suite 770 assertions / 55 test cases green. Sanitizer baseline complete: the only remaining finding is the pinned AUTHOR misaligned `%2d` store; see the findings log |
| 1 — Restructure legacy code into classes | Complete: `FieldReader` (cbcd7aa), `FieldWriter` (a6968fc), `RecordClassifier` (fd078d3), and the `PdbFormat` read (97303e4) and write (b274b94) sides are extracted and committed; `PDB` is a thin shell over them and the varargs `PDB::sscanf`/`PDB::sprintf` remain as compatibility wrappers. Table mutations were verified to fail the suite, and parse/write throughput on a 172k-record replay stays within 20% of the pre-refactor baseline |
| 2 — New modern record class (C++17) | Not started |
| 3 — Adapter shim and consumer migration | Not started |
| 4 — Cleanup, documentation, deprecation | Not started |

Tooling notes:

- `compile_commands.json` is generated with `make clean && bear -- make -j4 all`
  so clangd can provide per-edit diagnostics during the rewrite; the file is
  generated output and is git-ignored.
- The characterization suite is developed batch by batch (one small group of
  record fixtures per edit, each batch built and run before the next) so
  hand-derived golden lines are verified against actual behavior immediately.

Characterization findings log (updated as batches complete):

- ANISOU/SIGUIJ: verified the writer emits five blanks between the residue
  sequence number and the first `u[0]` field (columns 27–31); an initial
  four-space hand derivation failed and was corrected. Round-trip confirmed.
- `%D` conversions right-justify continuation numbers with zero rendered as
  blanks, and no separator follows before text fields; canonical heads
  therefore carry more leading blanks than wwPDB examples (AUTHOR, KEYWDS,
  MDLTYP, TITLE, EXPDTA, OBSLTE, SPRSDE, SPLIT, TURN, TVECT).
- `%Ns` string fields strip trailing blanks on read, but several writers
  right-justify (`%Ns`): HELIX/SHEET/TURN ids, CMPONT residue names, and the
  68-character AGRDES/SYMDES/MTXDES/CMPDES texts. Round-trips are only stable
  in the right-justified canonical form. Fields written with `%-Ns` (ATOM /
  HETATM / ANISOU names, SEQRES) round-trip stably after right-trimming.
- The v5 `%-38s`-style conversions are precision (maximum length), not width;
  v5 CNAME/COLOR/GFX COLOR/GFX FONT/GFX LABEL emit unpadded text.
- v5 `GFX LABEL` quoting degenerates: the reader captures the quotes into the
  text and the writer re-wraps them, so each parse/format cycle adds one pair
  of quotes. Pinned as a one-cycle assertion.
- v5 `USER  DISTANCE` uses `%8.3f` (not the v6 `%9.3f`); v5 `ANGLE` carries
  the extra `which` field.
- Defect (fixed): the ZDATOM write layout is 82 columns but `chars()`
  formatted into an 81-byte static buffer — a two-byte overflow — before the
  result was clamped and truncated to 80 characters. The record buffer now
  accommodates the widest legacy layout, so no write passes its end while the
  output stays clamped to the 80-character record length (dropping the last
  two charge columns); the Phase 2 modern class must keep records bounded by
  construction. Regression test: the full-width ZDATOM case in
  `test/pdb_records.cpp`.
- Defect (sanitizer-confirmed): AUTHOR parsing stores the `%2d` continuation
  through `&author.continuation`, a `char` at union offset 70 — a misaligned
  4-byte store (UB) that only stays inside the union because larger members
  follow. The modern class must use correctly typed fields.
- Quirk (pinned): a failed standard ATOM parse partially writes occupancy
  from its first `%6f` field before the ZDATOM fallback runs, so ZDATOM
  records surface a leaked occupancy value.
- Quirk (pinned): the UNKNOWN junk fallback calls `sprintf(junk, "%72s", buf)`
  from inside a member function, so unqualified lookup resolves it to the
  custom `PDB::sprintf`, not libc — `%72s` right-aligns short lines and
  renders overlong lines as 72 overflow `*` characters (see the overlong-line
  test). The table-driven parser must keep calling `PDB::sprintf` here.
- `isalpha()` returns true for parsed CA atoms: parsed names are
  right-trimmed to `" CA"`, which matches the 5-byte comparison. Only
  hand-built blank-padded names fail it. (An earlier draft of this document
  claimed the opposite; the characterization suite corrected it.)
- `MODEL` reads its number from columns 7–14 (`%8d`) but writes `%5d` at
  columns 11–15; both parse to the same value, only the writer's form
  round-trips.

Agreed strategy (decisions locked at planning time):

| Decision        | Choice                                                            |
| --------------- | ----------------------------------------------------------------- |
| Migration style | Phased, with an adapter shim; `libpdb::PDB` stays compiling until consumers are migrated, then is deprecated |
| C++ standard    | C++17 for new code (default flag bump; remains CLI-overridable)   |
| Record scope    | Parity first (1992 format + PDBRUN 5/6 + ZDATOM); new wwPDB v3.3+ records are a follow-up phase |

## 0. Ground rules (from `$libzdock-project-principles`)

- GNU Make stays the only build system; Catch2 and Eigen stay pinned submodules.
- Native code builds warning-free with `-Wall -Wextra -pedantic` on GCC, Clang,
  and Apple Clang.
- Upstream-derived files (`src/libpdb++/*`, `include/pdb++.h`, `doc/man/pdb++.3`)
  MUST keep their UCFS notices verbatim; they are also excluded from pre-commit
  rewriting. New project-owned files carry the `copyright.txt` header.
- Fixed-size record formatting MUST be bounded; tests force truncation and
  maximum-width paths.
- Every commit keeps `make clean && make -j4 all && make cpp-test &&
  make python-test` green. Safety-sensitive steps also run ASan/UBSan.
- `$git-commit-hygiene` governs staging/committing; commit messages must pass
  `tools/check_commit_message.py`.

## 1. Current-state inventory

### 1.1 Files and responsibilities

| File                               | Responsibility |
| ---------------------------------- | -------------- |
| `include/pdb++.h`                  | `libpdb::PDB`: `RecordType` enum, union of ~70 C structs, static PDBRUN version state, `memcmp` equality |
| `src/libpdb++/pdb_read.cpp`        | `PDB(const char*)`: giant `switch` over `RecordType`, drives `PDB::sscanf` with `read_format.i` / `pdbrun{5,6}_read.i` format strings; `UNKNOWN` fallback; ATOM retried as ZDATOM |
| `src/libpdb++/pdb_chars.cpp`       | `chars()`: giant `switch` driving `PDB::sprintf` with `write_format.i` / `pdbrun{5,6}_write.i`; writes into a **function-local static buffer**; trims trailing spaces |
| `src/libpdb++/pdb_sscanf.cpp`      | Custom fixed-width `sscanf` (varargs): `%d`/`%f`/`%s`/`%c`/`%<w> ` fields, defaults on short lines, `-1` on mismatch |
| `src/libpdb++/pdb_sprntf.cpp`      | Custom fixed-width `sprintf` (varargs): overflow char `*`, zero/space fill, unbounded writes into caller buffer |
| `src/libpdb++/pdb_type.cpp`        | `getType` (first-4-chars dispatch), `getGfxType`, `gfxChars`, `pdbrun{5,6}Type`; defines the static `pdbrunInputVersion`/`pdbrunOutputVersion` globals |
| `src/libpdb++/pdb++.cpp`           | `type()` (`memset(this)` + occupancy default 1.0 for ATOM), `isalpha()`, `byteCmp` |
| `src/libpdb++/pdbinput.cpp`        | `operator>>`: `getline` into 4×`BufLen` buffer, then parse |

### 1.2 Consumers inside this repository

- `zdock::PDB` (`include/PDB.hpp`, `src/pdb/PDB.cpp`): holds
  `std::shared_ptr<libpdb::PDB>` records; builds the Eigen coordinate matrix,
  models, filters, and lookups; touches `record->atom.*` and `record->model.num`.
- CLI tools (`createlig`, `createmultimer`, `pruning`, `constraints`,
  `centroids`, `zdsplit`, `zdunsplit`) consume `zdock::PDB` and stream records
  via `operator<<` → `chars()`. Their text output is regression-checked, so
  record text must stay byte-identical.
- Tests: `test/pdb.cpp` (wrapper behavior), `test/formatting.cpp` (bounded
  v5 USER record output).

### 1.3 Semantic quirks to characterize (each becomes a test)

1. `chars()` returns a pointer to a function-local static: last call wins;
   not reentrant; two records' text cannot be held simultaneously.
2. `PDB::sscanf` fixed-width semantics: end-of-line yields defaults
   (0 / 0.0 / empty string / spaces) instead of failure; `nmatch` counts every
   started field, including ignored fields; mismatch returns `-1`.
   Truncated records therefore parse with defaults rather than falling back.
3. The `0 > sscanf(...)` guard only detects hard mismatches (`-1`), never
   missing trailing fields.
4. Unknown records are stored as `sprintf(junk, "%72s", buf)` — right-aligned
   in a 72-char field.
5. `USER_PDBRUN` parsing mutates the **global** `pdbrunInputVersion`; output
   version is mutated as a side effect of `chars()` on USER_PDBRUN. Version
   leaks across files, instances, and threads.
6. `type()` byte-clears the whole object, then defaults ATOM occupancy to 1.0.
7. Equality is `memcmp`-based, including `rType` and padding.
8. Atom parsing falls back ATOM → ZDATOM format → UNKNOWN, setting `iszdatom`.
9. `operator>>` uses a 4×81-char buffer; overlong lines fail the stream and
   silently stop the consumer's read loop.
10. `%4s` fields keep leading spaces (e.g. atom name `" N  "`); consumers trim.
11. v5/v6 USER record layouts differ; several records format as empty (`count =
    0`) under the "wrong" output version.

### 1.4 Known coverage gaps

No per-record-type round-trip tests, no direct tests of the fixed-width
`sscanf` semantics, no tests for the PDBRUN version-state machine, no
malformed/overlong-line tests, and no UNKNOWN-fallback tests.

## 2. Phase 0 — Document and pin current behavior

**Goal: an executable specification of libpdb++ as it stands.**

1. Add `test/pdb_records.cpp`: table-driven characterization suite.
   - One canonical fixture line per `RecordType` (drawn from `test/data/`
     where possible, synthetic otherwise), plus v5 and v6 USER variants.
   - Round-trip assertions: `PDB(line).chars()` equals the canonical line
     (exact match, including trailing-space trimming).
   - Field assertions: parse each fixture and spot-check struct members against
     independently hand-derived column values (never restating production
     expressions).
   - Fallback assertions: truncated, mismatched, and overlong lines produce the
     documented quirk behavior from §1.3.
   - State assertions: pin the global PDBRUN version behavior (including its
     cross-instance leak) so the later removal is a deliberate, reviewed change.
2. Document semantics in this file (§1.3 is the seed) and cross-reference each
   record against its wwPDB v3.3 section; note the 4-char residue name and
   ZDATOM deviations from the spec.
3. Run the suite under ASan/UBSan to enumerate existing UB before touching
   anything: rebuild `CXXFLAGS="-O1 -g -fsanitize=address,undefined"` and run
   `test/test`.
4. Extend `test/formatting.cpp` truncation coverage from v5 USER records to
   all record types and both PDBRUN versions.

**Exit criteria:** characterization suite green in normal and sanitizer
builds; quirks documented; no production code changed.

## 3. Phase 1 — Restructure legacy code into cohesive classes (no behavior change)

Progress log (branch `pdbcpprewrite`):

- FieldReader extracted from `pdb_sscanf.cpp`; `PDB::sscanf` is now a thin
  format-string walker over the field primitives (commit cbcd7aa).
- FieldWriter extracted from `pdb_sprntf.cpp`, including the transient digit
  scratch; `PDB::sprintf` is now a thin format-string walker (commit a6968fc).
- RecordClassifier extracted from `pdb_type.cpp`; the PDBRUN input version is
  passed per call and the pinned globals stay on `PDB` (commit fd078d3).
- PdbFormat read side (complete): the `.i` format strings remain the single
  source of widths and literals, and per-record field lists (kind + offset)
  replace the varargs argument bindings. The seven v5 USER paths that use
  libc `sscanf` conversions (BGCOLOR, CNAME, COLOR, CHAIN, GFX COLOR,
  GFX FONT, GFX LABEL) keep explicit bindings because libc `va_list` calls
  cannot be table-driven. Tables were added one record family at a time with
  a build after each batch; a dev-time structural check verified that every
  field list pairs one-to-one with the value conversions of its format
  string before `pdb_read.cpp` switched over. The characterization suite
  passed unmodified.
- PdbFormat write side (complete): the write tables mirror the read tables,
  sharing field lists except where member types differ (AUTHOR's char
  continuation loads as one byte when formatting). The v5 USER paths that
  use libc `printf` conversions keep explicit `snprintf` bindings, null
  formats keep the count-zero suppressions, and the USER GFX BEGIN
  primitive-name substitution and the USER PDBRUN output-version side
  effect are preserved. `pdb_chars.cpp` now only selects the table-driven
  formatter and applies the pinned static-buffer, clamp, and trailing-blank
  trim.
- Table validation: mutating one read-side entry (v6 USER ANGLE pointed at
  the wrong field list) and one write-side entry (TVECT pointed at the TER
  list) each failed exactly one round-trip test case; both were reverted.
- Throughput: on a 172,200-record replay (2OOB ligand repeated 300x,
  `-O3 -std=c++14`, best of repeated runs), the table-driven parse is
  ~13% slower (125 to 143 ms) and the write path ~17% slower (45 to 52 ms)
  than the pre-refactor switch implementation — within the accepted 20%
  budget; the Phase 2 clean-room class recovers this by dropping the
  format-string walk entirely.
- Sanitizer baseline: with ASan/UBSan, the only remaining finding is the
  pinned AUTHOR misaligned `%2d` store; the ZDATOM record-buffer overrun is
  fixed by sizing the `chars()` record buffer for the widest legacy layout
  while keeping the 80-character output and the shared-static-buffer
  behavior unchanged.

The parse/format switches currently live inside the `PDB` constructor and
`chars()`. Extract them into internal, project-owned classes while keeping the
public `libpdb::PDB` API byte-for-byte compatible:

1. `src/libpdb++/PdbFormat.hpp/.cpp` (new, BSD-2 header): owns the format
   tables (`read_format.i`, `write_format.i`, `pdbrun*.i`) and per-record
   field descriptors as data instead of switch cases.
2. `src/libpdb++/FieldReader.hpp/.cpp`: encapsulates the fixed-width field
   extraction from `pdb_sscanf.cpp` (stateful cursor over the line; same
   defaults/nmatch semantics).
3. `src/libpdb++/FieldWriter.hpp/.cpp`: encapsulates bounded fixed-width
   formatting from `pdb_sprntf.cpp` (overflow `*`, fill, justification), always
   writing into a caller-owned `char(&)[BufLen]`.
4. `src/libpdb++/RecordClassifier.hpp/.cpp`: `getType`, `getGfxType`,
   `pdbrun{5,6}Type` dispatch as one cohesive class.
5. `libpdb::PDB` becomes a thin shell delegating to these classes. Keep the
   static-buffer `chars()` and the global version state for now — both are
   pinned by Phase 0 tests and removed only in Phase 2/3.
6. The varargs `PDB::sscanf`/`PDB::sprintf` remain as compatibility wrappers
   implemented on top of `FieldReader`/`FieldWriter`.

**Exit criteria:** clean build; characterization suite unchanged and green;
no new behavior; each extraction is a separate reviewable commit.

## 4. Phase 2 — New modern record class (C++17)

**Goal: a clean-room class with the same observable record text and parsed
values, no C-style string manipulation.**

New public header `include/Record.hpp`, namespace `zdock::pdb` (mirrors the
existing `zdock::PDB` home):

1. `enum class RecordType` — same record names as the legacy enum, plus
   `to_string(RecordType)`.
2. `struct Residue` and per-kind payloads (`AtomRecord`, `HetAtomRecord`,
   `ConnectRecord`, `Cryst1Record`, `ModelRecord`, `UserRecord`, …) with real
   types: `std::string` names, `std::optional` for spec-optional columns
   (element, charge, segment id, symmetry ops, ZDATOM extensions).
3. `class Record` holding `RecordType type_` and a
   `std::variant<Monostate, AtomRecord, …>` payload, with `isAtom()`,
   `asAtom()` accessors and `operator<<`.
   - Rationale: variant payloads make illegal field access a compile-time/
     runtime-typed error instead of a silent union misread. (Alternative
     considered: one aggregate with optionals — simpler but reintroduces
     "wrong-record field" bugs.)
4. Parsing:
   - `static std::optional<Record> parse(std::string_view line, const
     ParseContext&)` — no exceptions in the hot path; `std::nullopt` only for
     hard mismatches, matching legacy fallback behavior.
   - `ParseContext` carries the PDBRUN input version (replaces the global);
     parsing `USER PDBRUN` updates the context owned by the caller, never
     process-global state.
   - Field extraction via `std::string_view` slicing + `std::from_chars`;
     bounded `strtod` fallback where libc++ lacks floating-point
     `from_chars` (portability guard, isolated in one helper).
5. Formatting:
   - `std::string chars() const` — thread-safe by construction; identical
     field layout, overflow `*`, and trailing-space trim, validated against the
     Phase 0 golden lines for every record type and both PDBRUN versions.
   - Documented deviation (reviewed when reached): overlong unknown records are
     preserved verbatim instead of right-aligned junk, and `chars()` no longer
     aliases a static buffer. Golden tests updated deliberately, not silently.
6. Column layouts live in one data table per record type (column ranges +
   kind), replacing format strings entirely — no varargs, no format/argument
   mismatch class of bugs.
7. Unit tests for the new class in `test/record.cpp` mirroring the
   characterization table, plus adversarial cases (empty lines, NUL bytes,
   80-column overflow, negative coordinates at field width).

**Exit criteria:** new class builds warning-free; golden-line parity suite
green for every record; ASan/UBSan clean; legacy code untouched.

## 5. Phase 3 — Adapter shim and consumer migration

1. Add conversion shims: `zdock::pdb::Record::fromLegacy(const libpdb::PDB&)`
   and `toLegacy()` (temporary, tested both directions).
2. Switch `zdock::PDB` internals to `zdock::pdb::Record` while keeping its
   public API (records/attempts/matrix accessors) source-compatible for the
   tools; `Record` becomes the shared-pointer payload type.
3. Migrate `test/pdb.cpp` expectations to the new payload accessors; keep a
   bridge test asserting legacy↔modern round-trips.
4. Tools should require no changes (they only use `zdock::PDB` and
   `operator<<`); verify by running the full tool-level tests.

**Exit criteria:** all tools and tests build and pass with the new class in
the serving path; `libpdb::PDB` still compiles untouched for external users.

## 6. Phase 4 — Cleanup, documentation, deprecation

1. Mark `libpdb::PDB` deprecated in `include/PDB.hpp`-adjacent docs (upstream
   header itself stays untouched); README gains a "modern PDB record API"
   section pointing at `include/Record.hpp` and the wwPDB format reference.
2. Remove now-dead internal glue if the compiler no longer needs it; keep
   upstream files intact.
3. Full validation: `make clean && make -j4 all && make cpp-test &&
   make python-test && make python-check`, sanitizer build of the complete
   test suite, `git diff --check`, and `pre-commit run --all-files`.
4. Follow-up (out of scope here): extend `RecordType` and layouts to full
   wwPDB v3.3+ records (HETNAM, LINK, CISPEP, DBREF, …) using the same
   column-table mechanism.

## 7. Milestones and gates

| Phase | Deliverable | Gate |
| ----- | ----------- | ---- |
| 0 | Characterization suite + documented semantics + sanitizer baseline | Suite green; quirks documented; zero prod changes |
| 1 | `PdbFormat`/`FieldReader`/`FieldWriter`/`RecordClassifier` classes | Suite unchanged and green after every commit |
| 2 | `zdock::pdb::Record` (C++17) + golden parity tests | Byte-identical record text; sanitizers clean |
| 3 | `zdock::PDB` on new class + shims | Tools/tests pass unchanged |
| 4 | Deprecation notes, README, full validation | All repo checks green |

## 8. Risks and mitigations

| Risk | Mitigation |
| ---- | ---------- |
| Tool output must stay byte-identical (regression-checked downstream) | Phase 0 golden lines; parity tests in Phase 2; only reviewed deviations |
| `from_chars` double support varies by standard library | Isolated bounded fallback helper with a portability test |
| Removing global PDBRUN state changes multi-file behavior (currently leaks) | Pin current behavior in Phase 0; make the fix an explicit, tested change in Phase 2 (`ParseContext` per file) |
| Static-buffer `chars()` removal can expose latent double-use bugs | Shim keeps legacy signature; consumers already stream immediately |
| C++17 default bump vs macOS/FreeBSD toolchains | Flag stays overridable; C++17 is supported by all CI compilers; document minimum versions |
| Upstream license/notice rules on vendored code | New code lives in new project-owned files; upstream files untouched |
