---
name: git-commit-hygiene
description: Create, validate, and perform focused Git commits and pushes in libzdock, including staged-diff review, commit-message structure, hook installation, unrelated-change protection, and post-push verification. Use whenever staging, committing, tagging, pushing, or preparing commit text.
license: BSD-2-Clause
metadata:
  author: Arjan van der Velde, Weng Lab
  copyright: Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
  version: 1.0.0
---

# Git Commit Hygiene

## Normative Language

The keywords **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, and **MAY** are normative.

- **MUST** indicates an absolute requirement.
- **MUST NOT** indicates an absolute prohibition.
- **SHOULD** indicates a strong recommendation that may be departed from only for a documented reason.
- **SHOULD NOT** indicates a discouraged action that may be taken only for a documented reason.
- **MAY** indicates an optional action.

Interpret these keywords according to RFC 2119 and RFC 8174 when they appear in uppercase.

## Before Staging

- Run `git status --short --branch` before staging or committing.
- Stage only files relevant to the requested change.
- Leave unrelated, generated, untracked, or user-owned changes untouched.
- Agents MUST NOT amend, reset, discard, overwrite, or include user changes unless the user explicitly requests it.

## Before Committing

- Inspect the full staged diff with `git diff --cached`.
- Run `git diff --cached --check`.
- Construct the complete commit message in a file before invoking `git commit`.
- The message MUST describe the net staged change, not only the last action performed.
- Pre-commit and commit-msg hooks MUST be installed before committing.

## Commit Message Format

- The subject line MUST contain at most 79 characters.
- The subject MUST be followed by exactly one empty line.
- The body MUST be one paragraph containing at least two sentences.
- Every physical body line MUST contain at most 80 characters.
- The body MUST describe what changed. It MUST NOT describe failed attempts, omitted work, or the chronology used to reach the result.

## Required Validation

Validate the exact message file before committing:

```console
uv --directory python run python ../tools/check_commit_message.py /tmp/commit-message
```

If the Python environment is not yet available, `python3 tools/check_commit_message.py /tmp/commit-message` MAY be used because the checker has no third-party dependencies.

Use the validated file:

```console
git commit -F /tmp/commit-message
```

## Push And Tag Safety

- Verify the current branch, upstream, and new commit before pushing.
- Push only the branch or tag the user requested.
- Agents MUST NOT force-push unless the user explicitly authorizes rewriting that exact remote ref.
- Before creating a tag, verify the target commit and whether the tag already exists locally or remotely.
- After pushing, run `git status --short --branch` and confirm the local branch is synchronized.
