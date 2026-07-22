# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Validate repository commit message formatting rules."""

import re
import sys
from pathlib import Path


def validate_commit_message(message: str) -> list[str]:
    """Return formatting errors found in a complete commit message."""
    lines = message.splitlines()
    errors: list[str] = []

    subject = lines[0] if lines else ""
    if len(subject) > 79:
        errors.append("Commit subject exceeds 79 characters.")

    second_line = lines[1] if len(lines) > 1 else None
    if second_line != "":
        errors.append("Commit message second line must be empty.")

    body_lines = lines[2:]
    if not body_lines or "\n".join(body_lines) == "":
        errors.append("Commit message body is required.")
    elif any(line == "" for line in body_lines):
        errors.append("Commit message body must be one paragraph.")

    body_text = " ".join(body_lines)
    sentence_count = len(re.findall(r"[.!?](?:\s+|$)", body_text))
    if sentence_count < 2:
        errors.append("Commit message body must contain at least two sentences.")

    errors.extend(
        f"Commit body line {line_number} exceeds 80 characters."
        for line_number, line in enumerate(lines, start=1)
        if line_number > 2 and len(line) > 80
    )
    return errors


def main(argv: list[str]) -> int:
    """Validate the commit message file named by the command-line arguments."""
    if len(argv) != 1:
        print("Usage: check_commit_message.py <commit-message-file>", file=sys.stderr)
        return 2

    message_file = Path(argv[0])
    errors = validate_commit_message(message_file.read_text(encoding="utf-8"))
    for error in errors:
        print(error, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
