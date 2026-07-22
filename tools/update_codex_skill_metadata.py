# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Update Codex skill frontmatter with repository metadata."""

import sys
from pathlib import Path
from typing import Any

import yaml

LICENSE = "BSD-2-Clause"
AUTHOR = "Arjan van der Velde, Weng Lab"
COPYRIGHT = "Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab"
VERSION = "1.0.0"


def updated_skill_text(text: str) -> str:
    """Return skill Markdown with normalized frontmatter license metadata."""
    frontmatter, body = _split_frontmatter(text)
    loaded = yaml.safe_load(frontmatter)
    if not isinstance(loaded, dict):
        raise ValueError("skill frontmatter must be a YAML mapping")

    data: dict[str, Any] = dict(loaded)
    metadata = data.get("metadata", {})
    if not isinstance(metadata, dict):
        raise ValueError("skill metadata must be a YAML mapping when present")

    data["license"] = LICENSE
    data["metadata"] = {
        **metadata,
        "author": AUTHOR,
        "copyright": COPYRIGHT,
        "version": VERSION,
    }
    return _format_skill(data, body)


def _split_frontmatter(text: str) -> tuple[str, str]:
    """Split skill Markdown into YAML frontmatter and Markdown body."""
    if not text.startswith("---\n"):
        raise ValueError("skill must start with YAML frontmatter")
    marker = "\n---\n"
    end = text.find(marker, len("---\n"))
    if end == -1:
        raise ValueError("skill frontmatter must end with a closing marker")
    return text[len("---\n") : end], text[end + len(marker) :]


def _format_skill(frontmatter: dict[str, Any], body: str) -> str:
    """Format skill frontmatter followed by its original Markdown body."""
    ordered = _ordered_frontmatter(frontmatter)
    frontmatter_text = yaml.safe_dump(
        ordered,
        sort_keys=False,
        allow_unicode=False,
        width=1000,
    )
    return f"---\n{frontmatter_text}---\n{body}"


def _ordered_frontmatter(frontmatter: dict[str, Any]) -> dict[str, Any]:
    """Return frontmatter with stable top-level key order."""
    ordered: dict[str, Any] = {}
    for key in ("name", "description", "license", "metadata"):
        if key in frontmatter:
            ordered[key] = frontmatter[key]
    for key, value in frontmatter.items():
        if key not in ordered:
            ordered[key] = value
    return ordered


def discover_skill_files(root: Path) -> list[Path]:
    """Return all root-managed Codex skill files."""
    return sorted(root.glob(".codex/skills/*/SKILL.md"))


def main(argv: list[str]) -> int:
    """Update or check skill metadata for named or discovered skill files."""
    check_only = False
    arguments = list(argv)
    if arguments and arguments[0] == "--check":
        check_only = True
        arguments = arguments[1:]

    paths = [Path(argument) for argument in arguments] or discover_skill_files(
        Path.cwd()
    )
    if not paths:
        print("No Codex skill files found.", file=sys.stderr)
        return 1

    changed_paths: list[str] = []
    for path in paths:
        try:
            original = path.read_text(encoding="utf-8")
            updated = updated_skill_text(original)
            if updated != original:
                changed_paths.append(str(path))
                if not check_only:
                    path.write_text(updated, encoding="utf-8")
        except ValueError as exc:
            print(f"{path}: {exc}", file=sys.stderr)
            return 1

    if changed_paths:
        action = "Would update" if check_only else "Updated"
        print(f"{action} Codex skill metadata in:")
        for changed_path in changed_paths:
            print(f"  {changed_path}")
        if check_only:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
