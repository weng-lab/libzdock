#!/usr/bin/env python3
# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Shared constants and JSON coercion helpers for the parity tools."""

from __future__ import annotations

# The harness summary counter names, in output order.
COUNT_FIELDS = (
    "files",
    "lines",
    "type_mismatches",
    "text_mismatches",
    "field_mismatches",
    "state_mismatches",
    "rejections",
)


def as_int(value: object) -> int:
    """Coerce a decoded JSON value to int, rejecting anything else."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"expected a number, got {type(value).__name__}")
    return int(value)


def as_float(value: object) -> float:
    """Coerce a decoded JSON value to float, rejecting anything else."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"expected a number, got {type(value).__name__}")
    return float(value)


def as_str(value: object) -> str:
    """Coerce a decoded JSON value to str, rejecting anything else."""
    if not isinstance(value, str):
        raise ValueError(f"expected a string, got {type(value).__name__}")
    return value


def as_dict(value: object) -> dict[str, object]:
    """Coerce a decoded JSON value to a string-keyed object mapping."""
    if not isinstance(value, dict):
        raise ValueError(f"expected an object, got {type(value).__name__}")
    return value


def as_list(value: object) -> list[object]:
    """Coerce a decoded JSON value to a list."""
    if not isinstance(value, list):
        raise ValueError(f"expected an array, got {type(value).__name__}")
    return value
