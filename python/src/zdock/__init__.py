# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Public interfaces for parsing ZDOCK and M-ZDOCK output."""

from .zdock import ZDOCK, Prediction, Structure, ZDOCKError

__all__ = ["Prediction", "Structure", "ZDOCK", "ZDOCKError"]
