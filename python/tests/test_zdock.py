# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Tests for the installed ZDOCK parser package."""

import unittest
from pathlib import Path

from zdock import ZDOCK, ZDOCKError

DATA = Path(__file__).resolve().parents[2] / "test" / "data"


class ZDOCKTests(unittest.TestCase):
    """Verify parsing limits, metadata, and format validation."""

    zdock_file = DATA / "2OOB" / "zdock.out.pruned"

    def test_prediction_limit_is_exact(self) -> None:
        """Requesting n predictions returns exactly n predictions."""
        self.assertEqual(ZDOCK(self.zdock_file, n=3).npredictions, 3)

    def test_zero_prediction_limit_is_empty(self) -> None:
        """A zero limit still parses metadata while retaining no poses."""
        docking = ZDOCK(self.zdock_file, n=0)

        self.assertEqual(docking.npredictions, 0)
        self.assertEqual(docking.boxsize, 72)

    def test_negative_prediction_limit_is_rejected(self) -> None:
        """A negative limit cannot silently mean an unbounded parse."""
        with self.assertRaisesRegex(ValueError, "cannot be negative"):
            ZDOCK(self.zdock_file, n=-1)

    def test_symmetry_error_describes_supported_format(self) -> None:
        """Symmetry errors identify M-ZDOCK as the supported format."""
        with self.assertRaisesRegex(ZDOCKError, "only supported for M-ZDOCK"):
            _ = ZDOCK(self.zdock_file, n=1).symmetry

    def test_mzdock_metadata_is_parsed(self) -> None:
        """M-ZDOCK metadata exposes symmetry and its sole structure."""
        docking = ZDOCK(DATA / "ZDOCK" / "mzdock.out", n=1)

        self.assertTrue(docking.ismzdock)
        self.assertEqual(docking.symmetry, 24)
        self.assertEqual(docking.receptor.filename, "/tmp/structure.pdb")
        with self.assertRaisesRegex(ZDOCKError, "not supported"):
            _ = docking.ligand

    def test_invalid_and_mixed_predictions_are_rejected(self) -> None:
        """Malformed rows and format changes fail with parser errors."""
        cases = (
            "invalid_prediction.out",
            "mixed_after_mzdock.out",
            "mixed_after_zdock.out",
        )
        for filename in cases:
            with self.subTest(filename=filename):
                with self.assertRaises(ZDOCKError):
                    ZDOCK(DATA / "ZDOCK" / filename)


if __name__ == "__main__":
    unittest.main()
