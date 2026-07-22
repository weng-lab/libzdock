# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Integration tests for the compiled libzdock command-line tools."""

import subprocess
import tempfile
import unittest
from pathlib import Path

from zdock import ZDOCK

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "test" / "data"
BIN = ROOT / "bin"


class CLITests(unittest.TestCase):
    """Exercise CLI validation and output boundaries through built binaries."""

    def run_tool(self, tool: str, *arguments: str) -> subprocess.CompletedProcess[str]:
        """Run one built tool with a bounded timeout and captured text output.

        Args:
            tool: Executable name beneath the repository's bin directory.
            arguments: Command-line arguments passed to the executable.

        Returns:
            The completed process without raising for a nonzero exit status.
        """
        return subprocess.run(
            [str(BIN / tool), *arguments],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )

    def test_split_preserves_a_partial_final_chunk(self) -> None:
        """Splitting retains predictions from a non-full final chunk."""
        source = DATA / "2OOB" / "zdock.out.pruned"
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "chunk."
            result = self.run_tool("zdsplit", "-n", "100", "-p", str(prefix), str(source))

            self.assertEqual(result.returncode, 0, result.stderr)
            chunks = [Path(f"{prefix}{suffix}") for suffix in ("aaaa", "aaab", "aaac")]
            self.assertEqual(
                [ZDOCK(path).npredictions for path in chunks],
                [100, 100, 92],
            )
            self.assertEqual(
                sum(ZDOCK(path).npredictions for path in chunks),
                ZDOCK(source).npredictions,
            )

    def test_numeric_options_fail_without_terminating(self) -> None:
        """Malformed numeric options produce ordinary nonzero CLI exits."""
        cases = (
            ("createlig", "-n"),
            ("createmultimer", "-n"),
            ("centroids", "-n"),
            ("pruning", "-c"),
            ("zdsplit", "-n"),
        )
        for tool, option in cases:
            with self.subTest(tool=tool):
                result = self.run_tool(tool, option, "not-a-number")
                self.assertEqual(result.returncode, 1)
                self.assertIn("Invalid numeric option", result.stderr)
                trailing = self.run_tool(tool, option, "1x")
                self.assertEqual(trailing.returncode, 1)
                self.assertIn("Invalid numeric option", trailing.stderr)

    def test_split_rejects_zero_chunk_size(self) -> None:
        """Every requested output chunk must contain at least one pose."""
        source = DATA / "2OOB" / "zdock.out.pruned"
        result = self.run_tool("zdsplit", "-n", "0", str(source))

        self.assertEqual(result.returncode, 1)
        self.assertIn("Chunk size must be greater than zero", result.stderr)

    def test_multimer_rejects_unsupported_components(self) -> None:
        """Component selection cannot index outside the chain alphabet."""
        source = DATA / "ZDOCK" / "mzdock.out"
        result = self.run_tool("createmultimer", "-m", "-2", str(source))

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid component", result.stderr)

    def test_multimer_rejects_symmetry_above_chain_capacity(self) -> None:
        """Generating all components requires one chain ID per component."""
        source = DATA / "ZDOCK" / "symmetry_53.out"
        result = self.run_tool("createmultimer", str(source))

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid component", result.stderr)


if __name__ == "__main__":
    unittest.main()
