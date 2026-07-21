import os
import subprocess
import tempfile
import unittest

from zdock import ZDOCK


class CLITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
        cls.data = os.path.join(cls.root, "test", "data")
        cls.bin = os.path.join(cls.root, "bin")

    def run_tool(self, tool, *arguments):
        return subprocess.run(
            [os.path.join(self.bin, tool), *arguments],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_split_preserves_a_partial_final_chunk(self):
        # Invariant: splitting never drops predictions from a non-full last chunk.
        source = os.path.join(self.data, "2OOB", "zdock.out.pruned")
        with tempfile.TemporaryDirectory() as directory:
            prefix = os.path.join(directory, "chunk.")
            result = self.run_tool("zdsplit", "-n", "100", "-p", prefix, source)

            self.assertEqual(result.returncode, 0, result.stderr)
            chunks = [prefix + suffix for suffix in ("aaaa", "aaab", "aaac")]
            self.assertEqual([ZDOCK(path).npredictions for path in chunks],
                             [100, 100, 92])
            self.assertEqual(sum(ZDOCK(path).npredictions for path in chunks),
                             ZDOCK(source).npredictions)

    def test_numeric_options_fail_without_terminating(self):
        # Invariant: malformed numeric options produce ordinary nonzero CLI exits.
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

    def test_split_rejects_zero_chunk_size(self):
        # Invariant: every requested output chunk must contain at least one pose.
        source = os.path.join(self.data, "2OOB", "zdock.out.pruned")
        result = self.run_tool("zdsplit", "-n", "0", source)

        self.assertEqual(result.returncode, 1)
        self.assertIn("Chunk size must be greater than zero", result.stderr)

    def test_multimer_rejects_unsupported_components(self):
        # Invariant: component selection cannot index outside the chain alphabet.
        source = os.path.join(self.data, "ZDOCK", "mzdock.out")
        result = self.run_tool("createmultimer", "-m", "-2", source)

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid component", result.stderr)

    def test_multimer_rejects_symmetry_above_chain_capacity(self):
        # Invariant: generating all components requires one chain ID per component.
        source = os.path.join(self.data, "ZDOCK", "symmetry_53.out")
        result = self.run_tool("createmultimer", source)

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid component", result.stderr)


if __name__ == "__main__":
    unittest.main()
