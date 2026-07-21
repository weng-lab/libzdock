import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from zdock import ZDOCK


class ZDOCKTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.zdock_file = os.path.join(
            os.path.dirname(__file__), "..", "test", "data", "2OOB", "zdock.out.pruned"
        )

    def test_prediction_limit_is_exact(self):
        self.assertEqual(ZDOCK(self.zdock_file, n=3).npredictions, 3)

    def test_symmetry_error_describes_supported_format(self):
        with self.assertRaisesRegex(Exception, "only supported for M-ZDOCK"):
            ZDOCK(self.zdock_file, n=1).symmetry


if __name__ == "__main__":
    unittest.main()
