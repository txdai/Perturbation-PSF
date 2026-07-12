import unittest

import numpy as np
import pandas as pd

from perturbation_psf.audit import benjamini_hochberg, compare_replication_units


class AuditTests(unittest.TestCase):
    def test_bh_is_monotone_in_rank(self):
        q = benjamini_hochberg([0.001, 0.02, 0.2, 0.8])
        self.assertTrue(np.allclose(q, [0.004, 0.04, 0.26666667, 0.8]))

    def test_clustered_uncertainty_detects_source_reuse(self):
        rng = np.random.default_rng(7)
        rows = []
        for source in range(24):
            treat = int(source >= 12)
            shared = rng.normal(0, 1)
            for _ in range(100):
                rows.append(
                    {
                        "source_id": source,
                        "treat": treat,
                        "y": 0.15 * treat + shared + rng.normal(0, 0.25),
                    }
                )
        result = compare_replication_units(pd.DataFrame(rows), "y")
        self.assertGreater(result.se_source, result.se_naive * 3)
        self.assertEqual(result.n_sources, 24)


if __name__ == "__main__":
    unittest.main()
