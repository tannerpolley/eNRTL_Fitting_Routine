"""Unequal source sizes expose overwrite and average-of-averages defects."""
from types import SimpleNamespace
import unittest

import pandas as pd

from Profile_DeltaCp import equilibrium_metrics


class Block(dict):
    def __init__(self, name, values, kind):
        super().__init__({i: SimpleNamespace(**{
            "fug_phase_comp" if kind == "VLE" else "mole_frac_phase_comp_true":
            {("Liq", "CO2" if kind == "VLE" else "MEA"): value}
        }) for i, value in enumerate(values)})
        self.local_name = name


class MetricAggregation(unittest.TestCase):
    def test_sources_counts_order_and_units(self):
        blocks = []
        for kind, sources in (("VLE", (("z", [2000.]), ("a", [1000., 3000.]))),
                              ("ChEq", (("y", [.2]), ("b", [.1, .4])))):
            for name, values in sources:
                data = pd.DataFrame({"CO2_loading": [.2] * len(values),
                                     "temperature": [40] * len(values),
                                     "MEA_weight_fraction": [.3] * len(values),
                                     "CO2_pressure" if kind == "VLE" else "MEA":
                                     [1. if kind == "VLE" else .1] * len(values)})
                blocks.append((kind, data, Block(name, values, kind)))
        result, rows = equilibrium_metrics(blocks)
        reverse, reverse_rows = equilibrium_metrics(list(reversed(blocks)))
        self.assertEqual(result, reverse)
        pd.testing.assert_frame_equal(rows, reverse_rows)
        self.assertEqual(result["vle"]["n"], 3)
        self.assertEqual(result["vle"]["sources"]["a"]["n"], 2)
        self.assertEqual(result["vle"]["sources"]["z"]["n"], 1)
        self.assertEqual(result["vle"]["pressure_mape_fraction"], 1.)
        self.assertEqual(result["vle"]["pressure_mape_percent"], 100.)
        self.assertEqual(result["vle"]["unit"], "kPa")
        self.assertEqual(result["speciation"]["n"], 3)
        self.assertAlmostEqual(result["speciation"]["mae"], .4 / 3)
        self.assertEqual(result["speciation"]["unit"], "mol/mol true species")


if __name__ == "__main__":
    unittest.main()
