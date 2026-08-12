#!/usr/bin/env python3
"""Focused regression checks for squared coupling-order process suffixes."""

import pathlib
import sys
import unittest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import process_list  # noqa: E402


class CouplingOrderProcessSyntaxTest(unittest.TestCase):
    def setUp(self):
        self.saved_options = dict(process_list.options)
        process_list.options["include_resonance"] = False

    def tearDown(self):
        process_list.options.clear()
        process_list.options.update(self.saved_options)

    def parse(self, suffix=""):
        return process_list.ParseCollision("p p > z 2j" + suffix)

    def header(self):
        return process_list.WriteUniqueProcsIntoList(
            [("g", "g", "g", "g")], []
        )[0]

    def test_no_suffix_selects_automatic_global_maximum_as(self):
        request = self.parse()
        self.assertEqual(
            request["coupling_orders"],
            process_list.DefaultCouplingOrderSelection(),
        )
        self.assertEqual(self.header(), "4 1 3 0 -1 -1 -1 -1")

    def test_exact_orders_are_squared_orders(self):
        request = self.parse(" aS=2 aEW=1.5")
        self.assertEqual(
            request["coupling_orders"],
            {
                "mode": process_list.COUPLING_ORDER_MODE_EXPLICIT,
                "as_min2": 4,
                "as_max2": 4,
                "aew_min2": 3,
                "aew_max2": 3,
            },
        )
        self.assertEqual(self.header(), "4 1 3 1 4 4 3 3")

    def test_repeated_bounds_are_intersected(self):
        request = self.parse(" aS>=1 aS>=1.5 aS<=3 aS<=2 aEW>=0")
        self.assertEqual(
            request["coupling_orders"],
            {
                "mode": process_list.COUPLING_ORDER_MODE_EXPLICIT,
                "as_min2": 3,
                "as_max2": 4,
                "aew_min2": 0,
                "aew_max2": -1,
            },
        )

    def test_explicit_suffix_does_not_change_particle_request(self):
        unrestricted = self.parse()
        restricted = self.parse(" aS>=0 aEW>=0")
        for field in ("initial_state", "jet_count", "rest", "lep"):
            self.assertEqual(unrestricted[field], restricted[field])

    def test_invalid_restrictions_are_rejected(self):
        invalid = (
            " aS=-1",
            " aS=1.25",
            " aS>1",
            " AS=1",
            " aS>=2 aS<=1",
            " aEW=2 z",
        )
        for suffix in invalid:
            with self.subTest(suffix=suffix), self.assertRaises(ValueError):
                self.parse(suffix)


if __name__ == "__main__":
    unittest.main()
