#!/usr/bin/env python3
"""Regression checks for automatic, occurrence-tagged resonance histories."""

from __future__ import annotations

import subprocess
import sys
import unittest
import itertools
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import process_list  # noqa: E402


class ProcessListResonanceRegression(unittest.TestCase):
    def setUp(self):
        process_list._resonance_history_cache.clear()

    def histories(self, process):
        return set(process_list.DiscoverResonanceHistories(process))

    def test_old_resonance_switch_is_removed(self):
        help_text = subprocess.run(
            [sys.executable, str(ROOT / "process_list.py"), "--help"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        self.assertNotIn("--resonance", help_text)
        rejected = subprocess.run(
            [sys.executable, str(ROOT / "process_list.py"), "--resonance",
             "p p > e+ e-"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(rejected.returncode, 0)

    def test_explicit_leptons_use_the_complete_valid_order_set(self):
        process = ("u", "u~", "e+", "e-", "mu+", "mu-")
        brute_force = {
            order for order in itertools.permutations(range(len(process)))
            if process_list.ValidColorOrd(process, order)
        }
        self.assertEqual(
            set(process_list.GenerateValidColorOrders(process)),
            brute_force,
        )

    def test_two_and_four_lepton_single_resonances(self):
        two_lepton = self.histories(("u", "u~", "e+", "e-"))
        self.assertIn(((23, (2, 3)),), two_lepton)

        four_lepton = self.histories(
            ("u", "u~", "e+", "e-", "mu+", "mu-")
        )
        self.assertIn(((23, (2, 3, 4, 5)),), four_lepton)
        self.assertIn(((23, (2, 3)),), four_lepton)
        self.assertIn(((23, (4, 5)),), four_lepton)
        self.assertIn(((23, (2, 3)), (23, (4, 5))), four_lepton)

        charged_four_lepton = self.histories(
            ("d", "u~", "e+", "ve", "mu+", "mu-")
        )
        self.assertIn(((24, (2, 3, 4, 5)),), charged_four_lepton)

    def test_ww_and_zz_are_competing_histories(self):
        histories = self.histories(
            ("u", "u~", "e+", "e-", "ve", "ve~")
        )
        self.assertIn(((23, (2, 3)), (23, (4, 5))), histories)
        self.assertIn(((24, (2, 4)), (-24, (3, 5))), histories)
        # The individual W/Z maps remain present alongside both pair maps.
        self.assertIn(((23, (2, 3)),), histories)
        self.assertIn(((24, (2, 4)),), histories)

    def test_higgs_requires_a_reachable_production_tree(self):
        vbf = self.histories(
            ("u~", "d~", "u", "d", "e+", "e-", "mu+", "mu-")
        )
        self.assertIn(((25, (4, 5, 6, 7)),), vbf)
        self.assertIn(
            ((23, (4, 5)), (23, (6, 7)), (25, (4, 5, 6, 7))),
            vbf,
        )
        no_higgs_tree = self.histories(
            ("g", "g", "e+", "e-", "mu+", "mu-")
        )
        self.assertFalse(any(
            resonance[0] == 25
            for history in no_higgs_tree
            for resonance in history
        ))

    def test_top_supports_explicit_and_leptonic_w_daughters(self):
        explicit = self.histories(
            ("u", "u~", "b", "w+", "b~", "w-")
        )
        self.assertIn(((6, (2, 3)),), explicit)
        self.assertIn(((-6, (4, 5)),), explicit)

        leptonic = self.histories(
            ("u", "u~", "b", "e+", "ve", "b~", "e-", "ve~")
        )
        self.assertIn(((24, (3, 4)), (6, (2, 3, 4))), leptonic)
        self.assertIn(((-24, (6, 7)), (-6, (5, 6, 7))), leptonic)
        self.assertIn(((6, (2, 3, 4)), (-6, (5, 6, 7))), leptonic)
        self.assertIn((
            (24, (3, 4)), (-24, (6, 7)),
            (6, (2, 3, 4)), (-6, (5, 6, 7)),
        ), leptonic)

        hadronic = self.histories(
            ("u", "u~", "b", "u", "d~", "b~", "d", "u~")
        )
        self.assertFalse(any(
            abs(resonance[0]) == 6
            for history in hadronic
            for resonance in history
        ))

    def test_six_leptons_have_three_disjoint_poles(self):
        histories = self.histories((
            "u", "u~", "e+", "e-", "mu+", "mu-", "ta+", "ta-"
        ))
        self.assertIn((
            (23, (2, 3)), (23, (4, 5)), (23, (6, 7))
        ), histories)

    def test_version_three_serializes_provenance_and_resonances(self):
        saved_options = dict(process_list.options)
        saved_provenance = process_list.process_provenance
        try:
            process_list.options.update({
                "flavour_scheme": 5,
                "include_3qqbar_processes": False,
                "include_cc_processes": False,
            })
            process_list.process_provenance = "p p > e+ e-"
            unique = process_list.WriteUniqueProcsIntoList([
                ("u", "u~", "e+", "e-")
            ])
            self.assertEqual(unique[0], "4 1 3")
            self.assertEqual(unique[1], "# process: p p > e+ e-")
            self.assertIn("resonance_discovery=automatic", unique[2])

            process = ("u", "u~", "e+", "e-")
            process_list.phase_space_orders = process_list.ProcessProcess(process)
            process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
            process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
            output = process_list.WriteAllProcsIntoList()
            resonance_headers = []
            for index, line in enumerate(output):
                fields = line.split()
                if len(fields) == 8 and fields[-1] == "1":
                    resonance_headers.append((index, fields))
            self.assertTrue(resonance_headers)
            for index, _ in resonance_headers:
                self.assertEqual(output[index + 1].split(), ["23", "2", "3", "4"])
        finally:
            process_list.options.clear()
            process_list.options.update(saved_options)
            process_list.process_provenance = saved_provenance


if __name__ == "__main__":
    unittest.main()
