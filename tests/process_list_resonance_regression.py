#!/usr/bin/env python3
"""Regression checks for automatic, occurrence-tagged resonance histories."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import process_list  # noqa: E402


class ProcessListResonanceRegression(unittest.TestCase):
    def setUp(self):
        process_list.SwitchFlavourScheme(5)
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

    def test_explicit_leptons_are_kept_in_physical_pair_blocks(self):
        process = ("u", "u~", "e+", "e-", "mu+", "mu-")
        orders = set(process_list.GenerateValidColorOrders(process))
        self.assertEqual(len(orders), 2)
        for order in orders:
            self.assertEqual(order[:2], (0, 1))
            for first, second in (order[2:4], order[4:6]):
                first_pdg = int(process_list.pdgs[process[first]])
                second_pdg = int(process_list.pdgs[process[second]])
                self.assertTrue(process_list._lepton_pair_currents(
                    first_pdg, second_pdg
                ))

        competing = ("u", "u~", "e+", "e-", "ve", "ve~")
        competing_orders = process_list.GenerateValidColorOrders(competing)
        matchings = {
            frozenset((frozenset(order[2:4]), frozenset(order[4:6])))
            for order in competing_orders
        }
        self.assertEqual(matchings, {
            frozenset((frozenset((2, 3)), frozenset((4, 5)))),
            frozenset((frozenset((2, 4)), frozenset((3, 5)))),
        })
        self.assertEqual(len(competing_orders), 4)

    def test_identical_leptons_keep_distinct_occurrence_pairings(self):
        process = ("u", "u~", "e+", "e-", "e+", "e-")
        pairings = process_list._lepton_pairings(process)
        matching_labels = {
            frozenset(frozenset(labels) for labels, _ in pairing)
            for pairing in pairings
        }
        self.assertEqual(matching_labels, {
            frozenset((frozenset((2, 3)), frozenset((4, 5)))),
            frozenset((frozenset((2, 5)), frozenset((3, 4)))),
        })

        orders = process_list.GenerateValidColorOrders(process)
        self.assertEqual(len(orders), 4)
        self.assertEqual({
            frozenset((frozenset(order[2:4]), frozenset(order[4:6])))
            for order in orders
        }, matching_labels)

        histories = self.histories(process)
        self.assertIn(((23, (2, 3)), (23, (4, 5))), histories)
        self.assertIn(((23, (2, 5)), (23, (3, 4))), histories)

    def test_unpaired_or_flavour_mismatched_leptons_are_rejected(self):
        saved_options = dict(process_list.options)
        try:
            process_list.options.update({
                "include_3qqbar_processes": False,
                "include_cc_processes": True,
            })
            self.assertFalse(process_list.ValidProc(("u", "u~", "ve")))
            self.assertFalse(process_list.ValidProc(
                ("u", "u~", "ve", "vm~")
            ))
            self.assertTrue(process_list.ValidProc(
                ("u", "u~", "ve", "ve~")
            ))
        finally:
            process_list.options.clear()
            process_list.options.update(saved_options)

    def test_lepton_pair_is_not_split_between_quark_lines(self):
        process = ("u", "d", "u~", "d~", "e+", "e-")
        base_order = (0, 2, 1, 3, 4, 5)
        orders = process_list.SingletMultiChannelOrders(
            process, base_order
        )
        self.assertEqual(len(orders), 2)
        for order, _ in orders:
            self.assertEqual(abs(order.index(4) - order.index(5)), 1)

    def test_two_and_four_leptons_use_complete_vector_pairings(self):
        two_process = ("u", "u~", "e+", "e-")
        two_pairing = process_list._lepton_pairings(two_process)
        self.assertEqual(two_pairing, ((((2, 3), (22, 23)),),))
        two_lepton = self.histories(two_process)
        self.assertIn((), two_lepton)
        self.assertIn(((23, (2, 3)),), two_lepton)
        self.assertFalse(any(
            resonance[0] == 22
            for history in two_lepton
            for resonance in history
        ))

        four_process = ("u", "u~", "e+", "e-", "mu+", "mu-")
        four_pairing = process_list._lepton_pairings(four_process)
        self.assertEqual(four_pairing, (
            (((2, 3), (22, 23)), ((4, 5), (22, 23))),
        ))
        four_lepton = self.histories(four_process)
        self.assertIn((), four_lepton)
        self.assertIn(((23, (2, 3)), (23, (4, 5))), four_lepton)
        self.assertTrue(all(
            not history or self.complete_lepton_pair_coverage(
                history, (2, 3, 4, 5)
            )
            for history in four_lepton
        ))

        charged_process = ("d", "u~", "e+", "ve", "mu+", "mu-")
        charged_pairing = process_list._lepton_pairings(charged_process)
        self.assertEqual(charged_pairing, (
            (((2, 3), (24,)), ((4, 5), (22, 23))),
        ))
        charged_four_lepton = self.histories(charged_process)
        self.assertIn((), charged_four_lepton)
        self.assertIn(((24, (2, 3)), (23, (4, 5))),
                      charged_four_lepton)
        self.assertFalse(any(
            history and not self.complete_lepton_pair_coverage(
                history, (2, 3, 4, 5)
            )
            for history in charged_four_lepton
        ))

    def test_ww_and_zz_are_competing_histories(self):
        histories = self.histories(
            ("u", "u~", "e+", "e-", "ve", "ve~")
        )
        self.assertIn(((23, (2, 3)), (23, (4, 5))), histories)
        self.assertIn(((24, (2, 4)), (-24, (3, 5))), histories)
        self.assertIn((), histories)
        self.assertFalse(any(len(history) == 1 for history in histories
                             if history))

    def test_higgs_requires_a_reachable_production_tree(self):
        vbf = self.histories(
            ("u~", "d~", "u", "d", "e+", "e-", "mu+", "mu-")
        )
        self.assertIn(
            ((23, (4, 5)), (23, (6, 7)), (25, (4, 5, 6, 7))),
            vbf,
        )
        self.assertTrue(all(
            (23, (4, 5)) in history and (23, (6, 7)) in history
            for history in vbf
            if any(resonance[0] == 25 for resonance in history)
        ))
        no_higgs_tree = self.histories(
            ("g", "g", "e+", "e-", "mu+", "mu-")
        )
        self.assertFalse(any(
            resonance[0] == 25
            for history in no_higgs_tree
            for resonance in history
        ))

    def test_higgs_yukawa_currents_follow_the_flavour_scheme(self):
        bottom_process = (
            "g", "g", "b", "b~", "e+", "e-", "ve", "ve~"
        )
        charm_process = (
            "g", "g", "c", "c~", "e+", "e-", "ve", "ve~"
        )

        process_list.SwitchFlavourScheme(5)
        self.assertFalse(any(
            resonance[0] == 25
            for history in self.histories(bottom_process)
            for resonance in history
        ))

        process_list.SwitchFlavourScheme(4)
        self.assertTrue(any(
            resonance[0] == 25
            for history in self.histories(bottom_process)
            for resonance in history
        ))
        self.assertFalse(any(
            resonance[0] == 25
            for history in self.histories(charm_process)
            for resonance in history
        ))

        process_list.SwitchFlavourScheme(3)
        self.assertTrue(any(
            resonance[0] == 25
            for history in self.histories(charm_process)
            for resonance in history
        ))

        for scheme in range(1, 6):
            process_list.SwitchFlavourScheme(scheme)
            higgs_quark_flavours = {
                abs(particle)
                for vertex in process_list.THREE_POINT_VERTICES
                if 25 in vertex
                for particle in vertex
                if 1 <= abs(particle) <= 6
            }
            self.assertEqual(
                higgs_quark_flavours,
                set(range(scheme + 1, 7)),
            )

    def test_top_supports_explicit_and_leptonic_w_daughters(self):
        explicit = self.histories(
            ("u", "u~", "b", "w+", "b~", "w-")
        )
        self.assertIn(((6, (2, 3)),), explicit)
        self.assertIn(((-6, (4, 5)),), explicit)

        leptonic = self.histories(
            ("u", "u~", "b", "e+", "ve", "b~", "e-", "ve~")
        )
        self.assertIn((
            (24, (3, 4)), (-24, (6, 7)), (6, (2, 3, 4)),
        ), leptonic)
        self.assertIn((
            (24, (3, 4)), (-24, (6, 7)), (-6, (5, 6, 7)),
        ), leptonic)
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

    def test_version_four_serializes_provenance_and_resonances(self):
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
            self.assertEqual(unique[0], "4 1 4")
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
            serialized_pdgs = set()
            for index, _ in resonance_headers:
                record = output[index + 1].split()
                self.assertEqual(record[1:], ["2", "3", "4"])
                serialized_pdgs.add(int(record[0]))
            self.assertEqual(serialized_pdgs, {23})
        finally:
            process_list.options.clear()
            process_list.options.update(saved_options)
            process_list.process_provenance = saved_provenance

    @staticmethod
    def complete_lepton_pair_coverage(history, lepton_labels):
        covered = []
        lepton_labels = set(lepton_labels)
        for resonance, labels in history:
            if resonance not in (22, 23, 24, -24) or len(labels) != 2:
                continue
            if set(labels) <= lepton_labels:
                covered.extend(labels)
        return sorted(covered) == sorted(lepton_labels)


if __name__ == "__main__":
    unittest.main()
