#!/usr/bin/env python3
"""Regression checks for occurrence-labelled, diagram-derived PS trees."""

from __future__ import annotations

import subprocess
import sys
import unittest
from functools import lru_cache
from itertools import combinations_with_replacement, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import process_list  # noqa: E402


class ProcessListResonanceRegression(unittest.TestCase):
    def setUp(self):
        process_list.SwitchFlavourScheme(5)
        process_list.diagram_decay_products = None
        process_list._resonance_history_cache.clear()
        process_list._diagram_channel_cache.clear()
        process_list._diagram_topology_cache.clear()
        process_list._physical_diagram_cache.clear()

    def histories(self, process):
        return set(process_list.DiscoverResonanceHistories(process))

    def topologies(self, process):
        return set(process_list.DiscoverDiagramTopologies(process))

    def channels(self, process):
        return set(process_list.DiscoverDiagramChannels(process))

    @staticmethod
    def nodes_by_mask(topology):
        return {node.mask: node for node in topology}

    @staticmethod
    def brute_force_physical_tree(external_pdgs, catalogue):
        """Small independent labelled-tree enumerator used as an oracle."""

        physical_pdgs = {
            pdg for vertex in catalogue.vertices for pdg in vertex
        }

        @lru_cache(maxsize=None)
        def partitions(mask, count):
            labels = [index for index in range(mask.bit_length())
                      if mask & (1 << index)]
            found = set()
            for assignment in product(range(count), repeat=len(labels)):
                groups = [0] * count
                for label, group in zip(labels, assignment):
                    groups[group] |= 1 << label
                if all(groups):
                    found.add(tuple(sorted(groups)))
            return tuple(sorted(found))

        @lru_cache(maxsize=None)
        def currents(mask):
            if mask & (mask - 1) == 0:
                return frozenset((external_pdgs[mask.bit_length() - 1],))
            found = set()
            for child_count in (2, 3):
                for child_masks in partitions(mask, child_count):
                    for children in product(*(currents(item)
                                              for item in child_masks)):
                        for parent in physical_pdgs:
                            if catalogue.contains((
                                    process_list._anti_pdg(parent), *children
                            )):
                                found.add(parent)
            return frozenset(found)

        full_mask = (1 << len(external_pdgs)) - 1
        return any(
            process_list._anti_pdg(external_pdgs[root]) in
            currents(full_mask ^ (1 << root))
            for root in range(len(external_pdgs))
        )

    def test_physical_catalogue_matches_brute_force_and_forbidden_rules(self):
        catalogue = process_list.SM_CATALOGUE
        allowed = (
            (21, 21, 21),
            (21, 21, 21, 21),
            (2, -2, -11, 11),
            (24, -24, 22, 22),
            (24, 1, -2),
        )
        forbidden = (
            (23, 23, 23),
            (24, 1, -4),
            (23, 11, -13),
        )
        for external in (*allowed, *forbidden):
            expected = self.brute_force_physical_tree(external, catalogue)
            self.assertEqual(
                process_list.has_physical_tree(external, catalogue), expected,
                external,
            )
        self.assertTrue(all(self.brute_force_physical_tree(
            external, catalogue) for external in allowed))
        self.assertFalse(any(self.brute_force_physical_tree(
            external, catalogue) for external in forbidden))

        self.assertTrue(catalogue.contains((24, 1, -2)))
        self.assertFalse(catalogue.contains((24, 1, -4)))
        self.assertFalse(catalogue.contains((23, 23, 23, 23)))
        self.assertTrue(catalogue.contains((24, -24, 22, 23)))
        self.assertFalse(any(
            set(vertex) & process_list.AUXILIARY_PDGS
            for vertex in catalogue.vertices
        ))

        # Exhaust every three- and four-leg multiset in a representative SM
        # basis.  The basis contains both charges, active and massive quarks,
        # two lepton generations, all physical vectors, and the Higgs; flavour
        # and charge permutations outside it are catalogue symmetries.
        basis = (
            -24, 24, 21, 22, 23, 25,
            -6, 6, -5, 5, -2, 2, -1, 1,
            -14, 14, -13, 13, -12, 12, -11, 11,
        )
        for multiplicity in (3, 4):
            for external in combinations_with_replacement(
                    basis, multiplicity):
                self.assertEqual(
                    process_list.has_physical_tree(external, catalogue),
                    self.brute_force_physical_tree(external, catalogue),
                    external,
                )
        for auxiliary in process_list.AUXILIARY_PDGS:
            self.assertFalse(process_list.has_physical_tree(
                (auxiliary, 21, 21), catalogue
            ))

    def test_diagram_trees_capture_lepton_line_radiation(self):
        process = ("u", "u~", "e+", "e-", "ve", "ve~")
        topologies = self.topologies(process)
        neutrino_pair = (1 << 4) | (1 << 5)
        electron_radiation = (1 << 3) | neutrino_pair
        positron_radiation = (1 << 2) | neutrino_pair

        self.assertTrue(any(
            self.nodes_by_mask(topology).get(neutrino_pair) is not None and
            self.nodes_by_mask(topology)[neutrino_pair].pdg == 23 and
            self.nodes_by_mask(topology).get(electron_radiation) is not None and
            self.nodes_by_mask(topology)[electron_radiation].pdg == 11 and
            frozenset((
                self.nodes_by_mask(topology)[electron_radiation].left,
                self.nodes_by_mask(topology)[electron_radiation].right,
            )) == frozenset(((1 << 3), neutrino_pair))
            for topology in topologies
        ))
        self.assertTrue(any(
            self.nodes_by_mask(topology).get(positron_radiation) is not None and
            self.nodes_by_mask(topology)[positron_radiation].pdg == -11
            for topology in topologies
        ))

    def test_complete_diagrams_keep_distinct_production_spines(self):
        process = ("u", "u~", "e+", "e-", "ve", "ve~")
        channels = self.channels(process)
        topologies = {channel.topology for channel in channels}
        self.assertEqual(len(channels), 19)
        self.assertEqual(len(topologies), 17)
        self.assertTrue(all(len(channel.topology) >= 2
                            for channel in channels))

        final_mask = sum(1 << label for label in range(2, len(process)))
        self.assertTrue(all(
            process_list.validate_topology(channel.topology, final_mask)
            for channel in channels
        ))
        self.assertTrue(all(
            process_list.validate_tree_vertices(
                channel.vertices, process_list.SM_CATALOGUE
            )
            for channel in channels
        ))

        ee = (1 << 2) | (1 << 3)
        nunu = (1 << 4) | (1 << 5)
        neutral_current = next(
            topology for topology in topologies
            if {(node.pdg, node.mask) for node in topology} ==
            {(23, ee), (23, nunu)}
        )
        self.assertEqual({
            channel.order for channel in channels
            if channel.topology == neutral_current
        }, {
            (0, 2, 3, 4, 5, 1),
            (0, 4, 5, 2, 3, 1),
        })

        # The complete diagram maps, rather than arbitrary recursion views,
        # are the specialized channels attached to the four safety maps.
        process_list.phase_space_orders = process_list.ProcessProcess(process)
        process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
        process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
        self.assertEqual(len(process_list.all_keys_sorted), 4 + len(channels))
        specialized = {
            (process_list.channel_topologies[key],
             process_list.PhaseSpaceOrderFromKey(key))
            for key in process_list.all_keys_sorted
            if process_list.channel_topologies[key]
        }
        self.assertEqual(specialized, {
            (channel.topology, channel.order) for channel in channels
        })
        self.assertTrue(all(
            len(row[2]) == len(process_list.all_keys_sorted)
            for key in process_list.all_keys_sorted
            for row in process_list.phase_space_orders[key]
        ))

    def test_qcd_attachments_reuse_the_colour_aware_production_maps(self):
        process_list.diagram_decay_products = ("e+", "e-")
        process = ("u", "u~", "g", "e+", "e-")
        channels = self.channels(process)
        self.assertEqual(len(channels), 2)
        self.assertEqual({channel.order for channel in channels}, {
            (0, 2, 3, 4, 1),
        })
        self.assertEqual({
            tuple((node.pdg, node.mask) for node in channel.topology)
            for channel in channels
        }, {
            ((22, (1 << 3) | (1 << 4)),),
            ((23, (1 << 3) | (1 << 4)),),
        })

        # The low-level enumerator can still distinguish which incoming line
        # emitted the QCD branch. The process-list layer deliberately merges
        # those views because its existing colour maps already cover them.
        external_pdgs = tuple(
            int(process_list.pdgs[particle]) for particle in process
        )
        raw_channels = process_list.discover_diagram_channels(
            external_pdgs,
            process_list.CURRENT_COMBINATIONS,
            process_list._anti_pdg,
            mapped_final_mask=(1 << 3) | (1 << 4),
        )
        self.assertEqual({channel.order for channel in raw_channels}, {
            (0, 2, 3, 4, 1),
            (0, 3, 4, 2, 1),
        })

        process_list.phase_space_orders = process_list.ProcessProcess(process)
        process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
        ordinary_order = process_list.PhaseSpaceOrderFromKey(
            process_list.all_keys_sorted[0]
        )
        process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
        self.assertEqual(len(process_list.all_keys_sorted), 3)
        self.assertEqual({
            process_list.PhaseSpaceOrderFromKey(key)
            for key in process_list.all_keys_sorted
        }, {ordinary_order})

    def test_qcd_projection_keeps_electroweak_descendants_only(self):
        process = ("u", "u~", "b", "b~")
        external_pdgs = tuple(int(process_list.pdgs[item])
                              for item in process)
        channels = process_list.discover_diagram_channels(
            external_pdgs,
            catalogue=process_list.SM_CATALOGUE,
        )
        mapped_pdgs = {
            node.pdg for channel in channels for node in channel.topology
        }
        self.assertIn(22, mapped_pdgs)
        self.assertIn(23, mapped_pdgs)
        self.assertNotIn(21, mapped_pdgs)

        radiative = ("u", "u~", "b", "b~", "e+", "e-")
        radiative_pdgs = tuple(int(process_list.pdgs[item])
                               for item in radiative)
        ee_mask = (1 << 4) | (1 << 5)
        radiative_channels = process_list.discover_diagram_channels(
            radiative_pdgs,
            catalogue=process_list.SM_CATALOGUE,
            mapped_final_mask=sum(1 << label
                                  for label in range(2, len(radiative))),
        )
        self.assertTrue(any(
            any(node.mask == ee_mask and node.pdg in (22, 23)
                for node in channel.topology)
            for channel in radiative_channels
        ))
        self.assertFalse(any(
            node.pdg == 21
            for channel in radiative_channels
            for node in channel.topology
        ))

    def test_density_canonicalization_uses_every_numerical_map_field(self):
        mask = (1 << 2) | (1 << 3)
        left = 1 << 2
        right = 1 << 3
        order = (0, 2, 3, 1)
        w_plus = process_list.DiagramChannel((process_list.TopologyNode(
            24, mask, left, right, process_list.BREIT_WIGNER, 24
        ),), order)
        w_minus = process_list.DiagramChannel((process_list.TopologyNode(
            -24, mask, left, right, process_list.BREIT_WIGNER, 24
        ),), order)
        z_map = process_list.DiagramChannel((process_list.TopologyNode(
            23, mask, left, right, process_list.BREIT_WIGNER, 23
        ),), order)
        reversed_order = process_list.DiagramChannel((
            process_list.TopologyNode(
                24, mask, left, right, process_list.BREIT_WIGNER, 24
            ),
        ), (0, 3, 2, 1))
        canonical = process_list.canonicalize_channels_by_density((
            w_plus, w_minus, z_map, reversed_order,
        ))
        self.assertEqual(len(canonical), 3)
        collapsed = next(channel for channel in canonical
                         if channel.order == order and
                         channel.topology[0].parameter == 24)
        self.assertEqual(collapsed.multiplicity, 2)

    def test_diagram_trees_keep_competing_ww_zz_and_gamma_currents(self):
        process = ("u", "u~", "e+", "e-", "ve", "ve~")
        topologies = self.topologies(process)
        ee = (1 << 2) | (1 << 3)
        nunu = (1 << 4) | (1 << 5)
        wp = (1 << 2) | (1 << 4)
        wm = (1 << 3) | (1 << 5)

        self.assertTrue(any(
            {(node.pdg, node.mask) for node in topology} ==
            {(23, ee), (23, nunu)}
            for topology in topologies
        ))
        self.assertTrue(any(
            {(node.pdg, node.mask) for node in topology} ==
            {(24, wp), (-24, wm)}
            for topology in topologies
        ))
        self.assertTrue(any(
            (22, ee) in {(node.pdg, node.mask) for node in topology}
            for topology in topologies
        ))

    def test_diagram_discovery_is_flavour_and_process_general(self):
        electron = self.topologies(
            ("u", "u~", "e+", "e-", "ve", "ve~")
        )
        muon = self.topologies(
            ("d", "d~", "mu+", "mu-", "vm", "vm~")
        )
        electron_shapes = {
            tuple((node.mask, node.left, node.right, abs(node.pdg) in (11, 12))
                  for node in topology)
            for topology in electron
        }
        muon_shapes = {
            tuple((node.mask, node.left, node.right, abs(node.pdg) in (13, 14))
                  for node in topology)
            for topology in muon
        }
        self.assertEqual(electron_shapes, muon_shapes)

        top_decay = self.topologies(
            ("u", "u~", "b", "e+", "ve", "b~", "e-", "ve~")
        )
        self.assertTrue(any(
            (24, (1 << 3) | (1 << 4)) in
            {(node.pdg, node.mask) for node in topology} and
            (6, (1 << 2) | (1 << 3) | (1 << 4)) in
            {(node.pdg, node.mask) for node in topology}
            for topology in top_decay
        ))

    def test_physical_contacts_use_flat_scaffolds_but_qcd_stays_gen23(self):
        def raw_topologies(process):
            external_pdgs = tuple(
                int(process_list.pdgs[particle]) for particle in process
            )
            return process_list.discover_diagram_topologies(
                external_pdgs,
                catalogue=process_list.SM_CATALOGUE,
            )

        higgs_vector = raw_topologies(("u", "u~", "h", "h", "z"))
        self.assertTrue(any(
            node.pdg == 0 and node.kind == process_list.FLAT_CONTACT
            for topology in higgs_vector
            for node in topology
        ))
        vector_contact = raw_topologies(
            ("u", "u~", "w+", "w-", "a", "a")
        )
        self.assertTrue(any(
            node.pdg == 0 and node.kind == process_list.FLAT_CONTACT
            for topology in vector_contact
            for node in topology
        ))
        self.assertFalse(any(
            node.pdg in process_list.AUXILIARY_PDGS
            for topologies in (higgs_vector, vector_contact)
            for topology in topologies
            for node in topology
        ))
        self.assertEqual(self.topologies(("g", "g", "g", "g")), set())

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

    def test_version_six_serializes_four_compact_catalogues(self):
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
            self.assertEqual(unique[0], "4 1 6")
            self.assertEqual(unique[1], "# process: p p > e+ e-")
            self.assertIn("channel_discovery=diagrams", unique[2])

            process = ("u", "u~", "e+", "e-")
            process_list.phase_space_orders = process_list.ProcessProcess(process)
            process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
            process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
            output = process_list.WriteAllProcsIntoList()
            self.assertTrue(any(line.startswith("PERMUTATIONS ")
                                for line in output))
            self.assertTrue(any(line.startswith("PHASE_MAPS ")
                                for line in output))
            self.assertTrue(any(line.startswith("PARTNER_SETS ")
                                for line in output))
            self.assertTrue(any(line.startswith("INTEGRATION_FAMILIES ")
                                for line in output))
            self.assertEqual(output[-1], "END_PROCESSES")

            node_records = [line.split() for line in output
                            if line.startswith("N ")]
            self.assertEqual({int(record[1]) for record in node_records},
                             {22, 23})
            self.assertTrue(all(len(record) == 7 for record in node_records))
            self.assertFalse(any(
                int(record[1]) in process_list.AUXILIARY_PDGS
                for record in node_records
            ))

            catalogues = process_list.BuildIntegrationCatalogues()
            self.assertEqual(len(catalogues.families), 1)
            self.assertEqual(len(catalogues.partner_sets), 1)
            self.assertEqual(len(catalogues.permutations), 1)
            self.assertEqual(len(catalogues.maps), 3)
            self.assertEqual(catalogues.partner_sets[0], tuple(sorted(
                set(catalogues.partner_sets[0])
            )))
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
