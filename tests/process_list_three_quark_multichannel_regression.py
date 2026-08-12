#!/usr/bin/env python3
"""Regression checks for compact 0--3-line phase-space metadata."""

from __future__ import annotations

import sys
import unittest
from collections import Counter, defaultdict
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import process_list  # noqa: E402


def finalized_rows(proc):
    """Build and finalize the process-list rows for one concrete subprocess."""
    process_list.phase_space_orders = process_list.ProcessProcess(proc)
    process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
    process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
    process_list.CheckConsistency()

    rows = []
    for channel, key in enumerate(process_list.all_keys_sorted):
        for row in process_list.phase_space_orders[key]:
            rows.append((channel, row))
    return rows


def finalized_request(process_string, flavour_scheme, include_three_lines):
    """Finalize one inclusive request without invoking the command-line tool."""
    saved_options = dict(process_list.options)
    try:
        process_list.SwitchFlavourScheme(flavour_scheme)
        process_list.options.update({
            "include_3qqbar_processes": include_three_lines,
            "include_cc_processes": False,
            "include_resonance": False,
            "serial": True,
        })
        request = process_list.ParseCollision(process_string)
        unique = process_list.GenerateAllUniqueProcs(request)
        concrete = process_list.GenerateAllProcs(unique, request)
        process_list.phase_space_orders = process_list.CombineResults(
            [process_list.ProcessProcess(proc) for proc in concrete]
        )
        process_list.all_keys_sorted = sorted(process_list.phase_space_orders)
        process_list.DetermineMultiChannelPartnersAndSymmetryFactor()
        process_list.CheckConsistency()

        rows = [
            row
            for key in process_list.all_keys_sorted
            for row in process_list.phase_space_orders[key]
        ]
        three_line_keys = [
            key for key in process_list.all_keys_sorted
            if len(key[1]) > 1 and key[1][1] == 3
        ]
        return {
            "groups": len(process_list.all_keys_sorted),
            "rows": len(rows),
            "headers": len({
                process_list.PhaseSpaceOrderFromKey(key)
                for key in process_list.all_keys_sorted
            }),
            "three_line_groups": len(three_line_keys),
            "three_line_rows": sum(
                len(process_list.phase_space_orders[key])
                for key in three_line_keys
            ),
            "rank_one_groups": sum(key[2] == 1 for key in three_line_keys),
            "rank_one_rows": sum(
                len(process_list.phase_space_orders[key])
                for key in three_line_keys if key[2] == 1
            ),
        }
    finally:
        process_list.SwitchFlavourScheme(5)
        process_list.options.clear()
        process_list.options.update(saved_options)


def rotated_at_zero(order):
    zero = order.index(0)
    return tuple(order[zero:] + order[:zero])


class ThreeQuarkLineMultiChannelRegression(unittest.TestCase):
    def assert_permutation_metadata(self, rows):
        """Check base/target maps and their alignment with partner channels."""
        rows_by_channel = defaultdict(list)
        for channel, row in rows:
            rows_by_channel[channel].append(row)

        for channel, row in rows:
            self.assertEqual(len(row), 6)
            proc, order, partners, factor, permutation, partner_maps = row
            del factor
            nlegs = len(proc)
            self.assertEqual(tuple(sorted(set(partners))), partners)
            self.assertIn(channel, partners)
            self.assertEqual(sorted(permutation), list(range(nlegs)))
            self.assertEqual(len(partner_maps), len(partners))
            self.assertTrue(all(
                sorted(partner_map) == list(range(nlegs))
                for partner_map in partner_maps
            ))
            self.assertEqual(
                partner_maps[partners.index(channel)],
                permutation,
            )

            key = process_list.all_keys_sorted[channel]
            canonical_header = process_list.PhaseSpaceOrderFromKey(key)
            natural_header = rotated_at_zero(order)
            self.assertEqual(
                tuple(permutation[label] for label in canonical_header),
                natural_header,
            )
            if sum(p in process_list.quarks for p in proc) < 3:
                self.assertEqual(permutation, tuple(range(nlegs)))

            # Every partner map is serialized in exactly the same order as its
            # channel id.  A component has one fixed process but may have a
            # different colour order and label map in each partner channel.
            for partner, partner_map in zip(partners, partner_maps):
                candidates = [
                    other
                    for other in rows_by_channel[partner]
                    if other[0] == proc
                    and other[2] == partners
                    and other[4] == partner_map
                ]
                self.assertTrue(
                    candidates,
                    msg=f"missing partner row in channel {partner}",
                )

            if sum(p in process_list.quarks for p in proc) == 3:
                blocks = process_list.ThreeQuarkLineBlocks(proc, order)
                raw_swapped_order = tuple(blocks[0] + blocks[2] + blocks[1])
                swapped_proc, swapped_order = process_list.OrderProcPerm(
                    proc, raw_swapped_order
                )
                swapped_channels = process_list.process_order_to_index[
                    (swapped_proc, swapped_order)
                ]
                overlap = set(swapped_channels) & set(partners)
                if swapped_proc == proc and swapped_order == raw_swapped_order:
                    self.assertTrue(overlap)
                else:
                    # Equality only after a final-leg relabelling is suitable
                    # for adaptive-grid sharing, not pointwise MIS.
                    self.assertFalse(overlap - {channel})

    def test_canonical_header_keeps_the_fixed_target(self):
        proc = ("d~", "u~", "d", "u", "s", "s~", "z")
        raw = process_list.ProcessProcess(proc)
        self.assertEqual(sum(map(len, raw.values())), 36)
        for rows in raw.values():
            for target_proc, target_order, _ in rows:
                canonical, permutation, base_proc, base_order = \
                    process_list.CanonicalPhaseMap(target_proc, target_order)
                self.assertEqual(target_proc, proc)
                self.assertEqual(sorted(permutation), list(range(len(proc))))
                self.assertEqual(permutation[:2], (0, 1))
                self.assertEqual(permutation[proc.index("z")], proc.index("z"))
                self.assertEqual(
                    tuple(permutation[label] for label in canonical),
                    rotated_at_zero(target_order),
                )
                self.assertEqual(
                    base_proc,
                    tuple(target_proc[permutation[i]] for i in range(len(proc))),
                )
                self.assertEqual(rotated_at_zero(base_order), canonical)

    def test_distinct_three_line_maps_are_genuine_partners(self):
        rows = finalized_rows(("d~", "u~", "d", "u", "s", "s~"))
        self.assertEqual(len(process_list.all_keys_sorted), 12)
        self.assertEqual(len(rows), 12)
        self.assertEqual(
            Counter((len(row[2]), row[3]) for _, row in rows),
            Counter({(2, 1.0): 12}),
        )
        self.assert_permutation_metadata(rows)

    def test_q3z_retains_six_single_coefficient_maps(self):
        rows = finalized_rows(
            ("d~", "u~", "d", "u", "s", "s~", "z")
        )
        self.assertEqual(len(process_list.all_keys_sorted), 36)
        self.assertEqual(len(rows), 36)
        self.assertEqual(
            {(len(row[2]), row[3]) for _, row in rows},
            {(6, 1.0)},
        )
        self.assertEqual(len({row[4] for _, row in rows}), 12)
        self.assertEqual({key[2] for key in process_list.all_keys_sorted}, {0})
        self.assertEqual(
            Counter(key[1] for key in process_list.all_keys_sorted),
            Counter({
                ("flow", 3, (0, 1, 2)): 6,
                ("flow", 3, (0, 2, 1)): 6,
                ("flow", 3, (1, 0, 2)): 6,
                ("flow", 3, (1, 2, 0)): 6,
                ("flow", 3, (2, 0, 1)): 6,
                ("flow", 3, (2, 1, 0)): 6,
            }),
        )
        self.assert_permutation_metadata(rows)

    def test_component_local_source_rank_prevents_map_loss(self):
        rows = finalized_rows(("d~", "d", "u", "u~", "s", "s~"))
        self.assertEqual(len(rows), 12)
        self.assertEqual(len(process_list.all_keys_sorted), 8)
        self.assertEqual(
            Counter(key[2] for key in process_list.all_keys_sorted),
            Counter({0: 6, 1: 2}),
        )
        self.assertEqual({(len(row[2]), row[3]) for _, row in rows}, {(2, 1.0)})
        self.assert_permutation_metadata(rows)

    def test_identical_flavours_preserve_normalization(self):
        for proc in (
            ("d~", "d", "u", "u~", "u", "u~"),
            ("u~", "u", "u", "u~", "u", "u~"),
        ):
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(process_list.all_keys_sorted), 3)
                self.assertEqual(
                    Counter((len(row[2]), row[3]) for _, row in rows),
                    Counter({(1, 2.0): 3}),
                )
                # The first two targets are exchanged only after relabelling
                # the identical final-state pairs.  They may share a grid,
                # but must not enter one another's pointwise MIS sum.
                self.assertEqual(rows[0][1][2], (rows[0][0],))
                self.assertEqual(rows[1][1][2], (rows[1][0],))
                self.assert_permutation_metadata(rows)

    def test_singlets_compose_with_both_string_orders(self):
        cases = (
            (("d~", "u~", "d", "u", "s", "s~", "a"), 36, 6),
            (("d~", "u~", "d", "u", "s", "s~", "a", "z"), 144, 24),
        )
        for proc, expected_rows, expected_partners in cases:
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(rows), expected_rows)
                self.assertEqual(len(process_list.all_keys_sorted), expected_rows)
                self.assertEqual(
                    {(len(row[2]), row[3]) for _, row in rows},
                    {(expected_partners, 1.0)},
                )
                self.assert_permutation_metadata(rows)

    def test_gluon_maps_share_only_canonical_topologies(self):
        cases = (
            (("d~", "u~", "d", "u", "s", "s~", "g"), 24, 36, 2, 1.0),
            (("d~", "u~", "d", "u", "s", "s~", "g", "g"), 30, 72, 2, 2.0),
            (("d~", "d", "u", "u~", "s", "s~", "g", "z"), 48, 108, 6, 1.0),
        )
        for proc, groups, nrows, partners, factor in cases:
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(process_list.all_keys_sorted), groups)
                self.assertEqual(len(rows), nrows)
                metadata = Counter((len(row[2]), row[3]) for _, row in rows)
                if proc[-2:] == ("g", "g"):
                    self.assertEqual(
                        metadata,
                        Counter({(2, 2.0): 60, (1, 1.0): 12}),
                    )
                else:
                    self.assertEqual(set(metadata), {(partners, factor)})
                self.assert_permutation_metadata(rows)

    def test_inclusive_z_four_jet_channel_count_is_compact(self):
        without_three = finalized_request("p p > z 4j", 5, False)
        self.assertEqual(without_three, {
            "groups": 90,
            "rows": 4200,
            "headers": 30,
            "three_line_groups": 0,
            "three_line_rows": 0,
            "rank_one_groups": 0,
            "rank_one_rows": 0,
        })

        with_three = finalized_request("p p > z 4j", 3, True)
        self.assertEqual(with_three, {
            "groups": 318,
            "rows": 4230,
            "headers": 30,
            "three_line_groups": 228,
            "three_line_rows": 2790,
            "rank_one_groups": 18,
            "rank_one_rows": 108,
        })

    def test_inclusive_four_jet_generation_enumerates_multisets(self):
        saved_options = dict(process_list.options)
        try:
            process_list.SwitchFlavourScheme(5)
            process_list.options.update({
                "include_3qqbar_processes": True,
                "include_cc_processes": False,
                "include_resonance": False,
                "serial": True,
            })
            request = process_list.ParseCollision("p p > 4j")
            with mock.patch.object(
                process_list,
                "ValidProc",
                wraps=process_list.ValidProc,
            ) as valid_proc:
                unique = process_list.GenerateAllUniqueProcs(request)

            # There are only 1 + 5^2 + 15^2 + 35^2 possible flavour
            # multisets through three quark lines.  The old ordered-product
            # construction called ValidProc 11^6 times for this request.
            self.assertLessEqual(valid_proc.call_count, 1476)
            self.assertEqual(
                Counter(
                    sum(part in process_list.quarks for part in proc)
                    for proc in unique
                ),
                Counter({3: 35, 2: 15, 1: 5, 0: 1}),
            )
        finally:
            process_list.SwitchFlavourScheme(5)
            process_list.options.clear()
            process_list.options.update(saved_options)

    def test_zero_one_and_two_line_metadata_remain_compact(self):
        cases = (
            (("g", "g", "g", "g", "g", "g"), 5, 24.0),
            (("u~", "g", "u", "g", "g", "g"), 4, 6.0),
            (("d~", "u~", "d", "u", "g", "g"), 6, 2.0),
            (("u~", "u~", "u", "u", "g", "g"), 3, 4.0),
        )
        for proc, groups, factor in cases:
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(process_list.all_keys_sorted), groups)
                self.assertEqual(len(rows), groups)
                self.assertEqual(
                    {(len(row[2]), row[3]) for _, row in rows},
                    {(1, factor)},
                )
                self.assert_permutation_metadata(rows)

    def test_serialization_appends_current_and_partner_maps(self):
        proc = ("d~", "u~", "d", "u", "s", "s~", "z")
        finalized_rows(proc)
        process_list.options["include_resonance"] = False
        output = process_list.WriteAllProcsIntoList([])
        self.assertEqual(output[0], "36")

        position = 2
        for group in range(36):
            header = tuple(int(value) for value in output[position].split())
            self.assertEqual(header[0], group + 1)
            nprocesses = header[1]
            self.assertEqual(len(header), 3 + len(proc))
            position += 1
            for _ in range(nprocesses):
                fields = output[position].split()
                npartners = int(fields[0])
                offset = 1 + npartners + 2 * len(proc)
                float(fields[offset])
                permutation = tuple(
                    int(value) for value in
                    fields[offset + 1:offset + 1 + len(proc)]
                )
                flattened_partners = tuple(
                    int(value) for value in fields[offset + 1 + len(proc):]
                )
                self.assertEqual(sorted(permutation), list(range(1, len(proc) + 1)))
                self.assertEqual(len(flattened_partners), len(proc) * npartners)
                position += 1
            position += 3

    def test_process_file_header_marks_permutation_metadata(self):
        saved_options = dict(process_list.options)
        try:
            process_list.options["include_resonance"] = False
            header = process_list.WriteUniqueProcsIntoList(
                [("g", "g", "g", "g")], []
            )[0]
            self.assertEqual(header, "4 1 3 0 -1 -1 -1 -1")
        finally:
            process_list.options.clear()
            process_list.options.update(saved_options)

    def test_four_field_conversion_uses_identity_maps(self):
        proc = ("g", "g", "g", "g")
        line = process_list.ConvertProcToString(
            (proc, (0, 1, 2, 3), (0,), 2.0)
        ).split()
        # One partner, process, order, factor, one current map, one partner map.
        self.assertEqual(len(line), 1 + 1 + 4 + 4 + 1 + 4 + 4)
        self.assertEqual(tuple(map(int, line[-8:])), (1, 2, 3, 4) * 2)


if __name__ == "__main__":
    unittest.main()
