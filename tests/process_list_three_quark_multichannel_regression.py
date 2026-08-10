#!/usr/bin/env python3
"""Check three-quark-line phase-space multichannel metadata."""

from __future__ import annotations

import sys
import unittest
from collections import Counter, defaultdict
from math import isclose
from pathlib import Path


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


class ThreeQuarkLineMultiChannelRegression(unittest.TestCase):
    def assert_valid_three_line_partners(self, rows):
        rows_by_channel = defaultdict(list)
        for channel, row in rows:
            rows_by_channel[channel].append(row)

        rows_by_coefficient = defaultdict(Counter)
        factors_by_coefficient = defaultdict(lambda: defaultdict(float))
        for channel, (proc, order, partners, factor) in rows:
            self.assertEqual(tuple(sorted(set(partners))), partners)
            self.assertIn(channel, partners)
            signature = process_list.ThreeQuarkLineBlockSignature(proc, order)
            coefficient = (partners, proc, signature)
            rows_by_coefficient[coefficient][channel] += 1
            factors_by_coefficient[coefficient][channel] += factor
            swapped_proc, swapped_order = process_list.SwapThreeQuarkLineOrder(
                proc, order
            )
            self.assertEqual(swapped_proc, proc)
            swapped_channel = process_list.process_order_to_index[
                (swapped_proc, swapped_order)
            ]
            same_coefficient = (
                process_list.ThreeQuarkLineBlockSignature(proc, order)
                == process_list.ThreeQuarkLineBlockSignature(
                    swapped_proc, swapped_order
                )
            )
            if same_coefficient:
                self.assertIn(swapped_channel, partners)
            for partner in partners:
                self.assertTrue(
                    any(
                        partner_proc == proc
                        and process_list.ThreeQuarkLineBlockSignature(
                            partner_proc, partner_order
                        )
                        == signature
                        for partner_proc, partner_order, _, _ in rows_by_channel[
                            partner
                        ]
                    )
                )
        for coefficient, channel_counts in rows_by_coefficient.items():
            partners = coefficient[0]
            self.assertEqual(set(channel_counts), set(partners))
            channel_factors = factors_by_coefficient[coefficient]
            reference = next(iter(channel_factors.values()))
            self.assertTrue(
                all(
                    isclose(value, reference, rel_tol=1.0e-12)
                    for value in channel_factors.values()
                )
            )

    def test_distinct_scattering_uses_both_orderings(self):
        rows = finalized_rows(("d~", "u~", "d", "u", "s", "s~"))
        self.assertEqual(len(rows), 12)
        self.assertEqual(
            Counter((len(row[2]), row[3]) for _, row in rows),
            Counter({(2, 1.0): 12}),
        )
        self.assertEqual(len(process_list.all_keys_sorted), 12)
        self.assert_valid_three_line_partners(rows)

    def test_distinguishable_final_labels_are_not_reordered(self):
        rows = finalized_rows(("d~", "d", "u", "u~", "s", "s~"))
        self.assertEqual(
            Counter((len(row[2]), row[3]) for _, row in rows),
            Counter({(2, 1.0): 12}),
        )
        self.assertEqual(len(process_list.all_keys_sorted), 12)
        self.assert_valid_three_line_partners(rows)

    def test_identical_flavours_preserve_orbit_weights(self):
        for proc in (
            ("d~", "d", "u", "u~", "u", "u~"),
            ("u~", "u", "u", "u~", "u", "u~"),
        ):
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(
                    Counter((len(row[2]), row[3]) for _, row in rows),
                    Counter({(1, 2.0): 3}),
                )
                self.assert_valid_three_line_partners(rows)

    def test_qq_initial_crossing_uses_partner_maps(self):
        cases = (
            (
                ("u~", "u~", "u", "u", "d", "d~"),
                6,
                Counter({(2, 2.0): 4, (1, 1.0): 2}),
            ),
            (
                ("u~", "u~", "u", "u", "d", "d~", "a"),
                18,
                Counter({(6, 2.0): 12, (3, 1.0): 6}),
            ),
        )
        for proc, expected_rows, expected_metadata in cases:
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(rows), expected_rows)
                self.assertEqual(
                    Counter(
                        (len(row[2]), round(row[3], 12)) for _, row in rows
                    ),
                    expected_metadata,
                )
                self.assert_valid_three_line_partners(rows)

    def test_identical_flavours_compose_with_singlets(self):
        rows = finalized_rows(("u~", "u", "u", "u~", "u", "u~", "a"))
        self.assertEqual(
            Counter((len(row[2]), round(row[3], 12)) for _, row in rows),
            Counter({(3, 2.0): 9}),
        )
        self.assert_valid_three_line_partners(rows)

    def test_singlet_channels_compose_with_string_orderings(self):
        cases = (
            (("d~", "u~", "d", "u", "s", "s~", "a"), 6, 36),
            (("d~", "u~", "d", "u", "s", "s~", "a", "z"), 24, 144),
        )
        for proc, expected_partners, expected_rows in cases:
            with self.subTest(proc=proc):
                rows = finalized_rows(proc)
                self.assertEqual(len(rows), expected_rows)
                self.assertEqual(
                    {(len(row[2]), round(row[3], 12)) for _, row in rows},
                    {(expected_partners, 1.0)},
                )
                self.assert_valid_three_line_partners(rows)

    def test_gluons_move_with_their_open_strings(self):
        rows = finalized_rows(("d~", "u~", "d", "u", "s", "s~", "g"))
        self.assertEqual(len(rows), 36)
        self.assertEqual({(len(row[2]), row[3]) for _, row in rows}, {(2, 1.0)})
        self.assert_valid_three_line_partners(rows)

    def test_identical_gluons_are_paired_only_when_strings_are_unchanged(self):
        rows = finalized_rows(
            ("d~", "u~", "d", "u", "s", "s~", "g", "g")
        )
        self.assertEqual(len(rows), 72)
        self.assertEqual(
            Counter((len(row[2]), round(row[3], 12)) for _, row in rows),
            Counter({(2, 2.0): 60, (1, 1.0): 12}),
        )
        self.assert_valid_three_line_partners(rows)

    def test_gluon_singlet_maps_keep_a_common_external_labelling(self):
        rows = finalized_rows(
            ("d~", "d", "u", "u~", "s", "s~", "g", "z")
        )
        metadata = Counter((len(row[2]), round(row[3], 12)) for _, row in rows)
        self.assertEqual(metadata, Counter({(6, 1.0): 108}))
        self.assertEqual(len(process_list.all_keys_sorted), 108)
        self.assert_valid_three_line_partners(rows)

    def test_two_line_process_is_unchanged(self):
        rows = finalized_rows(("d~", "u~", "d", "u", "g", "g"))
        self.assertEqual(len(rows), 6)
        self.assertEqual(
            {(len(row[2]), row[3]) for _, row in rows}, {(1, 2.0)}
        )


if __name__ == "__main__":
    unittest.main()
