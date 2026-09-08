#!/usr/bin/env python3
"""Unit tests for selecting whole query-target pairs for plotting."""
from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

# The selection helper has no SW dependency, but plot_sw imports the package at
# module load time. Keep this unit test runnable in the lightweight test env.
if "ginfinity_sw" not in sys.modules:
    fake_sw = types.ModuleType("ginfinity_sw")
    fake_sw.ScoringParameters = object
    fake_sw.align = object()
    fake_sw.similarity_matrix = object()
    fake_sw.transform_scores = object()
    sys.modules["ginfinity_sw"] = fake_sw

from plot_r4rna import selected_pair_rows as r4_selected_pair_rows  # noqa: E402
from plot_rnartistcore import selected_pair_rows as rnartist_selected_pair_rows  # noqa: E402
from plot_sw import selected_pair_rows as sw_selected_pair_rows  # noqa: E402


class TestPlotPairSelection(unittest.TestCase):
    def test_limit_uses_merged_pair_score_and_retains_deduplicated_hsps(self) -> None:
        rows = [
            {"query_id": "q1", "target_id": "low", "cluster_id": "low", "score": "5", "query_start": "0", "query_end": "10", "target_start": "0", "target_end": "10"},
            {"query_id": "q1", "target_id": "high", "cluster_id": "high-1", "score": "8", "query_start": "20", "query_end": "30", "target_start": "20", "target_end": "30"},
            {"query_id": "q1", "target_id": "high", "cluster_id": "high-2", "score": "7", "query_start": "40", "query_end": "50", "target_start": "40", "target_end": "50"},
            {"query_id": "q1", "target_id": "high", "cluster_id": "duplicate", "score": "6", "query_start": "20", "query_end": "30", "target_start": "20", "target_end": "30"},
            {"query_id": "q1", "target_id": "middle", "cluster_id": "middle", "score": "10", "query_start": "60", "query_end": "70", "target_start": "60", "target_end": "70"},
        ]
        expected = [(1, "high-1"), (2, "high-2")]

        for selector in (
            rnartist_selected_pair_rows,
            r4_selected_pair_rows,
            sw_selected_pair_rows,
        ):
            selected = selector(rows, 1)
            self.assertEqual(
                [(index, row["cluster_id"]) for index, row in selected],
                expected,
            )


if __name__ == "__main__":
    raise SystemExit(unittest.main())
