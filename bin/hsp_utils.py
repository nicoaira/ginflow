"""Shared HSP de-duplication and merged-pair plot selection helpers."""
from __future__ import annotations

import math


Row = dict[str, str]
IndexedRow = tuple[int, Row]


def parse_float(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def hsp_score(row: Row) -> float | None:
    """Return the finite score used to rank one raw alignment row."""
    for field in ("score", "total_score"):
        value = parse_float(row.get(field), math.nan)
        if math.isfinite(value):
            return value
    return None


def _span(row: Row, axis: str) -> tuple[int, int] | None:
    start = parse_int(row.get(f"{axis}_start"))
    end = parse_int(row.get(f"{axis}_end"))
    if end <= start:
        return None
    return start, end


def _overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] < right[1] and right[0] < left[1]


def spans_overlap(left: Row, right: Row) -> bool:
    """Return whether two HSPs share query or target coordinates.

    Pair-level E-values are calibrated from disjoint local HSPs. Sharing either
    sequence span means the two alignments cannot both contribute independent
    evidence, even if their other span is different.
    """
    for axis in ("query", "target"):
        left_span = _span(left, axis)
        right_span = _span(right, axis)
        if left_span is not None and right_span is not None and _overlap(left_span, right_span):
            return True
    return False


def _signature(row: Row) -> tuple[str, ...]:
    """Capture the alignment identity independently of its seed cluster."""
    return tuple(
        str(row.get(field, ""))
        for field in (
            "query_start",
            "query_end",
            "target_start",
            "target_end",
            "query_aligned",
            "target_aligned",
        )
    )


def _priority(item: IndexedRow) -> tuple[object, ...]:
    index, row = item
    return (
        -(hsp_score(row) or 0.0),
        -parse_int(row.get("aligned_columns")),
        -parse_int(row.get("match_count")),
        parse_int(row.get("query_start")),
        parse_int(row.get("target_start")),
        parse_int(row.get("query_end")),
        parse_int(row.get("target_end")),
        str(row.get("cluster_id", "")),
        index,
    )


def deduplicate_indexed_hsps(items: list[IndexedRow]) -> list[IndexedRow]:
    """Keep the highest-scoring non-overlapping HSPs deterministically.

    Clusters have disjoint seed membership, but their padded SW crops can
    overlap and converge on the same local traceback. Greedy score order keeps
    the strongest evidence and prevents duplicate or overlapping alignments
    from being counted as independent HSPs.
    """
    kept: list[IndexedRow] = []
    seen: set[tuple[str, ...]] = set()
    for item in sorted(items, key=_priority):
        _index, row = item
        signature = _signature(row)
        if signature in seen:
            continue
        if any(spans_overlap(row, existing_row) for _kept_index, existing_row in kept):
            continue
        kept.append(item)
        seen.add(signature)
    return sorted(kept, key=lambda item: item[0])


def deduplicate_hsps(rows: list[Row]) -> list[Row]:
    """Return raw HSP rows with duplicate/overlapping evidence removed."""
    return [row for _index, row in deduplicate_indexed_hsps(list(enumerate(rows)))]


def deduplicate_hsps_by_pair(rows: list[Row]) -> list[Row]:
    """Deduplicate rows independently for each query-target pair."""
    groups: dict[tuple[str, str], list[IndexedRow]] = {}
    for index, row in enumerate(rows):
        key = (row.get("query_id", ""), row.get("target_id", ""))
        groups.setdefault(key, []).append((index, row))
    kept = [item for items in groups.values() for item in deduplicate_indexed_hsps(items)]
    return [row for _index, row in sorted(kept, key=lambda item: item[0])]


def select_merged_pair_rows(rows: list[Row], max_pairs: int) -> list[IndexedRow]:
    """Select HSPs from the highest-scoring merged query-target pairs.

    Plot processes receive raw cluster alignments so they can draw every HSP.
    The pair order is therefore reconstructed from the same deduplicated HSP
    scores used by MERGE_ALIGNMENTS. Within one query, total score and database
    E-value have the same ordering.
    """
    groups: dict[tuple[str, str], list[IndexedRow]] = {}
    for index, row in enumerate(rows):
        key = (row.get("query_id", ""), row.get("target_id", ""))
        groups.setdefault(key, []).append((index, row))

    ranked: list[tuple[float, float, str, str, int, list[IndexedRow]]] = []
    for (query_id, target_id), items in groups.items():
        kept = deduplicate_indexed_hsps(items)
        if not kept:
            continue
        scores = [hsp_score(row) or 0.0 for _index, row in kept]
        ranked.append(
            (
                sum(scores),
                max(scores),
                query_id,
                target_id,
                min(index for index, _row in kept),
                kept,
            )
        )

    ranked.sort(key=lambda item: (-item[0], -item[1], item[2], item[3], item[4]))
    selected: list[IndexedRow] = []
    for _total, _maximum, _query_id, _target_id, _first_index, items in ranked[: max(0, int(max_pairs))]:
        selected.extend(items)
    return selected
