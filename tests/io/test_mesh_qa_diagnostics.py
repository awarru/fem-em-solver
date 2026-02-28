"""Unit-level checks for mesh QA diagnostic formatting helpers."""

from __future__ import annotations

from fem_em_solver.io.mesh_qa import (
    format_cell_tag_summary,
    format_expected_tag_counts,
    print_required_tag_failure_summary,
)


def test_format_expected_tag_counts_reports_missing_and_present_tags() -> None:
    required = {1: "coil_1", 2: "coil_2", 3: "phantom"}
    counts = {1: 8, 3: 12, 4: 400}

    rendered = format_expected_tag_counts(counts, required)

    assert "coil_1(tag=1): expected>=1, actual=8 [OK]" in rendered
    assert "coil_2(tag=2): expected>=1, actual=0 [MISSING]" in rendered
    assert "phantom(tag=3): expected>=1, actual=12 [OK]" in rendered


def test_print_required_tag_failure_summary_includes_expected_and_observed(capsys) -> None:
    required = {1: "coil_1", 2: "coil_2", 3: "phantom"}
    counts = {1: 8, 3: 12, 4: 400}

    print_required_tag_failure_summary(counts, required, comm=None, prefix="[mesh-qa] ")
    out = capsys.readouterr().out

    assert "[mesh-qa] required-tag expected vs actual:" in out
    assert "coil_2(tag=2): expected>=1, actual=0 [MISSING]" in out
    assert "[mesh-qa] observed-tag summary:" in out
    assert "coil_1=8" in out
    assert "phantom=12" in out
    assert "tag_4=400" in out


def test_format_cell_tag_summary_empty_counts_is_explicit() -> None:
    assert format_cell_tag_summary({}, tag_names={1: "coil_1"}) == "(no tagged cells)"
