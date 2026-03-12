"""Tests for quick-look phantom metrics reporting helpers (chunk E3)."""

from __future__ import annotations

from pathlib import Path

from fem_em_solver.post.quicklook import (
    build_phantom_quicklook_report,
    format_phantom_quicklook_report,
    write_phantom_quicklook_report,
)


def _sample_metrics() -> dict:
    return {
        "sampling": {
            "requested_cells": 20,
            "sampling_cells": 16,
            "valid_sample_cells": 16,
            "boundary_adjacent_cells_dropped": 4,
            "invalid_samples_dropped_e": 0,
            "invalid_samples_dropped_b": 0,
        },
        "E_magnitude": {
            "count": 16,
            "min": 1.0,
            "max": 3.0,
            "mean": 2.0,
        },
        "B_magnitude": {
            "count": 16,
            "min": 0.2,
            "max": 0.8,
            "mean": 0.4,
        },
        "consistency": {
            "e_to_b_mean_ratio": 5.0,
            "e_to_b_max_ratio": 6.5,
            "mean_balance_rel_diff": 0.12,
            "warnings": [],
        },
    }


def test_build_phantom_quicklook_report_ok_status_without_warning_flags():
    report = build_phantom_quicklook_report(_sample_metrics())

    assert report["status"] == "OK"
    assert report["warning_flags"] == []
    assert report["consistency_warnings"] == []
    assert report["summary"]["e_mean"] == 2.0
    assert report["summary"]["b_mean"] == 0.4


def test_build_phantom_quicklook_report_warn_status_when_consistency_warnings_present():
    metrics = _sample_metrics()
    metrics["consistency"]["warnings"] = ["Large E/B imbalance detected"]

    report = build_phantom_quicklook_report(metrics)

    assert report["status"] == "WARN"
    assert report["warning_flags"] == []
    assert report["consistency_warnings"] == ["Large E/B imbalance detected"]


def test_format_and_write_phantom_quicklook_report_produce_expected_artifacts(tmp_path: Path):
    report = build_phantom_quicklook_report(_sample_metrics())
    text = format_phantom_quicklook_report(report)

    assert "quick-look phantom metrics:" in text
    assert "status: OK" in text
    assert "warning flags: none" in text

    outputs = write_phantom_quicklook_report(
        report,
        output_dir=tmp_path,
        basename="mri_coil_phantom",
        write_json=True,
        write_markdown=True,
    )

    json_path = Path(outputs["json"])
    md_path = Path(outputs["markdown"])

    assert json_path.exists()
    assert md_path.exists()

    md_text = md_path.read_text(encoding="utf-8")
    assert "quick-look phantom metrics:" in md_text
    assert "status: OK" in md_text
