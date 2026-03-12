"""Quick-look reporting helpers for MRI phantom E/B diagnostics."""

from __future__ import annotations

from pathlib import Path
import json
import math
from typing import Any


def _is_finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def build_phantom_quicklook_report(metrics: dict[str, Any]) -> dict[str, Any]:
    """Build a compact, warning-oriented quick-look report from phantom metrics."""
    e_stats = metrics["E_magnitude"]
    b_stats = metrics["B_magnitude"]
    consistency = metrics.get("consistency", {})
    sampling = metrics.get("sampling", {})

    warning_flags: list[str] = []

    for label, stats in (("E", e_stats), ("B", b_stats)):
        if int(stats.get("count", 0)) <= 0:
            warning_flags.append(f"{label}_NO_SAMPLES")
        for key in ("min", "max", "mean"):
            if not _is_finite_number(stats.get(key)):
                warning_flags.append(f"{label}_{key.upper()}_NON_FINITE")

    valid_cells = int(sampling.get("valid_sample_cells", 0))
    sampling_cells = int(sampling.get("sampling_cells", 0))
    if sampling_cells > 0 and valid_cells < sampling_cells:
        warning_flags.append("SAMPLE_COVERAGE_PARTIAL")

    for key, label in (
        ("e_to_b_mean_ratio", "E_TO_B_MEAN_RATIO_NON_FINITE"),
        ("e_to_b_max_ratio", "E_TO_B_MAX_RATIO_NON_FINITE"),
        ("mean_balance_rel_diff", "MEAN_BALANCE_REL_DIFF_NON_FINITE"),
    ):
        if key in consistency and not _is_finite_number(consistency[key]):
            warning_flags.append(label)

    consistency_warnings = list(consistency.get("warnings", []))
    status = "WARN" if warning_flags or consistency_warnings else "OK"

    return {
        "status": status,
        "warning_flags": warning_flags,
        "consistency_warnings": consistency_warnings,
        "summary": {
            "e_min": float(e_stats["min"]),
            "e_max": float(e_stats["max"]),
            "e_mean": float(e_stats["mean"]),
            "b_min": float(b_stats["min"]),
            "b_max": float(b_stats["max"]),
            "b_mean": float(b_stats["mean"]),
            "e_to_b_mean_ratio": float(consistency.get("e_to_b_mean_ratio", float("nan"))),
            "e_to_b_max_ratio": float(consistency.get("e_to_b_max_ratio", float("nan"))),
            "mean_balance_rel_diff": float(consistency.get("mean_balance_rel_diff", float("nan"))),
        },
        "sampling": {
            "requested_cells": int(sampling.get("requested_cells", 0)),
            "sampling_cells": sampling_cells,
            "valid_sample_cells": valid_cells,
            "boundary_adjacent_cells_dropped": int(sampling.get("boundary_adjacent_cells_dropped", 0)),
            "invalid_samples_dropped_e": int(sampling.get("invalid_samples_dropped_e", 0)),
            "invalid_samples_dropped_b": int(sampling.get("invalid_samples_dropped_b", 0)),
        },
    }


def format_phantom_quicklook_report(report: dict[str, Any]) -> str:
    """Format a deterministic console/markdown quick-look report."""
    s = report["summary"]
    sampling = report["sampling"]
    warning_flags = report.get("warning_flags", [])
    consistency_warnings = report.get("consistency_warnings", [])

    lines = [
        "quick-look phantom metrics:",
        f"  status: {report['status']}",
        (
            "  |E| min/max/mean: "
            f"{s['e_min']:.6e} / {s['e_max']:.6e} / {s['e_mean']:.6e}"
        ),
        (
            "  |B| min/max/mean: "
            f"{s['b_min']:.6e} / {s['b_max']:.6e} / {s['b_mean']:.6e}"
        ),
        (
            "  consistency ratios (mean/max): "
            f"{s['e_to_b_mean_ratio']:.6e} / {s['e_to_b_max_ratio']:.6e}"
        ),
        f"  mean-balance relative diff: {s['mean_balance_rel_diff']:.6f}",
        (
            "  sampling coverage (valid/sampling/requested): "
            f"{sampling['valid_sample_cells']}/{sampling['sampling_cells']}/{sampling['requested_cells']}"
        ),
        (
            "  dropped samples (boundary-adjacent, invalid E, invalid B): "
            f"{sampling['boundary_adjacent_cells_dropped']}, "
            f"{sampling['invalid_samples_dropped_e']}, "
            f"{sampling['invalid_samples_dropped_b']}"
        ),
    ]

    if warning_flags:
        lines.append("  warning flags:")
        lines.extend(f"    - {flag}" for flag in warning_flags)
    else:
        lines.append("  warning flags: none")

    if consistency_warnings:
        lines.append("  consistency warnings:")
        lines.extend(f"    - {warning}" for warning in consistency_warnings)
    else:
        lines.append("  consistency warnings: none")

    return "\n".join(lines)


def write_phantom_quicklook_report(
    report: dict[str, Any],
    *,
    output_dir: str | Path,
    basename: str,
    write_markdown: bool = True,
    write_json: bool = True,
) -> dict[str, str | None]:
    """Write optional markdown/json quick-look report artifacts."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, str | None] = {
        "json": None,
        "markdown": None,
    }

    if write_json:
        json_path = output_dir / f"{basename}_quicklook.json"
        json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        paths["json"] = str(json_path)

    if write_markdown:
        md_path = output_dir / f"{basename}_quicklook.md"
        md_path.write_text(format_phantom_quicklook_report(report) + "\n", encoding="utf-8")
        paths["markdown"] = str(md_path)

    return paths
