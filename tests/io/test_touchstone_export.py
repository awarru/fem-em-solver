"""Touchstone export/load tests (chunk D5)."""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.ports import PortDefinition
from fem_em_solver.ports.sparameters import SParameterSweepResult
from fem_em_solver.ports.touchstone import export_touchstone, load_touchstone


def test_touchstone_export_and_roundtrip_loader(tmp_path):
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw", z0_ohm=50.0),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw", z0_ohm=50.0),
    ]

    sweep_results = [
        SParameterSweepResult(
            frequency_hz=100.0e6,
            port_ids=("P1", "P2"),
            s_matrix=np.array(
                [
                    [0.10 + 0.01j, 0.20 - 0.02j],
                    [0.30 + 0.03j, 0.40 - 0.04j],
                ],
                dtype=np.complex128,
            ),
            excitation_results={},
        ),
        SParameterSweepResult(
            frequency_hz=110.0e6,
            port_ids=("P1", "P2"),
            s_matrix=np.array(
                [
                    [0.11 + 0.015j, 0.21 - 0.025j],
                    [0.31 + 0.035j, 0.41 - 0.045j],
                ],
                dtype=np.complex128,
            ),
            excitation_results={},
        ),
    ]

    touchstone_path = export_touchstone(
        sweep_results,
        ports,
        output_dir=tmp_path,
        write_csv_companion=True,
    )

    assert touchstone_path.exists()
    assert touchstone_path.suffix.lower() == ".s2p"
    assert "100.000MHz" in touchstone_path.name
    assert "110.000MHz" in touchstone_path.name
    assert "Z0_50.00" in touchstone_path.name

    csv_path = touchstone_path.with_suffix(".csv")
    assert csv_path.exists()

    touchstone_text = touchstone_path.read_text(encoding="utf-8")
    assert "! generated_utc:" in touchstone_text
    assert "! port_order: P1,P2" in touchstone_text
    assert "! frequency_points_hz: 1.000000000e+08,1.100000000e+08" in touchstone_text
    assert "! z0_ohm: 50.000000" in touchstone_text
    assert "# Hz S RI R 50.000000" in touchstone_text

    frequencies_hz, s_matrices, port_ids = load_touchstone(touchstone_path)

    assert port_ids == ("P1", "P2")
    np.testing.assert_allclose(frequencies_hz, np.array([100.0e6, 110.0e6]))
    np.testing.assert_allclose(s_matrices[0], sweep_results[0].s_matrix)
    np.testing.assert_allclose(s_matrices[1], sweep_results[1].s_matrix)


def test_touchstone_loader_rejects_frequency_metadata_mismatch(tmp_path):
    bad_touchstone_path = tmp_path / "bad_freq.s2p"
    bad_touchstone_path.write_text(
        "\n".join(
            [
                "! FEM-EM Solver Touchstone export",
                "! generated_utc: 2026-03-09T17:00:00Z",
                "! port_order: P1,P2",
                "! frequency_points_hz: 1.000000000e+08,1.200000000e+08",
                "! z0_ohm: 50.000000",
                "# Hz S RI R 50.000000",
                "1.000000000e+08 1.0e-01 0.0e+00 2.0e-01 0.0e+00 3.0e-01 0.0e+00 4.0e-01 0.0e+00",
                "1.100000000e+08 1.1e-01 0.0e+00 2.1e-01 0.0e+00 3.1e-01 0.0e+00 4.1e-01 0.0e+00",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="frequency_points_hz metadata does not match numeric data rows"):
        load_touchstone(bad_touchstone_path)


def test_touchstone_loader_rejects_z0_metadata_mismatch(tmp_path):
    bad_touchstone_path = tmp_path / "bad_z0.s2p"
    bad_touchstone_path.write_text(
        "\n".join(
            [
                "! FEM-EM Solver Touchstone export",
                "! generated_utc: 2026-03-09T17:00:00Z",
                "! port_order: P1,P2",
                "! frequency_points_hz: 1.000000000e+08,1.100000000e+08",
                "! z0_ohm: 75.000000",
                "# Hz S RI R 50.000000",
                "1.000000000e+08 1.0e-01 0.0e+00 2.0e-01 0.0e+00 3.0e-01 0.0e+00 4.0e-01 0.0e+00",
                "1.100000000e+08 1.1e-01 0.0e+00 2.1e-01 0.0e+00 3.1e-01 0.0e+00 4.1e-01 0.0e+00",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="z0_ohm metadata does not match Touchstone option line"):
        load_touchstone(bad_touchstone_path)
