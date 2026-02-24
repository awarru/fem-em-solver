"""Touchstone export/load tests (chunk E5)."""

from __future__ import annotations

import numpy as np

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

    frequencies_hz, s_matrices, port_ids = load_touchstone(touchstone_path)

    assert port_ids == ("P1", "P2")
    np.testing.assert_allclose(frequencies_hz, np.array([100.0e6, 110.0e6]))
    np.testing.assert_allclose(s_matrices[0], sweep_results[0].s_matrix)
    np.testing.assert_allclose(s_matrices[1], sweep_results[1].s_matrix)
