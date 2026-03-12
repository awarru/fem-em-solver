"""Example: MRI-like two-coil + gelled saline phantom field workflow.

This example demonstrates an end-to-end coarse workflow:
1. Build coil + phantom + air mesh
2. Apply phantom material properties
3. Solve magnetostatic B and frequency-domain proxy E fields
4. Print phantom-region field metrics
5. Export ParaView-ready outputs
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core import (
    HomogeneousMaterial,
    MagnetostaticProblem,
    MagnetostaticSolver,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.io.paraview_utils import write_combined_paraview_output
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.post import (
    build_phantom_quicklook_report,
    compute_phantom_eb_metrics_and_export,
    evaluate_vector_field_parallel,
    format_phantom_quicklook_report,
    write_phantom_quicklook_report,
)
from fem_em_solver.utils.constants import MU_0


def _print_tag_summary(cell_tags, comm: MPI.Intracomm) -> None:
    """Print deterministic global cell counts for core domain tags."""
    if cell_tags is None:
        if comm.rank == 0:
            print("[mesh] WARNING: no cell tags found")
        return

    tag_names = {1: "coil_1", 2: "coil_2", 3: "phantom", 4: "air"}
    if comm.rank == 0:
        print("[mesh] cell-tag summary:")

    for tag, name in tag_names.items():
        local_count = int(np.count_nonzero(cell_tags.values == tag))
        global_count = comm.allreduce(local_count, op=MPI.SUM)
        if comm.rank == 0:
            print(f"  tag {tag} ({name}): {global_count} cells")


RESOLUTION_PRESETS = {
    "coarse": 0.02,
    "medium": 0.015,
    "fine": 0.01,
}

SCENARIO_PRESETS = {
    "benchmark-lite": {
        "resolution_preset": "fine",
        "centerline_sample_count": 17,
        "ksp_max_it": 450,
        "frequency_probe_scale_factors": [0.995, 1.0, 1.005],
    },
    "debug": {
        "resolution_preset": "coarse",
        "centerline_sample_count": 5,
        "ksp_max_it": 180,
        "frequency_probe_scale_factors": [1.0],
    },
    "dev": {
        "resolution_preset": "medium",
        "centerline_sample_count": 9,
        "ksp_max_it": 300,
        "frequency_probe_scale_factors": [0.999, 1.0, 1.001],
    },
}


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MRI-like coil+phantom field workflow and export diagnostics.",
    )
    parser.add_argument(
        "--frequency-hz",
        type=float,
        default=127.74e6,
        help="Drive frequency in Hz (default: 127.74e6).",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(SCENARIO_PRESETS),
        default="debug",
        help="Scenario preset controlling mesh/detail/sweep diagnostics.",
    )
    parser.add_argument(
        "--resolution-preset",
        choices=sorted(RESOLUTION_PRESETS),
        default=None,
        help="Optional mesh resolution override for the selected scenario preset.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("paraview_output"),
        help="Directory for ParaView and phantom metric exports.",
    )
    return parser.parse_args(argv)


def _resolve_git_commit_hash() -> str:
    """Return current git commit hash, or 'unknown' when unavailable."""
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "HEAD"], text=True)
            .strip()
        )
    except Exception:
        return "unknown"


def _resolve_scenario(args: argparse.Namespace) -> dict[str, object]:
    """Build effective scenario config from preset + optional CLI overrides."""
    preset_config = dict(SCENARIO_PRESETS[args.preset])

    effective_resolution_preset = (
        args.resolution_preset
        if args.resolution_preset is not None
        else str(preset_config["resolution_preset"])
    )

    if effective_resolution_preset not in RESOLUTION_PRESETS:
        raise ValueError(
            f"Unsupported resolution preset: {effective_resolution_preset!r}. "
            f"Expected one of {sorted(RESOLUTION_PRESETS)}"
        )

    frequency_probe_scale_factors = [
        float(scale) for scale in preset_config["frequency_probe_scale_factors"]
    ]
    frequency_probe_scale_factors = sorted(frequency_probe_scale_factors)

    return {
        "name": args.preset,
        "resolution_preset": effective_resolution_preset,
        "resolution_m": float(RESOLUTION_PRESETS[effective_resolution_preset]),
        "centerline_sample_count": int(preset_config["centerline_sample_count"]),
        "ksp_max_it": int(preset_config["ksp_max_it"]),
        "frequency_probe_scale_factors": frequency_probe_scale_factors,
    }


def _write_output_manifest(
    *,
    output_dir: Path,
    args: argparse.Namespace,
    scenario: dict[str, object],
    current_a: float,
    mesh_params: dict[str, float],
    written_files: dict[str, str],
    metrics: dict,
    quicklook_artifacts: dict[str, str | None] | None = None,
) -> Path:
    """Write a reproducible run manifest capturing config + produced artifacts."""
    artifact_paths = [Path(path) for path in written_files.values()]
    phantom_metrics_json = output_dir / "mri_coil_phantom_phantom_metrics.json"
    phantom_e_csv = output_dir / "mri_coil_phantom_phantom_E_samples.csv"
    phantom_b_csv = output_dir / "mri_coil_phantom_phantom_B_samples.csv"
    artifact_paths.extend([phantom_metrics_json, phantom_e_csv, phantom_b_csv])

    if quicklook_artifacts is not None:
        for maybe_path in quicklook_artifacts.values():
            if maybe_path is not None:
                artifact_paths.append(Path(maybe_path))

    manifest_payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _resolve_git_commit_hash(),
        "parameters": {
            "frequency_hz": float(args.frequency_hz),
            "scenario_preset": str(scenario["name"]),
            "resolution_preset": str(scenario["resolution_preset"]),
            "resolution_m": float(scenario["resolution_m"]),
            "frequency_probe_hz": [
                float(args.frequency_hz) * float(scale)
                for scale in scenario["frequency_probe_scale_factors"]
            ],
            "output_dir": str(output_dir),
            "current_per_coil_a": float(current_a),
            "mesh": mesh_params,
        },
        "artifacts": [
            {
                "name": path.name,
                "path": str(path),
                "exists": path.exists(),
            }
            for path in sorted(artifact_paths, key=lambda p: p.name)
        ],
        "key_metrics": {
            "phantom_metrics_json": str(phantom_metrics_json),
            "e_magnitude_mean": float(metrics["E_magnitude"]["mean"]),
            "b_magnitude_mean": float(metrics["B_magnitude"]["mean"]),
            "e_to_b_mean_ratio": float(metrics["consistency"]["e_to_b_mean_ratio"]),
            "consistency_warning_count": int(len(metrics["consistency"]["warnings"])),
        },
    }

    manifest_path = output_dir / "mri_coil_phantom_manifest.json"
    manifest_path.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n")
    return manifest_path


def main(argv: list[str] | None = None):
    """Run MRI coil + phantom field workflow."""
    comm = MPI.COMM_WORLD
    args = _parse_args(argv)

    print("=" * 72)
    print("Example: MRI coil + gelled saline phantom fields")
    print("=" * 72)

    # Coarse VPS-safe geometry/mesh settings.
    current_a = 1.0
    frequency_hz = args.frequency_hz
    scenario = _resolve_scenario(args)
    resolution = float(scenario["resolution_m"])

    if frequency_hz <= 0.0:
        raise ValueError(f"--frequency-hz must be positive; received {frequency_hz}")

    probe_frequencies = [
        frequency_hz * float(scale)
        for scale in scenario["frequency_probe_scale_factors"]
    ]

    print("\nParameters:")
    print(f"  Current per driven coil region: {current_a:.3f} A")
    print(f"  Frequency: {frequency_hz:.6e} Hz")
    print(f"  Scenario preset: {scenario['name']}")
    print(f"  Mesh resolution preset: {scenario['resolution_preset']}")
    print(f"  Mesh resolution: {resolution:.4f} m")
    print(f"  Centerline sample count: {scenario['centerline_sample_count']}")
    print(f"  Solver ksp_max_it: {scenario['ksp_max_it']}")
    print(
        "  Frequency probe points (Hz): "
        + ", ".join(f"{freq:.6e}" for freq in probe_frequencies)
    )
    print(f"  Output directory: {args.output_dir}")

    mesh_params = {
        "coil_major_radius": 0.08,
        "coil_minor_radius": 0.01,
        "coil_separation": 0.08,
        "phantom_radius": 0.04,
        "phantom_height": 0.10,
        "air_padding": 0.04,
        "resolution": resolution,
    }

    print("\nGenerating coil + phantom mesh...")
    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        **mesh_params,
        comm=comm,
    )

    num_cells = mesh.topology.index_map(mesh.topology.dim).size_global
    num_vertices = mesh.topology.index_map(0).size_global
    if comm.rank == 0:
        print(f"  Mesh created: {num_cells} cells, {num_vertices} vertices")
    _print_tag_summary(cell_tags, comm)

    # Shared current source for both coil tags.
    coil_cross_section = np.pi * (0.01 ** 2)
    j_magnitude = current_a / coil_cross_section

    def current_density(_x):
        return ufl.as_vector([0.0, 0.0, j_magnitude])

    # Magnetostatics for B field.
    print("\nSolving magnetostatic field (A, B)...")
    mag_problem = MagnetostaticProblem(mesh=mesh, cell_tags=cell_tags, facet_tags=facet_tags, mu=MU_0)
    mag_solver = MagnetostaticSolver(mag_problem, degree=1)
    a_field = mag_solver.solve(current_density=current_density, subdomain_ids=[1, 2], gauge_penalty=1e-3)
    b_field = mag_solver.compute_b_field()
    if comm.rank == 0:
        print("  Magnetostatic solve complete")

    # Frequency-domain proxy E field with phantom material override.
    print("\nSolving frequency-domain proxy E field...")
    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0)
    phantom = GelledSalinePhantomMaterial(
        sigma=0.72,
        epsilon_r=76.5,
        frequency_hz=frequency_hz,
        mu_r=1.0,
    )

    th_problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=frequency_hz,
        material=background,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        phantom_material=phantom,
        phantom_tag=3,
        solver_petsc_options={
            "ksp_type": "gmres",
            "pc_type": "jacobi",
            "ksp_rtol": 1e-8,
            "ksp_max_it": int(scenario["ksp_max_it"]),
        },
        collect_solver_diagnostics=True,
    )
    th_solver = TimeHarmonicSolver(th_problem, degree=1)
    th_fields = th_solver.solve(current_density=current_density, subdomain_ids=[1, 2], gauge_penalty=1e-3)
    e_field = th_fields.e_imag
    if comm.rank == 0:
        print("  Time-harmonic solve complete (using E_imag as reported E field)")
        if th_fields.solve_diagnostics is not None:
            d = th_fields.solve_diagnostics
            print("  solve health diagnostics:")
            print(
                "    "
                f"ksp={d.ksp_type}, pc={d.pc_type}, converged={d.converged} "
                f"(reason={d.converged_reason}), iterations={d.iterations}"
            )
            print(
                "    "
                f"residual_norm={d.residual_norm:.6e}, residual_trend={d.residual_trend}, "
                f"history_samples={len(d.residual_history)}"
            )

    # Phantom metrics + phantom-only CSV/JSON exports.
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = compute_phantom_eb_metrics_and_export(
        e_field,
        b_field,
        cell_tags,
        phantom_tag=3,
        output_dir=output_dir,
        basename="mri_coil_phantom",
        comm=comm,
    )

    quicklook = build_phantom_quicklook_report(metrics)
    quicklook_artifacts = {"json": None, "markdown": None}
    if comm.rank == 0:
        quicklook_artifacts = write_phantom_quicklook_report(
            quicklook,
            output_dir=output_dir,
            basename="mri_coil_phantom",
            write_markdown=True,
            write_json=True,
        )

    # Interpolate to Lagrange for XDMF compatibility, then export combined outputs.
    v_lagrange = fem.functionspace(mesh, ("Lagrange", 1, (3,)))

    a_lagrange = fem.Function(v_lagrange, name="A")
    a_lagrange.interpolate(a_field)

    b_lagrange = fem.Function(v_lagrange, name="B")
    b_lagrange.interpolate(b_field)

    e_lagrange = fem.Function(v_lagrange, name="E")
    e_lagrange.interpolate(e_field)

    written_files = write_combined_paraview_output(
        output_dir=output_dir,
        basename="mri_coil_phantom_fields",
        mesh=mesh,
        cell_tags=cell_tags,
        fields={
            "A": (a_field, a_lagrange),
            "B": (b_field, b_lagrange),
            "E": (e_field, e_lagrange),
        },
        comm=comm,
    )

    manifest_path = None
    if comm.rank == 0:
        manifest_path = _write_output_manifest(
            output_dir=output_dir,
            args=args,
            scenario=scenario,
            current_a=current_a,
            mesh_params=mesh_params,
            written_files=written_files,
            metrics=metrics,
            quicklook_artifacts=quicklook_artifacts,
        )

    # Centerline diagnostics through phantom.
    z_line = np.linspace(-0.045, 0.045, int(scenario["centerline_sample_count"]))
    sample_points = np.zeros((len(z_line), 3), dtype=np.float64)
    sample_points[:, 2] = z_line
    e_samples, e_valid = evaluate_vector_field_parallel(e_lagrange, sample_points, comm=comm)
    b_samples, b_valid = evaluate_vector_field_parallel(b_lagrange, sample_points, comm=comm)

    if comm.rank == 0:
        print("\n" + format_phantom_quicklook_report(quicklook))
        print("\nPhantom diagnostics:")
        print(
            "  |E| min/max/mean: "
            f"{metrics['E_magnitude']['min']:.6e} / "
            f"{metrics['E_magnitude']['max']:.6e} / "
            f"{metrics['E_magnitude']['mean']:.6e}"
        )
        print(
            "  |B| min/max/mean: "
            f"{metrics['B_magnitude']['min']:.6e} / "
            f"{metrics['B_magnitude']['max']:.6e} / "
            f"{metrics['B_magnitude']['mean']:.6e}"
        )
        print(
            f"  centerline samples valid (E/B): "
            f"{np.count_nonzero(e_valid)}/{len(e_valid)} / {np.count_nonzero(b_valid)}/{len(b_valid)}"
        )
        consistency = metrics["consistency"]
        print("  consistency diagnostics:")
        print(
            f"    |E|/|B| mean ratio: {consistency['e_to_b_mean_ratio']:.6e} "
            f"(max ratio: {consistency['e_to_b_max_ratio']:.6e})"
        )
        print(
            f"    span ratios (E/B): {consistency['e_span_ratio']:.6f} / {consistency['b_span_ratio']:.6f}"
        )
        print(f"    mean-balance relative diff: {consistency['mean_balance_rel_diff']:.6f}")
        if consistency["warnings"]:
            print("    warnings:")
            for warning in consistency["warnings"]:
                print(f"      - {warning}")
        else:
            print("    warnings: none")
        print("\nCenterline sample magnitudes (z, |E|, |B|):")
        for i, z in enumerate(z_line):
            e_mag = np.linalg.norm(e_samples[i])
            b_mag = np.linalg.norm(b_samples[i])
            print(f"  z={z:+.4f} m -> |E|={e_mag:.6e}, |B|={b_mag:.6e}")

        print("\nParaView output:")
        for name, path in written_files.items():
            print(f"  {name}: {Path(path).name}")
        print("  phantom metrics json: mri_coil_phantom_phantom_metrics.json")
        print("  phantom E csv: mri_coil_phantom_phantom_E_samples.csv")
        print("  phantom B csv: mri_coil_phantom_phantom_B_samples.csv")
        if quicklook_artifacts["json"] is not None:
            print(f"  quick-look json: {Path(quicklook_artifacts['json']).name}")
        if quicklook_artifacts["markdown"] is not None:
            print(f"  quick-look markdown: {Path(quicklook_artifacts['markdown']).name}")
        if manifest_path is not None:
            print(f"  manifest json: {manifest_path.name}")

    print("\n" + "=" * 72)
    print("Example completed")
    print("=" * 72)


if __name__ == "__main__":
    main()
