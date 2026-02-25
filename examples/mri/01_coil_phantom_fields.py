"""Example: MRI-like two-coil + gelled saline phantom field workflow.

This example demonstrates an end-to-end coarse workflow:
1. Build coil + phantom + air mesh
2. Apply phantom material properties
3. Solve magnetostatic B and frequency-domain proxy E fields
4. Print phantom-region field metrics
5. Export ParaView-ready outputs
"""

from __future__ import annotations

from pathlib import Path
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
    compute_phantom_eb_metrics_and_export,
    evaluate_vector_field_parallel,
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


def main():
    """Run MRI coil + phantom field workflow."""
    comm = MPI.COMM_WORLD

    print("=" * 72)
    print("Example: MRI coil + gelled saline phantom fields")
    print("=" * 72)

    # Coarse VPS-safe geometry/mesh settings.
    current_a = 1.0
    frequency_hz = 127.74e6
    resolution = 0.02

    print("\nParameters:")
    print(f"  Current per driven coil region: {current_a:.3f} A")
    print(f"  Frequency: {frequency_hz:.6e} Hz")
    print(f"  Mesh resolution: {resolution:.4f} m")

    print("\nGenerating coil + phantom mesh...")
    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=resolution,
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
    )
    th_solver = TimeHarmonicSolver(th_problem, degree=1)
    th_fields = th_solver.solve(current_density=current_density, subdomain_ids=[1, 2], gauge_penalty=1e-3)
    e_field = th_fields.e_imag
    if comm.rank == 0:
        print("  Time-harmonic solve complete (using E_imag as reported E field)")

    # Phantom metrics + phantom-only CSV/JSON exports.
    output_dir = Path("paraview_output")
    output_dir.mkdir(exist_ok=True)

    metrics = compute_phantom_eb_metrics_and_export(
        e_field,
        b_field,
        cell_tags,
        phantom_tag=3,
        output_dir=output_dir,
        basename="mri_coil_phantom",
        comm=comm,
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

    # Centerline diagnostics through phantom.
    z_line = np.linspace(-0.045, 0.045, 9)
    sample_points = np.zeros((len(z_line), 3), dtype=np.float64)
    sample_points[:, 2] = z_line
    e_samples, e_valid = evaluate_vector_field_parallel(e_lagrange, sample_points, comm=comm)
    b_samples, b_valid = evaluate_vector_field_parallel(b_lagrange, sample_points, comm=comm)

    if comm.rank == 0:
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

    print("\n" + "=" * 72)
    print("Example completed")
    print("=" * 72)


if __name__ == "__main__":
    main()
