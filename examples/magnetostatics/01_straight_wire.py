"""
Example: Magnetic field of a straight wire.

This example demonstrates the magnetostatic solver by computing
the magnetic field around a current-carrying wire and comparing
with the analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
from pathlib import Path

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0
from fem_em_solver.post import evaluate_vector_field_parallel
from fem_em_solver.io.paraview_utils import write_combined_paraview_output

# Import dolfinx I/O for ParaView output
from dolfinx import io, fem


def main():
    """Run straight wire example."""
    comm = MPI.COMM_WORLD
    
    print("=" * 60)
    print("Example: Magnetic field of straight wire")
    print("=" * 60)
    
    # Problem parameters (balanced for accuracy and speed)
    current = 1.0              # Current [A]
    wire_length = 0.3          # Wire length [m]
    domain_radius = 0.04       # Domain radius [m] (4 cm)
    resolution = 0.01          # Mesh resolution [m] (coarse, cron-safe runtime)
    wire_radius = 0.003       # Wire radius [m] (1.5 mm)
    
    print(f"\nParameters:")
    print(f"  Current: {current} A")
    print(f"  Wire length: {wire_length} m")
    print(f"  Wire radius: {wire_radius} m")
    print(f"  Domain radius: {domain_radius} m")
    print(f"  Mesh resolution: {resolution} m")
    
    # Generate mesh
    print("\nGenerating mesh...")
    mesh, cell_tags, facet_tags = MeshGenerator.straight_wire_domain(
        wire_length=wire_length,
        wire_radius=wire_radius,
        domain_radius=domain_radius,
        resolution=resolution,
        comm=comm
    )

    # Diagnostic: Check mesh properties
    num_cells = mesh.topology.index_map(3).size_global
    num_vertices = mesh.topology.index_map(0).size_global
    print(f"  Mesh created: {num_cells} cells, {num_vertices} vertices")

    # Check cell tags
    if cell_tags is not None:
        unique_tags = np.unique(cell_tags.values)
        print(f"  Cell tags found: {unique_tags}")
        for tag in unique_tags:
            count = np.sum(cell_tags.values == tag)
            if tag == 1:
                print(f"    Tag {tag} (wire): {count} cells")
            elif tag == 2:
                print(f"    Tag {tag} (air/domain): {count} cells")
    else:
        print("  WARNING: No cell tags found!")
    
    # Set up problem with cell tags for subdomain integration
    problem = MagnetostaticProblem(mesh=mesh, cell_tags=cell_tags, mu=MU_0)
    
    # Create solver
    solver = MagnetostaticSolver(problem, degree=1)
    
    # Define current density in wire
    wire_area = np.pi * wire_radius**2
    J_magnitude = current / wire_area
    
    import ufl
    def current_density(x):
        """Uniform current density in z-direction."""
        return ufl.as_vector([0.0, 0.0, J_magnitude])
    
    # Solve with current restricted to wire subdomain (tag=1)
    print("\nSolving magnetostatic problem...")
    A = solver.solve(current_density=current_density, subdomain_id=1)
    print("  Solution computed (current restricted to wire volume)!")
    
    # Compute B-field
    print("\nComputing B-field...")
    B = solver.compute_b_field()

    # Interpolate B to Lagrange space for evaluation
    # (DG functions are discontinuous and need proper cell indices for eval)
    print("  Interpolating B to Lagrange space for evaluation...")
    V_lag = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
    B_lag = fem.Function(V_lag, name="B")
    B_lag.interpolate(B)

    # Evaluate along x-axis from wire edge to domain boundary
    n_points = 30
    r_eval = np.linspace(wire_radius, domain_radius * 0.95, n_points)  # From wire edge to near boundary

    print(f"\nEvaluating B-field along x-axis:")
    print(f"  From r = {r_eval[0]*1000:.2f} mm (wire edge) to r = {r_eval[-1]*100:.2f} cm")
    print(f"  Wire radius: {wire_radius*1000:.2f} mm, Domain radius: {domain_radius*100:.2f} cm")
    print(f"  Number of evaluation points: {n_points}")

    points = np.zeros((n_points, 3))
    points[:, 0] = r_eval  # x positions (radial distance from wire center)
    points[:, 1] = 0.0     # y = 0
    points[:, 2] = 0.0     # z = 0 (middle of wire)

    print(f"\n  DEBUG - Evaluation points:")
    print(f"  First 3 points: {points[:3]}")
    print(f"  Last 3 points: {points[-3:]}")

    # Evaluate points robustly across MPI partitions
    B_num, valid_mask = evaluate_vector_field_parallel(B_lag, points, comm=comm)

    invalid_count = np.count_nonzero(~valid_mask)
    print(f"\n  DEBUG - Valid evaluation points: {np.count_nonzero(valid_mask)}/{n_points}")
    if invalid_count > 0:
        print(f"  WARNING: {invalid_count} points were outside mesh partitions")
    print(f"\n  DEBUG - Numerical B-field evaluation:")
    print(f"  B_num shape: {B_num.shape}, dtype: {B_num.dtype}")
    print(f"  First point: B = ({B_num[0,0]:.6e}, {B_num[0,1]:.6e}, {B_num[0,2]:.6e}) T")
    print(f"  Last point:  B = ({B_num[-1,0]:.6e}, {B_num[-1,1]:.6e}, {B_num[-1,2]:.6e}) T")

    # Compute magnitude
    B_num_mag = np.linalg.norm(B_num, axis=1)
    print(f"  |B| at first point: {B_num_mag[0]:.6e} T")
    print(f"  |B| at last point:  {B_num_mag[-1]:.6e} T")
    print(f"  Ratio (should be ~{r_eval[-1]/r_eval[0]:.2f}): {B_num_mag[0]/B_num_mag[-1]:.2f}")

    # Analytical solution (wire at origin, matching mesh)
    wire_position_analytical = np.array([0.0, 0.0])  # Must match mesh generation at (x=0, y=0)
    B_ana = AnalyticalSolutions.straight_wire_magnetic_field(points, current, wire_position_analytical)

    print(f"\n  DEBUG - Analytical B-field:")
    print(f"  Wire position: ({wire_position_analytical[0]}, {wire_position_analytical[1]})")
    print(f"  First point: B = ({B_ana[0,0]:.6e}, {B_ana[0,1]:.6e}, {B_ana[0,2]:.6e}) T")
    print(f"  Last point:  B = ({B_ana[-1,0]:.6e}, {B_ana[-1,1]:.6e}, {B_ana[-1,2]:.6e}) T")

    B_ana_mag = np.linalg.norm(B_ana, axis=1)
    print(f"  |B| at first point: {B_ana_mag[0]:.6e} T")
    print(f"  |B| at last point:  {B_ana_mag[-1]:.6e} T")
    print(f"  Ratio (expected ~{r_eval[-1]/r_eval[0]:.2f}): {B_ana_mag[0]/B_ana_mag[-1]:.2f}")

    # Expected field from formula: B = μ₀I/(2πr)
    mu_0 = 4 * np.pi * 1e-7
    B_expected_first = mu_0 * current / (2 * np.pi * r_eval[0])
    B_expected_last = mu_0 * current / (2 * np.pi * r_eval[-1])
    print(f"\n  Expected from μ₀I/2πr formula:")
    print(f"  At r={r_eval[0]*1000:.2f}mm: {B_expected_first:.6e} T")
    print(f"  At r={r_eval[-1]*1000:.2f}mm: {B_expected_last:.6e} T")

    # Debug: Component analysis - field at (x,0,0) should be (0, By, 0)
    print(f"\n  Component verification at (x, 0, 0):")
    print(f"  Expected: Bx≈0, By=μ₀I/2πx, Bz≈0")
    print(f"  {'r [mm]':>10} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12} {'|B|':>12}")
    for i in [0, n_points//2, n_points-1]:
        print(f"  {r_eval[i]*1000:10.2f} {B_num[i,0]:12.6e} {B_num[i,1]:12.6e} {B_num[i,2]:12.6e} {B_num_mag[i]:12.6e}")

    # Debug: Print first and last few values to verify 1/r decay
    print(f"\n  Field magnitude verification (should decay as 1/r):")
    print(f"  {'r [mm]':>10} {'|B_num| [T]':>15} {'|B_ana| [T]':>15} {'Ratio':>10}")
    for i in [0, 1, 2, n_points//2, n_points-3, n_points-2, n_points-1]:
        ratio = B_num_mag[i] / B_ana_mag[i] if B_ana_mag[i] > 0 else 0
        print(f"  {r_eval[i]*1000:10.2f} {B_num_mag[i]:15.6e} {B_ana_mag[i]:15.6e} {ratio:10.4f}")

    # Verify 1/r scaling: B(r1)/B(r2) should equal r2/r1
    r1, r2 = r_eval[0], r_eval[-1]
    expected_ratio = r2 / r1
    actual_num_ratio = B_num_mag[0] / B_num_mag[-1]
    actual_ana_ratio = B_ana_mag[0] / B_ana_mag[-1]
    print(f"\n  1/r decay check:")
    print(f"  Expected B(r_min)/B(r_max) = {expected_ratio:.2f}")
    print(f"  Numerical: {actual_num_ratio:.2f}")
    print(f"  Analytical: {actual_ana_ratio:.2f}")
    
    # Error metrics
    rel_error = ErrorMetrics.l2_relative_error(B_num_mag, B_ana_mag)
    max_error = ErrorMetrics.max_relative_error(B_num_mag, B_ana_mag)
    
    print(f"\nResults:")
    print(f"  Max B-field (numerical): {np.max(B_num_mag):.6e} T")
    print(f"  Max B-field (analytical): {np.max(B_ana_mag):.6e} T")
    print(f"  Relative L2 error: {rel_error:.4%}")
    print(f"  Max relative error: {max_error:.4%}")
    
    # Compute magnetic energy
    energy = solver.compute_magnetic_energy()
    print(f"\n  Magnetic energy: {energy:.6e} J")

    # =========================================================================
    # Save results for ParaView visualization
    # =========================================================================
    print("\n" + "=" * 60)
    print("Saving results for ParaView visualization...")
    print("=" * 60)

    # Create output directory
    output_dir = Path("paraview_output")
    output_dir.mkdir(exist_ok=True)

    # Standardized XDMF exports (individual + combined tag/field output)
    print("\n  Writing XDMF files (traditional + combined)...")

    # V_lag and B_lag already created earlier for evaluation
    # Interpolate A to Lagrange space
    A_lag = fem.Function(V_lag, name="A")
    A_lag.interpolate(A)

    try:
        written_files = write_combined_paraview_output(
            output_dir=output_dir,
            basename="straight_wire",
            mesh=mesh,
            cell_tags=cell_tags,
            fields={
                "A": (A, A_lag),
                "B": (B, B_lag),
            },
            comm=comm,
        )
        if comm.rank == 0:
            print("    ✓ Standardized XDMF export complete")
            if "combined" in written_files:
                print(f"    ✓ Combined file saved to {written_files['combined'].name}")
                print("    This file has BOTH cell tags and fields on the same grid!")
    except Exception as e:
        import traceback
        print(f"    ⚠ XDMF export failed: {e}")
        print("\n  Full traceback:")
        traceback.print_exc()

    # Method 3: VTX format (modern, supports higher-order elements)
    # This uses ADIOS2 and creates .bp directory with data
    print("\n  Writing VTX files (modern ADIOS2 format)...")
    try:
        # VTXWriter for the vector potential A (Nedelec H(curl) space)
        vtx_A = io.VTXWriter(comm, output_dir / "straight_wire_A.bp", [A], engine="BP4")
        vtx_A.write(0.0)  # time = 0.0 for static problem
        vtx_A.close()
        print("    ✓ Vector potential A saved to straight_wire_A.bp/")

        # VTXWriter for the magnetic field B (DG space)
        vtx_B = io.VTXWriter(comm, output_dir / "straight_wire_B.bp", [B], engine="BP4")
        vtx_B.write(0.0)
        vtx_B.close()
        print("    ✓ Magnetic field B saved to straight_wire_B.bp/")
    except Exception as e:
        print(f"    ⚠ VTX output failed (ADIOS2 may not be available): {e}")
        print("    Note: XDMF files were still created and can be used instead")

    print("\n" + "=" * 60)
    print("ParaView Instructions:")
    print("=" * 60)
    print("\n  RECOMMENDED: Use the combined file!")
    print("    File -> Open -> straight_wire_combined.xdmf")
    print("    - Has BOTH cell tags AND B-field on same grid")
    print("    - Cell tags available in Threshold filter!")
    print("    - B field available for Glyph, Stream Tracer, etc.")
    print("\n  Filtering workflow:")
    print("    1. Apply Threshold filter:")
    print("       - Scalars: 'Cell tags'")
    print("       - Min: 2, Max: 2 (removes wire cells)")
    print("    2. Apply Glyph filter to thresholded data:")
    print("       - Orientation/Scale: B")
    print("    3. Enjoy clutter-free visualization!")
    print("\n  Alternative: Individual files")
    print("    - straight_wire_A.xdmf (vector potential)")
    print("    - straight_wire_B.xdmf (magnetic field)")
    print("    - straight_wire_mesh.xdmf (mesh with tags)")
    print("\n  For VTX files (if available):")
    print("    - Open straight_wire_B.bp/ directory")
    print("=" * 60)
    
    # Plot results
    if comm.rank == 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # B-field magnitude comparison
        ax = axes[0]
        ax.semilogy(r_eval * 1000, B_num_mag, 'b-o', label='Numerical (FEM)', markersize=4)
        ax.semilogy(r_eval * 1000, B_ana_mag, 'r--', label=f'Analytical (μ₀I/2πr)', linewidth=2)
        ax.set_xlabel('Radial Distance from Wire Center [mm]')
        ax.set_ylabel('|B| [T] (log scale)')
        ax.set_title(f'Magnetic Field Magnitude vs Distance\n(Should decay as 1/r)')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        # Relative error
        ax = axes[1]
        rel_err_pointwise = np.abs(B_num_mag - B_ana_mag) / B_ana_mag
        ax.semilogy(r_eval * 1000, rel_err_pointwise * 100, 'g-s', markersize=4)
        ax.set_xlabel('Distance from wire [mm]')
        ax.set_ylabel('Relative Error [%]')
        ax.set_title('Pointwise Relative Error')
        ax.axhline(y=rel_error * 100, color='r', linestyle='--', 
                   label=f'L2 error: {rel_error:.2%}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('straight_wire_validation.png', dpi=150, bbox_inches='tight')
        print(f"\n  Plot saved to: straight_wire_validation.png")
        plt.close()
    
    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)
    
    return solver, B


if __name__ == "__main__":
    main()
