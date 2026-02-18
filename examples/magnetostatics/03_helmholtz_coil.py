"""
Example: Magnetic field of a Helmholtz coil.

This example demonstrates the magnetostatic solver for a Helmholtz coil
(two loops separated by one radius) and shows the highly uniform field
in the central region.
"""

import numpy as np
from mpi4py import MPI
from pathlib import Path

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0

# Import dolfinx I/O for ParaView output
from dolfinx import io, fem


def main():
    """Run Helmholtz coil example."""
    comm = MPI.COMM_WORLD
    
    print("=" * 60)
    print("Example: Helmholtz coil magnetic field")
    print("=" * 60)
    
    # Problem parameters (optimized for faster meshing)
    current = 1.0              # Current in each loop [A]
    loop_radius = 0.05         # Loop radius [m] (5 cm)
    wire_radius = 0.002        # Wire cross-section [m] (2 mm, increased for faster meshing)
    domain_radius = 0.12       # Domain radius [m] (reduced from 0.15)
    resolution = 0.008         # Mesh resolution [m] (coarser: 8mm vs 5mm)
    
    separation = loop_radius   # Helmholtz condition
    
    print(f"\nParameters:")
    print(f"  Current: {current} A (in each loop)")
    print(f"  Loop radius: {loop_radius} m ({loop_radius*100:.1f} cm)")
    print(f"  Separation: {separation} m (Helmholtz condition)")
    print(f"  Wire radius: {wire_radius} m ({wire_radius*1000:.1f} mm)")
    print(f"  Loops at z = ±{separation/2*100:.1f} cm")
    
    # Generate mesh
    print("\nGenerating mesh...")
    mesh, cell_tags, facet_tags = MeshGenerator.helmholtz_coil_domain(
        loop_radius=loop_radius,
        wire_radius=wire_radius,
        domain_radius=domain_radius,
        resolution=resolution,
        comm=comm
    )
    print(f"  Mesh created: {mesh.topology.index_map(3).size_global} cells")
    
    # Set up problem
    problem = MagnetostaticProblem(
        mesh=mesh, 
        cell_tags=cell_tags,
        mu=MU_0
    )
    
    # Create solver
    solver = MagnetostaticSolver(problem, degree=1)
    
    # Define current density in wires
    wire_cross_section = np.pi * wire_radius**2
    J_magnitude = current / wire_cross_section
    
    import ufl
    def current_density(x):
        """Uniform current density."""
        return ufl.as_vector([0.0, 0.0, J_magnitude])
    
    # Solve with current restricted to wire subdomain (tag=1)
    print("\nSolving magnetostatic problem...")
    A = solver.solve(current_density=current_density, subdomain_id=1)
    print("  Solution computed!")
    
    # Compute B-field
    print("\nComputing B-field...")
    B = solver.compute_b_field()
    
    # Evaluate along z-axis
    n_points = 40  # Reduced from 50 for faster evaluation
    z_eval = np.linspace(-0.08, 0.08, n_points)  # Adjusted for smaller domain
    
    points = np.zeros((n_points, 3))
    points[:, 2] = z_eval  # z positions along axis
    
    B_num = B.eval(points, np.arange(n_points))
    B_num_z = B_num[:, 2]  # z-component only
    
    # Analytical solution
    B_ana_z = AnalyticalSolutions.helmholtz_coil_field_on_axis(
        z_eval, current, loop_radius
    )
    
    # Error metrics
    rel_error = ErrorMetrics.l2_relative_error(B_num_z, B_ana_z)
    max_error = ErrorMetrics.max_relative_error(B_num_z, B_ana_z)
    
    print(f"\nResults:")
    print(f"  B_z at center (numerical): {B_num_z[n_points//2]:.6e} T")
    print(f"  B_z at center (analytical): {B_ana_z[n_points//2]:.6e} T")
    
    # Expected B_z at center
    B_center_expected = (4/5)**(3/2) * MU_0 * current / loop_radius
    print(f"  Expected B_z(0) = (4/5)^(3/2) * μ₀I/R = {B_center_expected:.6e} T")
    
    print(f"\n  Relative L2 error: {rel_error:.4%}")
    print(f"  Max relative error: {max_error:.4%}")
    
    # Compute magnetic energy
    energy = solver.compute_magnetic_energy()
    print(f"\n  Magnetic energy: {energy:.6e} J")
    
    # Test uniformity in central region
    central_idx = (np.abs(z_eval) <= 0.005)  # Central 1 cm
    if np.any(central_idx):
        B_central = B_num_z[central_idx]
        mean_B = np.mean(B_central)
        std_B = np.std(B_central)
        cv = std_B / abs(mean_B)
        print(f"\n  Uniformity in central 1 cm:")
        print(f"    Mean B_z: {mean_B:.6e} T")
        print(f"    Std B_z: {std_B:.6e} T")
        print(f"    Coefficient of variation: {cv:.4%}")

    # =========================================================================
    # Save results for ParaView visualization
    # =========================================================================
    print("\n" + "=" * 60)
    print("Saving results for ParaView visualization...")
    print("=" * 60)

    # Create output directory
    output_dir = Path("paraview_output")
    output_dir.mkdir(exist_ok=True)

    # Method 1: XDMF format (traditional, widely compatible)
    # Note: XDMF only supports Lagrange elements, so we need to interpolate
    print("\n  Writing XDMF files...")

    # Save the mesh with cell tags for visualization
    with io.XDMFFile(comm, output_dir / "helmholtz_coil_mesh.xdmf", "w") as xdmf:
        xdmf.write_mesh(mesh)
        # Also write cell tags to visualize regions
        if cell_tags is not None:
            mesh.topology.create_connectivity(3, 3)
            xdmf.write_meshtags(cell_tags, mesh.geometry)
        print("    ✓ Mesh saved to helmholtz_coil_mesh.xdmf (with cell tags)")

    # Create Lagrange function space for visualization
    V_lag = fem.functionspace(mesh, ("Lagrange", 1, (3,)))

    # Interpolate A to Lagrange space
    A_lag = fem.Function(V_lag, name="A")
    A_lag.interpolate(A)

    with io.XDMFFile(comm, output_dir / "helmholtz_coil_A.xdmf", "w") as xdmf:
        xdmf.write_mesh(mesh)
        # Write cell tags for filtering in ParaView
        if cell_tags is not None:
            xdmf.write_meshtags(cell_tags, mesh.geometry)
        xdmf.write_function(A_lag)
        print("    ✓ Vector potential A saved to helmholtz_coil_A.xdmf (with cell tags)")

    # Interpolate B to Lagrange space
    B_lag = fem.Function(V_lag, name="B")
    B_lag.interpolate(B)

    with io.XDMFFile(comm, output_dir / "helmholtz_coil_B.xdmf", "w") as xdmf:
        xdmf.write_mesh(mesh)
        # Write cell tags so we can filter in ParaView
        if cell_tags is not None:
            xdmf.write_meshtags(cell_tags, mesh.geometry)
        xdmf.write_function(B_lag)
        print("    ✓ Magnetic field B saved to helmholtz_coil_B.xdmf (with cell tags for filtering)")

    # Method 2: VTX format (modern, supports higher-order elements)
    print("\n  Writing VTX files...")
    try:
        vtx_A = io.VTXWriter(comm, output_dir / "helmholtz_coil_A.bp", [A], engine="BP4")
        vtx_A.write(0.0)
        vtx_A.close()
        print("    ✓ Vector potential A saved to helmholtz_coil_A.bp/")

        vtx_B = io.VTXWriter(comm, output_dir / "helmholtz_coil_B.bp", [B], engine="BP4")
        vtx_B.write(0.0)
        vtx_B.close()
        print("    ✓ Magnetic field B saved to helmholtz_coil_B.bp/")
    except Exception as e:
        print(f"    ⚠ VTX output failed (ADIOS2 may not be available): {e}")

    print("\n  ✓ ParaView files saved to paraview_output/")
    print("    Open helmholtz_coil_B.xdmf in ParaView to visualize")
    print("    Tip: The highly uniform field in the center is the key feature!")

    # Save results (text format)
    if comm.rank == 0:
        print("\n  Saving results...")
        data = np.column_stack([z_eval, B_num_z, B_ana_z, 
                                 np.abs(B_num_z - B_ana_z)])
        np.savetxt('helmholtz_coil_results.txt', data, 
                   header='z[m] Bz_num[T] Bz_ana[T] error[T]',
                   fmt='%.6e')
        print("  Results saved to: helmholtz_coil_results.txt")
        
        # Print some values
        print("\n  Sample values:")
        print(f"    {'z [m]':<12} {'B_z num [T]':<15} {'B_z ana [T]':<15} {'Error':<10}")
        print("    " + "-" * 52)
        for i in [0, n_points//4, n_points//2, 3*n_points//4, n_points-1]:
            err_pct = 100 * abs(B_num_z[i] - B_ana_z[i]) / abs(B_ana_z[i])
            print(f"    {z_eval[i]:<12.4f} {B_num_z[i]:<15.6e} "
                  f"{B_ana_z[i]:<15.6e} {err_pct:<10.2f}%")
    
    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)
    
    return solver, B


if __name__ == "__main__":
    main()
