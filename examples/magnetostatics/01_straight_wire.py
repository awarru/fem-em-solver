"""
Example: Magnetic field of a straight wire.

This example demonstrates the magnetostatic solver by computing
the magnetic field around a current-carrying wire and comparing
with the analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


def main():
    """Run straight wire example."""
    comm = MPI.COMM_WORLD
    
    print("=" * 60)
    print("Example: Magnetic field of straight wire")
    print("=" * 60)
    
    # Problem parameters
    current = 1.0              # Current [A]
    wire_length = 0.5          # Wire length [m]
    domain_radius = 0.05       # Domain radius [m]
    resolution = 0.005         # Mesh resolution [m]
    wire_radius = 0.001        # Wire radius [m]
    
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
    print(f"  Mesh created: {mesh.topology.index_map(3).size_global} cells")
    
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
    
    # Evaluate along x-axis
    n_points = 20
    r_eval = np.linspace(0.002, domain_radius * 0.8, n_points)
    
    points = np.zeros((n_points, 3))
    points[:, 0] = r_eval  # x positions
    points[:, 2] = 0.0     # z = middle of wire
    
    B_num = B.eval(points, np.arange(n_points))
    B_num_mag = np.linalg.norm(B_num, axis=1)
    
    # Analytical solution
    B_ana = AnalyticalSolutions.straight_wire_magnetic_field(points, current)
    B_ana_mag = np.linalg.norm(B_ana, axis=1)
    
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
    
    # Plot results
    if comm.rank == 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # B-field magnitude comparison
        ax = axes[0]
        ax.semilogy(r_eval * 1000, B_num_mag, 'b-o', label='Numerical', markersize=4)
        ax.semilogy(r_eval * 1000, B_ana_mag, 'r--', label='Analytical', linewidth=2)
        ax.set_xlabel('Distance from wire [mm]')
        ax.set_ylabel('|B| [T]')
        ax.set_title('Magnetic Field Magnitude')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
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
