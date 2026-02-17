"""
Convergence study tests for magnetostatic solver.

These tests verify that the numerical solution converges to the analytical
solution as the mesh is refined (h-refinement) or polynomial degree is
increased (p-refinement).
"""

import pytest
import numpy as np
from mpi4py import MPI

# Import solver components
from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


class TestConvergence:
    """Convergence study tests for magnetostatic solver."""
    
    def test_h_refinement_straight_wire(self):
        """
        Test h-convergence (mesh refinement) for straight wire problem.
        
        As mesh size h decreases, error should decrease at expected rate.
        For linear Nedelec elements, expect rate ~1.0.
        """
        comm = MPI.COMM_WORLD
        
        # Problem parameters
        current = 1.0
        wire_length = 0.5
        wire_radius = 0.001
        
        # Mesh resolutions to test (reduced for speed)
        resolutions = [0.02, 0.01]
        errors = []
        
        # Evaluation points (along x-axis, away from wire)
        n_points = 10
        r_eval = np.linspace(0.005, 0.04, n_points)
        points = np.zeros((n_points, 3))
        points[:, 0] = r_eval
        points[:, 2] = 0.0
        
        # Analytical solution
        B_ana = AnalyticalSolutions.straight_wire_magnetic_field(points, current)
        B_ana_mag = np.linalg.norm(B_ana, axis=1)
        
        import ufl
        
        for res in resolutions:
            if comm.rank == 0:
                print(f"\n  Testing resolution h = {res}...")
            
            # Generate mesh
            mesh, cell_tags, _ = MeshGenerator.straight_wire_domain(
                wire_length=wire_length,
                wire_radius=wire_radius,
                domain_radius=0.05,
                resolution=res,
                comm=comm
            )
            
            # Solve
            problem = MagnetostaticProblem(mesh=mesh, cell_tags=cell_tags, mu=MU_0)
            solver = MagnetostaticSolver(problem, degree=1)
            
            wire_area = np.pi * wire_radius**2
            J_magnitude = current / wire_area
            
            def current_density(x):
                return ufl.as_vector([0.0, 0.0, J_magnitude])
            
            solver.solve(current_density=current_density, subdomain_id=1)
            B = solver.compute_b_field()
            
            # Evaluate
            B_num = B.eval(points, np.arange(n_points))
            B_num_mag = np.linalg.norm(B_num, axis=1)
            
            # Error
            rel_error = ErrorMetrics.l2_relative_error(B_num_mag, B_ana_mag)
            errors.append(rel_error)
            
            if comm.rank == 0:
                print(f"    Relative L2 error: {rel_error:.4%}")
        
        # Compute convergence rate
        # Fit log(error) = log(C) - rate * log(h)
        log_h = np.log(resolutions)
        log_err = np.log(errors)
        
        # Linear regression
        n = len(resolutions)
        rate = -np.sum((log_h - np.mean(log_h)) * (log_err - np.mean(log_err))) / np.sum((log_h - np.mean(log_h))**2)
        
        if comm.rank == 0:
            print(f"\n  Convergence Results:")
            print(f"    Resolutions (h): {resolutions}")
            print(f"    Errors: {[f'{e:.4%}' for e in errors]}")
            print(f"    Convergence rate: {rate:.2f}")
            print(f"    Expected rate for linear elements: ~1.0")
        
        # Assert reasonable convergence
        assert rate > 0.5, f"Convergence rate {rate:.2f} too low (expected > 0.5)"
    
    def test_p_refinement_straight_wire(self):
        """
        Test p-convergence (polynomial degree) for straight wire problem.
        
        As polynomial degree increases, error should decrease.
        
        Steps:
        1. Fixed mesh resolution (e.g., 0.01)
        2. Loop over degrees: [1, 2, 3]
        3. For each degree:
           - Create solver with degree=N
           - Solve and compute error
        4. Assert error decreases with higher degree
        
        TODO: Implement this test
        """
        pytest.skip("Not yet implemented - Chunk 7")
    
    def test_convergence_data_export(self):
        """
        Export convergence data for visualization.
        
        Save h vs error and degree vs error to files in results/convergence/
        
        TODO: Implement data export
        """
        pytest.skip("Not yet implemented - Chunk 8")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
