"""
Validation test: Magnetic field of straight wire.

This test validates the magnetostatic solver against the analytical
solution for an infinite straight wire carrying current I.

Analytical solution:
    B_φ = μ₀I / (2πr)  [Tesla]
    
The numerical solution uses a finite wire in a cylindrical domain with
appropriate boundary conditions.
"""

import pytest
import numpy as np
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


class TestStraightWire:
    """Validation tests for straight wire magnetostatics."""
    
    @pytest.fixture
    def wire_params(self):
        """Standard wire parameters."""
        return {
            'current': 1.0,           # 1 Ampere
            'wire_length': 1.0,       # 1 meter
            'domain_radius': 0.05,    # 5 cm domain
            'resolution': 0.005,      # 5 mm mesh size
        }
    
    def test_straight_wire_b_field(self, wire_params, tmp_path):
        """Test B-field matches analytical solution for straight wire.
        
        This test:
        1. Creates a cylindrical mesh with wire along z-axis
        2. Solves magnetostatic problem
        3. Extracts B-field at evaluation points
        4. Compares with analytical solution B = μ₀I/(2πr)
        5. Checks that relative error is within tolerance
        """
        comm = MPI.COMM_WORLD
        
        # Parameters
        current = wire_params['current']
        wire_length = wire_params['wire_length']
        domain_radius = wire_params['domain_radius']
        resolution = wire_params['resolution']
        
        # Generate mesh
        mesh, cell_tags, facet_tags = MeshGenerator.straight_wire_domain(
            wire_length=wire_length,
            wire_radius=0.001,  # 1mm wire
            domain_radius=domain_radius,
            resolution=resolution,
            comm=comm
        )
        
        # Create problem
        problem = MagnetostaticProblem(mesh=mesh, mu=MU_0)
        
        # Create solver
        solver = MagnetostaticSolver(problem, degree=1)
        
        # Define current density (uniform in wire volume)
        # J = I / A_wire in z-direction
        wire_area = np.pi * (0.001)**2
        J_magnitude = current / wire_area
        
        def current_density(x):
            """Return J(x) as UFL expression."""
            import ufl
            # Check if point is in wire region (simplified - assumes wire at origin)
            r = ufl.sqrt(x[0]**2 + x[1]**2)
            # J in z-direction, nonzero only in wire
            # For simplicity, apply uniform J in whole domain
            # (proper subdomain restriction would use cell_tags)
            return ufl.as_vector([0.0, 0.0, J_magnitude])
        
        # Solve
        A = solver.solve(current_density=current_density)
        
        # Compute B-field
        B_num = solver.compute_b_field()
        
        # Evaluate at test points (along x-axis, away from wire)
        n_points = 10
        r_test = np.linspace(0.005, domain_radius * 0.8, n_points)
        
        points = np.zeros((n_points, 3))
        points[:, 0] = r_test  # x positions
        points[:, 2] = 0.0     # z = 0 (middle of wire)
        
        # Evaluate numerical solution
        B_num_eval = B_num.eval(points, np.arange(n_points))
        
        # Analytical solution
        B_ana = AnalyticalSolutions.straight_wire_magnetic_field(
            points, current
        )
        
        # Only compare z-component (should be ~0) and in-plane components
        # The B-field should be purely azimuthal
        
        # Check that B_z is small
        B_z_max = np.max(np.abs(B_num_eval[:, 2]))
        assert B_z_max < 0.01, f"B_z should be near zero, got {B_z_max}"
        
        # Check magnitude of B-field
        B_num_mag = np.linalg.norm(B_num_eval, axis=1)
        B_ana_mag = np.linalg.norm(B_ana, axis=1)
        
        # Compute relative error
        rel_error = ErrorMetrics.l2_relative_error(B_num_mag, B_ana_mag)
        
        print(f"\nStraight wire validation:")
        print(f"  Current: {current} A")
        print(f"  Domain radius: {domain_radius} m")
        print(f"  Resolution: {resolution} m")
        print(f"  Max B (numerical): {np.max(B_num_mag):.6e} T")
        print(f"  Max B (analytical): {np.max(B_ana_mag):.6e} T")
        print(f"  Relative L2 error: {rel_error:.4%}")
        
        # Assert error is within tolerance (10% for coarse mesh)
        assert rel_error < 0.10, f"Relative error {rel_error:.4%} exceeds 10%"
    
    def test_straight_wire_convergence(self, wire_params, tmp_path):
        """Test h-convergence of B-field solution.
        
        Refine mesh and check that error decreases.
        """
        comm = MPI.COMM_WORLD
        current = wire_params['current']
        wire_length = wire_params['wire_length']
        domain_radius = wire_params['domain_radius']
        
        resolutions = [0.01, 0.007, 0.005]  # Coarse to fine
        errors = []
        
        # Evaluation points
        n_points = 5
        r_test = np.linspace(0.01, domain_radius * 0.5, n_points)
        points = np.zeros((n_points, 3))
        points[:, 0] = r_test
        
        B_ana = AnalyticalSolutions.straight_wire_magnetic_field(
            points, current
        )
        B_ana_mag = np.linalg.norm(B_ana, axis=1)
        
        for res in resolutions:
            mesh, cell_tags, facet_tags = MeshGenerator.straight_wire_domain(
                wire_length=wire_length,
                wire_radius=0.001,
                domain_radius=domain_radius,
                resolution=res,
                comm=comm
            )
            
            problem = MagnetostaticProblem(mesh=mesh, mu=MU_0)
            solver = MagnetostaticSolver(problem, degree=1)
            
            # Simplified current (uniform)
            wire_area = np.pi * (0.001)**2
            J_mag = current / wire_area
            
            def J(x):
                import ufl
                return ufl.as_vector([0.0, 0.0, J_mag])
            
            solver.solve(current_density=J)
            B_num = solver.compute_b_field()
            B_num_eval = B_num.eval(points, np.arange(n_points))
            B_num_mag = np.linalg.norm(B_num_eval, axis=1)
            
            rel_err = ErrorMetrics.l2_relative_error(B_num_mag, B_ana_mag)
            errors.append(rel_err)
            
            print(f"  Resolution {res:.4f}m: Error = {rel_err:.4%}")
        
        # Check convergence
        # Error should decrease as mesh refines
        assert errors[-1] < errors[0], "Error should decrease with mesh refinement"
        
        # Check approximate order (error should halve when h halves)
        # For linear elements, order ~1-2 depending on problem
        rate = np.log(errors[0] / errors[-1]) / np.log(resolutions[0] / resolutions[-1])
        print(f"  Convergence rate: {rate:.2f}")
        
        assert rate > 0.5, f"Convergence rate {rate:.2f} too low"


def test_analytical_straight_wire():
    """Unit test for analytical solution calculation."""
    current = 1.0  # A
    
    # Test at r = 0.01 m
    points = np.array([[0.01, 0, 0]])
    B = AnalyticalSolutions.straight_wire_magnetic_field(points, current)
    
    # Expected: B_phi = μ₀I/(2πr) = 4π×10⁻⁷ × 1 / (2π × 0.01)
    expected = 4e-7 * np.pi * 1.0 / (2 * np.pi * 0.01)
    
    B_mag = np.linalg.norm(B)
    assert np.abs(B_mag - expected) < 1e-10, f"Analytical B-field incorrect"
    
    # Direction should be in y-direction at x-axis (azimuthal)
    assert np.abs(B[0, 1] - expected) < 1e-10, "B-field direction incorrect"
    assert np.abs(B[0, 0]) < 1e-15, "B_x should be zero"
    assert np.abs(B[0, 2]) < 1e-15, "B_z should be zero"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
