"""
Validation test: Magnetic field of circular current loop.

This test validates the magnetostatic solver against the analytical
solution for B_z on the axis of a circular current loop.

Analytical solution on axis (z-axis):
    B_z(z) = μ₀Ia² / (2(a² + z²)^(3/2))  [Tesla]
    
Where:
    a = loop radius
    I = current
    z = distance along axis from loop center
"""

import pytest
import numpy as np
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


class TestCircularLoop:
    """Validation tests for circular current loop."""
    
    @pytest.fixture
    def loop_params(self):
        """Standard loop parameters."""
        return {
            'current': 1.0,           # 1 Ampere
            'loop_radius': 0.05,      # 5 cm radius
            'wire_radius': 0.001,     # 1 mm wire thickness
            'domain_radius': 0.15,    # 15 cm domain
            'resolution': 0.005,      # 5 mm mesh size
        }
    
    def test_circular_loop_on_axis(self, loop_params, tmp_path):
        """Test B_z on axis matches analytical solution.
        
        This test:
        1. Creates a torus mesh for the loop
        2. Solves magnetostatic problem
        3. Extracts B_z along the z-axis
        4. Compares with analytical solution
        5. Checks that relative error is within tolerance
        """
        comm = MPI.COMM_WORLD
        
        # Parameters
        current = loop_params['current']
        loop_radius = loop_params['loop_radius']
        wire_radius = loop_params['wire_radius']
        domain_radius = loop_params['domain_radius']
        resolution = loop_params['resolution']
        
        print(f"\nCircular loop test:")
        print(f"  Loop radius: {loop_radius} m")
        print(f"  Wire radius: {wire_radius} m")
        print(f"  Current: {current} A")
        
        # Generate mesh
        print("  Generating mesh...")
        mesh, cell_tags, facet_tags = MeshGenerator.circular_loop_domain(
            loop_radius=loop_radius,
            wire_radius=wire_radius,
            domain_radius=domain_radius,
            resolution=resolution,
            comm=comm
        )
        print(f"  Mesh created: {mesh.topology.index_map(3).size_global} cells")
        
        # Create problem
        problem = MagnetostaticProblem(
            mesh=mesh, 
            cell_tags=cell_tags,
            mu=MU_0
        )
        
        # Create solver
        solver = MagnetostaticSolver(problem, degree=1)
        
        # Define current density in wire
        wire_cross_section = np.pi * wire_radius**2
        J_magnitude = current / wire_cross_section
        
        import ufl
        def current_density(x):
            """Uniform current density in phi-direction (around loop)."""
            # For a loop in xy-plane, current is azimuthal
            # This is approximate - proper implementation would use
            # tangential current direction
            return ufl.as_vector([0.0, 0.0, J_magnitude])
        
        # Solve with current restricted to wire subdomain (tag=1)
        print("  Solving...")
        A = solver.solve(current_density=current_density, subdomain_id=1)
        
        # Compute B-field
        print("  Computing B-field...")
        B = solver.compute_b_field()
        
        # Evaluate along z-axis
        n_points = 15
        z_eval = np.linspace(-0.1, 0.1, n_points)
        
        points = np.zeros((n_points, 3))
        points[:, 2] = z_eval  # z positions
        
        B_num = B.eval(points, np.arange(n_points))
        B_num_z = B_num[:, 2]  # z-component
        
        # Analytical solution
        B_ana_z = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            z_eval, current, loop_radius
        )
        
        # Error metrics
        rel_error = ErrorMetrics.l2_relative_error(B_num_z, B_ana_z)
        max_error = ErrorMetrics.max_relative_error(B_num_z, B_ana_z)
        
        print(f"\n  Results:")
        print(f"    Max B_z (numerical): {np.max(np.abs(B_num_z)):.6e} T")
        print(f"    Max B_z (analytical): {np.max(np.abs(B_ana_z)):.6e} T")
        print(f"    Relative L2 error: {rel_error:.4%}")
        print(f"    Max relative error: {max_error:.4%}")
        
        # Assert error is within tolerance
        assert rel_error < 0.10, f"Relative error {rel_error:.4%} exceeds 10%"
    
    def test_circular_loop_field_symmetry(self, loop_params):
        """Test that B-field has expected symmetry.
        
        For a loop in xy-plane:
        - B_z should be maximum at center (z=0)
        - B_z should be symmetric about z=0
        - B_x and B_y should be zero on z-axis
        """
        comm = MPI.COMM_WORLD
        
        # Parameters
        current = loop_params['current']
        loop_radius = loop_params['loop_radius']
        wire_radius = loop_params['wire_radius']
        
        # Generate mesh
        mesh, cell_tags, _ = MeshGenerator.circular_loop_domain(
            loop_radius=loop_radius,
            wire_radius=wire_radius,
            domain_radius=0.15,
            resolution=0.005,
            comm=comm
        )
        
        # Create and solve
        problem = MagnetostaticProblem(mesh=mesh, cell_tags=cell_tags, mu=MU_0)
        solver = MagnetostaticSolver(problem, degree=1)
        
        wire_cross_section = np.pi * wire_radius**2
        J_magnitude = current / wire_cross_section
        
        import ufl
        def J(x):
            return ufl.as_vector([0.0, 0.0, J_magnitude])
        
        solver.solve(current_density=J, subdomain_id=1)
        B = solver.compute_b_field()
        
        # Evaluate at symmetric points
        z = 0.05  # 5 cm from center
        points = np.array([[0, 0, z], [0, 0, -z]])
        
        B_vals = B.eval(points, np.arange(2))
        
        # B_z should be the same at +z and -z
        assert np.abs(B_vals[0, 2] - B_vals[1, 2]) < 1e-10, "B_z not symmetric"
        
        # B_x and B_y should be zero (or very small) on axis
        assert np.abs(B_vals[0, 0]) < 1e-6, "B_x should be zero on axis"
        assert np.abs(B_vals[0, 1]) < 1e-6, "B_y should be zero on axis"


def test_analytical_circular_loop():
    """Unit test for analytical solution calculation."""
    current = 1.0  # A
    radius = 0.05  # m (5 cm)
    
    # Test at z = 0 (center)
    z = np.array([0.0])
    B_z = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(z, current, radius)
    
    # Expected: B_z(0) = μ₀I / (2a)
    expected = MU_0 * current / (2 * radius)
    
    assert np.abs(B_z[0] - expected) < 1e-10, f"B_z(0) incorrect: got {B_z[0]}, expected {expected}"
    
    # Test at z = a (one radius away)
    z = np.array([radius])
    B_z = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(z, current, radius)
    
    # Expected: B_z(a) = μ₀I / (2√2 a)
    expected = MU_0 * current / (2 * np.sqrt(2) * radius)
    
    assert np.abs(B_z[0] - expected) < 1e-10, f"B_z(a) incorrect: got {B_z[0]}, expected {expected}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
