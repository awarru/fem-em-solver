"""
Validation test: Helmholtz coil field uniformity.

A Helmholtz coil consists of two identical circular loops separated by
their radius. This configuration creates a highly uniform magnetic field
in the central region between the loops.

Analytical solution on axis:
    B_z(z) = B1(z) + B2(z)
    
Where B1 and B2 are the fields from each loop:
    B_loop(z) = μ₀IR² / (2(R² + z²)^(3/2))
    
For Helmholtz coil with loops at z = ±R/2:
    B_z(0) = (4/5)^(3/2) * μ₀I / R ≈ 0.7155 * μ₀I / R

The field is maximally uniform at the center.
"""

import pytest
import numpy as np
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticSolver, MagnetostaticProblem
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


class TestHelmholtzCoil:
    """Validation tests for Helmholtz coil."""
    
    @pytest.fixture
    def coil_params(self):
        """Standard Helmholtz coil parameters."""
        return {
            'current': 1.0,           # 1 Ampere in each loop
            'loop_radius': 0.05,      # 5 cm radius
            'wire_radius': 0.001,     # 1 mm wire thickness
            'domain_radius': 0.15,    # 15 cm domain
            'resolution': 0.005,      # 5 mm mesh size
        }
    
    def test_helmholtz_field_on_axis(self, coil_params, tmp_path):
        """Test B_z on axis matches analytical solution.
        
        This test verifies that the numerical solution matches the
        analytical Helmholtz coil field along the z-axis.
        """
        comm = MPI.COMM_WORLD
        
        # Parameters
        current = coil_params['current']
        loop_radius = coil_params['loop_radius']
        wire_radius = coil_params['wire_radius']
        domain_radius = coil_params['domain_radius']
        resolution = coil_params['resolution']
        
        print(f"\nHelmholtz coil test:")
        print(f"  Loop radius: {loop_radius} m")
        print(f"  Separation: {loop_radius} m (Helmholtz condition)")
        print(f"  Wire radius: {wire_radius} m")
        print(f"  Current: {current} A")
        
        # Generate mesh
        print("  Generating mesh...")
        mesh, cell_tags, facet_tags = MeshGenerator.helmholtz_coil_domain(
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
        
        # Define current density in wires
        wire_cross_section = np.pi * wire_radius**2
        J_magnitude = current / wire_cross_section
        
        import ufl
        def current_density(x):
            """Uniform current density in z-direction."""
            return ufl.as_vector([0.0, 0.0, J_magnitude])
        
        # Solve with current restricted to wire subdomain (tag=1)
        print("  Solving...")
        A = solver.solve(current_density=current_density, subdomain_id=1)
        
        # Compute B-field
        print("  Computing B-field...")
        B = solver.compute_b_field()
        
        # Evaluate along z-axis
        n_points = 21
        z_eval = np.linspace(-0.1, 0.1, n_points)
        
        points = np.zeros((n_points, 3))
        points[:, 2] = z_eval  # z positions
        
        B_num = B.eval(points, np.arange(n_points))
        B_num_z = B_num[:, 2]  # z-component
        
        # Analytical solution
        B_ana_z = AnalyticalSolutions.helmholtz_coil_field_on_axis(
            z_eval, current, loop_radius
        )
        
        # Error metrics
        rel_error = ErrorMetrics.l2_relative_error(B_num_z, B_ana_z)
        max_error = ErrorMetrics.max_relative_error(B_num_z, B_ana_z)
        
        print(f"\n  Results:")
        print(f"    B_z at center (numerical): {B_num_z[n_points//2]:.6e} T")
        print(f"    B_z at center (analytical): {B_ana_z[n_points//2]:.6e} T")
        print(f"    Expected B_z(0) = 0.7155 * μ₀I/R = {0.7155 * MU_0 * current / loop_radius:.6e} T")
        print(f"    Relative L2 error: {rel_error:.4%}")
        print(f"    Max relative error: {max_error:.4%}")
        
        # Assert error is within tolerance
        assert rel_error < 0.10, f"Relative error {rel_error:.4%} exceeds 10%"
    
    def test_helmholtz_field_uniformity(self, coil_params):
        """Test field uniformity in central region.
        
        The Helmholtz coil is designed to produce a uniform field in the
central region. We check that the coefficient of variation (CV) of
        B_z is less than 1% in the central 20% of the separation.
        """
        comm = MPI.COMM_WORLD
        
        # Parameters
        current = coil_params['current']
        loop_radius = coil_params['loop_radius']
        wire_radius = coil_params['wire_radius']
        
        # Generate mesh
        mesh, cell_tags, _ = MeshGenerator.helmholtz_coil_domain(
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
        
        # Evaluate in central region: z ∈ [-0.1R, 0.1R]
        # This is the central 20% of the coil separation
        central_range = 0.1 * loop_radius
        n_points = 11
        z_eval = np.linspace(-central_range, central_range, n_points)
        
        points = np.zeros((n_points, 3))
        points[:, 2] = z_eval
        
        B_vals = B.eval(points, np.arange(n_points))
        B_z = B_vals[:, 2]
        
        # Compute coefficient of variation (CV = std / mean)
        mean_Bz = np.mean(B_z)
        std_Bz = np.std(B_z)
        cv = std_Bz / abs(mean_Bz)
        
        print(f"\n  Uniformity test:")
        print(f"    Central region: z ∈ [{-central_range*100:.2f}, {central_range*100:.2f}] cm")
        print(f"    Mean B_z: {mean_Bz:.6e} T")
        print(f"    Std B_z: {std_Bz:.6e} T")
        print(f"    Coefficient of variation: {cv:.4%}")
        
        # Helmholtz coil should have CV < 1% in central region
        assert cv < 0.01, f"Field non-uniform: CV = {cv:.4%} (should be < 1%)"
    
    def test_helmholtz_symmetry(self, coil_params):
        """Test that field is symmetric about z=0."""
        comm = MPI.COMM_WORLD
        
        current = coil_params['current']
        loop_radius = coil_params['loop_radius']
        wire_radius = coil_params['wire_radius']
        
        mesh, cell_tags, _ = MeshGenerator.helmholtz_coil_domain(
            loop_radius=loop_radius,
            wire_radius=wire_radius,
            domain_radius=0.15,
            resolution=0.005,
            comm=comm
        )
        
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
        z = 0.03
        points = np.array([[0, 0, z], [0, 0, -z]])
        
        B_vals = B.eval(points, np.arange(2))
        
        # B_z should be the same at +z and -z
        assert np.abs(B_vals[0, 2] - B_vals[1, 2]) < 1e-6, "B_z not symmetric about z=0"
        
        # B_x and B_y should be zero on axis
        assert np.abs(B_vals[0, 0]) < 1e-6, "B_x should be zero on axis"
        assert np.abs(B_vals[0, 1]) < 1e-6, "B_y should be zero on axis"


def test_analytical_helmholtz():
    """Unit test for Helmholtz coil analytical solution."""
    current = 1.0  # A
    radius = 0.05  # m (5 cm)
    
    # Test at center (z=0)
    z = np.array([0.0])
    B_z = AnalyticalSolutions.helmholtz_coil_field_on_axis(z, current, radius)
    
    # Expected: B_z(0) = (4/5)^(3/2) * μ₀I / R
    expected = (4/5)**(3/2) * MU_0 * current / radius
    
    assert np.abs(B_z[0] - expected) < 1e-10, f"B_z(0) incorrect: got {B_z[0]}, expected {expected}"
    
    # Test symmetry: B_z(z) = B_z(-z)
    z_pos = np.array([0.02])
    z_neg = np.array([-0.02])
    
    B_pos = AnalyticalSolutions.helmholtz_coil_field_on_axis(z_pos, current, radius)
    B_neg = AnalyticalSolutions.helmholtz_coil_field_on_axis(z_neg, current, radius)
    
    assert np.abs(B_pos[0] - B_neg[0]) < 1e-15, "Analytical solution not symmetric"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
