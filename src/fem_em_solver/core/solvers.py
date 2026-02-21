"""
Magnetostatic solver using magnetic vector potential formulation.

The governing equation is:
    ∇ × (μ⁻¹ ∇ × A) = J

Where:
    A = magnetic vector potential [Wb/m]
    B = ∇ × A = magnetic flux density [T]
    H = μ⁻¹B = magnetic field intensity [A/m]
    J = current density [A/m²]
    μ = permeability [H/m]

Weak form:
    ∫ μ⁻¹ (∇ × A) · (∇ × v) dx = ∫ J · v dx

Boundary conditions:
    - Dirichlet: A × n = g (tangential A specified)
    - Natural: (∇ × A) × n = 0 (magnetic insulation)
"""

from typing import Optional, Callable, Union, List, Sequence
import numpy as np
from dataclasses import dataclass

import dolfinx
import ufl
from dolfinx import fem, mesh, io
from dolfinx.fem.petsc import LinearProblem
from ufl import curl, inner, dx, TrialFunction, TestFunction
from mpi4py import MPI
import gmsh

from ..utils.constants import MU_0


@dataclass
class MagnetostaticProblem:
    """Container for magnetostatic problem parameters."""
    mesh: dolfinx.mesh.Mesh
    cell_tags: Optional[dolfinx.mesh.MeshTags] = None
    facet_tags: Optional[dolfinx.mesh.MeshTags] = None
    mu: Union[float, Callable] = MU_0
    

class MagnetostaticSolver:
    """Solver for magnetostatic problems using vector potential formulation.
    
    This solver computes the magnetic vector potential A, from which
    B = ∇ × A and H = μ⁻¹B can be derived.
    
    Parameters
    ----------
    problem : MagnetostaticProblem
        Problem definition including mesh and material properties
    degree : int, optional
        Polynomial degree of Nedelec elements (default: 1)
    """
    
    def __init__(self, problem: MagnetostaticProblem, degree: int = 1):
        self.problem = problem
        self.degree = degree
        self.mesh = problem.mesh
        self.mu = problem.mu
        
        # Create function space (H(curl) - Nedelec elements)
        self.V = fem.functionspace(self.mesh, ("N1curl", degree))
        
        # Solution field
        self.A = fem.Function(self.V, name="A")
        self._solved = False
        
    def solve(self, current_density: Optional[Callable] = None, 
              bc_functions: Optional[List] = None,
              subdomain_id: Optional[int] = None,
              subdomain_ids: Optional[Sequence[int]] = None,
              gauge_penalty: float = 1e-3) -> fem.Function:
        """Solve the magnetostatic problem.
        
        Parameters
        ----------
        current_density : callable, optional
            Function returning J(x) for any point x.
            If None, assumes J = 0 (no sources).
        bc_functions : list, optional
            List of Dirichlet BC functions
        subdomain_id : int, optional
            If provided, restricts current density integration to cells
            with this tag (requires cell_tags in problem).
        subdomain_ids : Sequence[int], optional
            If provided, restricts current density integration to the union
            of these cell tags (requires cell_tags in problem).
            Mutually exclusive with ``subdomain_id``.
        gauge_penalty : float, optional
            Small regularization term added to remove the nullspace in
            pure curl-curl problems (default: 1e-3).
            
        Returns
        -------
        fem.Function
            Magnetic vector potential A
        """
        # Define trial and test functions
        A_trial = TrialFunction(self.V)
        v = TestFunction(self.V)
        
        # Permeability (could be spatially varying)
        if callable(self.mu):
            mu = self.mu(ufl.SpatialCoordinate(self.mesh))
        else:
            mu = fem.Constant(self.mesh, self.mu)
        mu_inv = 1.0 / mu
        
        # Bilinear form: a(A, v) = ∫ μ⁻¹ (∇ × A) · (∇ × v) dx
        # Add tiny gauge regularization to remove nullspace.
        gauge = fem.Constant(self.mesh, gauge_penalty)
        a = inner(mu_inv * curl(A_trial), curl(v)) * dx + gauge * inner(A_trial, v) * dx
        
        # Linear form: L(v) = ∫ J · v dx
        if subdomain_id is not None and subdomain_ids is not None:
            raise ValueError("Use either subdomain_id or subdomain_ids, not both")

        if subdomain_ids is None and subdomain_id is not None:
            subdomain_ids = [subdomain_id]

        if current_density is not None:
            x = ufl.SpatialCoordinate(self.mesh)
            J = current_density(x)
        else:
            J = fem.Constant(self.mesh, np.zeros(3))

        # If subdomain ids are provided, restrict integration to their union
        if subdomain_ids is not None:
            if self.problem.cell_tags is None:
                raise ValueError("subdomain_id(s) requested but problem.cell_tags is None")

            L = 0
            for tag in subdomain_ids:
                dx_sub = ufl.Measure(
                    "dx",
                    domain=self.mesh,
                    subdomain_data=self.problem.cell_tags,
                    subdomain_id=int(tag),
                )
                L += inner(J, v) * dx_sub
        else:
            # Integrate over whole domain
            L = inner(J, v) * dx
        
        # Apply boundary conditions
        bcs = []
        if bc_functions is not None:
            bcs = bc_functions
            
        # Solve
        problem = LinearProblem(
            a, L, 
            bcs=bcs,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
        )
        self.A = problem.solve()
        self._solved = True
        
        return self.A
    
    def compute_b_field(self) -> fem.Function:
        """Compute magnetic flux density B = ∇ × A.
        
        Returns
        -------
        fem.Function
            B-field in DG space (discontinuous Galerkin)
        """
        if not self._solved:
            raise RuntimeError("Must call solve() before computing B-field")
        
        # B = curl(A) needs to be interpolated to appropriate space
        # Use DG space for B (discontinuous, vector-valued)
        DG = fem.functionspace(self.mesh, ("DG", self.degree, (3,)))
        B = fem.Function(DG, name="B")
        
        # Project curl(A) onto DG space
        B_expr = fem.Expression(curl(self.A), DG.element.interpolation_points())
        B.interpolate(B_expr)
        
        return B
    
    def compute_h_field(self) -> fem.Function:
        """Compute magnetic field intensity H = B/μ.
        
        Returns
        -------
        fem.Function
            H-field in DG space
        """
        B = self.compute_b_field()
        
        DG = fem.functionspace(self.mesh, ("DG", self.degree, (3,)))
        H = fem.Function(DG, name="H")
        
        if callable(self.mu):
            x = ufl.SpatialCoordinate(self.mesh)
            mu = self.mu(x)
        else:
            mu = fem.Constant(self.mesh, self.mu)
        
        H_expr = fem.Expression(B / mu, DG.element.interpolation_points())
        H.interpolate(H_expr)
        
        return H
    
    def compute_magnetic_energy(self) -> float:
        """Compute total magnetic energy in the domain.
        
        W = ½ ∫ B · H dx = ½ ∫ μ⁻¹ |∇ × A|² dx
        
        Returns
        -------
        float
            Magnetic energy [Joules]
        """
        if not self._solved:
            raise RuntimeError("Must call solve() before computing energy")
        
        if callable(self.mu):
            x = ufl.SpatialCoordinate(self.mesh)
            mu_inv = 1.0 / self.mu(x)
        else:
            mu_inv = 1.0 / self.mu
        
        energy_expr = 0.5 * inner(mu_inv * curl(self.A), curl(self.A)) * dx
        energy = fem.assemble_scalar(fem.form(energy_expr))
        
        return energy
    
    def evaluate_at_points(self, points: np.ndarray, field: str = "A") -> np.ndarray:
        """Evaluate field at specific points.
        
        Parameters
        ----------
        points : np.ndarray
            Array of shape (n_points, 3) with coordinates
        field : str
            Field to evaluate: "A", "B", or "H"
            
        Returns
        -------
        np.ndarray
            Field values at points, shape (n_points, 3)
        """
        if field == "A":
            f = self.A
        elif field == "B":
            f = self.compute_b_field()
        elif field == "H":
            f = self.compute_h_field()
        else:
            raise ValueError(f"Unknown field: {field}")
        
        # Use dolfinx interpolation for evaluation
        values = f.eval(points, np.arange(len(points)))
        return values
    
    def save_to_vtk(self, filename: str, fields: Optional[List[str]] = None):
        """Save solution to VTK file for visualization.
        
        Parameters
        ----------
        filename : str
            Output filename (.pvd or .vtu)
        fields : list, optional
            List of fields to save: ["A", "B", "H"]
        """
        if fields is None:
            fields = ["A"]
        
        # Write to file
        with io.VTKFile(self.mesh.comm, filename, "w") as vtk:
            if "A" in fields:
                vtk.write_function(self.A)
            if "B" in fields:
                vtk.write_function(self.compute_b_field())
            if "H" in fields:
                vtk.write_function(self.compute_h_field())
