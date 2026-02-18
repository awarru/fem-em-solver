#!/usr/bin/env python3
"""Test cylindrical domain mesh generation."""

from fem_em_solver.io.mesh import MeshGenerator
from mpi4py import MPI

mesh, ct, ft = MeshGenerator.cylindrical_domain(
    inner_radius=0.01, outer_radius=0.1, length=0.2, resolution=0.02,
    comm=MPI.COMM_WORLD
)

print(f"Mesh cells: {mesh.topology.index_map(3).size_global}")
print(f"Inner cells: {(ct.values == 1).sum()}")
print(f"Outer cells: {(ct.values == 2).sum()}")
print("SUCCESS: Cylindrical domain mesh generated")