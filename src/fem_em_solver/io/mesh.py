"""Mesh generation utilities for EM simulations."""

from typing import Optional, Tuple, List
from pathlib import Path
import numpy as np
import gmsh
from mpi4py import MPI
import dolfinx
from dolfinx.io import gmshio


class MeshGenerator:
    """Generate meshes for common geometries using Gmsh."""
    
    @staticmethod
    def straight_wire_domain(
        wire_length: float = 1.0,
        wire_radius: float = 0.001,
        domain_radius: float = 0.1,
        resolution: float = 0.005,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for straight wire in cylindrical domain.
        
        Creates a cylindrical domain with a thin wire along the z-axis.
        Useful for validating against the analytical solution for an
        infinite straight wire.
        
        Parameters
        ----------
        wire_length : float
            Length of wire [m]
        wire_radius : float
            Radius of wire [m] (should be small for thin wire approximation)
        domain_radius : float
            Radius of surrounding cylindrical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire, surrounding)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("straight_wire")
            
            # Create wire (cylinder along z-axis)
            wire_tag = gmsh.model.occ.addCylinder(
                0, 0, -wire_length/2,  # center of bottom face
                0, 0, wire_length,      # axis direction and height
                wire_radius
            )
            
            # Create surrounding domain (hollow cylinder)
            domain_tag = gmsh.model.occ.addCylinder(
                0, 0, -wire_length/2,
                0, 0, wire_length,
                domain_radius
            )
            
            # Cut wire out of domain to create separate volumes
            # Actually, we want both as separate volumes, so we fragment
            ov, ovv = gmsh.model.occ.fragment(
                [(3, domain_tag)],
                [(3, wire_tag)]
            )
            gmsh.model.occ.synchronize()
            
            # Get volumes
            volumes = gmsh.model.getEntities(dim=3)
            wire_volume = None
            domain_volume = None
            
            # Tag volumes based on their bounding box
            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                # Check if this is the wire (small radius)
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                r_max = np.sqrt(max(x_max**2, y_max**2))
                if r_max < 2 * wire_radius:
                    wire_volume = vol[1]
                else:
                    domain_volume = vol[1]
            
            # Add physical groups
            if wire_volume:
                gmsh.model.addPhysicalGroup(3, [wire_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire")
            
            if domain_volume:
                gmsh.model.addPhysicalGroup(3, [domain_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "domain")
            
            # Tag boundaries
            surfaces = gmsh.model.getEntities(dim=2)
            boundary_surfaces = []
            wire_surfaces = []
            
            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                r_max = np.sqrt(max(x_max**2, y_max**2))
                
                # Cylindrical boundary of domain
                if abs(r_max - domain_radius) < resolution:
                    boundary_surfaces.append(surf[1])
                # Wire surface
                elif r_max < 2 * wire_radius:
                    wire_surfaces.append(surf[1])
            
            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "boundary")
            
            if wire_surfaces:
                gmsh.model.addPhysicalGroup(2, wire_surfaces, tag=2)
                gmsh.model.setPhysicalName(2, 2, "wire_surface")
            
            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            
            # Generate mesh
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

            # Save for debugging (best effort; do not fail mesh generation if unwritable)
            try:
                debug_dir = Path("paraview_output")
                debug_dir.mkdir(exist_ok=True)
                debug_mesh_path = debug_dir / "straight_wire.msh"
                gmsh.write(str(debug_mesh_path))
                if rank == 0:
                    print(f"  Debug: Mesh saved to {debug_mesh_path} (open in Gmsh to inspect)")
            except Exception as e:
                if rank == 0:
                    print(f"  Debug: Could not write debug mesh file ({e})")
            
        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )
        
        if comm.rank == rank:
            gmsh.finalize()
        
        return mesh, cell_tags, facet_tags
    
    @staticmethod
    def circular_loop_domain(
        loop_radius: float = 0.05,
        wire_radius: float = 0.001,
        domain_radius: float = 0.15,
        resolution: float = 0.005,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for circular current loop in spherical domain.
        
        Creates a torus (ring) for the wire surrounded by a spherical air domain.
        The loop lies in the xy-plane centered at origin.
        
        Parameters
        ----------
        loop_radius : float
            Major radius of loop (distance from center to wire center) [m]
        wire_radius : float
            Minor radius of wire cross-section [m]
        domain_radius : float
            Radius of surrounding spherical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire=1, air=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("circular_loop")
            
            # Create wire as torus in xy-plane
            # addTorus(x, y, z, major_radius, minor_radius)
            wire_tag = gmsh.model.occ.addTorus(0, 0, 0, loop_radius, wire_radius)
            
            # Create surrounding spherical domain
            domain_tag = gmsh.model.occ.addSphere(0, 0, 0, domain_radius)
            
            # Fragment to get separate volumes
            ov, ovv = gmsh.model.occ.fragment(
                [(3, domain_tag)],
                [(3, wire_tag)]
            )
            gmsh.model.occ.synchronize()
            
            # Get volumes and tag them
            volumes = gmsh.model.getEntities(dim=3)
            wire_volume = None
            air_volume = None
            
            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                
                # Compute bounding box dimensions
                dx = x_max - x_min
                dy = y_max - y_min
                dz = z_max - z_min
                
                # Wire has smaller bounding box
                max_dim = max(dx, dy, dz)
                if max_dim < 4 * loop_radius:  # Wire is smaller
                    wire_volume = vol[1]
                else:
                    air_volume = vol[1]
            
            # Add physical groups
            if wire_volume:
                gmsh.model.addPhysicalGroup(3, [wire_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire")
            
            if air_volume:
                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")
            
            # Tag boundaries
            surfaces = gmsh.model.getEntities(dim=2)
            boundary_surfaces = []
            
            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                
                # Check if on outer spherical boundary
                r_max = np.sqrt(x_max**2 + y_max**2 + z_max**2)
                if abs(r_max - domain_radius) < resolution:
                    boundary_surfaces.append(surf[1])
            
            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            
            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            
            # Generate mesh
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")
            
        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )
        
        if comm.rank == rank:
            gmsh.finalize()
        
        return mesh, cell_tags, facet_tags
    
    @staticmethod
    def helmholtz_coil_domain(
        loop_radius: float = 0.05,
        wire_radius: float = 0.002,  # Increased from 0.001 for simpler mesh
        domain_radius: float = 0.12,  # Reduced from 0.15
        resolution: float = 0.008,    # Coarser mesh
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for Helmholtz coil in spherical domain.
        
        A Helmholtz coil consists of two identical circular loops separated
        by a distance equal to their radius. This configuration creates a
        highly uniform magnetic field in the central region.
        
        The loops are positioned at z = -R/2 and z = +R/2, where R is the
        loop radius.
        
        Parameters
        ----------
        loop_radius : float
            Major radius of each loop (distance from center to wire center) [m]
        wire_radius : float
            Minor radius of wire cross-section [m]
        domain_radius : float
            Radius of surrounding spherical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire=1, air=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("helmholtz_coil")
            
            # Helmholtz condition: separation = loop radius
            separation = loop_radius
            z1 = -separation / 2
            z2 = separation / 2
            
            # Create two wire tori at z = -R/2 and z = +R/2
            # addTorus(x, y, z, major_radius, minor_radius)
            wire1_tag = gmsh.model.occ.addTorus(0, 0, z1, loop_radius, wire_radius)
            wire2_tag = gmsh.model.occ.addTorus(0, 0, z2, loop_radius, wire_radius)
            
            # Create surrounding spherical domain
            domain_tag = gmsh.model.occ.addSphere(0, 0, 0, domain_radius)
            
            # Fragment to get separate volumes
            # First fragment domain with wire1
            ov1, ovv1 = gmsh.model.occ.fragment(
                [(3, domain_tag)],
                [(3, wire1_tag)]
            )
            gmsh.model.occ.synchronize()
            
            # Then fragment result with wire2
            # Get the air volume from previous fragmentation
            volumes_after_first = gmsh.model.getEntities(dim=3)
            air_tag = None
            wire1_volume = None
            for vol in volumes_after_first:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                max_dim = max(x_max - x_min, y_max - y_min, z_max - z_min)
                # Wire is smaller than air domain
                if max_dim < 4 * loop_radius:
                    wire1_volume = vol[1]
                else:
                    air_tag = vol[1]
            
            # Fragment air with second wire
            if air_tag:
                ov2, ovv2 = gmsh.model.occ.fragment(
                    [(3, air_tag)],
                    [(3, wire2_tag)]
                )
            gmsh.model.occ.synchronize()
            
            # Get final volumes and tag them
            volumes = gmsh.model.getEntities(dim=3)
            wire_volumes = []
            air_volume = None
            
            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                
                # Compute bounding box dimensions
                dx = x_max - x_min
                dy = y_max - y_min
                dz = z_max - z_min
                max_dim = max(dx, dy, dz)
                
                # Wire volumes are smaller and have z-extent ~ 2*wire_radius
                if max_dim < 4 * loop_radius and dz < 4 * wire_radius:
                    wire_volumes.append(vol[1])
                else:
                    air_volume = vol[1]
            
            # Add physical groups
            if wire_volumes:
                # Tag both wires as "wire" (tag=1)
                gmsh.model.addPhysicalGroup(3, wire_volumes, tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire")
            
            if air_volume:
                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")
            
            # Tag outer boundary
            surfaces = gmsh.model.getEntities(dim=2)
            boundary_surfaces = []
            
            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                
                # Check if on outer spherical boundary
                r_max = np.sqrt(x_max**2 + y_max**2 + z_max**2)
                if abs(r_max - domain_radius) < resolution:
                    boundary_surfaces.append(surf[1])
            
            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            
            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            
            # Generate mesh
            gmsh.model.mesh.generate(3)
            # Skip optimization for Helmholtz coil to avoid Gmsh hangs
            # gmsh.model.mesh.optimize("Netgen")
            
        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )
        
        if comm.rank == rank:
            gmsh.finalize()
        
        return mesh, cell_tags, facet_tags
    
    @staticmethod
    def rectangular_domain(
        bounds: Tuple[float, float, float, float, float, float],
        resolution: float,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> dolfinx.mesh.Mesh:
        """Generate simple rectangular domain.
        
        Parameters
        ----------
        bounds : tuple
            (xmin, xmax, ymin, ymax, zmin, zmax)
        resolution : float
            Mesh size
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("rectangular")
            
            xmin, xmax, ymin, ymax, zmin, zmax = bounds
            
            # Create box
            box = gmsh.model.occ.addBox(xmin, ymin, zmin, 
                                         xmax-xmin, ymax-ymin, zmax-zmin)
            gmsh.model.occ.synchronize()
            
            gmsh.model.addPhysicalGroup(3, [box], tag=1)
            
            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
        
        mesh, _, _ = gmshio.model_to_mesh(gmsh.model, comm, rank, gdim=3)
        
        if comm.rank == rank:
            gmsh.finalize()
        
        return mesh
    
    @staticmethod
    def create_simple_box(
        L: float = 1.0,
        n: int = 10,
        comm: MPI.Intracomm = MPI.COMM_WORLD
    ) -> dolfinx.mesh.Mesh:
        """Create simple box mesh using dolfinx built-in generator.
        
        Parameters
        ----------
        L : float
            Box half-length (domain is [-L, L]³)
        n : int
            Number of cells in each direction
        comm : MPI.Intracomm
            MPI communicator
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            Box mesh
        """
        from dolfinx.mesh import create_box, CellType
        
        domain = [(-L, -L, -L), (L, L, L)]
        mesh = create_box(comm, domain, [n, n, n], CellType.tetrahedron)
        
        return mesh
    
    @staticmethod
    def cylindrical_domain(
        inner_radius: float = 0.01,
        outer_radius: float = 0.1,
        length: float = 0.2,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with cylindrical inner volume inside cylindrical domain.
        
        Creates two concentric cylinders along the z-axis. Useful for practicing
        multi-volume meshing and for problems with cylindrical symmetry.
        
        Parameters
        ----------
        inner_radius : float
            Radius of inner cylinder [m]
        outer_radius : float
            Radius of outer cylinder [m]
        length : float
            Length of cylinders along z-axis [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)
            
        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (inner=1, outer=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("cylindrical_domain")
            
            # Create inner cylinder along z-axis
            inner_tag = gmsh.model.occ.addCylinder(
                0, 0, -length/2,  # center of bottom face
                0, 0, length,      # axis direction and height
                inner_radius
            )
            
            # Create outer cylinder
            outer_tag = gmsh.model.occ.addCylinder(
                0, 0, -length/2,
                0, 0, length,
                outer_radius
            )
            
            # Fragment to create separate volumes (inner and outer region)
            ov, ovv = gmsh.model.occ.fragment(
                [(3, outer_tag)],
                [(3, inner_tag)]
            )
            gmsh.model.occ.synchronize()
            
            # Get volumes and tag them
            volumes = gmsh.model.getEntities(dim=3)
            inner_volume = None
            outer_volume = None
            
            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                
                # Check if this is the inner cylinder by radius
                r_max = np.sqrt(max(x_max**2, y_max**2))
                if r_max < (inner_radius + outer_radius) / 2:
                    inner_volume = vol[1]
                else:
                    outer_volume = vol[1]
            
            # Add physical groups
            if inner_volume:
                gmsh.model.addPhysicalGroup(3, [inner_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "inner")
            
            if outer_volume:
                gmsh.model.addPhysicalGroup(3, [outer_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "outer")
            
            # Tag boundaries
            surfaces = gmsh.model.getEntities(dim=2)
            outer_boundary_surfaces = []
            inner_boundary_surfaces = []
            
            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                r_max = np.sqrt(max(x_max**2, y_max**2))
                
                # Outer cylindrical boundary
                if abs(r_max - outer_radius) < resolution:
                    outer_boundary_surfaces.append(surf[1])
                # Inner cylinder surface
                elif abs(r_max - inner_radius) < resolution:
                    inner_boundary_surfaces.append(surf[1])
            
            if outer_boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            
            if inner_boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, inner_boundary_surfaces, tag=2)
                gmsh.model.setPhysicalName(2, 2, "inner_boundary")
            
            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            
            # Generate mesh
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")
            
        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )
        
        if comm.rank == rank:
            gmsh.finalize()
        
        return mesh, cell_tags, facet_tags

    @staticmethod
    def two_cylinder_domain(
        separation: float = 0.05,
        radius: float = 0.01,
        length: float = 0.1,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with two side-by-side cylinders inside a box domain.

        This intentionally avoids boolean operations/fragmentation and keeps
        all three volumes explicitly tagged.

        Parameters
        ----------
        separation : float
            Center-to-center distance between the two cylinders [m]
        radius : float
            Radius of each cylinder [m]
        length : float
            Cylinder length along z-axis [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (cylinder_1=1, cylinder_2=2, domain=3)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("two_cylinder_domain")

            x_offset = separation / 2

            cylinder_1 = gmsh.model.occ.addCylinder(
                -x_offset, 0, -length / 2,
                0, 0, length,
                radius
            )
            cylinder_2 = gmsh.model.occ.addCylinder(
                x_offset, 0, -length / 2,
                0, 0, length,
                radius
            )

            box_half_x = x_offset + 2.0 * radius
            box_half_y = 2.0 * radius
            domain = gmsh.model.occ.addBox(
                -box_half_x,
                -box_half_y,
                -length / 2,
                2.0 * box_half_x,
                2.0 * box_half_y,
                length,
            )

            gmsh.model.occ.synchronize()

            gmsh.model.addPhysicalGroup(3, [cylinder_1], tag=1)
            gmsh.model.setPhysicalName(3, 1, "cylinder_1")

            gmsh.model.addPhysicalGroup(3, [cylinder_2], tag=2)
            gmsh.model.setPhysicalName(3, 2, "cylinder_2")

            gmsh.model.addPhysicalGroup(3, [domain], tag=3)
            gmsh.model.setPhysicalName(3, 3, "domain")

            # Tag all box boundary surfaces as outer boundary
            boundary_surfaces = []
            x_tol = box_half_x + resolution
            y_tol = box_half_y + resolution
            z_tol = (length / 2) + resolution

            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                if (
                    abs(abs(x_min) - box_half_x) < x_tol or abs(abs(x_max) - box_half_x) < x_tol
                    or abs(abs(y_min) - box_half_y) < y_tol or abs(abs(y_max) - box_half_y) < y_tol
                    or abs(abs(z_min) - (length / 2)) < z_tol or abs(abs(z_max) - (length / 2)) < z_tol
                ):
                    boundary_surfaces.append(surf)

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def two_torus_domain(
        separation: float = 0.05,
        major_radius: float = 0.02,
        minor_radius: float = 0.005,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with two tori inside a box domain.

        Uses non-fragmenting geometry construction: two separate torus volumes
        plus one enclosing domain volume, each explicitly tagged.
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("two_torus_domain")

            z_offset = separation / 2

            wire_1 = gmsh.model.occ.addTorus(0, 0, -z_offset, major_radius, minor_radius)
            wire_2 = gmsh.model.occ.addTorus(0, 0, z_offset, major_radius, minor_radius)

            radial_extent = major_radius + minor_radius
            box_half_x = radial_extent + 2.0 * minor_radius
            box_half_y = radial_extent + 2.0 * minor_radius
            box_half_z = z_offset + minor_radius + 2.0 * minor_radius

            domain = gmsh.model.occ.addBox(
                -box_half_x,
                -box_half_y,
                -box_half_z,
                2.0 * box_half_x,
                2.0 * box_half_y,
                2.0 * box_half_z,
            )

            gmsh.model.occ.synchronize()

            gmsh.model.addPhysicalGroup(3, [wire_1], tag=1)
            gmsh.model.setPhysicalName(3, 1, "wire_1")

            gmsh.model.addPhysicalGroup(3, [wire_2], tag=2)
            gmsh.model.setPhysicalName(3, 2, "wire_2")

            gmsh.model.addPhysicalGroup(3, [domain], tag=3)
            gmsh.model.setPhysicalName(3, 3, "domain")

            boundary_surfaces = []
            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                if (
                    abs(abs(x_min) - box_half_x) < resolution or abs(abs(x_max) - box_half_x) < resolution
                    or abs(abs(y_min) - box_half_y) < resolution or abs(abs(y_max) - box_half_y) < resolution
                    or abs(abs(z_min) - box_half_z) < resolution or abs(abs(z_max) - box_half_z) < resolution
                ):
                    boundary_surfaces.append(surf)

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def coil_phantom_domain(
        coil_major_radius: float = 0.08,
        coil_minor_radius: float = 0.01,
        coil_separation: float = 0.08,
        phantom_radius: float = 0.04,
        phantom_height: float = 0.10,
        air_padding: float = 0.04,
        resolution: float = 0.015,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate a coarse two-coil + cylindrical phantom + air mesh.

        Cell tags:
        - 1: coil_1
        - 2: coil_2
        - 3: phantom
        - 4: air
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("coil_phantom_domain")

            z_offset = coil_separation / 2
            coil_1 = gmsh.model.occ.addTorus(0, 0, -z_offset, coil_major_radius, coil_minor_radius)
            coil_2 = gmsh.model.occ.addTorus(0, 0, z_offset, coil_major_radius, coil_minor_radius)
            phantom = gmsh.model.occ.addCylinder(
                0, 0, -phantom_height / 2,
                0, 0, phantom_height,
                phantom_radius,
            )

            radial_extent = max(coil_major_radius + coil_minor_radius, phantom_radius)
            z_extent = max(z_offset + coil_minor_radius, phantom_height / 2)
            air = gmsh.model.occ.addBox(
                -(radial_extent + air_padding),
                -(radial_extent + air_padding),
                -(z_extent + air_padding),
                2 * (radial_extent + air_padding),
                2 * (radial_extent + air_padding),
                2 * (z_extent + air_padding),
            )

            gmsh.model.occ.fragment(
                [(3, air)],
                [(3, coil_1), (3, coil_2), (3, phantom)],
            )
            gmsh.model.occ.synchronize()

            volumes = gmsh.model.getEntities(dim=3)
            masses = {tag: gmsh.model.occ.getMass(3, tag) for _, tag in volumes}
            air_tag = max(masses, key=masses.get)

            remaining = [tag for _, tag in volumes if tag != air_tag]
            z_centers = {
                tag: gmsh.model.occ.getCenterOfMass(3, tag)[2]
                for tag in remaining
            }

            coil_1_tag = min(remaining, key=lambda tag: z_centers[tag])
            coil_2_tag = max(remaining, key=lambda tag: z_centers[tag])
            phantom_tag = [tag for tag in remaining if tag not in (coil_1_tag, coil_2_tag)][0]

            gmsh.model.addPhysicalGroup(3, [coil_1_tag], tag=1)
            gmsh.model.setPhysicalName(3, 1, "coil_1")
            gmsh.model.addPhysicalGroup(3, [coil_2_tag], tag=2)
            gmsh.model.setPhysicalName(3, 2, "coil_2")
            gmsh.model.addPhysicalGroup(3, [phantom_tag], tag=3)
            gmsh.model.setPhysicalName(3, 3, "phantom")
            gmsh.model.addPhysicalGroup(3, [air_tag], tag=4)
            gmsh.model.setPhysicalName(3, 4, "air")

            outer_boundary_surfaces = []
            box_half_x = radial_extent + air_padding
            box_half_z = z_extent + air_padding
            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                if (
                    abs(abs(x_min) - box_half_x) < resolution or abs(abs(x_max) - box_half_x) < resolution
                    or abs(abs(y_min) - box_half_x) < resolution or abs(abs(y_max) - box_half_x) < resolution
                    or abs(abs(z_min) - box_half_z) < resolution or abs(abs(z_max) - box_half_z) < resolution
                ):
                    outer_boundary_surfaces.append(surf)

            if outer_boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def birdcage_port_domain(
        n_legs: int = 4,
        ring_radius: float = 0.07,
        leg_radius: float = 0.006,
        leg_height: float = 0.14,
        ring_minor_radius: float = 0.004,
        phantom_radius: float = 0.03,
        phantom_height: float = 0.08,
        port_box_size: Tuple[float, float, float] = (0.010, 0.008, 0.010),
        air_padding: float = 0.03,
        resolution: float = 0.015,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate a coarse birdcage-like geometry fixture with explicit port tags.

        Cell tags:
        - 1: conductor (rings + legs)
        - 2: air
        - 3: phantom
        - 101..(100+n_legs): per-port regions between adjacent legs
        """
        if n_legs < 3:
            raise ValueError("n_legs must be >= 3 for a birdcage-like fixture")

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("birdcage_port_domain")

            # Simplified conductor scaffold: two rings plus vertical legs.
            z_ring_offset = 0.5 * leg_height - 1.5 * ring_minor_radius
            top_ring = gmsh.model.occ.addTorus(0, 0, z_ring_offset, ring_radius, ring_minor_radius)
            bottom_ring = gmsh.model.occ.addTorus(0, 0, -z_ring_offset, ring_radius, ring_minor_radius)

            leg_tags: List[int] = []
            theta = np.linspace(0.0, 2.0 * np.pi, n_legs, endpoint=False)
            for angle in theta:
                x = ring_radius * np.cos(angle)
                y = ring_radius * np.sin(angle)
                leg = gmsh.model.occ.addCylinder(
                    x,
                    y,
                    -0.5 * leg_height,
                    0.0,
                    0.0,
                    leg_height,
                    leg_radius,
                )
                leg_tags.append(leg)

            phantom_tag = gmsh.model.occ.addCylinder(
                0.0,
                0.0,
                -0.5 * phantom_height,
                0.0,
                0.0,
                phantom_height,
                phantom_radius,
            )

            port_dx, port_dy, port_dz = port_box_size
            port_radius = ring_radius + 0.5 * port_dy
            port_tags: List[int] = []
            for idx, angle in enumerate(theta):
                next_angle = theta[(idx + 1) % n_legs]
                midpoint_angle = np.arctan2(
                    np.sin(angle) + np.sin(next_angle),
                    np.cos(angle) + np.cos(next_angle),
                )
                cx = port_radius * np.cos(midpoint_angle)
                cy = port_radius * np.sin(midpoint_angle)
                port = gmsh.model.occ.addBox(
                    cx - 0.5 * port_dx,
                    cy - 0.5 * port_dy,
                    -0.5 * port_dz,
                    port_dx,
                    port_dy,
                    port_dz,
                )
                port_tags.append(port)

            radial_extent = ring_radius + max(leg_radius, ring_minor_radius) + port_dy + air_padding
            z_extent = max(0.5 * leg_height, 0.5 * phantom_height) + air_padding
            air_tag = gmsh.model.occ.addBox(
                -radial_extent,
                -radial_extent,
                -z_extent,
                2.0 * radial_extent,
                2.0 * radial_extent,
                2.0 * z_extent,
            )

            conductor_tags = [top_ring, bottom_ring] + leg_tags
            cut_result, _ = gmsh.model.occ.cut(
                [(3, air_tag)],
                [(3, tag) for tag in conductor_tags + [phantom_tag] + port_tags],
                removeObject=True,
                removeTool=False,
            )
            gmsh.model.occ.synchronize()

            if not cut_result:
                raise RuntimeError("Failed to create air volume for birdcage_port_domain")
            air_volume_tags = [tag for dim, tag in cut_result if dim == 3]
            if not air_volume_tags:
                raise RuntimeError("No 3D air volume returned by cut operation")
            air_tag = air_volume_tags[0]

            gmsh.model.addPhysicalGroup(3, conductor_tags, tag=1)
            gmsh.model.setPhysicalName(3, 1, "conductor")
            gmsh.model.addPhysicalGroup(3, [air_tag], tag=2)
            gmsh.model.setPhysicalName(3, 2, "air")
            gmsh.model.addPhysicalGroup(3, [phantom_tag], tag=3)
            gmsh.model.setPhysicalName(3, 3, "phantom")

            for idx, port_tag in enumerate(port_tags, start=1):
                physical_tag = 100 + idx
                gmsh.model.addPhysicalGroup(3, [port_tag], tag=physical_tag)
                gmsh.model.setPhysicalName(3, physical_tag, f"port_P{idx}")

            outer_boundary_surfaces = []
            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                if (
                    abs(abs(x_min) - radial_extent) < resolution
                    or abs(abs(x_max) - radial_extent) < resolution
                    or abs(abs(y_min) - radial_extent) < resolution
                    or abs(abs(y_max) - radial_extent) < resolution
                    or abs(abs(z_min) - z_extent) < resolution
                    or abs(abs(z_max) - z_extent) < resolution
                ):
                    outer_boundary_surfaces.append(surf)

            if outer_boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, comm, rank, gdim=3
        )

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags
