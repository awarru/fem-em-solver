"""Mesh generation utilities for EM simulations."""

from typing import Optional, Tuple, List
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
            
            # Save for debugging (optional)
            # gmsh.write("straight_wire.msh")
            
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
        wire_radius: float = 0.001,
        domain_radius: float = 0.15,
        resolution: float = 0.005,
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
            gmsh.model.mesh.optimize("Netgen")
            
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
