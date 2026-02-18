"""
ParaView output utilities for FEM-EM solver.

This module provides custom XDMF writers that properly combine
cell tags and field data on the same grid for ParaView visualization.
"""

import h5py
import numpy as np
from pathlib import Path
from mpi4py import MPI


def write_xdmf_with_tags(filename, mesh, cell_tags, functions, comm=MPI.COMM_WORLD):
    """
    Write a combined XDMF file with cell tags and field data on the same grid.

    This function creates a properly structured XDMF file where cell tags (cell data)
    and field functions (point data) coexist on a single grid, making them both
    available for filters like Threshold in ParaView.

    Parameters
    ----------
    filename : str or Path
        Output filename (without extension, .xdmf and .h5 will be added)
    mesh : dolfinx.mesh.Mesh
        The mesh to write
    cell_tags : dolfinx.mesh.MeshTags
        Cell tags to write (e.g., wire vs air domain)
    functions : dict
        Dictionary of {name: function} pairs to write
        Functions should already be interpolated to Lagrange space if needed
    comm : MPI.Comm
        MPI communicator (default: MPI.COMM_WORLD)

    Example
    -------
    >>> from fem_em_solver.io.paraview_utils import write_xdmf_with_tags
    >>> write_xdmf_with_tags(
    ...     "straight_wire_combined",
    ...     mesh,
    ...     cell_tags,
    ...     {"B": B_lag, "A": A_lag},
    ...     comm=comm
    ... )
    """
    filename = Path(filename)
    h5_file = filename.with_suffix('.h5')
    xdmf_file = filename.with_suffix('.xdmf')

    rank = comm.rank
    size = comm.size

    # Get local mesh data
    tdim = mesh.topology.dim
    num_cells_local = mesh.topology.index_map(tdim).size_local

    # Get geometry and topology locally
    geom_dofmap = mesh.geometry.dofmap
    nodes_per_cell = tdim + 1  # tetrahedra=4, triangles=3
    local_topology = geom_dofmap.reshape(-1, nodes_per_cell)[:num_cells_local, :]
    local_geometry = mesh.geometry.x

    # Get local cell tags
    if cell_tags is not None:
        local_cell_tag_values = cell_tags.values[:num_cells_local]
    else:
        local_cell_tag_values = None

    # Get local function values (we'll need to gather these too)
    local_function_data = {}
    for name, func in functions.items():
        local_function_data[name] = func.x.array.copy()

    if rank == 0:
        print(f"    DEBUG: Running with {size} MPI rank(s)")
        print(f"    DEBUG: Local data - geometry: {local_geometry.shape}, topology: {local_topology.shape}")

    # Gather all data to rank 0
    if rank == 0:
        print(f"    DEBUG: Gathering geometry from all ranks...")
    all_geometry = comm.gather(local_geometry, root=0)

    if rank == 0:
        print(f"    DEBUG: Gathering topology from all ranks...")
    all_topology = comm.gather(local_topology, root=0)

    if rank == 0:
        print(f"    DEBUG: Gathering cell tags from all ranks...")
    all_cell_tags = comm.gather(local_cell_tag_values, root=0)

    if rank == 0:
        print(f"    DEBUG: Gathering function data from all ranks...")
    all_function_data = comm.gather(local_function_data, root=0)

    # Combine gathered data on rank 0
    if rank == 0:
        print(f"    DEBUG: All data gathered successfully!")
        print(f"    DEBUG: Processing {size} rank(s) to build global mesh...")

        # Build global geometry and vertex mapping
        # We need to identify unique vertices and renumber globally
        # Use spatial hashing for O(N) instead of O(N²) performance
        tolerance = 1e-10
        vertex_hash = {}  # Maps rounded coords to list of (exact_coords, global_idx)
        global_geometry_list = []
        global_topology_list = []
        global_cell_tags_list = []

        # Track function values per vertex
        vertex_function_values = {name: {} for name in functions.keys()}

        def hash_coords(coords, tol=tolerance):
            """Create a hash key for spatial lookup."""
            return tuple(np.round(coords / tol).astype(np.int64))

        for rank_id in range(size):
            rank_geometry = all_geometry[rank_id]
            rank_topology = all_topology[rank_id]
            rank_cell_tags = all_cell_tags[rank_id]
            rank_functions = all_function_data[rank_id]

            print(f"    DEBUG: Processing rank {rank_id}/{size-1} - {len(rank_geometry)} vertices, {len(rank_topology)} cells")

            # Build local-to-global vertex mapping for this rank
            local_to_global = {}
            num_new_vertices = 0
            num_duplicate_vertices = 0

            for local_idx in range(len(rank_geometry)):
                if local_idx % 1000 == 0 and local_idx > 0:
                    print(f"    DEBUG:   Rank {rank_id} - processed {local_idx}/{len(rank_geometry)} vertices (new: {num_new_vertices}, duplicates: {num_duplicate_vertices})")

                vertex_coords = rank_geometry[local_idx]
                hash_key = hash_coords(vertex_coords, tolerance)

                # Check if we've seen this vertex before
                found = False
                if hash_key in vertex_hash:
                    # Check all vertices in this hash bucket
                    for exact_coords, global_idx in vertex_hash[hash_key]:
                        if np.allclose(vertex_coords, exact_coords, rtol=1e-10, atol=1e-12):
                            local_to_global[local_idx] = global_idx
                            found = True
                            num_duplicate_vertices += 1
                            # Update function values (use the one from this rank)
                            for name in functions.keys():
                                try:
                                    value_shape = functions[name].function_space.element.value_shape
                                    if len(value_shape) > 0 and value_shape[0] > 1:
                                        # Vector function
                                        num_components = value_shape[0]
                                        vertex_function_values[name][global_idx] = rank_functions[name][local_idx*num_components:(local_idx+1)*num_components]
                                    else:
                                        # Scalar function
                                        vertex_function_values[name][global_idx] = rank_functions[name][local_idx]
                                except:
                                    # Fallback - assume vector with 3 components
                                    vertex_function_values[name][global_idx] = rank_functions[name][local_idx*3:(local_idx+1)*3]
                            break

                if not found:
                    # New unique vertex
                    global_idx = len(global_geometry_list)
                    num_new_vertices += 1

                    # Add to hash map
                    if hash_key not in vertex_hash:
                        vertex_hash[hash_key] = []
                    vertex_hash[hash_key].append((vertex_coords.copy(), global_idx))

                    local_to_global[local_idx] = global_idx
                    global_geometry_list.append(vertex_coords)

                    # Store function values for this vertex
                    for name in functions.keys():
                        try:
                            value_shape = functions[name].function_space.element.value_shape
                            if len(value_shape) > 0 and value_shape[0] > 1:
                                # Vector function
                                num_components = value_shape[0]
                                vertex_function_values[name][global_idx] = rank_functions[name][local_idx*num_components:(local_idx+1)*num_components]
                            else:
                                # Scalar function
                                vertex_function_values[name][global_idx] = rank_functions[name][local_idx]
                        except:
                            # Fallback - assume vector with 3 components
                            vertex_function_values[name][global_idx] = rank_functions[name][local_idx*3:(local_idx+1)*3]

            print(f"    DEBUG:   Rank {rank_id} complete - {num_new_vertices} new vertices, {num_duplicate_vertices} duplicates")

            # Renumber topology using local_to_global mapping
            print(f"    DEBUG:   Renumbering {len(rank_topology)} cells...")
            renumbered_topology = np.array([[local_to_global[local_vertex] for local_vertex in cell]
                                             for cell in rank_topology], dtype=np.int32)
            global_topology_list.append(renumbered_topology)
            print(f"    DEBUG:   Renumbering complete")

            # Add cell tags
            if rank_cell_tags is not None:
                global_cell_tags_list.extend(rank_cell_tags)

        print(f"    DEBUG: Converting to numpy arrays...")
        # Convert to numpy arrays
        geometry = np.array(global_geometry_list, dtype=np.float64)
        topology = np.vstack(global_topology_list).astype(np.int32)

        print(f"    DEBUG: Reconstructing function arrays in vertex order...")
        # Reconstruct function arrays in vertex order
        function_arrays = {}
        for name in functions.keys():
            print(f"    DEBUG:   Processing function '{name}'...")
            sorted_values = [vertex_function_values[name][i] for i in range(len(geometry))]
            function_arrays[name] = np.array(sorted_values, dtype=np.float64)
            print(f"    DEBUG:   Function '{name}' array shape: {function_arrays[name].shape}")

        print(f"    DEBUG: Global mesh constructed successfully!")
        print(f"    DEBUG: Global geometry shape: {geometry.shape}")
        print(f"    DEBUG: Global topology shape: {topology.shape}")

        if len(global_cell_tags_list) > 0:
            cell_tag_values = np.array(global_cell_tags_list, dtype=np.int32)
            print(f"    DEBUG: Global cell tags shape: {cell_tag_values.shape}")
        else:
            cell_tag_values = None

        num_cells = topology.shape[0]
        num_nodes = geometry.shape[0]
        gdim = geometry.shape[1]

        # Cell type mapping (FEniCS to XDMF)
        cell_type_map = {
            2: ("Triangle", 3),
            3: ("Tetrahedron", 4),
        }
        cell_type_name, nodes_per_cell_xdmf = cell_type_map[mesh.topology.dim]

        # Write HDF5 file
        print(f"    DEBUG: ========== Starting HDF5 file write ==========")
        print(f"    DEBUG: Writing to {h5_file}")
        print(f"    DEBUG: Geometry shape={geometry.shape}, dtype={geometry.dtype}")
        print(f"    DEBUG: Topology shape={topology.shape}, dtype={topology.dtype}")

        print(f"    DEBUG: Creating HDF5 file...")
        with h5py.File(h5_file, 'w') as h5:
            # Write mesh (ensure correct dtypes)
            print("    DEBUG: Writing geometry...")
            h5.create_dataset("Mesh/geometry", data=geometry.astype(np.float64))
            print("    DEBUG: Writing topology...")
            h5.create_dataset("Mesh/topology", data=topology.astype(np.int32))

            # Write cell tags
            if cell_tag_values is not None:
                print(f"    DEBUG: Cell tags shape={cell_tag_values.shape}, dtype={cell_tag_values.dtype}")
                h5.create_dataset("CellTags/values", data=cell_tag_values.astype(np.int32))

            # Write functions
            for name in functions.keys():
                print(f"    DEBUG: Processing function '{name}'...")
                values = function_arrays[name]
                print(f"    DEBUG:   Values shape={values.shape}, dtype={values.dtype}")

                # Ensure 2D shape
                if len(values.shape) == 1:
                    values = values.reshape(-1, 1)
                    print(f"    DEBUG:   Reshaped to: {values.shape}")

                print(f"    DEBUG:   Writing dataset Functions/{name}...")
                h5.create_dataset(f"Functions/{name}", data=values.astype(np.float64))
                print(f"    DEBUG:   Successfully wrote {name}")

        print(f"    DEBUG: HDF5 file write complete!")

        # Write XDMF file
        print(f"    DEBUG: Writing XDMF file to {xdmf_file}...")
        with open(xdmf_file, 'w') as xdmf:
            xdmf.write('<?xml version="1.0"?>\n')
            xdmf.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
            xdmf.write('<Xdmf Version="3.0">\n')
            xdmf.write('  <Domain>\n')
            xdmf.write('    <Grid Name="mesh" GridType="Uniform">\n')

            # Topology
            xdmf.write(f'      <Topology TopologyType="{cell_type_name}" ')
            xdmf.write(f'NumberOfElements="{num_cells}">\n')
            xdmf.write(f'        <DataItem Dimensions="{num_cells} {nodes_per_cell_xdmf}" ')
            xdmf.write('NumberType="Int" Format="HDF">\n')
            xdmf.write(f'          {h5_file.name}:/Mesh/topology\n')
            xdmf.write('        </DataItem>\n')
            xdmf.write('      </Topology>\n')

            # Geometry
            xdmf.write(f'      <Geometry GeometryType="XYZ">\n')
            xdmf.write(f'        <DataItem Dimensions="{num_nodes} {gdim}" ')
            xdmf.write('Format="HDF">\n')
            xdmf.write(f'          {h5_file.name}:/Mesh/geometry\n')
            xdmf.write('        </DataItem>\n')
            xdmf.write('      </Geometry>\n')

            # Cell tags (cell data)
            if cell_tag_values is not None:
                xdmf.write('      <Attribute Name="Cell tags" AttributeType="Scalar" ')
                xdmf.write('Center="Cell">\n')
                xdmf.write(f'        <DataItem Dimensions="{num_cells} 1" ')
                xdmf.write('Format="HDF">\n')
                xdmf.write(f'          {h5_file.name}:/CellTags/values\n')
                xdmf.write('        </DataItem>\n')
                xdmf.write('      </Attribute>\n')

            # Functions (point data)
            for name, func in functions.items():
                # Determine if scalar or vector
                try:
                    value_shape = func.function_space.element.value_shape
                    if len(value_shape) == 0 or value_shape[0] == 1:
                        attr_type = "Scalar"
                        dimensions = f"{num_nodes} 1"
                    else:
                        attr_type = "Vector"
                        dimensions = f"{num_nodes} {value_shape[0]}"
                except:
                    # Fallback - assume vector
                    attr_type = "Vector"
                    dimensions = f"{num_nodes} {gdim}"

                xdmf.write(f'      <Attribute Name="{name}" ')
                xdmf.write(f'AttributeType="{attr_type}" Center="Node">\n')
                xdmf.write(f'        <DataItem Dimensions="{dimensions}" ')
                xdmf.write('Format="HDF">\n')
                xdmf.write(f'          {h5_file.name}:/Functions/{name}\n')
                xdmf.write('        </DataItem>\n')
                xdmf.write('      </Attribute>\n')

            xdmf.write('    </Grid>\n')
            xdmf.write('  </Domain>\n')
            xdmf.write('</Xdmf>\n')

        print(f"    DEBUG: XDMF file write complete!")
        print(f"    DEBUG: ========== File writing finished ==========")

    print(f"    DEBUG: Rank {rank} - waiting at barrier...")
    comm.Barrier()
    print(f"    DEBUG: Rank {rank} - passed barrier")

    if rank == 0:
        return xdmf_file, h5_file
    else:
        return None, None


def write_combined_paraview_output(
    output_dir,
    basename,
    mesh,
    cell_tags,
    fields,
    comm=MPI.COMM_WORLD
):
    """
    Convenience function to write combined ParaView output.

    Creates both:
    1. Individual XDMF files for each field (standard DOLFINx output)
    2. A combined XDMF file with all fields and cell tags on one grid

    Parameters
    ----------
    output_dir : Path
        Output directory
    basename : str
        Base name for files (e.g., "straight_wire")
    mesh : dolfinx.mesh.Mesh
        The mesh
    cell_tags : dolfinx.mesh.MeshTags
        Cell tags
    fields : dict
        Dictionary of {name: (original_func, lagrange_func)} pairs
    comm : MPI.Comm
        MPI communicator

    Returns
    -------
    dict
        Dictionary of written files
    """
    from dolfinx import io

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    written_files = {}

    # Write individual files (standard approach)
    for name, (original_func, lagrange_func) in fields.items():
        xdmf_path = output_dir / f"{basename}_{name}.xdmf"
        with io.XDMFFile(comm, xdmf_path, "w") as xdmf:
            xdmf.write_mesh(mesh)
            if cell_tags is not None:
                mesh.topology.create_connectivity(3, 3)
                xdmf.write_meshtags(cell_tags, mesh.geometry)
            xdmf.write_function(lagrange_func)
        written_files[name] = xdmf_path

    # Write combined file (custom - all on one grid)
    lagrange_funcs = {name: lag_func for name, (_, lag_func) in fields.items()}
    combined_path = output_dir / f"{basename}_combined"
    xdmf_file, h5_file = write_xdmf_with_tags(
        combined_path,
        mesh,
        cell_tags,
        lagrange_funcs,
        comm=comm
    )

    if xdmf_file:
        written_files['combined'] = xdmf_file

    return written_files
