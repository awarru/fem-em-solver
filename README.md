# FEM Electromagnetics Solver

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![FEniCSX](https://img.shields.io/badge/FEniCSX-0.7+-green.svg)](https://fenicsproject.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A finite element method (FEM) solver for electromagnetic simulations of MRI coils with gelled saline phantoms, built on [FEniCSX/DolfinX](https://fenicsproject.org/).

## Features

- **Magnetostatics**: Magnetic vector potential formulation for static problems
- **Time-Harmonic**: Full Maxwell equations for frequency-domain analysis
- **Coil Models**: Loop coils, birdcage coils, TEM coils
- **Material Models**: Complex permittivity, dispersion models, gelled saline phantoms
- **MRI-Focused**: B1+ mapping, SAR calculation, coil loading analysis
- **Validation**: Benchmarked against analytic solutions and Ansys HFSS

## Quick Start

### Installation

#### Option 1: Using Docker (Recommended)

```bash
# Build and start the container
cd docker
docker compose up -d

# Enter the container
docker compose exec fem-em-solver bash
```

**Prerequisites:**
- Docker installed and running
- Your user in the `docker` group (to avoid sudo)

**To add your user to docker group:**
```bash
sudo usermod -aG docker $USER
# IMPORTANT: Log out and log back in for changes to take effect
```

**Test if Docker works:**
```bash
docker ps  # Should show containers without sudo
```

**If you get permission errors** after adding to docker group, you need to log out and back in, or run:
```bash
newgrp docker  # Applies group change to current shell (temporary)
```

**Note:** Modern Docker uses `docker compose` (space). Older versions use `docker-compose` (hyphen).

#### Option 2: Conda Environment

```bash
conda create -n femem python=3.11
conda activate femem
conda install -c conda-forge fenics-dolfinx gmsh pyvista
pip install -e ".[dev,docs]"
```

### Run Examples (without entering Docker)

```bash
# List available examples
./run_examples.sh --list

# Run one example by number
./run_examples.sh --example 1

# Run multiple examples and set MPI ranks
./run_examples.sh --example 1,3 --nproc 4

# Run all examples
./run_examples.sh --example all --nproc 2
```

The script automatically targets `docker/docker-compose.yml` and runs each selected example as:
`docker compose exec fem-em-solver ... mpiexec -n <nproc> python3 <example>`

### Python API Example

```python
from fem_em_solver import MagnetostaticSolver
from fem_em_solver.coils import CircularLoop

# Create a circular loop coil
coil = CircularLoop(radius=0.05, current=1.0, position=(0, 0, 0))

# Set up solver
solver = MagnetostaticSolver(mesh_resolution=0.005)
solver.add_coil(coil)

# Solve
A = solver.solve()
B = solver.compute_b_field(A)

# Visualize
solver.plot_field(B, component='z')
```

## Project Structure

```
fem-em-solver/
├── src/fem_em_solver/     # Main package
│   ├── core/              # FEM formulations and solvers
│   ├── coils/             # Coil geometry definitions
│   ├── materials/         # Material property models
│   ├── post/              # Post-processing and analysis
│   └── io/                # Mesh and data I/O
├── examples/              # Tutorial notebooks and scripts
├── tests/                 # Test suite
├── docs/                  # Documentation
└── meshes/                # Pre-generated mesh files
```

## Development Phases

1. **Phase 0**: Infrastructure & setup
2. **Phase 1**: Magnetostatics foundation
3. **Phase 2**: Time-harmonic Maxwell equations
4. **Phase 3**: Material models & phantoms
5. **Phase 4**: Coil modeling
6. **Phase 5**: Full MRI system integration
7. **Phase 6**: Advanced features

See [PROJECT_PLAN.md](PROJECT_PLAN.md) for detailed roadmap.

## Validation

The solver is validated against:
- Analytic solutions (wire, loop, Helmholtz coil)
- Commercial software (Ansys HFSS)
- Literature data and measured results

See `docs/validation/` for validation reports.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{fem_em_solver,
  author = {Awarru},
  title = {FEM Electromagnetics Solver for MRI},
  url = {https://github.com/awarru/fem-em-solver},
  year = {2026}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [FEniCS Project](https://fenicsproject.org/) for the excellent FEM framework
- [Gmsh](https://gmsh.info/) for mesh generation
- MRI research community for validation data and phantom specifications
