# Project: FEM Electromagnetics Solver for MRI Coil Simulation

## Overview
A finite element method (FEM) solver built on FEniCSX/DolfinX for electromagnetic simulations of MRI coils with gelled saline phantoms. The solver will handle magnetostatics and time-harmonic Maxwell equations with support for complex, frequency-dependent material properties.

## Motivation
- Commercial EM solvers (Ansys HFSS, CST) are expensive and black-box
- Open-source alternatives (Elmer, OpenFOAM EM) lack MRI-specific features
- Academic codes often use legacy frameworks (legacy FEniCS, custom C++ solvers)
- Need a modern, Python-based solver for research and educational purposes

---

## Project Phases

### Phase 0: Infrastructure & Setup (Week 1-2)
**Goal**: Establish development environment, CI/CD, and basic repository structure

**Deliverables**:
- [ ] Repository structure with proper Python packaging
- [ ] Docker/environment setup with FEniCSX 0.7+ 
- [ ] Continuous integration (GitHub Actions) for testing
- [ ] Documentation infrastructure (Sphinx/MkDocs)
- [ ] Unit test framework with pytest
- [ ] Basic mesh generation utilities

**Key Technologies**:
- dolfinx 0.7.x (latest stable)
- gmsh (mesh generation)
- pyvista (visualization)
- pytest + pytest-xdist (testing)
- Docker/Docker Compose

---

### Phase 1: Magnetostatics Foundation (Week 3-6)
**Goal**: Implement magnetostatic solver with basic test cases

**Deliverables**:
- [ ] Magnetostatic formulation (B-formulation: ∇×(μ⁻¹∇×A) = J)
- [ ] Dirichlet/Neumann boundary conditions
- [ ] Current source implementation (thin wire approximation)
- [ ] **Validation Case 1**: Magnetic field of straight wire (analytic solution)
- [ ] **Validation Case 2**: Magnetic field of circular loop (analytic solution)
- [ ] **Validation Case 3**: Helmholtz coil (analytic solution)
- [ ] Post-processing: B-field, H-field visualization
- [ ] Convergence studies (h-refinement, p-refinement)

**Mathematical Details**:
```
Weak form: ∫ μ⁻¹(∇×A)·(∇×v) dx = ∫ J·v dx
Where A is magnetic vector potential, J is current density
B = ∇×A, H = μ⁻¹B
```

**Validation Metrics**:
- L2 error vs analytic solution
- H(curl) semi-norm error
- Field uniformity in Helmholtz coil

---

### Phase 2: Time-Harmonic Maxwell's Equations (Week 7-10)
**Goal**: Extend to time-harmonic solver for frequency-domain analysis

**Deliverables**:
- [ ] Full Maxwell formulation: ∇×(μ⁻¹∇×E) - ω²εE = -iωJ
- [ ] Complex material properties (ε = ε' - jε'', μ = μ' - jμ'')
- [ ] Perfect Electric Conductor (PEC) boundary conditions
- [ ] Absorbing Boundary Conditions (ABC) or PML (Perfectly Matched Layer)
- [ ] **Validation Case 4**: Plane wave in free space
- [ ] **Validation Case 5**: Isotropic antenna/radiating dipole
- [ ] **Validation Case 6**: Cylindrical waveguide (TM/TE modes)
- [ ] S-parameter extraction (1-port, 2-port)
- [ ] Impedance calculation

**Mathematical Details**:
```
Time-harmonic E-field formulation:
∇×(μᵣ⁻¹∇×E) - k₀²εᵣE = 0
Where k₀ = ω√(ε₀μ₀), ω = 2πf

Lossy materials: εᵣ = ε' - j(σ/ωε₀)
```

**HFSS Comparison**:
- Export same geometry to HFSS
- Compare E-field magnitude/phase
- Compare H-field patterns
- Compare impedance (Z11)
- Compare quality factor (Q)

---

### Phase 3: Material Models & Phantoms (Week 11-14)
**Goal**: Implement realistic material properties for MRI applications

**Deliverables**:
- [ ] Debye/Drude dispersion models for dielectrics
- [ ] Cole-Cole model for biological tissues
- [ ] Gelled saline phantom model (known conductivity, permittivity)
- [ ] Temperature-dependent conductivity
- [ ] **Validation Case 7**: Sphere phantom in uniform B-field
- [ ] **Validation Case 8**: Cylindrical phantom with birdcage coil
- [ ] SAR (Specific Absorption Rate) calculation
- [ ] Power deposition analysis

**Material Properties (Gelled Saline)**:
```
Typical values at 128 MHz (3T MRI):
- Conductivity σ: 0.6-0.9 S/m (tunable with salt concentration)
- Relative permittivity εᵣ: 78-80
- Phantom diameter: 16-20 cm (head), 30-40 cm (body)
- Gel concentration: ~1% agarose
```

**Phantom Geometries**:
- Spherical phantom (MRI QA standards)
- Cylindrical phantom (head/body)
- Anthropomorphic head model
- Multi-compartment phantoms

---

### Phase 4: Coil Modeling (Week 15-18)
**Goal**: Implement various MRI coil types

**Deliverables**:
- [ ] Lumped element modeling (capacitors, inductors)
- [ ] Tuning and matching circuits
- [ ] **Coil Type 1**: Circular loop coil
- [ ] **Coil Type 2**: Figure-8 coil
- [ ] **Coil Type 3**: Birdcage coil (low-pass, high-pass, band-pass)
- [ ] **Coil Type 4**: TEM coil
- [ ] **Coil Type 5**: Surface array coil
- [ ] Multi-channel excitation (phase/amplitude control)
- [ ] B1+ field mapping
- [ ] B1- field (receive sensitivity)

**Coil Details**:
```
Birdcage Coil:
- N equally spaced legs
- End rings with capacitors
- Resonant frequency: f₀ = 1/(2π√(LC))
- Tuning: adjust capacitors
- Matching: match impedance to 50Ω

TEM Coil:
- Transverse electromagnetic mode
- High-pass transmission line structure
- Better for high-field (≥7T)
```

---

### Phase 5: Full MRI System Integration (Week 19-22)
**Goal**: Complete MRI coil + phantom simulations

**Deliverables**:
- [ ] Realistic 3T head coil (16-rung birdcage)
- [ ] Realistic 3T body coil
- [ ] Realistic 7T head coil (TEM or dipole array)
- [ ] Phantom loading effects (frequency shift, Q degradation)
- [ ] B1+ homogeneity analysis
- [ ] SAR hotspot identification
- [ ] Safety compliance (IEC/FDA limits)
- [ ] Publication-quality visualizations

**Quantitative Metrics**:
- B1+ efficiency (μT/√W)
- B1+ homogeneity (coefficient of variation)
- SAR10g (max local SAR averaged over 10g)
- Whole-body SAR
- Power efficiency

---

### Phase 6: Advanced Features (Week 23-26)
**Goal**: Professional-grade features

**Deliverables**:
- [ ] Parallel computing (MPI) for large meshes
- [ ] Adaptive mesh refinement
- [ ] Parameter sweeps (frequency, geometry)
- [ ] Sensitivity analysis
- [ ] Optimization framework (coil geometry)
- [ ] GPU acceleration (if feasible with FEniCSX)
- [ ] Mesh import from DICOM/STEP files
- [ ] Complete documentation with tutorials
- [ ] Publication of validation results

---

## Repository Structure

```
fem-em-solver/
├── README.md                      # Project overview, quickstart
├── LICENSE                        # Open source license (MIT/GPL)
├── pyproject.toml                 # Python package config
├── requirements.txt               # Dependencies
├── docker/
│   ├── Dockerfile                 # FEniCSX development image
│   └── docker-compose.yml         # Easy container setup
│
├── docs/                          # Documentation
│   ├── index.md                   # Main docs
│   ├── theory/                    # Mathematical theory
│   │   ├── magnetostatics.md
│   │   ├── time_harmonic.md
│   │   └── boundary_conditions.md
│   ├── tutorials/                 # Step-by-step tutorials
│   │   ├── 01_magnetostatics.md
│   │   ├── 02_coil_simulation.md
│   │   └── 03_phantom_analysis.md
│   └── validation/                # Validation reports
│       ├── magnetostatics.md
│       ├── helmholtz_coil.md
│       └── birdcage_comparison.md
│
├── src/
│   └── fem_em_solver/             # Main package
│       ├── __init__.py
│       ├── core/                  # Core FEM functionality
│       │   ├── __init__.py
│       │   ├── forms.py           # Weak forms
│       │   ├── solvers.py         # Solver classes
│       │   └── function_spaces.py # H(curl), H(div), etc.
│       │
│       ├── materials/             # Material models
│       │   ├── __init__.py
│       │   ├── dispersive.py      # Debye/Drude models
│       │   ├── biological.py      # Cole-Cole for tissue
│       │   └── phantoms.py        # Gelled saline models
│       │
│       ├── coils/                 # Coil geometries
│       │   ├── __init__.py
│       │   ├── loop.py            # Circular loops
│       │   ├── birdcage.py        # Birdcage coils
│       │   ├── tem.py             # TEM coils
│       │   └── circuits.py        # Lumped elements
│       │
│       ├── io/                    # Input/output
│       │   ├── __init__.py
│       │   ├── mesh.py            # Mesh utilities
│       │   ├── gmsh_interface.py  # Gmsh integration
│       │   ├── hfss_export.py     # Export for HFSS comparison
│       │   └── vtk_export.py      # Paraview output
│       │
│       ├── post/                  # Post-processing
│       │   ├── __init__.py
│       │   ├── fields.py          # Field calculations
│       │   ├── sar.py             # SAR computation
│       │   ├── impedance.py       # S/Z/Y parameters
│       │   └── visualization.py   # Plotting with pyvista
│       │
│       └── utils/                 # Utilities
│           ├── __init__.py
│           ├── constants.py       # Physical constants
│           ├── units.py           # Unit conversions
│           └── convergence.py     # Error analysis
│
├── tests/                         # Test suite
│   ├── __init__.py
│   ├── conftest.py                # Pytest fixtures
│   ├── unit/                      # Unit tests
│   │   ├── test_forms.py
│   │   └── test_materials.py
│   ├── integration/               # Integration tests
│   │   ├── test_magnetostatics.py
│   │   └── test_time_harmonic.py
│   └── validation/                # Validation tests
│       ├── test_wire_field.py
│       ├── test_helmholtz.py
│       └── test_phantom.py
│
├── examples/                      # Example notebooks/scripts
│   ├── magnetostatics/
│   │   ├── 01_straight_wire.py
│   │   ├── 02_circular_loop.py
│   │   └── 03_helmholtz_coil.py
│   ├── time_harmonic/
│   │   ├── 01_dipole_antenna.py
│   │   └── 02_waveguide.py
│   ├── coils/
│   │   ├── 01_loop_coil.py
│   │   ├── 02_birdcage_coil.py
│   │   └── 03_tem_coil.py
│   └── mri/
│       ├── 01_phantom_only.py
│       ├── 02_birdcage_loaded.py
│       └── 03_full_head_model.py
│
├── meshes/                        # Mesh files
│   ├── simple/
│   ├── coils/
│   └── phantoms/
│
├── data/                          # Reference data
│   ├── hfss_comparisons/
│   ├── analytical_solutions/
│   └── literature_data/
│
└── scripts/                       # Utility scripts
    ├── run_tests.sh
    ├── build_docs.sh
    └── generate_meshes.py
```

---

## Key Technical Considerations

### 1. Function Spaces
- **H(curl)** for electric field E and magnetic vector potential A
- **H(div)** for magnetic flux density B
- **L2** for scalar potential (if using A-V formulation)
- Nedelec elements (1st and 2nd order)

### 2. Boundary Conditions
- **PEC**: n×E = 0 (perfect conductor)
- **PMC**: n·B = 0 (magnetic symmetry)
- **ABC**: Absorbing boundary for radiation
- **PML**: Perfectly matched layer (complex coordinate stretching)
- **Port**: Waveguide ports for S-parameter extraction

### 3. Material Handling
- Complex permittivity: ε = ε₀(ε' - jε'')
- Complex permeability: μ = μ₀(μ' - jμ'')
- Anisotropic materials (tensor ε, μ)
- Frequency-dependent dispersion

### 4. Linear Solvers
- Direct: MUMPS (for small/medium problems)
- Iterative: GMRES + ILU preconditioner
- Multigrid for large problems
- Complex linear systems (for time-harmonic)

### 5. Post-Processing
- Field interpolation on regular grids
- SAR = σ|E|²/(2ρ) [W/kg]
- Current density J = σE
- Poynting vector S = ½Re(E×H*)

---

## Development Workflow

1. **Branch per feature**: `feature/magnetostatics`, `feature/coils/birdcage`
2. **Test-driven**: Write test → implement → validate
3. **HFSS comparison**: Export identical geometry, compare quantitatively
4. **Documentation**: Docstrings + notebooks + theory docs
5. **CI/CD**: Run tests on every commit, build docs on merge

---

## Success Criteria

### Minimum Viable Product (End of Phase 2):
- [ ] Solver runs Helmholtz coil simulation
- [ ] Results match analytic solution to <5%
- [ ] Results match HFSS to <10%

### Target (End of Phase 4):
- [ ] Birdcage coil with phantom simulation
- [ ] B1+ field matches literature/measured data
- [ ] SAR patterns are physically reasonable

### Stretch Goal (End of Phase 6):
- [ ] Multi-channel coil optimization
- [ ] Publication in relevant journal
- [ ] Community adoption (GitHub stars, issues, PRs)

---

## Resources

### Documentation
- [FEniCSX Tutorial](https://jorgensd.github.io/dolfinx-tutorial/)
- [FEniCSX API Docs](https://docs.fenicsproject.org/dolfinx/v0.7.0/)
- [IEEE Standards for Phantoms](https://standards.ieee.org/)
- [HFSS Theory Manual](https://ansyshelp.ansys.com/)

### Similar Projects
- [Elmer FEM](https://www.elmerfem.org/) (has EM module)
- [OpenEMS](https://openems.de/) (FDTD, not FEM)
- [scikit-rf](https://scikit-rf.readthedocs.io/) (RF network analysis)

---

Ready to proceed with repository creation and Phase 0 implementation?
