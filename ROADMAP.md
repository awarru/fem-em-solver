# ROADMAP.md — FEM EM Solver Development Roadmap

This file tracks the next actionable chunks of work. Each chunk is a single
session's worth of effort. When a chunk is completed, mark it ✅ and move on
to the next one. The heartbeat system reads this file to know what to do next.

---

## Phase 1: Magnetostatics Foundation

### ✅ Chunk 0: Initial solver structure
- Implemented MagnetostaticSolver with A-formulation (H(curl) Nedelec elements)
- Created straight wire mesh generation via Gmsh
- Added analytical solutions module (wire, loop, Helmholtz)
- Implemented error metrics (L2, max relative)
- Created validation test structure and example script
- **Committed:** `9086ecd` (2026-02-17)

### ⬜ Chunk 1: Fix current density restriction to wire subdomain
**What:** Current density J is applied to the whole mesh domain, not just the wire volume.
This causes incorrect results because J should be zero outside the wire.

**How to fix:**
1. Open `src/fem_em_solver/core/solvers.py`
2. Modify `solve()` to accept `cell_tags` and a `subdomain_id` parameter
3. Use `ufl.Measure("dx")` with subdomain data from `cell_tags` to restrict the integral:
   ```python
   dx_wire = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags, subdomain_id=wire_tag)
   L = inner(J, v) * dx_wire
   ```
4. The bilinear form (left side) still integrates over the whole domain
5. Update `examples/magnetostatics/01_straight_wire.py` to pass cell_tags
6. Run: `cd ~/Development/fem-em-solver && PYTHONPATH=src python3 examples/magnetostatics/01_straight_wire.py`
7. Verify the B-field magnitude matches analytical solution within ~5%

**Success criteria:** Validation test passes with relative L2 error < 10%
**Commit message:** `"Phase 1: Restrict current density to wire subdomain using cell_tags"`

### ⬜ Chunk 2: Circular loop mesh and validation
**What:** Add mesh generation for a circular current loop and validate B_z on axis.

**How:**
1. Add `circular_loop_domain()` to `src/fem_em_solver/io/mesh.py`
   - Create a torus (ring) for the wire using `gmsh.model.occ.addTorus()`
   - Surround with a spherical or cylindrical air domain
   - Tag the torus volume as "wire" (tag=1) and air as "domain" (tag=2)
2. Create `tests/validation/test_circular_loop.py`
   - Set up loop: radius=0.05m, current=1A, wire cross-section radius=0.001m
   - Evaluate B_z along the z-axis from z=-0.1 to z=0.1
   - Compare with `AnalyticalSolutions.circular_loop_magnetic_field_on_axis()`
   - Assert relative L2 error < 10%
3. Create `examples/magnetostatics/02_circular_loop.py`
   - Plot B_z(z) numerical vs analytical
   - Save plot as `circular_loop_validation.png`

**Success criteria:** On-axis B_z matches analytical solution within 10%
**Commit message:** `"Phase 1: Circular loop mesh generation and on-axis validation"`

### ⬜ Chunk 3: Helmholtz coil validation
**What:** Two loops separated by one radius. Test field uniformity in center.

**How:**
1. Create `helmholtz_coil_domain()` in `src/fem_em_solver/io/mesh.py`
   - Two tori at z = -R/2 and z = +R/2 (where R is loop radius)
   - Same approach as circular loop but with two wire volumes
   - Both wires tagged as "wire" (tag=1)
2. Create `tests/validation/test_helmholtz.py`
   - Radius = 0.05m, current = 1A in both loops
   - Evaluate B_z along z-axis
   - Compare with `AnalyticalSolutions.helmholtz_coil_field_on_axis()`
   - Check uniformity: coefficient of variation of B_z in central 20% should be < 1%
3. Create `examples/magnetostatics/03_helmholtz_coil.py`

**Success criteria:** B_z uniform to <1% in central region, matches analytical <10%
**Commit message:** `"Phase 1: Helmholtz coil mesh and field uniformity validation"`

### ⬜ Chunk 4: Convergence study infrastructure
**What:** Prove the solver converges at the expected rate as mesh refines.

**How:**
1. Create `tests/validation/test_convergence.py`
2. For straight wire problem, run at resolutions [0.02, 0.01, 0.007, 0.005, 0.003]
3. Compute L2 error at each resolution
4. Fit log(error) vs log(h) to get convergence rate
5. For Nedelec order 1: expect rate ~1.0
6. For Nedelec order 2: expect rate ~2.0
7. Test both orders
8. Save convergence plot to `docs/validation/convergence_plot.png`

**Success criteria:** Convergence rate > 0.8 for order 1, > 1.5 for order 2
**Commit message:** `"Phase 1: h-refinement convergence study for magnetostatics"`

### ⬜ Chunk 5: Phase 1 documentation and cleanup
**What:** Write up results and clean code for Phase 1 completion.

**How:**
1. Write `docs/validation/magnetostatics.md` — validation report with:
   - Problem description and formulation
   - Mesh details
   - Results tables (error vs resolution)
   - Convergence plots
   - Comparison figures
2. Ensure all public methods have docstrings
3. Run full test suite: `pytest tests/ -v`
4. Ensure all tests pass
5. Update `PROJECT_PLAN.md` — mark Phase 1 deliverables as complete
6. Update `README.md` if needed

**Success criteria:** All tests pass, docs written, Phase 1 marked complete
**Commit message:** `"Phase 1 complete: magnetostatics validated against analytical solutions"`

---

## Phase 2: Time-Harmonic Maxwell (future — do not start until Phase 1 is done)

### ⬜ Chunk 6: Time-harmonic E-field formulation
### ⬜ Chunk 7: PEC boundary conditions and plane wave test
### ⬜ Chunk 8: Dipole antenna validation
### ⬜ Chunk 9: Cylindrical waveguide modes
### ⬜ Chunk 10: HFSS comparison infrastructure

---

## How to use this file

1. Find the first ⬜ chunk
2. Read its instructions carefully
3. Do exactly what it says
4. When done, change ⬜ to ✅ and add the commit hash and date
5. Commit and push
6. Stop — do not start the next chunk in the same session
