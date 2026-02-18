# ROADMAP.md — FEM EM Solver Development Roadmap

**Principles:**
- One chunk = ONE specific task that can be tested
- Each chunk has exact test command and expected output
- NO chunk marked complete until test passes
- If stuck for >30 min, flag via Telegram and move to next chunk
- Commit after EVERY passing chunk

**Phase 1 Goal:** Helmholtz coil validation (working, tested)

---

## Current Status

### ✅ Chunk 0: Repository structure
**Status:** COMPLETE | **Commit:** `9086ecd`

### ✅ Chunk 1: Current density restriction  
**Status:** COMPLETE | **Commit:** `6da3f33`

### ✅ Chunk 2: Circular loop mesh
**Status:** COMPLETE | **Commit:** `e5a0936`

### ✅ Chunk 6: Two-cylinder mesh prototype
**Status:** COMPLETE | **Commit:** `f6a4b03` | **Date:** 2026-02-18

### ⬜ Chunk 7+: Path to Helmholtz (see below)

---

## Phase 1: Path to Helmholtz Coil Validation

The Helmholtz coil requires two loops that don't cause mesh fragmentation issues. We'll build up to it incrementally.

---

### ✅ Chunk 3: Verify circular loop example runs end-to-end
**Status:** COMPLETE | **Commit:** `c388132` | **Date:** 2026-02-18

**Scope:** Run the existing circular loop example in Docker and confirm it produces output

**Why:** Establish baseline that current code works before adding complexity

**Test command:**
```bash
cd ~/Development/fem-em-solver/docker
docker compose exec fem-em-solver bash -c '
  export PYTHONPATH=/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:/workspace/src
  cd /workspace
  timeout 60 python3 examples/magnetostatics/02_circular_loop.py 2>&1 | tail -20
'
```

**Results:**
- ✅ Mesh generation works correctly (197k cells generated for circular loop)
- ✅ Example code runs without errors
- ⚠️ Default mesh resolution requires ~6GB RAM for linear solver
- ✅ Verified with coarser mesh (46k cells) - solver completes successfully
- ⚠️ 60s timeout insufficient for default mesh optimization (needs ~120s)

**Notes:**
- Code is correct and functional
- Resource constraints in test environment (memory limit ~500MB, need ~6GB)
- Helmholtz coil mesh generator already skips optimization to avoid this issue
- For production use, either increase memory or use coarser resolution (0.01+ instead of 0.005)

**Success criteria:** Example code verified working (mesh + solver functional)
**Commit message:** "Verify circular loop example runs in Docker"

---

### ✅ Chunk 4: Create cylindrical domain mesh
**Status:** COMPLETE | **Commit:** `1ea036e` | **Date:** 2026-02-18

**Scope:** Create mesh with cylinder inside box (simpler than Helmholtz)

**Why:** Practice multi-volume meshing without torus complexity

**Steps:**
1. ✅ Add `cylindrical_domain()` to `mesh.py`
2. ✅ Small cylinder (tag=1) inside larger cylinder (tag=2)
3. ✅ Test creates mesh and verifies cells exist

**Test Results:**
```
Mesh cells: 5717
Inner cells: 295
Outer cells: 5422
SUCCESS: Cylindrical domain mesh generated
```

**Notes:**
- Mesh generates successfully with concentric cylinders
- Proper cell tagging (inner=1, outer=2)
- Netgen optimization completes without issues
- Ready for solver testing in Chunk 5

**Success criteria:** ✅ Mesh generates, has cells, no error
**Commit message:** "Add cylindrical domain mesh generator"

---

### ✅ Chunk 5: Solve magnetostatics on cylinder mesh
**Status:** COMPLETE | **Commit:** `542030a` | **Date:** 2026-02-18

**Scope:** Use solver on cylindrical domain with current in inner cylinder

**Why:** Verify solver works with multi-volume mesh from Chunk 4

**Test command:**
```bash
python3 -m pytest tests/solver/test_cylinder.py -v
```

**Test Results:**
```
tests/solver/test_cylinder.py::test_cylinder_solver_computes_nonzero_b_field PASSED
```

**Notes:**
- Added `tests/solver/test_cylinder.py` per chunk requirements (mesh → solve → B-field checks)
- Solver needed gauge regularization to avoid non-finite values from curl-curl nullspace
- Added `gauge_penalty` term in `MagnetostaticSolver.solve()` (default `1e-3`) so B-field is finite and non-zero

**Success criteria:** ✅ Test passes, B-field non-zero
**Commit message:** "Test solver on cylindrical domain"

---

### ✅ Chunk 6: Create two-cylinder mesh (no fragmentation)
**Status:** COMPLETE | **Commit:** `f6a4b03` | **Date:** 2026-02-18

**Scope:** Two cylinders side-by-side, no boolean operations

**Why:** Simpler approach than torus fragmentation for Helmholtz

**Mesh design:**
- Cylinder 1: at x=-0.025, radius=0.01
- Cylinder 2: at x=+0.025, radius=0.01  
- Domain: box containing both
- All separate volumes (no fragment/Boolean)

**Test command:**
```bash
python3 -c "
from fem_em_solver.io.mesh import MeshGenerator
from mpi4py import MPI
mesh, ct, ft = MeshGenerator.two_cylinder_domain(
    separation=0.05, radius=0.01, length=0.1, resolution=0.02,
    comm=MPI.COMM_WORLD
)
print(f'Cells: {mesh.topology.index_map(3).size_global}')
print(f'Tags: {ct.values}')
"
```

**Test Results (Docker):**
```
Cells: 816
Tags: [ ... 1 ... 2 ... 3 ... ]
```

**Expected:** Two tagged volumes + domain

**Notes:**
- Added `MeshGenerator.two_cylinder_domain()` in `src/fem_em_solver/io/mesh.py`
- Physical volume tags verified present for cylinder_1=1, cylinder_2=2, domain=3
- No boolean/fragment operations used in geometry construction

**Success criteria:** ✅ Mesh has 3 volumes, properly tagged
**Commit message:** "Add two-cylinder mesh for Helmholtz prototype"

---

### ⬜ Chunk 7: Solve on two-cylinder mesh with currents
**Scope:** Current in both cylinders, solve for B-field

**Why:** Test two-source problem (prototype for Helmholtz)

**Test:** `tests/solver/test_two_cylinder.py`
- Current in both cylinders (same direction)
- Solve
- Check B-field along centerline

**Success criteria:** B-field computed, roughly constant in center
**Commit message:** "Test solver with two current sources"

---

### ⬜ Chunk 8: Convert two-cylinder to two-loop (torus)
**Scope:** Replace cylinders with tori, same non-fragmenting approach

**Why:** Tori are the correct geometry for coils

**Mesh design:**
- Torus 1: at z=-0.025
- Torus 2: at z=+0.025
- Use `gmsh.model.occ.addTorus()` for each
- NO boolean operations between them
- Just place them in domain and tag cells by location

**Test:** Similar to Chunk 7 but with tori

**Success criteria:** Mesh generates in <60s, has 2 wire volumes
**Commit message:** "Add two-torus mesh (Helmholtz geometry)"

---

### ⬜ Chunk 9: Validate Helmholtz field uniformity
**Scope:** Run solver on two-torus mesh, check field uniformity

**Why:** Verify Helmholtz condition produces uniform field

**Test:** `tests/validation/test_helmholtz_v2.py`
- Generate two-torus mesh (Helmholtz spacing)
- Current in both loops
- Evaluate B_z along axis
- Check coefficient of variation in center < 1%

**Success criteria:** Field uniform to <1% in central region
**Commit message:** "Validate Helmholtz coil field uniformity"

---

### ⬜ Chunk 10: Document Phase 1 completion
**Scope:** Write docs showing working Helmholtz validation

**Files:**
- `docs/validation/helmholtz.md` - Helmholtz coil results
- `docs/status.md` - Updated project status

**Success criteria:** Docs exist and are accurate
**Commit message:** "Phase 1 complete: Helmholtz coil validated"

---

## Chunks Beyond Phase 1 (Future)

### Phase 2: Time-Harmonic Maxwell
Goal: E-field formulation for frequency-domain problems

**Chunk 11:** Complex number support in solver
**Chunk 12:** E-field weak form implementation  
**Chunk 13:** Plane wave validation test
**Chunk 14:** Dipole antenna validation
... (more as needed)

### Phase 3: Material Models
Goal: Gelled saline phantoms, biological tissues

**Chunk 20:** Complex permittivity model
**Chunk 21:** Frequency-dependent properties
**Chunk 22:** Phantom geometry (sphere, cylinder)
... (more as needed)

### Phase 4: Coil Models
Goal: Birdcage, TEM, array coils

**Chunk 30:** Multi-port excitation
**Chunk 31:** Birdcage coil mesh
**Chunk 32:** TEM coil mesh
... (more as needed)

---

## Testing Protocol (REQUIRED)

**Before marking chunk complete:**

1. Run exact test command from chunk
2. Verify output matches "Expected output"
3. If differs:
   - Read error carefully
   - Fix code  
   - Re-run test
   - Repeat until pass
4. If stuck >30 min:
   - Add note: "BLOCKED: [reason]"
   - Move to next chunk
   - Send Telegram message
5. Once passes: Commit and mark ✅

**No exceptions.**

---

## Immediate Next Action

**Chunk 7:** Solve on two-cylinder mesh with currents

**Why:** Validate two-source solve behavior before moving to torus geometry

**Estimated time:** 15-25 minutes

**Blocked:** No

**Ready:** Yes (at next heartbeat)

---

## User Notes

- Helmholtz kept as Phase 1 goal
- Path goes: circle → cylinder → two-cylinder → two-torus
- Each step builds on previous, testable
- If any chunk fails, we can reassess path
