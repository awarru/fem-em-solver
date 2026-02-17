# ROADMAP.md — FEM EM Solver Development Roadmap

**Principles:**
- One chunk = ONE specific task (e.g., "create and verify a cube mesh")
- Each chunk has exact test command and expected output
- NO chunk is marked complete until test passes
- If stuck for >30 min, flag via Telegram and move to next chunk
- Commit after EVERY passing chunk

---

## Current Status

### ✅ Chunk 0: Repository structure
**Status:** COMPLETE
**Commit:** `9086ecd`
**Test:** `git log --oneline -1` shows commit

### ✅ Chunk 1: Current density restriction  
**Status:** COMPLETE
**Commit:** `6da3f33`
**Test:** Code review only (solver API change)

### ✅ Chunk 2: Circular loop mesh
**Status:** COMPLETE  
**Commit:** `e5a0936`
**Test:** `python3 -c "from fem_em_solver.io.mesh import MeshGenerator; print('OK')"`

### ⬜ Chunk 3-6: KNOWN ISSUES (see below)

---

## KNOWN ISSUES (Flagged for Review)

### ⬜ ISSUE-1: Helmholtz coil mesh hangs
**Status:** CODE WRITTEN, DOES NOT WORK
**File:** `src/fem_em_solver/io/mesh.py::helmholtz_coil_domain()`
**Problem:** Two tori fragmentation creates 1.6M elements, dolfinx conversion hangs
**Tested:** 2026-02-17, times out after 120s
**Options:** 
1. Use box domain instead of sphere (simpler)
2. Don't fragment - use single wire volume
3. Remove Helmholtz from Phase 1 scope
**Recommendation:** Option 3 - Helmholtz not critical for MRI coil goals

### ⬜ ISSUE-2: Convergence test too slow  
**Status:** TEST SKELETON WRITTEN, SOLVER TOO SLOW
**File:** `tests/validation/test_convergence.py`
**Problem:** Straight wire mesh + solve takes >60s per resolution
**Tested:** 2026-02-17, times out on 2 resolutions
**Options:**
1. Use built-in dolfinx mesh (no Gmsh)
2. Coarser mesh + fewer resolutions  
3. Skip convergence tests for now
**Recommendation:** Option 2 - simplify and retest

**ACTION REQUIRED:** User to decide on Issues 1 and 2

---

## Phase 1 Revised: Working Magnetostatics (Small Chunks)

Goal: Have working, tested magnetostatic solver with basic validation.

### ⬜ Chunk 7: Create simple cube mesh with Gmsh
**Scope:** Use Gmsh Python API to create a cube and export to .msh file

**Why:** Learn Gmsh API, establish mesh testing pattern

**Steps:**
1. Create `tests/mesh/test_cube.py`
2. Write test that:
   ```python
   import gmsh
   gmsh.initialize()
   gmsh.model.add("cube")
   box = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
   gmsh.model.occ.synchronize()
   gmsh.model.mesh.generate(3)
   gmsh.write("/tmp/test_cube.msh")
   gmsh.finalize()
   ```
3. Assert file `/tmp/test_cube.msh` exists and has >100 lines

**Test command:**
```bash
cd ~/Development/fem-em-solver/docker
docker compose exec fem-em-solver bash -c '
  export PYTHONPATH=/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:/workspace/src
  cd /workspace
  python3 -m pytest tests/mesh/test_cube.py -v
'
```

**Expected output:**
```
tests/mesh/test_cube.py::test_gmsh_cube PASSED
```

**Success criteria:** Test passes, cube.msh file created
**Commit message:** "Add Gmsh cube mesh test"

---

### ⬜ Chunk 8: Convert cube mesh to dolfinx
**Scope:** Load .msh file into dolfinx mesh

**Why:** Verify Gmsh → dolfinx pipeline works

**Steps:**
1. Add to `tests/mesh/test_cube.py`:
   ```python
   from dolfinx.io import gmshio
   from mpi4py import MPI
   
   def test_cube_to_dolfinx():
       mesh, cell_tags, facet_tags = gmshio.read_from_msh(
           "/tmp/test_cube.msh", MPI.COMM_WORLD, 0, gdim=3
       )
       assert mesh.topology.index_map(3).size_global > 0
   ```

**Test command:** Same as Chunk 7

**Expected output:**
```
tests/mesh/test_cube.py::test_gmsh_cube PASSED
tests/mesh/test_cube.py::test_cube_to_dolfinx PASSED
```

**Success criteria:** Both tests pass, dolfinx mesh has cells
**Commit message:** "Add dolfinx mesh conversion test"

---

### ⬜ Chunk 9: Solve magnetostatics on cube mesh (no source)
**Scope:** Run solver on cube with J=0, verify it doesn't crash

**Why:** Test solver works with simple dolfinx mesh

**Steps:**
1. Create `tests/solver/test_cube_solver.py`
2. Use `MeshGenerator.create_simple_box()` to make mesh
3. Create MagnetostaticProblem with no current
4. Call solver.solve() 
5. Verify solver._solved is True

**Test command:**
```bash
python3 -m pytest tests/solver/test_cube_solver.py -v
```

**Expected output:**
```
tests/solver/test_cube_solver.py::test_solve_no_source PASSED
```

**Success criteria:** Solver completes without error
**Commit message:** "Test solver on simple cube mesh"

---

### ⬜ Chunk 10: Solve with uniform current on cube
**Scope:** Add constant J to cube, solve, check B-field is reasonable

**Steps:**
1. Create cube mesh
2. Define uniform J = [0, 0, 1]
3. Solve
4. Compute B-field
5. Check B-field is not all zeros (has some magnitude)

**Test command:**
```bash
python3 -m pytest tests/solver/test_cube_with_current.py -v
```

**Expected output:**
```
tests/solver/test_cube_with_current.py::test_uniform_current PASSED
```

**Success criteria:** B-field computed and has non-zero values
**Commit message:** "Test solver with uniform current"

---

### ⬜ Chunk 11: Document what works
**Scope:** Write docs showing verified working features

**Steps:**
1. Create `docs/status.md`
2. List working features:
   - Repository structure
   - Docker setup
   - Straight wire mesh generation (Gmsh)
   - Circular loop mesh generation (Gmsh)
   - Straight wire validation test
   - Circular loop validation test
   - Cube mesh (built-in and Gmsh)
   - Solver on simple meshes
3. List known issues:
   - Helmholtz coil mesh too complex
   - Convergence tests too slow
4. Add example commands that work

**Test:** Manual review

**Success criteria:** Document exists and is accurate
**Commit message:** "Add project status documentation"

---

### ⬜ Chunk 12: Clean up broken code
**Scope:** Remove or disable code that doesn't work

**Files to address:**
1. `examples/magnetostatics/03_helmholtz_coil.py` - move to `examples/broken/`
2. `tests/validation/test_convergence.py` - mark all tests as skip with reason
3. Add comments explaining why

**Test:**
```bash
python3 -m pytest tests/ -v 2>&1 | grep -E "(PASSED|FAILED|SKIPPED)"
```

**Expected:** No FAILED, broken tests show as SKIPPED

**Success criteria:** Clean test output, broken code isolated
**Commit message:** "Isolate broken code, document known issues"

---

### ⬜ Chunk 13: Final Phase 1 verification
**Scope:** Run everything, confirm working state

**Steps:**
1. Run all tests:
   ```bash
   python3 -m pytest tests/ -v --tb=short 2>&1 | tail -20
   ```
2. Run working examples:
   ```bash
   python3 examples/magnetostatics/01_straight_wire.py
   python3 examples/magnetostatics/02_circular_loop.py
   ```
3. Verify imports:
   ```bash
   python3 -c "import fem_em_solver; print('OK')"
   ```

**Success criteria:**
- All tests pass or skip cleanly (no failures)
- Examples run without error
- Package imports correctly

**Commit message:** "Phase 1 verified: working magnetostatics foundation"

---

## Phase 2 Planning (Post-Phase 1)

After Phase 1 is verified working:

### ⬜ Chunk 14: Design time-harmonic formulation
**Scope:** Write math documentation for E-field formulation

### ⬜ Chunk 15: Implement complex material properties
**Scope:** Add support for ε = ε' - jε''

### ⬜ Chunk 16: Add PEC boundary conditions
**Scope:** Implement n×E = 0 boundary

... (more chunks as needed)

---

## Testing Protocol (REQUIRED)

**Before marking ANY chunk complete:**

1. **Run exact test command from chunk**
2. **Verify output matches "Expected output"**
3. **If output differs:**
   - Read error carefully
   - Fix code
   - Re-run test
   - Repeat until pass
4. **If stuck for >30 min:**
   - DO NOT mark chunk complete
   - Add note to chunk: "BLOCKED: [reason]"
   - Move to next chunk
   - Send Telegram message flagging issue
5. **Once test passes:**
   - Commit
   - Mark chunk ✅ with commit hash
   - Move to next chunk

**NO EXCEPTIONS.** If a chunk takes multiple sessions, that's fine. Don't fake completion.

---

## Immediate Priority

**Next chunk to work on:** Chunk 7 (Gmsh cube mesh)

**Why:** Establishes working Gmsh → dolfinx pipeline for future mesh work

**Blocked on:** Nothing

**Estimated time:** 15-30 minutes

**Ready to start:** Yes

---

## User Decision Needed

Please reply with decisions on:

1. **Helmholtz coil:** Keep trying to fix, or remove from Phase 1?
2. **Convergence tests:** Simplify and retry, or skip for now?

Once decided, I'll update ROADMAP and proceed with Chunk 7.
