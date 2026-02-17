# ROADMAP.md — FEM EM Solver Development Roadmap

Each chunk is a single focused task that can be completed and tested in one session.
**CRITICAL:** Each chunk MUST be tested before marking complete. If it doesn't run, fix it.

---

## Phase 1: Magnetostatics Foundation

### ✅ Chunk 0: Initial solver structure
**Status:** Complete and tested in Docker
**Commit:** `9086ecd`

### ✅ Chunk 1: Fix current density restriction to wire subdomain  
**Status:** Complete
**Commit:** `6da3f33`

### ✅ Chunk 2: Circular loop mesh and validation
**Status:** Complete
**Commit:** `e5a0936`

### ⬜ Chunk 3: Helmholtz coil validation
**Status:** CODE WRITTEN, NOT WORKING
**Commit:** `030eb48`
**⚠️ ISSUE FOUND:** Mesh generation hangs/times out
**Problem:** Two tori fragmentation is too complex for Gmsh

---

## Current Priority: Fix Broken Chunks

### ⬜ Chunk 4: FIX Helmholtz coil mesh generation
**Scope:** Make Helmholtz coil actually run without hanging

**Problem Analysis:**
- Current implementation creates two tori and fragments them with domain
- Gmsh mesh.optimize("Netgen") hangs on complex topology
- Need simpler approach

**Fix Options:**
1. Remove mesh.optimize() call (may produce lower quality mesh)
2. Use boolean union for wires first, then fragment with domain
3. Simplify geometry (larger wire radius, coarser mesh)
4. Use different meshing algorithm

**Steps:**
1. Try Option 1 first (simplest):
   - Edit `src/fem_em_solver/io/mesh.py`
   - Comment out `gmsh.model.mesh.optimize("Netgen")` in helmholtz_coil_domain()
   
2. Test in Docker:
   ```bash
   cd ~/Development/fem-em-solver/docker
   docker compose exec fem-em-solver bash
   cd /workspace
   export PYTHONPATH=/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:/workspace/src
   timeout 60 python3 examples/magnetostatics/03_helmholtz_coil.py
   ```

3. If that doesn't work, try Option 2:
   - Create both tori
   - Use `gmsh.model.occ.fuse()` to merge them into one object
   - Then fragment with domain
   
4. Repeat test until it completes in <60 seconds

5. Once working, commit fix

**Success criteria:** Helmholtz example runs to completion in under 60 seconds
**Commit message:** "Fix Helmholtz coil mesh generation - prevent hang"

---

## Remaining Phase 1 Work (After Fix)

### ⬜ Chunk 5: VERIFY Helmholtz coil results
**Scope:** Run and check output is physically reasonable

**Steps:**
1. Run example and capture output
2. Check B_z at center matches analytical ~0.7155 * μ₀I/R
3. Check field uniformity in central region
4. If results look wrong, debug and fix

**Success criteria:** Numerical results match analytical within 10%

---

### ⬜ Chunk 6: Add convergence study test file (skeleton only)
**Scope:** Create empty test file with proper structure

**Steps:**
1. Create `tests/validation/test_convergence.py`
2. Add imports (dolfinx, fem_em_solver modules)
3. Create test class `TestConvergence` with empty methods:
   - `test_h_refinement_straight_wire()` 
   - `test_p_refinement_straight_wire()`
4. Add docstrings explaining what each test will do
5. Run `python3 -c "import tests.validation.test_convergence"` to verify no import errors
6. Commit

**Success criteria:** File exists, imports work, pytest can discover tests
**Commit message:** "Add convergence study test skeleton"

---

### ⬜ Chunk 7: Implement h-refinement convergence test
**Scope:** Test error vs mesh size for straight wire problem

**Steps:**
1. In `test_convergence.py`, implement `test_h_refinement_straight_wire()`:
   - Loop over resolutions: [0.02, 0.01, 0.007, 0.005]
   - For each resolution:
     a. Generate straight wire mesh
     b. Solve magnetostatic problem
     c. Compute B-field at sample points
     d. Calculate L2 error vs analytical
   - Store (h, error) pairs
   - Fit log(error) vs log(h) to get convergence rate
   
2. Assert rate > 0.8 (expected for linear elements)

3. Run test in Docker:
   ```bash
   pytest tests/validation/test_convergence.py::TestConvergence::test_h_refinement_straight_wire -v
   ```

4. Fix any errors until test passes

**Success criteria:** Test runs and reports convergence rate > 0.8
**Commit message:** "Implement h-refinement convergence test"

---

### ⬜ Chunk 8: Implement p-refinement convergence test  
**Scope:** Test error vs polynomial degree

**Steps:**
1. Implement `test_p_refinement_straight_wire()`:
   - Fixed mesh resolution (e.g., 0.01)
   - Loop over degrees: [1, 2, 3]
   - For each degree:
     a. Create solver with degree=N
     b. Solve and compute error
   - Expect error to decrease with higher degree
   
2. Assert that degree 2 error < degree 1 error

3. Run and verify in Docker

**Success criteria:** Higher degree gives lower error
**Commit message:** "Implement p-refinement convergence test"

---

### ⬜ Chunk 9: Document Phase 1 results
**Scope:** Write validation report

**Steps:**
1. Create `docs/validation/magnetostatics.md`
2. Document:
   - Problem formulations (weak form)
   - Mesh details for each geometry
   - Convergence study results table
   - Comparison with analytical solutions
3. Include example command outputs
4. Link to test files

**Success criteria:** Document is complete and readable
**Commit message:** "Add Phase 1 validation documentation"

---

### ⬜ Chunk 10: Final Phase 1 cleanup
**Scope:** Ensure everything works end-to-end

**Steps:**
1. Run all tests in Docker:
   ```bash
   pytest tests/ -v --tb=short
   ```
   
2. Verify all pass or document known failures

3. Run all examples:
   ```bash
   python3 examples/magnetostatics/01_straight_wire.py
   python3 examples/magnetostatics/02_circular_loop.py
   python3 examples/magnetostatics/03_helmholtz_coil.py
   ```

4. Check git status - everything committed?

5. Update PROJECT_PLAN.md to mark Phase 1 complete

**Success criteria:** All tests run, examples execute, docs complete
**Commit message:** "Phase 1 complete: magnetostatics validated"

---

## Testing Protocol (STRICT)

**Before marking ANY chunk complete:**

1. **Run the code in Docker:**
   ```bash
   cd ~/Development/fem-em-solver/docker
   docker compose exec fem-em-solver bash
   cd /workspace
   export PYTHONPATH=/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:/workspace/src
   python3 <test_or_example>
   ```

2. **Verify it produces expected output:**
   - No Python exceptions
   - No hangs/timeouts
   - Reasonable numerical values
   - Tests pass (if applicable)

3. **If it fails:**
   - Read error messages carefully
   - Fix the issue
   - Re-run
   - Repeat until it works

4. **Only then:** Commit and mark chunk complete

**DO NOT** mark chunks complete based on code review alone. They must actually run.

---

## Immediate Next Actions

1. **Chunk 4:** Fix Helmholtz coil mesh (it's broken)
2. **Chunk 5:** Verify Helmholtz results (after fix)
3. **Chunks 6-10:** Continue with convergence and documentation

**DO NOT** start Phase 2 until all Phase 1 chunks are verified working.
