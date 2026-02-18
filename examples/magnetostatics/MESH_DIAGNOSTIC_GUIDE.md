# Mesh Diagnostic Guide for Straight Wire Problem

## What the Mesh SHOULD Be

### Geometry Structure

The straight wire problem should have **TWO separate 3D volumetric regions**:

```
┌─────────────────────────────────┐
│                                 │  ← Outer boundary (cylinder surface)
│         AIR DOMAIN              │     r = domain_radius (4 cm)
│         (tag = 2)               │
│                                 │
│    ┌─────┐                      │
│    │WIRE │  ← Inner cylinder    │  ← Wire (current source)
│    │tag=1│     r = wire_radius  │     tag = 1
│    └─────┘     (1.5 mm)         │     Tetrahedral mesh
│                                 │
│                                 │  ← Air/surrounding domain
│                                 │     tag = 2
│                                 │     Tetrahedral mesh
└─────────────────────────────────┘

        ↑                 ↑
     z-axis         Annular region
                    (air domain)
```

### Volume Breakdown

1. **Wire (tag=1)**:
   - Small cylinder along z-axis
   - Radius: wire_radius (e.g., 1.5 mm)
   - Contains tetrahedral mesh elements
   - Current density J applied here
   - Should be ~1-5% of total mesh volume

2. **Air Domain (tag=2)**:
   - Annular region: outer cylinder minus inner wire
   - From wire_radius to domain_radius
   - Contains tetrahedral mesh elements
   - B-field computed here (no current)
   - Should be ~95-99% of total mesh volume

### What You Should See

**In console output:**
```
Cell tags found: [1 2]
  Tag 1 (wire): ~500 cells
  Tag 2 (air/domain): ~30000 cells
```

The air domain should have MANY more cells than the wire!

---

## How to Diagnose Issues

### Step 1: Check Console Output

Run the example and look for:
```
Mesh created: XXXXX cells, YYYYY vertices
Cell tags found: [1 2]
  Tag 1 (wire): ### cells
  Tag 2 (air/domain): ##### cells
```

**Red flags:**
- ❌ Only one tag found (not two)
- ❌ Wire has more cells than air
- ❌ Wire has comparable number of cells to air
- ❌ "WARNING: No cell tags found!"

**Good signs:**
- ✓ Two tags: [1 2]
- ✓ Air domain has 50-100x more cells than wire
- ✓ Total reasonable (~30,000-50,000 cells)

### Step 2: Inspect Gmsh File

The code now saves `straight_wire.msh` for inspection.

**Open in Gmsh:**
```bash
gmsh straight_wire.msh
```

**What to check:**
1. Tools → Options → Mesh → Surface faces: Enable
2. Tools → Options → Mesh → Volume faces: Enable
3. Tools → Visibility:
   - Check "Physical groups"
   - You should see "wire" and "domain"

**In Gmsh 3D view:**
- ✓ Should see two distinct volumes
- ✓ Wire should be thin cylinder in center
- ✓ Domain should surround wire
- ✓ Both should be filled with tetrahedra (not just surfaces)

### Step 3: Check Cell Tags in ParaView

**Load mesh with tags:**
```
Open: paraview_output/straight_wire_mesh.xdmf
Reader: Xdmf3ReaderT
Apply
```

**Visualize tags:**
1. In "Coloring" dropdown, select "cell_tags" or "meshtags"
2. You should see two colors:
   - One color for wire (center, thin)
   - Another color for air domain (surrounding)

**What you should see:**
- ✓ Thin cylindrical region in center (wire, tag=1)
- ✓ Large surrounding region (air, tag=2)
- ✓ Clear boundary between regions
- ✓ Both regions extend full length in z-direction

---

## Common Issues and Fixes

### Issue 1: No Cell Tags Found

**Symptom:** Console shows "WARNING: No cell tags found!"

**Cause:** Mesh generation failed to create physical groups

**Fix:** Check that gmsh fragment operation succeeded:
```python
# In mesh.py, add debugging:
print(f"Volumes after fragment: {volumes}")
print(f"Wire volume ID: {wire_volume}")
print(f"Domain volume ID: {domain_volume}")
```

### Issue 2: Only Surface Mesh (Not Volumetric)

**Symptom:** Mesh looks like hollow cylinders in ParaView

**Cause:** 3D mesh generation failed

**Fix:** Check gmsh output - should say "generate(3)" not "generate(2)"

### Issue 3: Wrong Volume Tagging

**Symptom:** All cells have same tag, or tags reversed

**Cause:** Bounding box detection not working

**Fix:** Improve volume identification in mesh.py:
```python
# Instead of just bounding box, check actual radius
vol_mass = gmsh.model.occ.getMass(vol[0], vol[1])
# Smaller volume = wire
if vol_mass < threshold:
    wire_volume = vol[1]
```

### Issue 4: B-field Only at Ends

**Symptom:** Field concentration at z = ±L/2 (ends of domain)

**Possible causes:**
1. **Boundary conditions:** Natural BC might be causing field exit at ends
2. **Domain too short:** Wire length should be >> domain radius
3. **Current not applied correctly:** Check subdomain_id=1 matches wire tag
4. **Gauge issue:** Coulomb gauge penalty might be too strong

**Fixes:**
- Increase wire_length (make it longer relative to radius)
- Check that current is applied to tag=1 (wire)
- Reduce gauge_penalty from 1e-3 to 1e-6

### Issue 5: Field Doesn't Look Cylindrical

**Symptom:** Field pattern not radially symmetric

**Possible causes:**
1. Mesh quality issues
2. Wire not centered at origin
3. Solver convergence problems

**Fixes:**
- Check mesh optimization worked
- Verify wire center = (0, 0, z)
- Increase mesh resolution
- Check solver residuals

---

## Expected Physical Behavior

### Correct B-field Pattern

For a straight wire along z-axis:

1. **Inside wire (r < wire_radius)**:
   - B increases linearly with r
   - B(r) = (μ₀J/2) × r for uniform J

2. **Outside wire (r > wire_radius)**:
   - B decreases as 1/r
   - B(r) = μ₀I/(2πr)
   - Pure azimuthal direction (θ)

3. **On z-axis (r=0)**:
   - B should be ZERO (no field on axis)

4. **At ends (z = ±L/2)**:
   - Field should gradually go to zero
   - Some "fringing" expected

### What You Should NOT See

- ❌ Strong field only at domain ends
- ❌ Radial field components (should be azimuthal only)
- ❌ Field along wire axis (should be perpendicular)
- ❌ Discontinuities in field
- ❌ Spiral field patterns

---

## Verification Checklist

Run through these checks:

**Mesh Structure:**
- [ ] Console shows two tags: [1 2]
- [ ] Air domain has 50-100x more cells than wire
- [ ] Gmsh file shows two distinct volumes
- [ ] ParaView mesh shows tagged regions

**Mesh Quality:**
- [ ] All elements are tetrahedra (3D)
- [ ] Both regions have volumetric mesh
- [ ] Wire is thin cylinder in center
- [ ] Domain surrounds wire completely

**Physics:**
- [ ] Current applied to tag=1 (wire) only
- [ ] B-field computed throughout domain
- [ ] Field is azimuthal (circular around wire)
- [ ] Field decreases with distance (1/r)

**Visualization:**
- [ ] Field lines form circles in ParaView
- [ ] Magnitude decreases from wire outward
- [ ] Cylindrical symmetry visible
- [ ] No field concentration at ends only

If all checks pass ✓, the mesh is correct!

---

## Debugging Script

Quick Python script to check mesh properties:

```python
import numpy as np
from mpi4py import MPI
from fem_em_solver.io.mesh import MeshGenerator

comm = MPI.COMM_WORLD

mesh, cell_tags, facet_tags = MeshGenerator.straight_wire_domain(
    wire_length=0.2,
    wire_radius=0.0015,
    domain_radius=0.04,
    resolution=0.006,
    comm=comm
)

print(f"Total cells: {mesh.topology.index_map(3).size_global}")
print(f"Total vertices: {mesh.topology.index_map(0).size_global}")

if cell_tags is not None:
    tags = cell_tags.values
    for tag in np.unique(tags):
        count = np.sum(tags == tag)
        percent = 100 * count / len(tags)
        print(f"Tag {tag}: {count} cells ({percent:.1f}%)")
else:
    print("ERROR: No cell tags!")
```

Expected output:
```
Total cells: ~35000
Total vertices: ~7000
Tag 1: ~200 cells (0.6%)
Tag 2: ~34800 cells (99.4%)
```

---

## If Problem Persists

If you still see issues after these checks:

1. **Share diagnostics:**
   - Console output (cell tag counts)
   - Screenshot from Gmsh
   - Screenshot from ParaView

2. **Check solver setup:**
   - Verify current_density applied to subdomain_id=1
   - Check gauge_penalty value
   - Review boundary conditions

3. **Try simpler problem:**
   - Test with coarser mesh first
   - Verify with single-domain problem
   - Compare with analytical solution quantitatively
