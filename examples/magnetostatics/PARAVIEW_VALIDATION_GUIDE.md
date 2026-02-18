# ParaView Validation Guide for Magnetostatic Examples

## Straight Wire: What You Should See

### Physical Expectations (Analytical Solution)

For a straight wire carrying current **I** along the z-axis, the magnetic field should have:

1. **Direction**: Purely azimuthal (circumferential) around the wire
   - Field circles around the wire
   - No component along the wire (B_z ≈ 0)
   - No radial component (B_r ≈ 0)
   - Uses right-hand rule: thumb = current direction, fingers = B direction

2. **Magnitude**: Decreases as **1/r**
   - B(r) = μ₀I/(2πr)
   - Strongest near wire surface
   - Decays with distance
   - At r = 1 cm with I = 1A: B ≈ 20 μT

3. **Symmetry**:
   - Cylindrical symmetry about wire axis
   - Same at all z-positions (for infinite wire approximation)
   - Same at all azimuthal angles θ

---

## How to Verify in ParaView

### Step 1: Load the Data

**RECOMMENDED:** Use the combined file which has cell tags and B-field on the same grid:
```
File → Open → paraview_output/straight_wire_combined.xdmf
Reader: Xdmf3ReaderT
Click "Apply"
```

**Alternative:** Use individual files:
```
File → Open → paraview_output/straight_wire_B.xdmf
Reader: Xdmf3ReaderT
Click "Apply"
```
(Note: Individual files require extra steps to access cell tags)

### Step 2: Check Field Magnitude Distribution

**What to do:**
1. In "Coloring" dropdown, select **"B"** (or the magnitude array)
2. Click "Rescale to Data Range" button
3. Change representation to "Surface with Edges" to see mesh

**What to look for:**
- ✓ **Strongest field (red)** near the wire surface
- ✓ **Weakest field (blue)** far from wire
- ✓ **Smooth gradient** from wire to boundary
- ✓ **Cylindrical pattern** - looks the same at all angles around wire
- ✗ **Avoid**: Irregular patches, discontinuities, noise

**Expected values** (for I=1A):
- Near wire (r ≈ 2mm): ~100-200 μT
- Mid-domain (r ≈ 2cm): ~10 μT
- Far field (r ≈ 4cm): ~5 μT

---

### Step 3: Visualize Field Lines (Stream Tracer)

**What to do:**
1. Filters → Common → **Stream Tracer**
2. Set **Vectors** = "B"
3. **Seed Type** = "Line Source"
   - Point1: [0.005, 0, 0.05]  (start near wire)
   - Point2: [0.005, 0, -0.05] (line along wire)
   - Resolution: 20
4. Click "Apply"
5. Optional: Add **Tube** filter to make lines thicker

**What to look for:**
- ✓ Field lines form **circles** around the wire
- ✓ Circles lie in planes **perpendicular to wire** (xy-planes)
- ✓ Lines are **smooth and continuous**
- ✓ **Denser spacing** near wire (stronger field)
- ✓ All circles **centered on wire axis**
- ✗ **Avoid**: Lines along wire, spirals, lines going radially

---

### Step 4: Check Field Vectors (Glyph Filter)

**What to do:**

**IMPORTANT: First filter out the wire cells!**
1. Before applying Glyph, apply **Threshold** filter:
   - Filters → Common → Threshold
   - Scalars: "meshtags" or "cell_tags"
   - Minimum: 2, Maximum: 2 (keeps only air domain)
   - Click Apply

2. Now apply **Glyph** to the thresholded data:
   - Filters → Common → Glyph
   - **Glyph Type**: Arrow
   - **Orientation Array**: B
   - **Scale Array**: B
   - **Scale Factor**: 0.0002 (adjust to see arrows)
   - **Glyph Mode**: Every Nth Point (N=5)
   - Click "Apply"

**Why?** Without threshold, glyphs in the tiny wire cells dominate and clutter the view!

**What to look for:**
- ✓ Arrows **tangent to circles** around wire
- ✓ Arrows point in **azimuthal direction** (θ-direction)
- ✓ **Longer arrows** near wire (stronger field)
- ✓ **Shorter arrows** far from wire
- ✓ All arrows **perpendicular to radial direction**
- ✗ **Avoid**: Arrows pointing toward/away from wire, arrows along wire axis

---

### Step 5: Cross-Section View (Slice Filter)

**What to do:**
1. Delete Glyph filter
2. Filters → Common → **Slice**
3. **Slice Type**: Plane
4. **Origin**: [0, 0, 0]
5. **Normal**: [0, 0, 1] (for xy-plane) or [1, 0, 0] (for yz-plane)
6. Click "Apply"
7. Color by "B" magnitude

**What to look for (xy-plane slice):**
- ✓ **Concentric circles** of color around wire
- ✓ Wire appears as small circle in center
- ✓ **Radially symmetric** pattern
- ✓ **Smooth gradient** from center to edge

**What to look for (yz-plane slice through wire):**
- ✓ **Vertical symmetry** (same above and below any z-level)
- ✓ Field **strongest at wire edges** (r = wire_radius)
- ✓ **Horizontal bands** of constant color (constant z)
- ✓ Pattern looks like two vertical "hot zones" at wire edges

---

### Step 6: Quantitative Check (Plot Over Line)

**What to do:**
1. Remove filters, back to original data
2. Filters → Data Analysis → **Plot Over Line**
3. **Line Parameters**:
   - Point1: [0.002, 0, 0] (near wire)
   - Point2: [0.04, 0, 0]  (to boundary)
4. **X Axis Parameters**: Use "arc_length" (this is radial distance)
5. Click "Apply"
6. In the line chart, select "B" arrays to plot

**What to look for:**
- ✓ B_X and B_Y should be **comparable magnitude** (not both zero)
- ✓ B_Z should be **nearly zero** (< 1% of total |B|)
- ✓ |B| should **decrease** with distance (1/r trend)
- ✓ Plot log(|B|) vs log(r) → should be **straight line** with slope ≈ -1

**Quick calculation check:**
At r = 0.01 m (1 cm), with I = 1 A:
- B_analytical = (4π×10⁻⁷ × 1) / (2π × 0.01) ≈ 20 μT
- B_numerical should be within 5-15% of this

---

## Common Issues and Fixes

### Issue: Field lines spiral along wire
**Problem**: Current density might be in wrong direction
**Fix**: Check that J = [0, 0, J_magnitude] (along z)

### Issue: Field too weak everywhere
**Problem**: Current density too low or mesh too coarse
**Fix**: Check wire_area calculation, increase mesh resolution

### Issue: Irregular field pattern
**Problem**: Mesh quality issues or solver didn't converge
**Fix**: Refine mesh, check solver output for warnings

### Issue: Field has z-component on axis
**Problem**: Geometry not perfectly symmetric or boundary conditions
**Fix**: Check that wire is centered at origin, domain is long enough

### Issue: Field doesn't decay properly
**Problem**: Domain too small (boundary effects)
**Fix**: Increase domain_radius

---

## Quick Validation Checklist

Run through this checklist in ParaView:

- [ ] Load data successfully with Xdmf3ReaderT
- [ ] Color by B magnitude shows cylindrical pattern
- [ ] Stream Tracer shows circular field lines
- [ ] Field lines perpendicular to wire
- [ ] Glyph arrows tangent to circles
- [ ] Slice (xy) shows concentric circles
- [ ] Plot Over Line shows 1/r decay
- [ ] B_z component is nearly zero on axis
- [ ] Field magnitude matches analytical (~20 μT at 1cm)
- [ ] No discontinuities or artifacts

If all checks pass ✓, your solution is **physically accurate**!

---

## Advanced: Comparing to Analytical Solution

### Export numerical data from ParaView:
1. Use Plot Over Line as above
2. File → Save Data → save as CSV
3. Import into Python/MATLAB

### Python comparison script:
```python
import numpy as np
import matplotlib.pyplot as plt

# Load ParaView data
data = np.loadtxt('lineplot.csv', delimiter=',', skiprows=1)
r = data[:, 0]  # radial distance
B_num = data[:, 1]  # B magnitude

# Analytical solution
I = 1.0
mu_0 = 4 * np.pi * 1e-7
B_ana = mu_0 * I / (2 * np.pi * r)

# Plot comparison
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.loglog(r*100, B_num*1e6, 'o', label='FEM')
plt.loglog(r*100, B_ana*1e6, '-', label='Analytical')
plt.xlabel('Distance from wire [cm]')
plt.ylabel('|B| [μT]')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
rel_error = np.abs(B_num - B_ana) / B_ana * 100
plt.semilogx(r*100, rel_error, 'o-')
plt.xlabel('Distance from wire [cm]')
plt.ylabel('Relative Error [%]')
plt.grid(True)
plt.tight_layout()
plt.show()
```

Expected error: < 10% for r > 2×wire_radius

---

## For Other Examples

### Circular Loop (02_circular_loop.py):
- **Expected**: Dipole field pattern
- **On axis**: Field maximum at center, decreases with |z|
- **Off axis**: Field curves through loop
- **Validation**: Compare B_z on axis to analytical solution

### Helmholtz Coil (03_helmholtz_coil.py):
- **Expected**: Nearly uniform field in center region
- **Key feature**: Flat B_z between the two coils
- **Validation**: Plot B_z along axis - should be plateau in middle
- **Quantify**: Coefficient of variation in central 1cm < 1%

---

## References

- Jackson, "Classical Electrodynamics", Chapter 5
- Griffiths, "Introduction to Electrodynamics", Chapter 5
- Biot-Savart Law
- Ampère's Circuital Law
