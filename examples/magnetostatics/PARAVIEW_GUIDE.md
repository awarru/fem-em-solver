# ParaView Visualization Guide

## Running the Example

```bash
cd /home/taz5297/Development/fem-em-solver
python examples/magnetostatics/01_straight_wire.py
```

This will create a `paraview_output/` directory with visualization files.

---

## Output Files

### XDMF Format (Traditional)
- `straight_wire_mesh.xdmf` + `.h5` - Mesh only
- `straight_wire_A.xdmf` + `.h5` - Vector potential A field
- `straight_wire_B.xdmf` + `.h5` - Magnetic field B

### VTX Format (Modern, requires ADIOS2)
- `straight_wire_A.bp/` - Vector potential A field
- `straight_wire_B.bp/` - Magnetic field B

---

## Opening in ParaView

### Method 1: XDMF Files (Recommended for beginners)

1. **Open ParaView**
2. **File → Open**
3. Navigate to `paraview_output/`
4. Select `straight_wire_B.xdmf` (the magnetic field)
5. In the dialog, choose **"Xdmf3ReaderT"** as the reader
6. Click **"Apply"** in the Properties panel on the left

You should now see the mesh loaded!

### Method 2: VTX Files (Modern)

1. **Open ParaView**
2. **File → Open**
3. Navigate to `paraview_output/`
4. Select the `straight_wire_B.bp` **directory**
5. ParaView should auto-detect the format
6. Click **"Apply"**

---

## Visualizing the Magnetic Field

### View Field as Colored Surface

1. After loading the data, in the **Properties** panel:
   - Under "Coloring", select **"B"** from the dropdown
   - This colors the mesh by magnetic field magnitude

2. Click the **"Rescale to Data Range"** button (⟳ icon) to adjust colors

### View Field Vectors (Glyphs)

1. With your data loaded, go to **Filters → Common → Glyph**
2. In Glyph properties:
   - **Glyph Type**: Arrow
   - **Orientation Array**: B
   - **Scale Array**: B
   - **Scale Factor**: Adjust to make arrows visible (try 0.001 to start)
   - **Glyph Mode**: All Points or Every Nth Point
3. Click **"Apply"**

This shows arrows indicating field direction and magnitude!

### View Field Lines (Stream Tracer)

1. With your data loaded, go to **Filters → Common → Stream Tracer**
2. In Stream Tracer properties:
   - **Vectors**: B
   - **Seed Type**: Point Cloud (or Line Source)
   - Adjust seed positions to interesting regions
3. Click **"Apply"**
4. Optionally add **Tube** filter to make lines thicker

This shows magnetic field lines flowing around the wire!

### View Cross-Section (Slice)

1. With your data loaded, go to **Filters → Common → Slice**
2. In Slice properties:
   - **Slice Type**: Plane
   - **Origin**: Center of your domain
   - **Normal**: [0, 0, 1] for XY plane, [0, 1, 0] for XZ plane, etc.
3. Click **"Apply"**
4. Color by **"B"** magnitude

This shows the field in a 2D slice through the domain!

---

## Viewing the Mesh

If you just want to see the mesh structure:

1. Open `straight_wire_mesh.xdmf`
2. In the toolbar, change representation from "Surface" to **"Wireframe"** or **"Surface with Edges"**

---

## Troubleshooting

### "Cannot find reader for XDMF file"
- Make sure you select **"Xdmf3ReaderT"** (note the "3" and "T")
- Not "XdmfReader" (old version)

### "Vectors not showing"
- Check that **Orientation Array** and **Scale Array** are both set to "B"
- Increase the **Scale Factor** if arrows are too small
- Try "Every Nth Point" with N=10 if there are too many glyphs

### "Field looks wrong"
- Click **"Rescale to Data Range"** button
- Check units - B field is in Tesla (very small values expected)
- For the straight wire, expect cylindrical symmetry

### "VTX files not created"
- This is normal if ADIOS2 is not installed
- Use XDMF files instead - they work just as well

### "ParaView crashes when opening file"
- Try opening only the mesh first: `straight_wire_mesh.xdmf`
- If that works, then open the field files
- Check ParaView version - need 5.10+ for best XDMF support

---

## Expected Physics for Straight Wire

For a straight wire carrying current I along the z-axis:

- **B-field pattern**: Circular/azimuthal around the wire
- **B-field magnitude**: Decreases as 1/r (where r = distance from wire)
- **Analytical**: B = μ₀I/(2πr) in the azimuthal direction
- **Symmetry**: Cylindrical symmetry about the wire axis

### What you should see:
- Field lines forming circles around the wire
- Strongest field near the wire surface
- Field decreasing with distance
- No field component along the wire axis (z-direction)

---

## Advanced Visualizations

### 1. Magnitude Contours
- **Filters → Common → Contour**
- **Contour By**: B_Magnitude (or magnitude of B vector)
- Creates iso-surfaces of constant field strength

### 2. Vector Field Animation
- Use **Calculator** filter to compute derived quantities
- Example: Compute |B| = sqrt(B_X^2 + B_Y^2 + B_Z^2)

### 3. Compare with Analytical Solution
- Export points from ParaView
- Compare in Python/MATLAB with analytical solution
- The example already does this comparison!

---

## File Format Details

### XDMF Format
- **Pros**:
  - Widely supported, works on all platforms
  - Human-readable XML descriptor
  - Efficient HDF5 binary data storage
- **Cons**:
  - Limited to 2nd order geometry
  - Two files per dataset (.xdmf + .h5)

### VTX Format
- **Pros**:
  - Modern ADIOS2 backend
  - Supports arbitrary order elements
  - Better parallel I/O performance
  - Single directory per dataset
- **Cons**:
  - Requires ADIOS2 installation
  - Newer format, less widespread

**Recommendation**: Use XDMF for maximum compatibility, VTX for performance.

---

## Additional Resources

**ParaView Documentation**:
- [Official ParaView Guide](https://docs.paraview.org/en/latest/)
- [ParaView Tutorial](https://www.paraview.org/Wiki/The_ParaView_Tutorial)

**FEniCSx + ParaView**:
- [FEniCSx Tutorial - ParaView Chapter](https://jsdokken.com/dolfinx-tutorial/chapter1/membrane_paraview.html)
- [DOLFINx I/O Documentation](https://docs.fenicsproject.org/dolfinx/main/python/generated/dolfinx.io.html)

**Understanding the Physics**:
- Magnetic field of current-carrying wires
- Biot-Savart law
- Magnetostatics and Maxwell's equations

---

## Questions?

If you encounter issues or want to visualize other quantities:
1. Check the ParaView built-in help (Help → Guide)
2. Try different filters and experiment!
3. Consult the FEniCSx documentation for other output options

Happy visualizing! 🎨🔬
