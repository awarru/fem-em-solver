# Combined XDMF Output for ParaView

## Problem Solved

When DOLFINx writes both meshtags (cell data) and functions (point data) to the same XDMF file, they end up on **separate grids**. This means ParaView can't use both simultaneously - you can't apply a Threshold filter based on cell tags to data that contains the B-field.

## Solution

We've created a custom XDMF writer (`fem_em_solver.io.paraview_utils.write_xdmf_with_tags`) that properly combines cell tags and field data on a **single unified grid**.

## What's New

### File Structure

After running the examples, you'll now get:

**Individual files (standard DOLFINx output):**
- `straight_wire_mesh.xdmf` - mesh with cell tags
- `straight_wire_A.xdmf` - vector potential A
- `straight_wire_B.xdmf` - magnetic field B

**NEW - Combined file (recommended):**
- `straight_wire_combined.xdmf` - **mesh + cell tags + A + B on one grid** ✨
- `straight_wire_combined.h5` - associated HDF5 data

## Usage in ParaView

### Open the Combined File

1. **File → Open** → `paraview_output/straight_wire_combined.xdmf`
2. **Reader**: Select `Xdmf3ReaderT`
3. **Properties panel**: Check both:
   - ☑ Cell Arrays: `Cell tags`
   - ☑ Point Arrays: `A`, `B`
4. Click **Apply**

### Now Cell Tags Work Everywhere!

**Apply Threshold Filter:**
1. **Filters → Common → Threshold**
2. **Scalars**: Select `Cell tags` (now it shows up!)
3. **Minimum**: 2, **Maximum**: 2
4. Click **Apply**

This removes all wire cells (tag=1), leaving only the air domain (tag=2).

**Apply Glyph Filter** (to the thresholded data):
1. **Filters → Common → Glyph**
2. **Orientation Array**: B
3. **Scale Array**: B
4. **Scale Factor**: 0.0002 (adjust as needed)
5. **Glyph Mode**: Every Nth Point (N=5)
6. Click **Apply**

Result: Clean B-field vectors without wire cell clutter!

### Color by Cell Tags

You can also directly color the mesh by cell tags:
- In the **Coloring** dropdown, select `Cell tags`
- Wire cells (value=1) show in one color
- Air domain (value=2) shows in another

## Technical Details

### What's in the Combined File?

The XDMF file structure:
```xml
<Grid Name="mesh" GridType="Uniform">
  <Topology ...>       <!-- Mesh elements -->
  <Geometry ...>       <!-- Node coordinates -->

  <!-- Cell data (one value per cell) -->
  <Attribute Name="Cell tags" Center="Cell">
    ...
  </Attribute>

  <!-- Point data (one value per node) -->
  <Attribute Name="B" Center="Node">
    ...
  </Attribute>

  <Attribute Name="A" Center="Node">
    ...
  </Attribute>
</Grid>
```

All attributes are on the **same grid**, making them simultaneously accessible.

### Comparison: Individual vs Combined

| Feature | Individual Files | Combined File |
|---------|------------------|---------------|
| Cell tags visible in Threshold | ❌ Separate grid | ✅ Same grid |
| Field data available | ✅ Yes | ✅ Yes |
| Number of files | 3-4 files | 1 file |
| Ease of use | Medium | **Easy** |
| Recommended | For separate visualization | **For filtering** |

## Implementation

The custom writer is in:
```python
from fem_em_solver.io.paraview_utils import write_xdmf_with_tags

write_xdmf_with_tags(
    "output_name",
    mesh,
    cell_tags,
    {"B": B_lag, "A": A_lag},  # Dictionary of functions
    comm=comm
)
```

It manually constructs the XDMF XML to ensure proper grid structure.

## Applied to All Examples

This feature is now integrated into:
- ✅ `01_straight_wire.py`
- TODO: `02_circular_loop.py`
- TODO: `03_helmholtz_coil.py`

## Validation

After running the straight wire example:

1. Open `straight_wire_combined.xdmf` in ParaView
2. Check Properties panel - you should see:
   - Cell Arrays: `Cell tags`
   - Point Arrays: `A`, `B`
3. Apply Threshold → Scalars dropdown should show `Cell tags`
4. Success! ✨

If `Cell tags` doesn't appear in the Threshold filter, you're using the wrong file.
