# FC3 Visual Explorer Implementation Notes

## Project Goal

FC3 Visual Explorer is a lightweight browser-based application for introducing, visualizing, and comparing third-order interatomic force constants (FC3) from ShengBTE/thirdorder-style workflows.

The project turns dense `FORCE_CONSTANTS_3RD` blocks into interactive geometry and compact analysis plots. It is intended to help users inspect FC3 block geometry, tensor-component distributions, block-norm distributions, and differences between multiple FC3 files.

## Default Example

The page opens with an embedded example dataset:

- FC3 display name: `FORCE_CONSTANTS_3RD_Example`
- Structure display name: `BASE.si_supper.scf.in_example`
- Browser payload: `app/data/fc3_045.json`

The embedded JSON contains the prebuilt structure and FC3 block data needed for immediate 3D visualization and analysis without requiring the user to upload local files.

## User Inputs

The app supports user-supplied files through browser file inputs:

- A ShengBTE/thirdorder-compatible `FORCE_CONSTANTS_3RD` file.
- A matching Quantum ESPRESSO structure file containing lattice vectors and atomic positions.

If only an FC3 file is uploaded, the app may reuse the currently loaded structure. For data from a different material, primitive cell, or supercell, the matching structure file should be uploaded as well.

## FC3 File Interpretation

Each `FORCE_CONSTANTS_3RD` block is interpreted as one three-body FC3 block with:

- `block_id`: write-order identifier.
- `Rj`, `Rk`: periodic translation vectors for the second and third atoms.
- `i`, `j`, `k`: primitive-cell atom labels.
- 27 Cartesian tensor components.

The real-space triangle used for visualization is built from:

```text
xi = ri
xj = rj + Rj
xk = rk + Rk
```

Here `ri`, `rj`, and `rk` are basis positions from the structure file. `Rj = 0` or `Rk = 0` means the corresponding atom is in the reference primitive-cell image; it does not mean the atom is at the coordinate origin.

## Derived Quantities

For each FC3 block, the app computes:

- Frobenius norm:

```text
FC3_block = ||Phi3_block||F
          = sqrt(sum_{alpha,beta,gamma} Phi_{alpha,beta,gamma}^2)
```

- Triangle edge lengths:

```text
edge_ij = |xj - xi|
edge_ik = |xk - xi|
edge_jk = |xj - xk|
```

- Triangle perimeter:

```text
perimeter = edge_ij + edge_ik + edge_jk
```

- Maximum edge length:

```text
max_edge_length = max(edge_ij, edge_ik, edge_jk)
```

- Tensor-component magnitudes:

```text
abs(Phi3_component)
```

Zero values and values below plotting thresholds are retained in counts and tables where applicable, but they are not drawn directly on logarithmic axes.

## Introduction Page

The Introduction page provides a compact explanation of:

- What FC3 represents physically.
- How two-body and three-body force responses differ.
- Why FC3 blocks can be represented as three-atom triangles.
- What the web app does and which analysis modules it provides.
- Which input files are required and how users should operate the app.

The Concept Illustration section uses static images from `app/assets/` to show two-body and three-body force-response sketches.

## 3D Visualization Page

The 3D page renders a VESTA-like FC3 inspection view:

- Atoms are drawn as shaded spheres.
- Bonds are drawn as neutral lines based on element-pair distance criteria.
- FC3 blocks are drawn as triangle edges without filled faces.
- FC3 block norm controls triangle color and line width.
- A colorbar indicates the norm-color mapping.
- Lattice boundary inputs control the displayed image range.
- View direction controls apply common crystallographic views.
- Zoom controls and mouse wheel change camera distance.
- Toggles control atoms, bonds, FC3 network, and cell boundary visibility.

The FC3 object tree is organized by:

```text
Element -> atom label -> FC3 blocks sorted by norm
```

Users can expand or collapse tree levels, hide elements or atom labels, hide individual blocks, and highlight selected FC3 blocks. Pick-atom mode identifies a concrete visible atom instance and emphasizes FC3 blocks involving that atom within the current norm range.

## FC3 Analysis Page

The Analysis page provides four plots:

1. `FC3_component count distribution`
2. `FC3_block count distribution`
3. `FC3_block vs perimeter`
4. `FC3_block vs max_edge_length`

The component plot uses the 27 Cartesian tensor components from each FC3 block. The app checks whether positive and negative component magnitudes are paired for the current dataset. If they are paired, only the positive-component distribution is plotted; otherwise, positive and negative magnitude distributions are plotted separately.

The block distribution uses the Frobenius norm of each FC3 block. The perimeter and max-edge plots relate block strength to the three-body triangle geometry.

The default plotting threshold is:

```text
1e-4 eV/Angstrom^3
```

Values below this threshold are omitted from log-scale plots but retained in notes and tables. Analysis plots use a white background, black axes, clear ticks, high-resolution canvas rendering, and ordinary numeric tick labels.

## FC3 Compare Page

The Compare page overlays distribution-level descriptors from multiple FC3 files. It does not require the compared files to have identical block counts.

The compare workflow includes:

- One shared structure input.
- Up to four FC3 file slots.
- Editable legend text for each slot.
- Editable colors for each slot.
- Load and unload controls for each slot.
- Optional batch loading into the compare slots.

The Compare page draws:

1. `FC3_component count distribution`
2. `FC3_block count distribution`
3. `FC3_block vs perimeter`
4. `FC3_block vs max_edge_length`

All compare plots use global axis ranges across the loaded files so that distributions are visually comparable.
