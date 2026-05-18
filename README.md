# FC3 Visual Explorer

FC3 Visual Explorer is a browser-based tool for visualizing, analyzing, and comparing third-order interatomic force constants from ShengBTE/thirdorder workflows.

## Online App

Use the web app here:

https://zxy-1104.github.io/FC3-Visual-Explorer/

## What It Does

FC3 Visual Explorer includes three main analysis modules:

- **3D Visualization**  
  Displays atoms, bonds, unit-cell boundaries, and FC3 three-body interaction triangles. FC3 block norms are mapped to color and line width, with controls for norm filtering, view direction, lattice boundary selection, atom picking, highlighting, and hiding FC3 objects.

- **FC3 Analysis**  
  Provides plots for FC3 tensor-component distributions, FC3 block-norm distributions, and the relationship between FC3 block strength and three-body geometry, including triangle perimeter and maximum edge length.

- **FC3 Compare**  
  Compares multiple FC3 files using consistent distribution and geometry-based plots, which is useful for checking differences caused by cutoff, supercell size, displacement amplitude, or other calculation settings.

## Input Files

To analyze your own data, provide:

- A ShengBTE/thirdorder-compatible `FORCE_CONSTANTS_3RD` file.
- A matching Quantum ESPRESSO structure file containing lattice vectors and atomic positions.

The structure file should include information such as:

```text
CELL_PARAMETERS
ATOMIC_POSITIONS
```

## Example Data

The web app includes an embedded default example that loads automatically when the page opens.

You can also upload your own `FORCE_CONSTANTS_3RD` and QE structure files to generate new 3D visualizations, analysis plots, and FC3 comparisons.

## Contact

For questions or suggestions, please contact:

zhaoxingyu@hnu.edu.cn
