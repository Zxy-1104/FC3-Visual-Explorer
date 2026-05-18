# FC3 Visual Explorer SPEC

## Purpose

Build and maintain a lightweight browser-based tool for introducing and visualizing third-order force constants (FC3) from ShengBTE/thirdorder workflows. The tool supports FC3 physical explanation, 3D three-body network inspection, and 1D/2D block-norm analysis for uncertainty-quantification studies of lattice thermal conductivity.

## Default Inputs

- `FORCE_CONSTANTS_3RD_0.45`: default third-order force-constant file.
- `BASE.si_supper.scf.in`: default Quantum ESPRESSO Si supercell structure.
- `app/data/fc3_045.json`: prebuilt default browser payload generated from the two files above.

## User-Supplied Inputs

The browser UI should allow users to upload:

- A ShengBTE/thirdorder `FORCE_CONSTANTS_3RD` file.
- Optionally, a Quantum ESPRESSO structure file containing `ATOMIC_POSITIONS crystal` and `CELL_PARAMETERS angstrom`.

If the user uploads only an FC3 file, the app may reuse the currently loaded structure. If the uploaded FC3 belongs to a different structure, the user should upload the matching structure file as well.

## Pages

### FC3 Introduction

Summarize the project context from `AGENTS.md` in a web-readable form:

- FC3 definition and physical meaning.
- Finite-difference reconstruction logic.
- `FORCE_CONSTANTS_3RD` block geometry and the meaning of `Rj/Rk`.
- Relationship between FC3, three-phonon scattering, phonon lifetime, and thermal conductivity.
- Clear warning that FC3 block norm is a consistency descriptor, not a direct predictor of `kappa`.

### 3D Visualization

Display a VESTA-like FC3 three-body network:

- Si atoms as shaded spheres.
- Si bonds as neutral gray lines.
- FC3 three-body interactions as triangle edges without filled faces.
- Edge color and width mapped to block Frobenius norm.
- Visible colorbar with colors distinct from atom bonds.
- Custom lattice-index boundary.
- Preset boundary buttons for common `1x1x1`, `2x2x2`, and `3x3x3` views.
- View-direction presets: free, `[100]`, `[010]`, `[001]`, `[110]`, `[111]`.
- Display toggles for atoms, bonds, FC3 network, and cell boundary.
- Mouse drag rotation.
- Mouse-wheel zoom and explicit zoom/reset controls.
- Object tree: element -> atom label -> FC3 block list sorted by norm.
- Clicking a block in the tree highlights the corresponding FC3 triangle.

### FC3 Analysis

Provide publication-style FC3 block-norm analysis:

- `FC3_component count distribution` using all 27 tensor components from every FC3 block. Check the positive/negative magnitude symmetry for the current dataset; plot only the positive-component distribution when the check passes, and retain zero counts in the note and table.
- `FC3_block count distribution` using the Frobenius norm of each FC3 block.
- `FC3_block vs perimeter`, where `FC3_block = ||Phi3_block||F` and `perimeter = |xj - xi| + |xk - xi| + |xj - xk|`. Plot perimeter on a linear axis and `FC3_block` on a log axis, retaining zero-perimeter blocks at `x = 0`.
- `FC3_block vs max_edge_length`, where `max_edge_length` is the longest side of the FC3 three-body triangle.
- Use `1e-4 eV/Å^3` as the default plotting threshold for component magnitudes and block norms. Values below the threshold are retained in notes/tables but omitted from log-scale plots.
- Provide scrollable data tables for the component distribution and the block-norm/perimeter grouping.
- White background, black axes, clear ticks, restrained color, high-resolution canvas rendering.
- Log axes should use decade major ticks and `2..9 x 10^n` minor ticks with ordinary numeric tick labels; do not add an extra `log10 scale` text label.

### FC3 Compare

Compare multiple `FORCE_CONSTANTS_3RD` files with a shared structure. This page focuses on distribution comparison, not block-by-block parity, so files with different block counts are allowed.

- Provide four FC3 compare slots. Each slot has a file input, load/unload controls, editable legend text, and an editable color. Slot 1 follows the current main dataset by default.
- Provide one shared Compare structure input; if unset, reuse the current main structure.
- Plot multi-file `FC3_component count distribution` curves using the same thresholded component-distribution mode as the 2D page. If every loaded file has paired positive/negative component magnitudes, show one positive-component comparison plot. If any loaded file is not paired, show separate positive and negative magnitude comparison plots.
- Plot multi-file `FC3_block count distribution`, `FC3_block vs perimeter`, and `FC3_block vs max_edge_length` using the same thresholded modes as the Analysis page.
- Use distinct colors and editable legends.
- Report block counts and tensor-component counts for each loaded slot.
- Use global axis ranges across all loaded FC3 files, not per-file ranges.

## Derived Quantities

For each FC3 block:

- Frobenius norm:

```text
||Phi3_block||F = sqrt(sum_{alpha,beta,gamma} Phi_{alpha,beta,gamma}^2)
```

- Real-space vertices:

```text
xi = ri
xj = rj + Rj
xk = rk + Rk
```

- Maximum pair distance:

```text
dmax = max(|xj - xi|, |xk - xi|, |xj - xk|)
```

- Triangle perimeter:

```text
perimeter = |xj - xi| + |xk - xi| + |xj - xk|
```

- Tensor-component distribution values:

```text
abs(Phi3_component)
```

Use only nonzero component values for log-scale histograms and report the zero count separately.

## Acceptance Criteria

- The app opens with the default Si 0.45 nm dataset already visible.
- Uploading a valid `FORCE_CONSTANTS_3RD` file rebuilds the 3D, 2D, and 1D views in the browser.
- Uploading a QE structure file updates the structure used for geometry reconstruction.
- Custom boundary controls update atoms, bonds, FC3 blocks, and cell boundary.
- 3D view supports rotation, wheel zoom, zoom buttons, and reset view.
- The FC3 object tree lists only atoms participating in FC3 blocks and sorts each atom's block list by norm.
- Selecting an object-tree block highlights the corresponding triangle clearly.
- 2D and 1D charts are crisp enough for screenshots and use a PRB-like visual style.
- Controls do not overlap the visualization on desktop or narrow windows.

## Prompt, AGENTS.md, And Skill Guidance

- Put task-specific goals, files, UI features, and acceptance criteria in prompts or this `SPEC.md`.
- Keep durable physics explanations, FC3 interpretation rules, and writing boundaries in `AGENTS.md`.
- Create a skill only after the FC3 parsing and visualization workflow is reused enough to stabilize into a general method.
- Use sub-agents only when tasks are independent, for example one worker for frontend interaction, one for parser validation, and one for documentation. This project is currently small enough to implement locally in one pass.
