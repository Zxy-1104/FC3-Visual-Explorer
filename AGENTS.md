# Project Writing Context: FC3, thirdorder, and ShengBTE UQ

This project studies uncertainty quantification (UQ) of lattice thermal conductivity calculated with MLIP/MatterSim-assisted ShengBTE workflows. The goal is not only to explain microscopic phonon physics, but to use MLIP to reduce DFT sampling cost, quantify how computational parameters affect thermal conductivity, and provide practical guidance for first-principles ShengBTE users.

The following notes preserve the key context from the FC3 discussion and should be used when writing, explaining, or extending this project.

## Core Scope

- Treat FC3 analysis as a supporting physical consistency check, not the main UQ objective.
- Main UQ objective: quantify uncertainty in `kappa` caused by computational settings such as supercell, atomic displacement, FC3 cutoff, and q-mesh.
- FC3-specific analysis is useful for explaining what changes when cutoff/supercell/displacement changes, but it should not be forced to directly fit thermal conductivity trends.
- Avoid overclaiming simple monotonic relationships between FC3 norm and thermal conductivity.

## FC3 Definition and Physical Meaning

Crystal potential energy near equilibrium is expanded as:

```latex
U = U_0
+ \frac{1}{2!}
\sum_{ij}\sum_{\alpha\beta}
\Phi_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta
+ \frac{1}{3!}
\sum_{ijk}\sum_{\alpha\beta\gamma}
\Phi_{ijk}^{\alpha\beta\gamma}
u_i^\alpha u_j^\beta u_k^\gamma
+ \cdots .
```

Second- and third-order IFCs are:

```latex
\Phi_{ij}^{\alpha\beta}
=
\frac{\partial^2 U}
{\partial u_i^\alpha \partial u_j^\beta},
\qquad
\Phi_{ijk}^{\alpha\beta\gamma}
=
\frac{\partial^3 U}
{\partial u_i^\alpha
\partial u_j^\beta
\partial u_k^\gamma}.
```

Because force is:

```latex
F_i^\alpha = -\frac{\partial U}{\partial u_i^\alpha},
```

FC3 can be obtained from force curvature:

```latex
\Phi_{ijk}^{\alpha\beta\gamma}
=
-
\frac{\partial^2 F_i^\alpha}
{\partial u_j^\beta \partial u_k^\gamma}.
```

Interpretation:

- FC2 describes the linear force-displacement response, like a harmonic spring constant.
- FC3 describes whether this effective spring constant changes with deformation.
- For repeated displacement indices,

```latex
\Phi_{ijj}^{\alpha\beta\beta}
=
-
\frac{\partial^2 F_i^\alpha}
{\partial (u_j^\beta)^2},
```

this is still a third-order force constant, not FC2. It measures the curvature of force as one displacement freedom changes.

## Finite Difference Logic

For one displacement degree of freedom:

```latex
\Phi_{ijj}^{\alpha\beta\beta}
\approx
-
\frac{
F_i^\alpha(+\Delta)
-2F_i^\alpha(0)
+F_i^\alpha(-\Delta)
}{\Delta^2}.
```

For two different displacement degrees of freedom:

```latex
\Phi_{ijk}^{\alpha\beta\gamma}
\approx
-
\frac{
F_i^\alpha(+\Delta,+\Delta)
-F_i^\alpha(+\Delta,-\Delta)
-F_i^\alpha(-\Delta,+\Delta)
+F_i^\alpha(-\Delta,-\Delta)
}
{4\Delta^2}.
```

Important explanation:

- The calculated force after displacement is the total force, not only a third-order force contribution.
- The finite-difference combination cancels undesired low-order terms and retains the target derivative.
- In mixed central difference, the operator keeps terms odd in both variables, e.g. `xy`, and cancels constants, single-variable linear terms, and single-variable quadratic terms.

For a monomial `x^m y^n`, the mixed-difference numerator contains:

```latex
\Delta^{m+n}
\left[
1-(-1)^n-(-1)^m+(-1)^{m+n}
\right].
```

Only terms with both `m` and `n` odd survive. The target `xy` survives; higher-order terms such as `x^3y` and `xy^3` survive as `O(Delta^2)` truncation errors.

## Finite-Difference Error Sources

There are two main error classes.

Numerical force noise:

```latex
\epsilon_{\mathrm{noise}}
\sim
\frac{\epsilon_F}{\Delta^2}.
```

High-order anharmonic truncation error:

```latex
\epsilon_{\mathrm{anharm}}
\sim
C\Delta^2.
```

Combined heuristic:

```latex
\epsilon_{\mathrm{total}}
\sim
\frac{\epsilon_F}{\Delta^2}
+C\Delta^2.
```

Writing guidance:

- If displacement is too small, force noise is amplified.
- If displacement is too large, fifth- and higher-order anharmonic terms contaminate the FC3 estimate.
- Atomic displacement controls finite-difference accuracy, not interaction range.

## Displacement Count Estimate for a 128-Atom Si Supercell

For a naive finite-displacement calculation with `N=128` atoms:

```latex
M = 3N = 384.
```

Assumptions for this estimate:

- No FC3 cutoff.
- No space-group symmetry.
- No translational invariance.
- No thirdorder irreducible triplet reduction.
- Only the mathematical equality of mixed derivative order is used, so `(q_a,q_b)` and `(q_b,q_a)` are the same unordered degree-of-freedom pair.

Different degree-of-freedom pairs:

```latex
\frac{M(M-1)}{2}
=
\frac{384\times383}{2}
=
73536.
```

Each requires four double-displacement structures, so:

```latex
73536\times4=294144.
```

Same degree-of-freedom pairs require two single-displacement structures each:

```latex
384\times2=768.
```

Including one equilibrium structure:

```latex
294144+768+1=294913.
```

This huge number is why `thirdorder` relies on symmetry, cutoff, and force-output reuse.

## Symmetry in FC3

### Exchange Symmetry

FC3 satisfies exchange symmetry of full atom-direction pairs:

```latex
\Phi_{ijk}^{\alpha\beta\gamma}
=
\Phi_{ikj}^{\alpha\gamma\beta}.
```

This means:

```latex
\frac{\partial^2 F_i^\alpha}
{\partial u_j^\beta \partial u_k^\gamma}
=
\frac{\partial^2 F_i^\alpha}
{\partial u_k^\gamma \partial u_j^\beta}.
```

Do not confuse this with swapping only physical directions. In general:

```latex
\frac{\partial^2 F_i^\alpha}
{\partial u_2^x \partial u_3^y}
\ne
\frac{\partial^2 F_i^\alpha}
{\partial u_2^y \partial u_3^x}.
```

### Space-Group Symmetry

For a space-group operation:

```latex
g = \{R|\mathbf{t}\},
```

FC3 transforms as:

```latex
\Phi_{i'j'k'}^{\alpha'\beta'\gamma'}
=
\sum_{\alpha\beta\gamma}
R_{\alpha'\alpha}
R_{\beta'\beta}
R_{\gamma'\gamma}
\Phi_{ijk}^{\alpha\beta\gamma}.
```

This maps symmetry-equivalent atom triplets and tensor components onto one another, reducing independent FC3 blocks and displacement structures.

### Translational Invariance

Overall translation of the crystal should not change the potential energy:

```latex
U(\{\mathbf{r}_i+\mathbf{a}\})
=
U(\{\mathbf{r}_i\}).
```

For FC3:

```latex
\sum_i
\Phi_{ijk}^{\alpha\beta\gamma}
=0.
```

This is an acoustic sum-rule-like constraint. It helps ensure that acoustic modes near Gamma behave correctly.

Gamma-region writing guidance:

- The Gamma region corresponds to long-wavelength acoustic phonons.
- These modes often have large group velocity and long lifetime, especially in high-thermal-conductivity crystals such as Si.
- If FC2 or FC3 violates translational invariance, acoustic branches and low-frequency scattering near Gamma can be polluted, strongly affecting `kappa`.

## thirdorder Internal Workflow

Use this implementation logic when explaining `thirdorder`.

Overall flow:

```text
unit cell
-> symmetry operations
-> supercell
-> distance table
-> irreducible FC3 set
-> displaced supercells
```

For Si, a primitive cell with 2 atoms and a `4x4x4` supercell gives:

```latex
2\times4\times4\times4 = 128
```

atoms.

Actual `thirdorder` logic:

1. Read unit cell and identify symmetry operations using spglib.
2. Generate supercell.
3. Compute minimum periodic distances between atoms.
4. Keep only triplets satisfying cutoff:

```latex
d_{ij}<r_{\mathrm{cut}},
\qquad
d_{ik}<r_{\mathrm{cut}},
\qquad
d_{jk}<r_{\mathrm{cut}}.
```

5. For each candidate triplet, scan six permutations of `(i,j,k)`.
6. Apply all space-group operations to each permuted triplet.
7. Normalize equivalent triplets by translation so that the first atom is assigned to the reference cell.
8. Construct `27 x 27` transformation matrices for the 27 Cartesian FC3 components:

```latex
\boldsymbol{\Phi}'
=
\mathbf{T}\boldsymbol{\Phi}.
```

9. If a symmetry maps the triplet back to itself, impose:

```latex
(\mathbf{T}-\mathbf{I})\boldsymbol{\Phi}=0.
```

10. Use linear algebra to identify independent tensor components and minimal atom-direction displacement pairs.
11. For each irreducible displacement pair, generate four structures:

```latex
(+\Delta,+\Delta),
\quad
(+\Delta,-\Delta),
\quad
(-\Delta,+\Delta),
\quad
(-\Delta,-\Delta).
```

12. In `reap`, read forces and reconstruct the irreducible FC3 by central difference.
13. Expand irreducible FC3 back to full `FORCE_CONSTANTS_3RD`.

Important implementation detail:

- One double-displacement force calculation returns forces on all atoms and all Cartesian directions.
- Therefore one calculation supplies many `i,alpha` responses for a fixed `(j,beta),(k,gamma)` pair.

## FORCE_CONSTANTS_3RD Format

Each block has the form:

```text
block_id
Rj_x Rj_y Rj_z
Rk_x Rk_y Rk_z
i j k
alpha beta gamma value
...
```

Interpretation:

```latex
\Phi_{ijk}^{\alpha\beta\gamma}
(\mathbf{R}_j,\mathbf{R}_k)
=
-
\frac{\partial^2 F_i^\alpha(0)}
{\partial u_j^\beta(\mathbf{R}_j)
\partial u_k^\gamma(\mathbf{R}_k)}.
```

Key points:

- `block_id` is only the write-order ID, not physical importance.
- `i,j,k` are primitive-cell atom labels, not absolute 128-atom supercell indices.
- `i` is assigned to the reference primitive cell, usually cell `(0,0,0)` by convention.
- `Rj` and `Rk` are periodic lattice translation vectors in Cartesian coordinates.
- The true coordinates are:

```latex
\mathbf{x}_i = \mathbf{r}_i,
\qquad
\mathbf{x}_j = \mathbf{r}_j+\mathbf{R}_j,
\qquad
\mathbf{x}_k = \mathbf{r}_k+\mathbf{R}_k.
```

For relative plotting, put the reference force atom at the origin:

```latex
\mathbf{x}_i = 0,
\qquad
\mathbf{x}_j = \mathbf{r}_j+\mathbf{R}_j-\mathbf{r}_i,
\qquad
\mathbf{x}_k = \mathbf{r}_k+\mathbf{R}_k-\mathbf{r}_i.
```

Important clarification:

- `Rj = 0 0 0` does not mean atom `j` is at coordinate zero.
- It means atom `j` is in the reference primitive cell image.
- Its real position still includes its primitive-cell basis coordinate `r_j`.

## Primitive Cell, Reference Cell, and Images

Use the following language to avoid confusion.

- Reference cell or base cell: the primitive cell chosen as `(0,0,0)` for bookkeeping.
- Periodic image: another copy of the same primitive cell shifted by integer lattice vectors.
- The reference cell is not necessarily the geometric center of the supercell.
- It is a convention. Periodicity makes all choices physically equivalent.

For the Si primitive cell used in this project:

```latex
\mathbf{a}_1=(0,\ 2.731994625,\ 2.731994625),
```

```latex
\mathbf{a}_2=(2.731994625,\ 0,\ 2.731994625),
```

```latex
\mathbf{a}_3=(2.731994625,\ 2.731994625,\ 0).
```

Therefore a Cartesian vector:

```text
0 2.731994625 2.731994625
```

equals the lattice vector `a1`, i.e. lattice image index `(1,0,0)`, not Cartesian `(0,1,1)`.

Finite-displacement amplitude is a different quantity. `thirdorder` usually displaces atoms by a small value such as:

```latex
\Delta \approx 0.01\ \mathrm{\AA}.
```

Long-range interaction is represented by moving atoms located farther away, not by increasing displacement amplitude.

## FC3 and Thermal Conductivity

ShengBTE projects real-space FC3 into phonon mode space to form three-phonon matrix elements:

```latex
V_3(\lambda,\lambda',\lambda'')
\propto
\sum_{ijk}
\sum_{\alpha\beta\gamma}
\frac{
\Phi_{ijk}^{\alpha\beta\gamma}
e_i^\alpha(\lambda)
e_j^\beta(\lambda')
e_k^\gamma(\lambda'')
}
{\sqrt{
M_iM_jM_k
\omega_\lambda
\omega_{\lambda'}
\omega_{\lambda''}
}}
e^{i(\cdots)}.
```

Three-phonon scattering approximately follows:

```latex
\Gamma_\lambda
=
\frac{1}{\tau_\lambda}
\propto
\sum_{\lambda'\lambda''}
\left|
V_3(\lambda,\lambda',\lambda'')
\right|^2
\delta(
\omega_\lambda
\pm
\omega_{\lambda'}
-
\omega_{\lambda''}
)
\mathcal{N}(T).
```

Pathway:

```text
FC3 -> V3 -> scattering rate -> lifetime -> thermal conductivity
```

If all FC3 terms were scaled up uniformly while other quantities stayed fixed, stronger FC3 would generally increase scattering and reduce `kappa`. In real cutoff/supercell/displacement tests, FC3 changes are not uniform, so `sum ||Phi3||` or block norm distribution cannot be treated as a direct predictor of thermal conductivity.

## FC3 Visualization Standards

Use the following visualization approach for FC3.

### Block Norm

For each block:

```latex
\left\|
\Phi_{ijk}^{(3)}
\right\|_F
=
\sqrt{
\sum_{\alpha\beta\gamma}
\left(
\Phi_{ijk}^{\alpha\beta\gamma}
\right)^2
}.
```

This is the preferred scalar strength descriptor for one FC3 block.

### Recommended Plots

1. `FC3 block norm distribution`
   - x-axis: `||Phi3_block||_F`, usually log scale.
   - Purpose: show strong/weak block distribution and added weak terms as cutoff increases.

2. `FC3 strength vs dmax shell`
   - Define:

```latex
d_{\max}
=
\max
\left(
|\mathbf{x}_j-\mathbf{x}_i|,
|\mathbf{x}_k-\mathbf{x}_i|,
|\mathbf{x}_j-\mathbf{x}_k|
\right).
```

   - y-axis:

```latex
\sum_{\mathrm{shell}}
\|\Phi^{(3)}_{\mathrm{block}}\|_F.
```

   - Purpose: show spatial decay of FC3.

3. `common block parity plot`
   - Compare common block norms between cutoffs, e.g. small cutoff vs large cutoff.
   - Use log-log parity line `y=x`.
   - Purpose: test whether common short-range terms are stable.

4. `Delta FC3 shell plot`
   - For common blocks:

```latex
\|\Delta \Phi^{(3)}_{\mathrm{block}}\|_F
=
\|\Phi^{(3)}_{\mathrm{large\,cutoff}}
-
\Phi^{(3)}_{\mathrm{small\,cutoff}}\|_F.
```

   - Purpose: identify which distance shell changes most.

5. `Interactive FC3 three-body network`
   - Build a `3x3x3` primitive-cell Si lattice.
   - Convert each Cartesian `Rj/Rk` to integer primitive-lattice image indices.
   - For each block, draw a triangle with vertices:

```latex
\mathbf{x}_i=\mathbf{r}_i,
\qquad
\mathbf{x}_j=\mathbf{r}_j+\mathbf{R}_j,
\qquad
\mathbf{x}_k=\mathbf{r}_k+\mathbf{R}_k.
```

   - Do not fill triangles unless explicitly requested.
   - Draw only three edges of each triangle.
   - Edge color maps to `||Phi3_block||_F`.
   - Edge width also maps to norm, but keep width within a readable range, e.g. roughly 1.2 to 5.4 px.
   - Output as Plotly HTML when interactive rotation is needed.

### Existing Local Outputs

In this project, relevant generated files include:

- `FC3_derivation_notes.txt`: full LaTeX-formatted FC3 derivation notes.
- `FC3_analyze/05_three_body_network_035`: 0.35 nm static FC3 network plots.
- `FC3_analyze/06_three_body_network_045/plot_fc3_3x3x3_interactive_html.py`: script for interactive 3D FC3 network.
- `FC3_analyze/06_three_body_network_045/fc3_three_body_network_045_3x3x3_interactive.html`: interactive Plotly HTML for 0.45 nm cutoff.
- `FC3_analyze/06_three_body_network_045/fc3_three_body_network_045_lattice_indices.csv`: conversion of Cartesian `Rj/Rk` into primitive-cell image indices.

## Writing Style for This Topic

- Present formulas in LaTeX display format, not inline text-only equations.
- Avoid vague wording such as "FC3 is proportional to kappa." It is not generally true.
- Clearly separate:
  - finite displacement amplitude `Delta`,
  - periodic image translation `Rj/Rk`,
  - actual interatomic distance,
  - primitive-cell basis coordinate `r_i`.
- When explaining `FORCE_CONSTANTS_3RD`, always state that `Rj/Rk = 0` means same reference primitive cell image, not zero atomic coordinate.
- When discussing FC3 visualizations, make clear that a block's triangle geometry is based on atom basis positions plus periodic translations.
- If producing a document, use `ctexart` and compile with XeLaTeX for Chinese text.
