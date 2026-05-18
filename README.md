# FC3 Visual Explorer

FC3 Visual Explorer is a browser-based tool for visualizing and analyzing
ShengBTE/thirdorder `FORCE_CONSTANTS_3RD` files. It converts third-order force
constant blocks into interactive three-body geometry and compact analysis plots.

## Features

- Interactive 3D FC3 three-body network visualization.
- Atom, bond, unit-cell, FC3 network, object-tree, highlight, hide, and pick-atom controls.
- FC3 component and block-norm distributions.
- FC3 block norm versus triangle perimeter and maximum edge length.
- Multi-file FC3 comparison for distribution-level analysis.
- Browser-side parsing for user-uploaded FC3 and Quantum ESPRESSO structure files.

## Online Deployment

This project is designed as a static web app.

### GitHub Pages

This repository includes a GitHub Actions workflow that publishes the `app`
directory to GitHub Pages.

After pushing to GitHub, enable Pages with:

- Repository `Settings`
- `Pages`
- Source: `GitHub Actions`

The expected public URL is:

```text
https://zxy-1104.github.io/FC3-Visual-Explorer/
```

### Cloudflare Pages

For Cloudflare Pages, connect this GitHub repository and use:

- Build command: leave empty
- Build output directory: `app`

The app uses relative asset paths, so the same repository can be deployed by
GitHub Pages, Cloudflare Pages, or any static file server.

## Run Locally

```bash
cd app
python -m http.server 8766
```

Then open:

```text
http://127.0.0.1:8766
```

## Input Files

The default example is embedded in `app/data/fc3_045.json`.

For user uploads, provide:

- A ShengBTE/thirdorder-compatible `FORCE_CONSTANTS_3RD` file.
- A matching Quantum ESPRESSO `pw.x` structure file containing lattice vectors
  and atomic positions.

Uploaded files are parsed in the browser. They are not sent to a server by this
static web app.

## Notes

FC3 block norms are useful for visual inspection and physical consistency
checks, but they should not be interpreted as direct predictors of lattice
thermal conductivity trends.
