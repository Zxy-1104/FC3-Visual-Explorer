from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent
FC3_PATH = ROOT / "FORCE_CONSTANTS_3RD_0.45"
UNIT_QE_PATH = ROOT.parents[0] / "job_080" / "fc3_work" / "si_unit_relaxed.scf.in"

sys.path.insert(0, str(ROOT))
from _utils.fc3_utils import read_fc3


def read_qe_unit_cell(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    cell = []
    frac = []
    species = []

    for idx, line in enumerate(lines):
        if line.strip().lower().startswith("cell_parameters"):
            for row in lines[idx + 1 : idx + 4]:
                cell.append([float(x) for x in row.split()[:3]])
        if line.strip().lower().startswith("atomic_positions"):
            j = idx + 1
            while j < len(lines):
                parts = lines[j].split()
                if len(parts) < 4:
                    break
                try:
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                except ValueError:
                    break
                species.append(parts[0])
                frac.append(coords)
                j += 1
            break

    if len(cell) != 3 or not frac:
        raise ValueError(f"Failed to parse QE unit cell from {path}")

    cell_arr = np.asarray(cell, dtype=float)
    frac_arr = np.asarray(frac, dtype=float)
    cart_arr = frac_arr @ cell_arr
    return cell_arr, cart_arr, species


def cart_shift_to_lattice_index(shift_cart: tuple[float, float, float], cell: np.ndarray) -> tuple[int, int, int]:
    coeff = np.linalg.solve(cell.T, np.asarray(shift_cart, dtype=float))
    rounded = np.rint(coeff).astype(int)
    if not np.allclose(coeff, rounded, atol=2e-5):
        raise ValueError(f"Shift {shift_cart} is not an integer lattice vector; fractional coeff={coeff}")
    return tuple(int(x) for x in rounded)


def build_supercell_points(cell: np.ndarray, atom_cart: np.ndarray, species: list[str], span: int = 1):
    points = []
    labels = []
    lattice_indices = []
    atom_ids = []
    for a in range(-span, span + 1):
        for b in range(-span, span + 1):
            for c in range(-span, span + 1):
                shift = np.array([a, b, c], dtype=float) @ cell
                for atom_idx, atom_pos in enumerate(atom_cart, start=1):
                    points.append(atom_pos + shift)
                    labels.append(f"{species[atom_idx - 1]}{atom_idx} cell=({a},{b},{c})")
                    lattice_indices.append((a, b, c))
                    atom_ids.append(atom_idx)
    return np.asarray(points, dtype=float), labels, lattice_indices, atom_ids


def fc3_block_geometry(block, cell: np.ndarray, atom_cart: np.ndarray):
    i, j, k = (x - 1 for x in block.atoms)
    rj_idx = cart_shift_to_lattice_index(block.rj, cell)
    rk_idx = cart_shift_to_lattice_index(block.rk, cell)

    xi = atom_cart[i]
    xj = atom_cart[j] + np.asarray(rj_idx, dtype=float) @ cell
    xk = atom_cart[k] + np.asarray(rk_idx, dtype=float) @ cell
    return xi, xj, xk, rj_idx, rk_idx


def linewidth_from_norm(value: float, vmin: float, vmax: float) -> float:
    if vmax <= vmin:
        return 3.0
    t = (value - vmin) / (vmax - vmin)
    return float(1.2 + 4.2 * np.sqrt(max(0.0, min(1.0, t))))


def make_discrete_color(value: float, vmin: float, vmax: float, nbin: int = 10) -> str:
    if vmax <= vmin:
        t = 0.5
    else:
        t = (value - vmin) / (vmax - vmin)
    idx = int(np.clip(np.floor(t * nbin), 0, nbin - 1))
    palette = [
        "#440154",
        "#482878",
        "#3e4989",
        "#31688e",
        "#26828e",
        "#1f9e89",
        "#35b779",
        "#6ece58",
        "#b5de2b",
        "#fde725",
    ]
    return palette[idx]


def main() -> None:
    cell, atom_cart, species = read_qe_unit_cell(UNIT_QE_PATH)
    blocks = read_fc3(FC3_PATH)
    norms = np.asarray([b.norm for b in blocks], dtype=float)
    vmin = float(norms.min())
    vmax = float(norms.max())

    atom_points, atom_labels, _, _ = build_supercell_points(cell, atom_cart, species, span=1)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=atom_points[:, 0],
            y=atom_points[:, 1],
            z=atom_points[:, 2],
            mode="markers",
            marker=dict(size=4.5, color="rgba(80,80,80,0.72)", line=dict(width=0.5, color="black")),
            text=atom_labels,
            hovertemplate="%{text}<br>x=%{x:.3f} Å<br>y=%{y:.3f} Å<br>z=%{z:.3f} Å<extra></extra>",
            name="3x3x3 Si atoms",
            showlegend=True,
        )
    )

    edge_rows = [
        "block_id,i,j,k,Rj_lattice_a,Rj_lattice_b,Rj_lattice_c,Rk_lattice_a,Rk_lattice_b,Rk_lattice_c,block_norm"
    ]
    for block in blocks:
        xi, xj, xk, rj_idx, rk_idx = fc3_block_geometry(block, cell, atom_cart)
        color = make_discrete_color(block.norm, vmin, vmax)
        width = linewidth_from_norm(block.norm, vmin, vmax)
        x = [xi[0], xj[0], xk[0], xi[0], None]
        y = [xi[1], xj[1], xk[1], xi[1], None]
        z = [xi[2], xj[2], xk[2], xi[2], None]
        hover = (
            f"block {block.block_id}<br>"
            f"(i,j,k)=({block.atoms[0]},{block.atoms[1]},{block.atoms[2]})<br>"
            f"Rj lattice={rj_idx}<br>"
            f"Rk lattice={rk_idx}<br>"
            f"||Phi3||F={block.norm:.6g}"
        )
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color=color, width=width),
                hoverinfo="text",
                text=[hover] * len(x),
                name=f"block {block.block_id}",
                showlegend=False,
            )
        )
        edge_rows.append(
            ",".join(
                [
                    str(block.block_id),
                    str(block.atoms[0]),
                    str(block.atoms[1]),
                    str(block.atoms[2]),
                    str(rj_idx[0]),
                    str(rj_idx[1]),
                    str(rj_idx[2]),
                    str(rk_idx[0]),
                    str(rk_idx[1]),
                    str(rk_idx[2]),
                    f"{block.norm:.10g}",
                ]
            )
        )

    # Invisible colorbar trace, independent from the per-edge traces.
    fig.add_trace(
        go.Scatter3d(
            x=[None],
            y=[None],
            z=[None],
            mode="markers",
            marker=dict(
                size=0,
                color=[vmin, vmax],
                colorscale="Viridis",
                cmin=vmin,
                cmax=vmax,
                colorbar=dict(title=r"||Phi3||F", thickness=18),
                showscale=True,
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    # Draw primitive lattice vectors from the reference cell for orientation.
    origin = np.zeros(3)
    for idx, vec in enumerate(cell, start=1):
        fig.add_trace(
            go.Scatter3d(
                x=[origin[0], vec[0]],
                y=[origin[1], vec[1]],
                z=[origin[2], vec[2]],
                mode="lines+text",
                line=dict(color="black", width=5),
                text=["", f"a{idx}"],
                textposition="top center",
                name=f"a{idx}",
                showlegend=False,
            )
        )

    fig.update_layout(
        title=(
            "FC3 three-body network, cutoff = 0.45 nm<br>"
            "3x3x3 primitive Si cells; triangle edges colored and weighted by block norm"
        ),
        scene=dict(
            xaxis_title="x (Å)",
            yaxis_title="y (Å)",
            zaxis_title="z (Å)",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, b=0, t=70),
        template="plotly_white",
    )

    html_path = OUT_DIR / "fc3_three_body_network_045_3x3x3_interactive.html"
    fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
    (OUT_DIR / "fc3_three_body_network_045_lattice_indices.csv").write_text(
        "\n".join(edge_rows) + "\n", encoding="utf-8"
    )
    print(f"Saved {html_path}")


if __name__ == "__main__":
    main()
