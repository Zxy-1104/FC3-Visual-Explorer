from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter


ROOT = Path(__file__).resolve().parents[1]
FC3_PATH = ROOT / "FORCE_CONSTANTS_3RD_0.45"
QE_PATH = ROOT / "BASE.si_supper.scf.in"
OUT_COMPONENT = ROOT / "preview_fc3_component_count_distribution.png"
OUT_BLOCK = ROOT / "preview_fc3_block_vs_perimeter.png"


def dot_vec(coeff, cell):
    return np.asarray(coeff, dtype=float) @ np.asarray(cell, dtype=float)


def read_qe_supercell(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    positions = []
    cell = None
    for idx, line in enumerate(lines):
        head = line.strip().lower()
        if head.startswith("atomic_positions"):
            row = idx + 1
            while row < len(lines):
                parts = lines[row].split()
                if len(parts) < 4:
                    break
                try:
                    frac = [float(parts[1]), float(parts[2]), float(parts[3])]
                except ValueError:
                    break
                positions.append({"species": parts[0], "frac": frac})
                row += 1
        if head.startswith("cell_parameters"):
            cell = [[float(x) for x in lines[idx + n].split()[:3]] for n in range(1, 4)]
    if cell is None or not positions:
        raise ValueError(f"Failed to parse QE structure from {path}")
    return np.asarray(cell, dtype=float), positions


def infer_primitive(supercell, positions):
    nat = len(positions)
    repeat = round((nat / 2) ** (1 / 3))
    if 2 * repeat**3 != nat:
        raise ValueError(f"Cannot infer two-atom cubic primitive repeat from nat={nat}")
    cell = supercell / repeat
    basis = []
    seen = set()
    for pos in positions:
        frac = []
        for value in pos["frac"]:
            raw = value * repeat
            frac.append(raw - math.floor(raw + 1e-9))
        key = tuple(round(x * 1_000_000) for x in frac)
        if key in seen:
            continue
        seen.add(key)
        basis.append({"species": pos["species"], "frac": frac, "cart": dot_vec(frac, cell)})
        if len(basis) == 2:
            break
    if len(basis) != 2:
        raise ValueError("Failed to infer two primitive basis atoms")
    return cell, basis


def lattice_index(shift_cart, cell):
    coeff = np.linalg.solve(cell.T, np.asarray(shift_cart, dtype=float))
    rounded = np.rint(coeff).astype(int)
    if not np.allclose(coeff, rounded, atol=2e-5):
        raise ValueError(f"Shift is not an integer lattice vector: {shift_cart} -> {coeff}")
    return rounded


def read_fc3(path: Path):
    lines = [line.strip() for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    nblocks = int(lines[0])
    blocks = []
    idx = 1
    for _ in range(nblocks):
        block_id = int(lines[idx])
        rj = np.array([float(x) for x in lines[idx + 1].split()[:3]], dtype=float)
        rk = np.array([float(x) for x in lines[idx + 2].split()[:3]], dtype=float)
        atoms = [int(x) for x in lines[idx + 3].split()[:3]]
        values = np.array([float(lines[idx + 4 + n].split()[3]) for n in range(27)], dtype=float)
        blocks.append({"id": block_id, "rj": rj, "rk": rk, "atoms": atoms, "values": values})
        idx += 31
    if len(blocks) != nblocks:
        raise ValueError(f"Expected {nblocks} blocks, parsed {len(blocks)}")
    return blocks


def build_block_rows(blocks, cell, basis):
    rows = []
    for block in blocks:
        atom_records = [basis[idx - 1] for idx in block["atoms"]]
        rj_idx = lattice_index(block["rj"], cell)
        rk_idx = lattice_index(block["rk"], cell)
        xi = atom_records[0]["cart"]
        xj = atom_records[1]["cart"] + dot_vec(rj_idx, cell)
        xk = atom_records[2]["cart"] + dot_vec(rk_idx, cell)
        perimeter = np.linalg.norm(xj - xi) + np.linalg.norm(xk - xi) + np.linalg.norm(xj - xk)
        norm = float(np.sqrt(np.sum(block["values"] ** 2)))
        rows.append({"id": block["id"], "perimeter": float(perimeter), "norm": norm, "values": block["values"]})
    return rows


def configure_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 15,
            "legend.fontsize": 11,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "lines.linewidth": 2.0,
        }
    )


def decimal_log_formatter(value, _pos=None):
    if value <= 0:
        return ""
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.0f}"
    if value >= 1:
        return f"{value:.0f}" if abs(value - round(value)) < 1e-10 else f"{value:g}"
    if value >= 0.1:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    if value >= 0.01:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    if value >= 0.001:
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return f"{value:.0e}"


def apply_closed_axes(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.4)
    ax.tick_params(which="major", length=7, width=1.3, top=False, right=False)
    ax.tick_params(which="minor", length=4, width=1.0, top=False, right=False)


def apply_log_ticks(ax, axis: str):
    locator = LogLocator(base=10.0, subs=(1.0,))
    minor_locator = LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    formatter = FuncFormatter(decimal_log_formatter)
    if axis == "x":
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_minor_locator(minor_locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_major_locator(locator)
        ax.yaxis.set_minor_locator(minor_locator)
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.set_minor_formatter(NullFormatter())


def minor_log_edges(values):
    positive = np.asarray(values, dtype=float)
    positive = positive[np.isfinite(positive) & (positive > 0)]
    if positive.size == 0:
        raise ValueError("No positive values for log bins")
    vmin = float(np.min(positive))
    vmax = float(np.max(positive))
    emin = math.floor(math.log10(vmin))
    emax = math.ceil(math.log10(vmax))
    edges = []
    for exponent in range(emin, emax + 1):
        for multiplier in range(1, 10):
            edges.append(multiplier * 10.0**exponent)
    edges.append(10.0**emax)
    edges = np.array(sorted(set(edges)), dtype=float)
    edges = edges[(edges >= vmin / 1.000001) & (edges <= vmax * 1.000001)]
    if edges[0] > vmin:
        edges = np.insert(edges, 0, vmin)
    if edges[-1] < vmax:
        edges = np.append(edges, vmax)
    return edges


def plot_count_curve(ax, values, label, color):
    edges = minor_log_edges(values)
    counts, edges = np.histogram(values, bins=edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    nonzero = counts > 0
    plotted_label = False
    nz_idx = np.flatnonzero(nonzero)
    if nz_idx.size == 0:
        return
    ax.plot(
        centers[nz_idx],
        counts[nz_idx],
        linestyle="None",
        marker="o",
        color=color,
        markersize=4.2,
        label=label,
    )
    plotted_label = True
    for left, right in zip(nz_idx[:-1], nz_idx[1:]):
        linestyle = "-" if right == left + 1 else "--"
        ax.plot(
            [centers[left], centers[right]],
            [counts[left], counts[right]],
            color=color,
            linestyle=linestyle,
            marker=None,
            label=None if plotted_label else label,
        )


def check_component_symmetry(values, rtol=1e-8, atol=1e-10):
    positive = np.sort(np.abs(values[values > 0]))
    negative = np.sort(np.abs(values[values < 0]))
    zero_count = int(np.sum(values == 0))
    result = {
        "positive": positive,
        "negative": negative,
        "zero_count": zero_count,
        "same_count": positive.size == negative.size,
        "symmetric": False,
        "max_abs_diff": None,
        "largest_diffs": [],
    }
    if positive.size == negative.size and positive.size > 0:
        diffs = np.abs(positive - negative)
        result["max_abs_diff"] = float(np.max(diffs))
        result["symmetric"] = bool(np.allclose(positive, negative, rtol=rtol, atol=atol))
        worst = np.argsort(diffs)[-5:][::-1]
        result["largest_diffs"] = [
            (int(idx), float(positive[idx]), float(negative[idx]), float(diffs[idx])) for idx in worst
        ]
    return result


def plot_component_distribution(values, symmetry):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    if symmetry["symmetric"]:
        positive = values[values > 0]
        plot_count_curve(ax, positive, "component", "#005eb8")
    else:
        positive = np.abs(values[values > 0])
        negative = np.abs(values[values < 0])
        plot_count_curve(ax, positive, "positive", "#005eb8")
        plot_count_curve(ax, negative, "negative abs", "#c62828")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("FC3_component (eV/Å³)")
    ax.set_ylabel("count")
    ax.set_title("FC3_component count distribution")
    ax.grid(True, which="major", color="#bdbdbd", linewidth=0.8, alpha=0.8)
    ax.grid(True, which="minor", color="#dddddd", linewidth=0.5, alpha=0.65)
    ax.legend(loc="upper right", frameon=False)
    apply_log_ticks(ax, "x")
    apply_log_ticks(ax, "y")
    apply_closed_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT_COMPONENT, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_block_vs_perimeter(rows):
    perimeters = np.array([row["perimeter"] for row in rows], dtype=float)
    norms = np.array([row["norm"] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.scatter(perimeters, norms, s=30, color="#c62828", alpha=0.85, label="FC3 blocks")

    ax.set_yscale("log")
    ax.set_xlabel("perimeter (Å)")
    ax.set_ylabel("FC3_block (eV/Å³)")
    ax.set_title("FC3_block vs perimeter")
    ax.grid(True, which="major", color="#bdbdbd", linewidth=0.8, alpha=0.8)
    ax.grid(True, which="minor", color="#dddddd", linewidth=0.5, alpha=0.65)
    ax.legend(loc="upper right", frameon=False)
    apply_log_ticks(ax, "y")
    apply_closed_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT_BLOCK, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    configure_style()
    supercell, positions = read_qe_supercell(QE_PATH)
    cell, basis = infer_primitive(supercell, positions)
    blocks = read_fc3(FC3_PATH)
    rows = build_block_rows(blocks, cell, basis)
    values = np.concatenate([row["values"] for row in rows])
    symmetry = check_component_symmetry(values)

    plot_component_distribution(values, symmetry)
    plot_block_vs_perimeter(rows)

    print(f"positive count: {symmetry['positive'].size}")
    print(f"negative count: {symmetry['negative'].size}")
    print(f"zero count: {symmetry['zero_count']}")
    print(f"same positive/negative count: {symmetry['same_count']}")
    print(f"absolute-value symmetric within rtol=1e-8, atol=1e-10: {symmetry['symmetric']}")
    if symmetry["max_abs_diff"] is not None:
        print(f"max abs diff between sorted abs positive/negative: {symmetry['max_abs_diff']:.12g}")
        print("largest abs-diff samples: index, positive_abs, negative_abs, diff")
        for sample in symmetry["largest_diffs"]:
            print(f"  {sample[0]}, {sample[1]:.12g}, {sample[2]:.12g}, {sample[3]:.12g}")
    zero_perimeter_norms = [row["norm"] for row in rows if row["perimeter"] == 0]
    if zero_perimeter_norms:
        print(
            "zero-perimeter blocks: "
            f"{len(zero_perimeter_norms)}, FC3_block range "
            f"{min(zero_perimeter_norms):.12g} to {max(zero_perimeter_norms):.12g}"
        )
    print(f"saved: {OUT_COMPONENT}")
    print(f"saved: {OUT_BLOCK}")


if __name__ == "__main__":
    main()
