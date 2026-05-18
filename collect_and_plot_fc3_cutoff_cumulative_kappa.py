import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors


def resolve_work_dir():
    script_dir = Path(__file__).resolve().parent
    candidates = [
        Path.cwd(),
        script_dir,
        script_dir.parent,
    ]
    for candidate in candidates:
        if any(candidate.glob("run_fc3_cut_*")):
            return candidate.resolve()
    return Path.cwd().resolve()


WORK_DIR = resolve_work_dir()
RESULT_DIR = WORK_DIR / "result"
PLOT_PNG = RESULT_DIR / "fc3_cutoff_cumulative_kappa_vs_mfp_colored_by_cutoff.png"


def parse_cutoff_value(folder_name: str):
    prefix = "run_fc3_cut_"
    if not folder_name.startswith(prefix):
        return None
    raw = folder_name[len(prefix):]
    try:
        return float(raw)
    except ValueError:
        return None


def read_kappa_xx_conv(conv_file: Path):
    with conv_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                temp = float(parts[0])
                if abs(temp - 300.0) < 1.0:
                    return float(parts[1])
            except ValueError:
                continue
    return None


def find_preferred_file(run_dir: Path, relative_candidates):
    for relative_path in relative_candidates:
        direct_path = run_dir / relative_path
        if direct_path.exists():
            return direct_path

    matched_paths = []
    for relative_path in relative_candidates:
        matched_paths.extend(run_dir.rglob(relative_path.name))

    if not matched_paths:
        return None

    matched_paths = sorted(
        set(matched_paths),
        key=lambda path: (
            0 if "T300K" in path.parts else 1,
            len(path.parts),
            str(path),
        ),
    )
    return matched_paths[0]


def collect_data():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    run_dirs = glob.glob(str(WORK_DIR / "run_fc3_cut_*"))
    parsed_dirs = []
    for path_str in run_dirs:
        path = Path(path_str)
        cutoff = parse_cutoff_value(path.name)
        if cutoff is None:
            continue
        parsed_dirs.append((cutoff, path))
    parsed_dirs.sort(key=lambda item: item[0])

    rows = []
    for actual_cutoff, run_dir in parsed_dirs:

        conv_file = find_preferred_file(
            run_dir,
            [
                Path("BTE.KappaTensorVsT_CONV"),
                Path("T300K") / "BTE.KappaTensorVsT_CONV",
            ],
        )
        cumulative_file = find_preferred_file(
            run_dir,
            [
                Path("T300K") / "BTE.cumulative_kappa_scalar",
                Path("BTE.cumulative_kappa_scalar"),
            ],
        )
        if conv_file is None or cumulative_file is None:
            continue

        if conv_file.parent != cumulative_file.parent:
            sibling_conv_file = cumulative_file.parent / "BTE.KappaTensorVsT_CONV"
            if sibling_conv_file.exists():
                conv_file = sibling_conv_file

        kappa_conv = read_kappa_xx_conv(conv_file)
        if kappa_conv is None:
            continue

        try:
            curve = np.loadtxt(cumulative_file, dtype=float)
        except Exception:
            continue

        curve = np.atleast_2d(curve)
        if curve.shape[1] < 2:
            continue

        mfp = curve[:, 0]
        cumulative_kappa = curve[:, -1]
        mask = np.isfinite(mfp) & np.isfinite(cumulative_kappa) & (mfp > 0.0)
        mfp = mfp[mask]
        cumulative_kappa = cumulative_kappa[mask]
        if mfp.size == 0:
            continue

        rows.append(
            {
                "cutoff": float(actual_cutoff),
                "kappa": float(kappa_conv),
                "mfp": mfp,
                "cum_kappa": cumulative_kappa,
            }
        )

    if not rows:
        raise RuntimeError(
            f"No valid run_fc3_cut_* cumulative kappa data were found under {WORK_DIR}."
        )
    return rows


def plot_curves(rows):
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "font.size": 12,
            "axes.labelsize": 14,
            "lines.linewidth": 1.8,
        }
    )

    rows = sorted(rows, key=lambda row: row["cutoff"])
    cutoff_values = np.array([row["cutoff"] for row in rows], dtype=float)
    cmap = cm.get_cmap("turbo")
    norm = colors.Normalize(vmin=float(np.min(cutoff_values)), vmax=float(np.max(cutoff_values)))

    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    for row in rows:
        ax.plot(
            row["mfp"],
            row["cum_kappa"],
            color=cmap(norm(row["cutoff"])),
            alpha=0.92,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Mean free path")
    ax.set_ylabel(r"Cumulative $\kappa$ (W m$^{-1}$ K$^{-1}$)")
    ax.set_title(r"Cumulative $\kappa$ Colored by FC3 Cutoff")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("FC3 cutoff distance")

    plt.tight_layout()
    plt.savefig(PLOT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    rows = collect_data()
    plot_curves(rows)
    plotted_cutoffs = [
        f"{row['cutoff']:.3f} (kappa={row['kappa']:.2f})"
        for row in sorted(rows, key=lambda row: row["cutoff"])
    ]
    print(f"Plotted cutoff groups: {plotted_cutoffs}")
    print(f"Saved: {PLOT_PNG}")


if __name__ == "__main__":
    main()
