import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
DFT_FC3_JSON = BASE_DIR / "BAs_FC3_Kappa_Data_DFT.json"
DFT_FC4_JSON = BASE_DIR / "BAs_FC4_Kappa_Data_DFT.json"
MATTERSIM_FC3_JSON = BASE_DIR / "BAs_FC3_Kappa_Data_Mattersim.json"
MATTERSIM_FC4_JSON = BASE_DIR / "BAs_FC4_Kappa_Data_Mattersim.json"
OUTPUT_PNG = BASE_DIR / "BAs_Kappa_Comparison.png"

COLOR_DFT = "#1f77b4"
COLOR_MATTERSIM = "#d62728"


def apply_prb_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "font.size": 12,
            "axes.labelsize": 13,
            "legend.fontsize": 10,
            "legend.frameon": False,
            "lines.linewidth": 1.8,
            "lines.markersize": 5.5,
        }
    )


def load_kappa_series(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    temperatures = sorted((float(key) for key in data.keys()))
    conv = [data[str(int(temp))]["CONV"]["kappa"] for temp in temperatures]
    rta = [data[str(int(temp))]["RTA"]["kappa"] for temp in temperatures]
    return np.array(temperatures), np.array(conv), np.array(rta)


def main():
    apply_prb_style()

    temp_dft_fc3, conv_dft_fc3, rta_dft_fc3 = load_kappa_series(DFT_FC3_JSON)
    temp_dft_fc4, conv_dft_fc4, rta_dft_fc4 = load_kappa_series(DFT_FC4_JSON)
    temp_ml_fc3, conv_ml_fc3, rta_ml_fc3 = load_kappa_series(MATTERSIM_FC3_JSON)
    temp_ml_fc4, conv_ml_fc4, rta_ml_fc4 = load_kappa_series(MATTERSIM_FC4_JSON)

    fig, ax = plt.subplots(figsize=(6.8, 4.8))

    ax.plot(temp_dft_fc3, conv_dft_fc3, color=COLOR_DFT, marker="s", linestyle="-", label="DFT FC3 CONV")
    ax.plot(
        temp_dft_fc3,
        rta_dft_fc3,
        color=COLOR_DFT,
        marker="s",
        linestyle="--",
        markerfacecolor="none",
        label="DFT FC3 RTA",
    )
    ax.plot(
        temp_dft_fc4,
        conv_dft_fc4,
        color="#9467bd",
        marker="D",
        linestyle="-",
        label="DFT FC4 CONV",
    )
    ax.plot(
        temp_dft_fc4,
        rta_dft_fc4,
        color="#9467bd",
        marker="D",
        linestyle="--",
        markerfacecolor="none",
        label="DFT FC4 RTA",
    )
    ax.plot(temp_ml_fc3, conv_ml_fc3, color=COLOR_MATTERSIM, marker="o", linestyle="-", label="MatterSim FC3 CONV")
    ax.plot(
        temp_ml_fc3,
        rta_ml_fc3,
        color=COLOR_MATTERSIM,
        marker="o",
        linestyle="--",
        markerfacecolor="none",
        label="MatterSim FC3 RTA",
    )
    ax.plot(
        temp_ml_fc4,
        conv_ml_fc4,
        color="#2ca02c",
        marker="^",
        linestyle="-",
        label="MatterSim FC4 CONV",
    )
    ax.plot(
        temp_ml_fc4,
        rta_ml_fc4,
        color="#2ca02c",
        marker="^",
        linestyle="--",
        markerfacecolor="none",
        label="MatterSim FC4 RTA",
    )

    ax.set_xlabel(r"Temperature $T$ (K)")
    ax.set_ylabel(r"Thermal Conductivity $\kappa$ (W m$^{-1}$ K$^{-1}$)")
    ax.set_xlim(80, 1020)
    ax.set_xticks(np.arange(100, 1001, 100))
    ax.legend(loc="upper right", ncol=1)

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plot to: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
