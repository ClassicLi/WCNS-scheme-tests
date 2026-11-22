#!/usr/bin/env python3
"""Plot ADR dispersion/dissipation curves from adr/results CSV files."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# CSV file stems grouped by the figures they belong to.
GROUPS = {
    "weno3": [
        "WENO3_LINEAR", "WENO3_JS", "WENO3_LLR",
        "BAEZA_WENO3", "TENO_MOD_SEMI_ENO", "WENO_C",
    ],
    "weno5": ["WENO5_LINEAR", "WENO5_JS", "WENO5_M", "WENO5_Z", "TENO_333"],
    "wcns3": [
        "LINEAR3_M1P1_HCS4", "JS3_HCS4", "Z3_HCS4", "QR_HCS4",
        "MOC322_HCS4", "MOC33_GA_HCS4",
    ],
    "wcns5": [
        "LINEAR5_M2P2_HCS4", "JS5_HCS4", "Z5_HCS4",
        "MOC53_LINEAR_HCS4", "MOZ_HCS4", "RBF_GA5_HCS4",
    ],
}

OPERATOR_BASES = ["LINEAR5_M2P2", "Z5", "MOC53_LINEAR"]
OPERATORS = ["WCNSE2", "WCNSE4", "WCNSE6", "HCS4"]
GROUPS["operators"] = [f"{b}_{o}" for b in OPERATOR_BASES for o in OPERATORS]


def select_files(group: str, csv_dir: Path):
    """Return the CSV files that belong to a group."""
    if group == "quick":
        return sorted(csv_dir.glob("*.csv"))
    return [csv_dir / f"{name}.csv" for name in GROUPS[group]
            if (csv_dir / f"{name}.csv").exists()]


def plot_group(group: str, dir_path: Path):
    """Write dispersion and dissipation figures for one group."""
    csv_dir = dir_path / "csv"
    fig_dir = dir_path / "figures"
    files = select_files(group, csv_dir)
    if not files:
        print(f"plot_adr: no CSV files for group {group} in {csv_dir}")
        return
    fig_dir.mkdir(parents=True, exist_ok=True)

    curves = []
    for path in files:
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        curves.append((path.stem.replace("_", " "), data))

    for kind, index, ylabel, suffix in [
        ("dispersion", 1, r"$Re(\Phi)$", "dispersion"),
        ("dissipation", 2, r"$-Im(\Phi)$", "dissipation"),
    ]:
        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        for label, data in curves:
            values = data[:, index] if kind == "dispersion" else -data[:, index]
            ax.plot(data[:, 0], values, linewidth=1.2, label=label)
        if kind == "dispersion":
            ax.plot([0, np.pi], [0, np.pi], "k--", linewidth=0.8, label="exact")
        else:
            ax.axhline(0.0, color="k", linestyle="--", linewidth=0.8)
        ax.set_xlabel(r"$\varphi$")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, np.pi)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{group}_{suffix}.png", dpi=150)
        plt.close(fig)

    print(f"plot_adr: wrote {group} figures to {fig_dir}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=None, help="results directory")
    parser.add_argument("--group", default="quick",
                        choices=sorted(GROUPS) + ["quick"])
    parser.add_argument("--quick", action="store_true",
                        help="alias for --group quick")
    args = parser.parse_args()

    group = "quick" if args.quick else args.group
    if args.dir is None:
        dir_path = Path("results/quick") if args.quick else Path("results")
    else:
        dir_path = Path(args.dir)
    plot_group(group, dir_path)


if __name__ == "__main__":
    main()
