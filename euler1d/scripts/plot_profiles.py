#!/usr/bin/env python3
"""Plot euler1d profiles against exact or reference solutions."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIELDS = [("rho", 1, r"$\rho$"), ("u", 2, r"$u$"), ("p", 3, r"$p$")]


def load(path):
    """Return the x,rho,u,p columns of one profile CSV."""
    return np.loadtxt(path, delimiter=",", skiprows=1)


def label(name):
    """Readable label for a <scheme>_<operator> suffix."""
    return name.replace("_", " ")


def reference_name(case):
    """Exact file for Sod/Lax, numerical reference for Shu-Osher."""
    if case == "shu_osher":
        return "reference_shu_osher.csv"
    return f"exact_{case}.csv"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results")
    args = parser.parse_args()

    results = Path(args.results_dir)
    csv_dir = results / "csv"
    fig_dir = results / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    for case in ("sod", "lax", "shu_osher"):
        ref_path = csv_dir / reference_name(case)
        if not ref_path.exists():
            continue
        ref = load(ref_path)
        numeric = sorted(p for p in csv_dir.glob(f"{case}_*.csv"))
        if not numeric:
            continue
        ref_label = "exact" if case != "shu_osher" else "reference"

        # Three-panel profile for each scheme/operator pair.
        for path in numeric:
            data = load(path)
            name = path.stem[len(case) + 1:]
            fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.4))
            for ax, (field, col, ylabel) in zip(axes, FIELDS):
                ax.plot(data[:, 0], data[:, col], label="numerical")
                ax.plot(ref[:, 0], ref[:, col], "k--", label=ref_label)
                ax.set_xlabel("x")
                ax.set_ylabel(ylabel)
                ax.set_title(field)
            axes[0].legend(fontsize=8)
            fig.suptitle(f"{case} {label(name)}")
            fig.tight_layout()
            fig.savefig(fig_dir / f"{case}_{name}.png", dpi=150)
            plt.close(fig)

        # Density comparison across schemes.
        fig, ax = plt.subplots(figsize=(8.0, 4.0))
        ax.plot(ref[:, 0], ref[:, 1], "k--", label=ref_label)
        for path in numeric:
            data = load(path)
            name = path.stem[len(case) + 1:]
            ax.plot(data[:, 0], data[:, 1], label=label(name))
        ax.set_xlabel("x")
        ax.set_ylabel(r"$\rho$")
        ax.set_title(f"{case} density")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{case}_comparison.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
