#!/usr/bin/env python3
"""Plot complex-wave snapshots against the initial profile."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Result groups shown in the final-period comparison figures.
GROUPS = {
    "order3_four_periods": ["JS3_HCS4", "Z3_HCS4", "QR_HCS4"],
    "order5_four_periods": ["JS5_HCS6", "Z5_HCS6"],
}


def load(path):
    """Return x and u columns from one snapshot CSV."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1]


def snapshot_index(path):
    """Numeric period index from a name like Z5_HCS6_t3.csv."""
    return int(path.stem.split("_t")[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results")
    args = parser.parse_args()

    results = Path(args.results_dir)
    csv_dir = results / "csv"
    fig_dir = results / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # One evolution figure per scheme: initial profile plus every period.
    for path in sorted(csv_dir.glob("*_t0.csv")):
        name = path.stem[:-3]
        files = sorted(csv_dir.glob(f"{name}_t*.csv"), key=snapshot_index)
        plt.figure(figsize=(8, 4))
        for snapshot in files:
            x, u = load(snapshot)
            index = snapshot_index(snapshot)
            label = "initial" if index == 0 else f"t={index}T"
            plt.plot(x, u, label=label)
        plt.xlabel("x")
        plt.ylabel("u")
        plt.title(name)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / f"{name}_evolution.png", dpi=150)
        plt.close()

    # One comparison figure per order group using the final snapshot.
    for group, names in GROUPS.items():
        names = [name for name in names if (csv_dir / f"{name}_t0.csv").exists()]
        if not names:
            continue
        plt.figure(figsize=(8, 4))
        x0, u0 = load(csv_dir / f"{names[0]}_t0.csv")
        plt.plot(x0, u0, "k--", label="initial")
        for name in names:
            files = sorted(csv_dir.glob(f"{name}_t*.csv"), key=snapshot_index)
            x, u = load(files[-1])
            plt.plot(x, u, label=name)
        plt.xlabel("x")
        plt.ylabel("u")
        plt.title(group.replace("_", " "))
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / f"{group}.png", dpi=150)
        plt.close()


if __name__ == "__main__":
    main()
