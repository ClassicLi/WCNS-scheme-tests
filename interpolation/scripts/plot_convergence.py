#!/usr/bin/env python3
"""Plot L1 convergence curves from results/convergence/*.csv."""

import argparse
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    n_values = []
    l1_values = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            n_values.append(float(row["N"]))
            l1_values.append(float(row["L1_error"]))
    return n_values, l1_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(root, "results", "convergence")
    files = sorted(glob.glob(os.path.join(outdir, "*.csv")))
    if not files:
        raise SystemExit("no CSV files in results/convergence")

    fig, ax = plt.subplots(figsize=(8, 5))
    for path in files:
        n_values, l1_values = read_csv(path)
        if not n_values:
            continue
        label = os.path.basename(path)[:-4]
        ax.loglog(n_values, l1_values, marker="o", linewidth=1, label=label)

        fig_one, ax_one = plt.subplots(figsize=(6, 4))
        ax_one.loglog(n_values, l1_values, marker="o")
        ax_one.set_xlabel("N")
        ax_one.set_ylabel("L1 error")
        ax_one.set_title(label)
        ax_one.grid(True, which="both", alpha=0.3)
        fig_one.tight_layout()
        fig_one.savefig(path[:-4] + "_L1.png", dpi=150)
        plt.close(fig_one)

    ax.set_xlabel("N")
    ax.set_ylabel("L1 error")
    ax.set_title("L1 convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(outdir, "L1_convergence.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
