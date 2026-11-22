#!/usr/bin/env python3
"""Plot beta/omega/tau/h_noscale diagnostics from the small Tecplot samples."""

import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tecplot(path):
    with open(path) as handle:
        lines = handle.readlines()

    n_points = None
    start = None
    for index, line in enumerate(lines):
        if line.strip().startswith("ZONE"):
            match = re.search(r"I\s*=\s*(\d+)", line)
            if match:
                n_points = int(match.group(1))
                start = index + 1
                break

    if n_points is None or start is None:
        raise ValueError("no ZONE I= in " + path)

    values = []
    for line in lines[start:]:
        values.extend(float(token) for token in line.split())

    n_columns = 20
    if len(values) < n_points * n_columns:
        raise ValueError("not enough data in " + path)

    return [values[k * n_points:(k + 1) * n_points] for k in range(n_columns)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(root, "results", "diagnostics")
    os.makedirs(outdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(root, "results", "samples", "*.dat")))
    if not files:
        raise SystemExit("no sample .dat files in results/samples")

    for path in files:
        columns = read_tecplot(path)
        weight_max = max(abs(value) for k in range(3) for value in columns[8 + k])
        if weight_max < 1.0e-30:
            print("skip", path, "(no nonlinear weights)")
            continue

        x_values = columns[1]
        tau_values = columns[11]
        h_values = columns[17]
        tau_max = max(abs(value) for value in tau_values)
        h_max = max(abs(value) for value in h_values)

        panels = ["beta", "omega"]
        if tau_max > 1.0e-30:
            panels.append("tau")
        if h_max > 1.0e-30:
            panels.append("h_noscale")

        fig, axes = plt.subplots(len(panels), 1, sharex=True,
                                 figsize=(8, 2.2 * len(panels)), squeeze=False)
        ax = axes[:, 0]

        for k in range(3):
            beta = [max(value, 1.0e-300) for value in columns[5 + k]]
            ax[0].semilogy(x_values, beta, label=f"beta_{k + 1}")
            ax[1].plot(x_values, columns[8 + k], label=f"omega_{k + 1}")

        ax[0].set_ylabel("beta (log)")
        ax[1].set_ylabel("omega")

        panel_index = 2
        if "tau" in panels:
            ax[panel_index].semilogy(x_values, [max(value, 1.0e-300) for value in tau_values],
                                     color="k", label="tau")
            ax[panel_index].set_ylabel("tau (log)")
            ax[panel_index].legend(fontsize=7)
            panel_index += 1
        if "h_noscale" in panels:
            ax[panel_index].plot(x_values, h_values, color="k", label="h_noscale")
            ax[panel_index].set_ylabel("h_noscale")
            ax[panel_index].legend(fontsize=7)

        ax[0].legend(fontsize=7)
        ax[1].legend(fontsize=7)
        ax[-1].set_xlabel("x")
        fig.suptitle(os.path.basename(path))
        fig.tight_layout()

        out = os.path.join(outdir, os.path.basename(path)[:-4] + "_beta_omega.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
