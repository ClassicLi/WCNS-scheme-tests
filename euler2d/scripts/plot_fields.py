#!/usr/bin/env python3
"""Plot a density field from an euler2d Tecplot ASCII output.

Usage: plot_fields.py input.dat output.png

The 2D Euler driver writes a point-packed zone whose variables are
``x``, ``z``, and ``rho``.  The script reads that zone and saves a single
density plot.  matplotlib is optional: when it is unavailable the script
prints a message and exits successfully without writing a figure.
"""

import argparse
import re
import sys
from pathlib import Path


def variable_names(lines):
    """Return the lower-case variable names from the VARIABLES record."""
    for line in lines:
        if line.strip().upper().startswith("VARIABLES"):
            names = re.findall(r'"([^"]*)"', line)
            return [name.strip().lower() for name in names]
    return []


def load_tecplot(path):
    """Return (values, i_count, j_count, column indices) for one zone."""
    lines = Path(path).read_text().splitlines()
    zone_index = None
    i_count = None
    j_count = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.upper().startswith("ZONE"):
            zone_index = index
            match_i = re.search(r"\bI\s*=\s*(\d+)", stripped, re.IGNORECASE)
            match_j = re.search(r"\bJ\s*=\s*(\d+)", stripped, re.IGNORECASE)
            if match_i:
                i_count = int(match_i.group(1))
            if match_j:
                j_count = int(match_j.group(1))
            break
    if zone_index is None:
        raise ValueError(f"{path}: no ZONE record found")
    if i_count is None or j_count is None:
        raise ValueError(f"{path}: ZONE needs explicit I and J extents")

    names = variable_names(lines)
    if not names:
        names = ["x", "z", "rho"]
    ncols = len(names)

    values = []
    for line in lines[zone_index + 1:]:
        tokens = line.split()
        if not tokens:
            continue
        values.extend(float(token.replace("D", "E").replace("d", "e"))
                      for token in tokens)

    expected = i_count * j_count * ncols
    if len(values) != expected:
        raise ValueError(
            f"{path}: expected {expected} values for a "
            f"{i_count}x{j_count} zone with {ncols} variables, "
            f"found {len(values)}")

    def column(name, fallback):
        """Return the column index of one variable."""
        return names.index(name) if name in names else fallback

    xcol = column("x", 0)
    zcol = column("z", 1)
    rhocol = column("rho", ncols - 1)
    return values, i_count, j_count, xcol, zcol, rhocol


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot the density field of an euler2d Tecplot output.")
    parser.add_argument("input", help="Tecplot ASCII input file")
    parser.add_argument("output", help="PNG output file")
    args = parser.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        print(f"matplotlib is not available; skipping plot ({exc})",
              file=sys.stderr)
        return 0

    try:
        values, i_count, j_count, xcol, zcol, rhocol = load_tecplot(args.input)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    # Point packing lists i fastest, so each row is one j line.
    data = np.asarray(values).reshape(j_count, i_count, -1)
    x = data[:, :, xcol]
    z = data[:, :, zcol]
    rho = data[:, :, rhocol]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    field = ax.pcolormesh(x, z, rho, shading="auto", cmap="viridis")
    fig.colorbar(field, ax=ax, label=r"$\rho$")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    ax.set_title(Path(args.input).stem)
    fig.tight_layout()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
