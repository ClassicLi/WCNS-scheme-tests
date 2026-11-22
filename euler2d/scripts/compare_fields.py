#!/usr/bin/env python3
"""Compare the density field of two same-grid Tecplot ASCII files.

Usage:
    python3 compare_fields.py ref.dat new.dat

Both files must contain a single point-packed ZONE with the same I/J extents.
The density variable is located by name ("rho", "den", or "density") in the
VARIABLES header; if it is not named, the third variable is used.  The report
prints max_abs_diff, max_rel_diff, and nan_count, where max_rel_diff uses the
reference density as the denominator (a true zero reference falls back to
1.0).  The script exits nonzero on a parse/grid error or when either field
contains a non-finite value.
"""

import math
import re
import sys

NUMBER_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
DENSITY_NAMES = {"rho", "den", "density"}


def fail(message):
    print("error: " + message, file=sys.stderr)
    raise SystemExit(1)


def parse_tecplot(path):
    """Return (variables, i_extent, j_extent, values)."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as stream:
            lines = stream.readlines()
    except OSError as exc:
        fail("cannot read %s: %s" % (path, exc))

    variables = None
    i_extent = None
    j_extent = None
    data_start = None

    for index, line in enumerate(lines):
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("VARIABLES"):
            variables = re.findall(r'"([^"]*)"', stripped)
        elif upper.startswith("ZONE"):
            match_i = re.search(r"(?<![A-Za-z])I\s*=\s*(\d+)", stripped)
            match_j = re.search(r"(?<![A-Za-z])J\s*=\s*(\d+)", stripped)
            if not match_i or not match_j:
                fail("ZONE line in %s has no I/J extents" % path)
            i_extent = int(match_i.group(1))
            j_extent = int(match_j.group(1))
            data_start = index + 1
            break

    if data_start is None:
        fail("%s has no ZONE header" % path)
    if not variables:
        fail("%s has no VARIABLES header" % path)
    if i_extent <= 0 or j_extent <= 0:
        fail("%s has a non-positive I/J extent" % path)

    values = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.upper().startswith("ZONE"):
            fail("%s contains more than one ZONE" % path)
        tokens = NUMBER_RE.findall(stripped)
        for token in tokens:
            values.append(float(token.replace("D", "E").replace("d", "e")))

    nvar = len(variables)
    expected = i_extent * j_extent * nvar
    if len(values) != expected:
        fail(
            "%s has %d data values, expected %d (%d points x %d variables)"
            % (path, len(values), expected, i_extent * j_extent, nvar)
        )

    return variables, i_extent, j_extent, values


def density_index(variables):
    for index, name in enumerate(variables):
        if name.strip().lower() in DENSITY_NAMES:
            return index
    if len(variables) >= 3:
        return 2
    fail("cannot locate a density variable in VARIABLES=%s" % variables)


def extract_density(variables, values, index):
    return values[index::len(variables)]


def main(argv):
    if len(argv) != 3:
        print("usage: compare_fields.py ref.dat new.dat", file=sys.stderr)
        return 2

    ref_path, new_path = argv[1], argv[2]
    ref_vars, ref_i, ref_j, ref_values = parse_tecplot(ref_path)
    new_vars, new_i, new_j, new_values = parse_tecplot(new_path)

    if (ref_i, ref_j) != (new_i, new_j):
        fail(
            "grid mismatch: reference has I=%d J=%d, new has I=%d J=%d"
            % (ref_i, ref_j, new_i, new_j)
        )

    ref_index = density_index(ref_vars)
    new_index = density_index(new_vars)
    ref_rho = extract_density(ref_vars, ref_values, ref_index)
    new_rho = extract_density(new_vars, new_values, new_index)

    nan_count = 0
    max_abs_diff = 0.0
    max_rel_diff = 0.0
    finite_pairs = 0

    for ref_value, new_value in zip(ref_rho, new_rho):
        if not (math.isfinite(ref_value) and math.isfinite(new_value)):
            nan_count += 1
            continue
        diff = abs(ref_value - new_value)
        denominator = abs(ref_value) if ref_value != 0.0 else 1.0
        rel = diff / denominator
        if diff > max_abs_diff:
            max_abs_diff = diff
        if rel > max_rel_diff:
            max_rel_diff = rel
        finite_pairs += 1

    print("ref_file=" + ref_path)
    print("new_file=" + new_path)
    print("ref_variables=" + ",".join(ref_vars))
    print("new_variables=" + ",".join(new_vars))
    print("points=%d" % (ref_i * ref_j))
    print("finite_pairs=%d" % finite_pairs)
    print("max_abs_diff=%.16e" % max_abs_diff)
    print("max_rel_diff=%.16e" % max_rel_diff)
    print("nan_count=%d" % nan_count)

    if nan_count > 0:
        print(
            "error: %d point(s) contain a non-finite density value" % nan_count,
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
