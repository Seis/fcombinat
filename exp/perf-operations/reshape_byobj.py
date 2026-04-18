#!/usr/bin/env python3
"""
reshape_byobj.py — Per-operation pivot: structures on x, (n, backend) as bar series.

Usage:
    python3 reshape_byobj.py <operation> <results_dir>

Output (stdout): tab-separated table
Rows: structures (cycles, distinct_partitions, ...).
Columns: n10_fcombinat, n10_maple, n20_fcombinat, n20_maple, ... (interleaved).
Missing data (no CSV) → nan.
"""

import csv
import math
import os
import sys

STRUCTURES = [
    ("cycles",              "Cycles"),
    ("distinct_partitions", "Dist. Part."),
    ("integer_partitions",  "Int. Part."),
    ("permutations",        "Permutations"),
    ("surjections",         "Surjections"),
]

BACKENDS  = ["fcombinat", "maple"]
N_VALUES  = [10, 20, 50, 100, 200]


def read_csv(path):
    result = {}
    if not os.path.exists(path):
        return result
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                result[int(row["n"])] = float(row["cycles"])
            except (KeyError, ValueError):
                pass
    return result


def fmt(v):
    return f"{v:.0f}" if not math.isnan(v) else "nan"


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <operation> <results_dir>", file=sys.stderr)
        sys.exit(1)

    operation   = sys.argv[1]
    results_dir = sys.argv[2]

    cols   = [(n, b) for n in N_VALUES for b in BACKENDS]
    header = "\t".join(f"n{n}_{b}" for n, b in cols)
    print(f"label\t{header}")

    for struct_id, display in STRUCTURES:
        data = {}
        for backend in BACKENDS:
            path = os.path.join(results_dir, f"{backend}_{struct_id}_{operation}.csv")
            data[backend] = read_csv(path)

        values = [fmt(data[b].get(n, float("nan"))) for n, b in cols]
        print(f"{display}\t" + "\t".join(values))


if __name__ == "__main__":
    main()
