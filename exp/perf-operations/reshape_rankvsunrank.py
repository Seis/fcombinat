#!/usr/bin/env python3
"""
reshape_rankvsunrank.py — Pivot rank+unrank fcombinat CSVs into a .dat for pgfplots.

Usage:
    python3 reshape_rankvsunrank.py <results_dir>

Output (stdout): tab-separated table
Rows: structures. Columns: n10_rank, n10_unrank, n20_rank, n20_unrank, ... (interleaved).
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

OPERATIONS = ["rank", "unrank"]
N_VALUES   = [10, 20, 50, 100, 200]


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
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <results_dir>", file=sys.stderr)
        sys.exit(1)

    results_dir = sys.argv[1]

    cols   = [(n, op) for n in N_VALUES for op in OPERATIONS]
    header = "\t".join(f"n{n}_{op}" for n, op in cols)
    print(f"label\t{header}")

    for struct_id, display in STRUCTURES:
        data = {}
        for op in OPERATIONS:
            path = os.path.join(results_dir, f"fcombinat_{struct_id}_{op}.csv")
            data[op] = read_csv(path)

        values = [fmt(data[op].get(n, float("nan"))) for n, op in cols]
        print(f"{display}\t" + "\t".join(values))


if __name__ == "__main__":
    main()
