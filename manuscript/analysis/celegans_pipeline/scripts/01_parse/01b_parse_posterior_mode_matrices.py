#!/usr/bin/env python3
"""
01b_parse_posterior_mode_matrices.py
=====================================

PURPOSE
-------
Parse the nine 7x7 posterior-mode G-matrix text files (one for each of the
four populations x two environments, plus a mutation-accumulation baseline).
Each cell is formatted as "value  [lower_95;upper_95]".

These are the point-estimate matrices that appear directly in Mallard et al.
(2023). They are used here only as a sanity check against the posterior
modes we compute from the MCMC chains later (script 01c) -- the chains are
the primary data.

INPUT FILES (data/raw/)
-----------------------
  A6140_Low_Salt_posterior_mode.txt
  A6140_High_Salt_posterior_mode.txt
  GA150_Low_Salt_posterior_mode.txt
  GA150_High_Salt_posterior_mode.txt
  GA250_Low_Salt_posterior_mode.txt
  GA250_High_Salt_posterior_mode.txt
  GA450_Low_Salt_posterior_mode.txt
  GA450_High_Salt_posterior_mode.txt
  MAlines_posterior_mode.txt             (mutational input baseline, NGM only)

OUTPUT
------
data/processed/G_matrices_posterior_mode_long.csv
    Tidy long-format table with all 9 * 49 = 441 cells.
    Columns: population, environment, row_trait, col_trait, row_idx, col_idx,
             post_mode, lower_95CI, upper_95CI
"""

import re
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

TRAITS = ["SF", "SB", "FS", "FB", "BS", "BF", "Size"]

FILES = {
    "A6140_Low_Salt_posterior_mode.txt":   ("A6140", "NGM"),
    "A6140_High_Salt_posterior_mode.txt":  ("A6140", "NaCl"),
    "GA150_Low_Salt_posterior_mode.txt":   ("GA150", "NGM"),
    "GA150_High_Salt_posterior_mode.txt":  ("GA150", "NaCl"),
    "GA250_Low_Salt_posterior_mode.txt":   ("GA250", "NGM"),
    "GA250_High_Salt_posterior_mode.txt":  ("GA250", "NaCl"),
    "GA450_Low_Salt_posterior_mode.txt":   ("GA450", "NGM"),
    "GA450_High_Salt_posterior_mode.txt":  ("GA450", "NaCl"),
    "MAlines_posterior_mode.txt":          ("MA",    "NGM"),
}

# Regex to capture cells formatted "value  [lower;upper]".
CELL = re.compile(r"(-?\d+\.?\d*)\s*\[(-?\d+\.?\d*);(-?\d+\.?\d*)\]")


def parse_matrix(path):
    """Parse one matrix file into a list of (row_idx, col_idx, mode, lo, hi) tuples."""
    rows = []
    with open(path) as f:
        for row_idx, line in enumerate(f):
            if not line.strip():
                continue
            for col_idx, (mode, lo, hi) in enumerate(CELL.findall(line)):
                rows.append({
                    "row_idx": row_idx,
                    "col_idx": col_idx,
                    "row_trait": TRAITS[row_idx],
                    "col_trait": TRAITS[col_idx],
                    "post_mode": float(mode),
                    "lower_95CI": float(lo),
                    "upper_95CI": float(hi),
                })
    return rows


# ------------------------------------------------------------
# PARSE ALL FILES
# ------------------------------------------------------------
all_rows = []
for filename, (pop, env) in FILES.items():
    path = RAW / filename
    cells = parse_matrix(path)
    if len(cells) != 49:
        raise RuntimeError(
            f"{filename} parsed to {len(cells)} cells; expected 49. "
            "Check that trait count matches TRAITS."
        )
    for c in cells:
        all_rows.append({"population": pop, "environment": env, **c})

# ------------------------------------------------------------
# WRITE OUTPUT
# ------------------------------------------------------------
out_csv = OUT / "G_matrices_posterior_mode_long.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "population", "environment",
        "row_trait", "col_trait", "row_idx", "col_idx",
        "post_mode", "lower_95CI", "upper_95CI",
    ])
    writer.writeheader()
    writer.writerows(all_rows)

print(f"Wrote {out_csv}")
print(f"  rows: {len(all_rows)}  (9 matrices x 49 cells)")

# ------------------------------------------------------------
# SANITY CHECK: print the A6140 high-salt diagonal and compare to the paper
# ------------------------------------------------------------
import pandas as pd
df = pd.read_csv(out_csv)
diag = df[(df.population == "A6140")
          & (df.environment == "NaCl")
          & (df.row_trait == df.col_trait)]
print("\nA6140 high-salt diagonal (variances):")
print(diag[["row_trait", "post_mode", "lower_95CI", "upper_95CI"]].to_string(index=False))
