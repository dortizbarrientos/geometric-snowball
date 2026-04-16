#!/usr/bin/env python3
"""
01a_parse_selection_tables.py
==============================

PURPOSE
-------
The file `Combined_Selection_Tables.txt` supplied by Mallard et al. contains
THREE stacked tables with different column structures. This script parses
all three and writes tidy CSV outputs.

SOURCE TABLE LAYOUT
-------------------
Table 1 (lines 1-7, 6 columns): Standardized selection gradient (beta)
    Columns: trait, Post.mode, lower_83, upper_83, lower_95, upper_95
    Rows: SF, SB, FS, FB, BS, BF, Size
    Biological meaning:
        Post.mode = posterior mode of the standardized genetic selection
        gradient beta = d(log W) / d(z), measured in the ancestral A6140
        population under high-salt conditions.

Table 2 (lines 8-15, 17 columns): Selection gradients for evolved populations
    The header row labels blocks of CIs for GA150, GA250, GA450, and a
    "Beta G" block. The final column is the posterior mode of "Beta G",
    which differs from Table 1's beta in that it is computed from the
    genetic regression of fitness on breeding values, rather than from
    the phenotypic regression.

Table 3 (lines 16-31, 7 columns): Selection differentials (s = cov(trait, w))
    Columns: trait, environment, Post.mode, lower_83, upper_83, lower_95, upper_95
    Rows: 7 traits times 2 environments (High Salt, Low Salt), plus the
    fitness intercept ("w") for High Salt.
    Biological meaning: these are raw selection differentials
    s_i = cov(trait_i, fitness), measured directly without any multivariate
    partitioning. They are related to beta by s = P beta where P is the
    phenotypic covariance matrix.

WHICH TABLE DO WE USE FOR THE CHANNEL-OPERATOR COMPUTATION?
-----------------------------------------------------------
We use Table 1 (standardized phenotypic selection gradient) because Lande's
equation dz_bar/dt = G beta uses this phenotypic beta. The rank-1 surrogate
for fitness curvature is then
    H_tilde = (beta beta^T) / ||beta||^2 + epsilon I
which preserves the direction of selection while keeping H_tilde positive
definite so that the symmetric operator S_tilde = G^{1/2} H_tilde G^{1/2}
is well defined.

OUTPUTS
-------
data/processed/beta_ancestor_highsalt_standardized.csv
    Table 1 reformatted: one row per trait, columns for post.mode and CIs.
data/processed/beta_G_evolved_populations.csv
    Table 2 reformatted: one row per trait, columns for each population's
    CIs plus the "Beta G" posterior mode.
data/processed/selection_differentials.csv
    Table 3 reformatted: one row per (trait, environment).
"""

import os
import csv
from pathlib import Path

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                        # celegans_pipeline/
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE = RAW / "Combined_Selection_Tables.txt"

# ------------------------------------------------------------
# TRAIT ORDER (consistent across all tables and scripts)
# ------------------------------------------------------------
TRAITS = ["SF", "SB", "FS", "FB", "BS", "BF", "Size"]

# ------------------------------------------------------------
# READ FILE AND LOCATE TABLE BOUNDARIES
# ------------------------------------------------------------
with open(SOURCE, "r") as f:
    lines = [ln.rstrip("\n") for ln in f]

# Table 1: lines 0-7 (header + 7 trait rows).
# Table 2: lines 8-15 (header + 7 trait rows, 17 columns).
# Table 3: lines 16-end (header + 14 or 15 trait-environment rows, 7 columns).
#
# We detect boundaries by column count rather than hard-coded line numbers
# so that the script remains robust if the source is lightly edited.

def col_count(ln):
    return ln.count("\t") + 1

table_ranges = []
start = 0
current_cols = col_count(lines[0])
for i, ln in enumerate(lines[1:], start=1):
    c = col_count(ln)
    if c != current_cols:
        table_ranges.append((start, i, current_cols))
        start = i
        current_cols = c
table_ranges.append((start, len(lines), current_cols))

print(f"Detected {len(table_ranges)} tables in {SOURCE.name}:")
for s, e, c in table_ranges:
    print(f"  lines {s}-{e-1}, {c} columns, {e-s} rows")

# ------------------------------------------------------------
# TABLE 1 — standardized selection gradient (ancestor, high salt)
# ------------------------------------------------------------
s, e, _ = table_ranges[0]
t1_header = lines[s].split("\t")
t1_rows = [ln.split("\t") for ln in lines[s + 1:e]]

out1 = OUT / "beta_ancestor_highsalt_standardized.csv"
with open(out1, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["trait", "post_mode", "lower_83", "upper_83", "lower_95", "upper_95"])
    for row in t1_rows:
        trait, vals = row[0], [float(x) for x in row[1:]]
        w.writerow([trait] + vals)
print(f"Wrote {out1}")

# ------------------------------------------------------------
# TABLE 2 — evolved populations plus "Beta G"
# ------------------------------------------------------------
s, e, _ = table_ranges[1]
t2_header = lines[s].split("\t")
t2_rows = [ln.split("\t") for ln in lines[s + 1:e]]

# Column blocks: 4 per pop (83%, 95% lower/upper), 4 for Beta G CIs, 1 for Beta G mode.
# Order in file: GA150, GA250, GA450, Beta G (from the header).
# There are 7 rows (7 traits), in the SAME order as Table 1 (SF, SB, FS, FB, BS, BF, Size).
out2 = OUT / "beta_G_evolved_populations.csv"
pop_blocks = {"GA150": range(0, 4), "GA250": range(4, 8), "GA450": range(8, 12)}
betaG_ci_block = range(12, 16)
betaG_mode_col = 16

with open(out2, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "trait",
        "GA150_lower_83", "GA150_upper_83", "GA150_lower_95", "GA150_upper_95",
        "GA250_lower_83", "GA250_upper_83", "GA250_lower_95", "GA250_upper_95",
        "GA450_lower_83", "GA450_upper_83", "GA450_lower_95", "GA450_upper_95",
        "BetaG_lower_83", "BetaG_upper_83", "BetaG_lower_95", "BetaG_upper_95",
        "BetaG_post_mode",
    ])
    for i, row in enumerate(t2_rows):
        trait = TRAITS[i]
        vals = [float(x) for x in row]
        w.writerow([trait] + vals)
print(f"Wrote {out2}")

# ------------------------------------------------------------
# TABLE 3 — selection differentials by environment
# ------------------------------------------------------------
s, e, _ = table_ranges[2]
t3_header = lines[s].split("\t")
t3_rows = [ln.split("\t") for ln in lines[s + 1:e]]

out3 = OUT / "selection_differentials.csv"
with open(out3, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["trait", "environment", "post_mode", "lower_83", "upper_83", "lower_95", "upper_95"])
    for row in t3_rows:
        trait = row[0].replace("*w", "")           # "SF*w" -> "SF"
        environment = row[1].replace(" ", "")       # "High Salt" -> "HighSalt"
        vals = [float(x) for x in row[2:]]
        w.writerow([trait, environment] + vals)
print(f"Wrote {out3}")

# ------------------------------------------------------------
# SANITY CHECK: print the beta we will use downstream
# ------------------------------------------------------------
print("\nBeta (standardized selection gradient, ancestor, high salt):")
import pandas as pd
df = pd.read_csv(out1)
df["trait"] = pd.Categorical(df["trait"], categories=TRAITS, ordered=True)
df = df.sort_values("trait").reset_index(drop=True)
print(df[["trait", "post_mode", "lower_95", "upper_95"]].to_string(index=False))
import numpy as np
beta = df["post_mode"].values
print(f"\n||beta|| = {np.linalg.norm(beta):.4f}")
print(f"beta_normalized:\n{beta / np.linalg.norm(beta)}")
