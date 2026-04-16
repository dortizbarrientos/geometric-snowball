#!/usr/bin/env python3
"""
02a_posterior_analysis.py
==========================

PURPOSE
-------
Main analysis. For each of the four C. elegans populations and each of the
two environments, propagate the full MCMC posterior of G through the
channel-operator framework. For every posterior draw we compute:

    - lambda_1, lambda_2                : two largest eigenvalues of G
    - ratio = lambda_1 / lambda_2       : channel anisotropy
    - delta = lambda_1 - lambda_2       : raw spectral gap
    - rotation_deg                      : angle between g_max and the
                                          leading eigenvector of
                                          S_tilde = G^{1/2} H_tilde G^{1/2}
                                          (pulled back to phenotype space)

INPUT
-----
    data/processed/G_chains_high_salt_long.csv      (from 01c)
    data/processed/G_chains_low_salt_long.csv       (from 01c)
    data/processed/beta_ancestor_highsalt_standardized.csv   (from 01a)

ASSUMPTIONS
-----------
- The same beta (ancestral A6140, high salt) is used for all populations
  and both environments. This is conservative: it assumes no divergence
  in fitness curvature among evolved replicates, so any channel rotation
  we detect must arise from divergence in G alone.

- The rank-1 surrogate H_tilde captures only the DIRECTION of selection,
  not its magnitude. Magnitude drops out of the rotation calculation
  because rotation is invariant to positive scaling of H_tilde.

OUTPUTS
-------
    data/processed/celegans_posterior_perdraw.csv
        Per-draw quantities for every (environment, population, draw).
        Columns: environment, population, draw, lambda_1, lambda_2,
                 ratio, delta, rotation_deg, was_PD

    results/tables/celegans_posterior_summary.csv
        Posterior median + 95% CI per (environment, population, quantity).
        Columns: environment, population, n_measurement, n_draws_used,
                 quantity, median, lo95, hi95
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from geometry import compute_channel_geometry, posterior_summary  # noqa: E402

# ----------------------------------------------------------------------
# CONSTANTS
# ----------------------------------------------------------------------
TRAITS = ["SF", "SB", "FS", "FB", "BS", "BF", "Size"]
P = len(TRAITS)

ENV_FILES = {
    "NaCl": PROCESSED / "G_chains_high_salt_long.csv",
    "NGM":  PROCESSED / "G_chains_low_salt_long.csv",
}

# Sample sizes (number of inbred lines assayed) from Mallard et al.
# Used only as metadata in the summary table; not used in calculation.
N_MEASUREMENT = {"A6140": 362, "GA150": 118, "GA250": 121, "GA450": 83}


# ----------------------------------------------------------------------
# LOAD BETA (standardized phenotypic selection gradient, A6140 high salt)
# ----------------------------------------------------------------------
beta_df = pd.read_csv(PROCESSED / "beta_ancestor_highsalt_standardized.csv")
beta_df["trait"] = pd.Categorical(beta_df["trait"], categories=TRAITS, ordered=True)
beta_df = beta_df.sort_values("trait").reset_index(drop=True)
beta = beta_df["post_mode"].to_numpy()
print("Loaded beta (ancestor, high salt, standardized):")
for t, b in zip(TRAITS, beta):
    print(f"  {t:6s} = {b:+.4f}")
print(f"  ||beta|| = {np.linalg.norm(beta):.4f}\n")


# ----------------------------------------------------------------------
# FUNCTION TO RESHAPE ONE DRAW'S 49 ROWS INTO A 7X7 MATRIX
# ----------------------------------------------------------------------
def rows_to_matrix(group):
    """
    group is a DataFrame containing the 49 rows for one (population, draw).
    Return the 7x7 G matrix.
    """
    M = np.zeros((P, P))
    for _, r in group.iterrows():
        M[int(r.row_idx), int(r.col_idx)] = r.value
    return (M + M.T) / 2     # symmetrize defensively


# ----------------------------------------------------------------------
# MAIN LOOP: process each environment, each population, every draw
# ----------------------------------------------------------------------
per_draw_rows = []
summary_rows = []

for env, chain_file in ENV_FILES.items():
    if not chain_file.exists():
        print(f"[skip] {env}: {chain_file} not found")
        continue

    print(f"Processing environment {env}...")
    print(f"  Loading {chain_file.name} ...")
    df = pd.read_csv(chain_file)
    print(f"  Loaded {len(df):,} rows")

    for pop in ["A6140", "GA150", "GA250", "GA450"]:
        pop_df = df[df.population == pop]
        n_draws = pop_df["draw"].nunique()
        print(f"    {pop}: {n_draws} draws")

        # For speed, pivot once to a (n_draws x 49) array
        pivot = (
            pop_df
            .pivot(index="draw", columns=["row_idx", "col_idx"], values="value")
            .sort_index()
        )
        # Column order is sorted by (row_idx, col_idx); index is 0..48 in C order
        arr = pivot.to_numpy()

        # Per-draw geometry
        lambda_1 = np.empty(n_draws); lambda_2 = np.empty(n_draws)
        ratio = np.empty(n_draws);    delta = np.empty(n_draws)
        rotation = np.empty(n_draws); was_PD = np.empty(n_draws, dtype=bool)

        for d in range(n_draws):
            G = arr[d].reshape(P, P)
            G = (G + G.T) / 2
            res = compute_channel_geometry(G, beta)
            lambda_1[d] = res["lambda_1"]
            lambda_2[d] = res["lambda_2"]
            ratio[d]    = res["ratio"]
            delta[d]    = res["delta"]
            rotation[d] = res["rotation_deg"]
            was_PD[d]   = res["was_PD"]

        draws_array = pivot.index.to_numpy()
        for d_i, draw_id in enumerate(draws_array):
            per_draw_rows.append(dict(
                environment=env, population=pop, draw=int(draw_id),
                lambda_1=lambda_1[d_i], lambda_2=lambda_2[d_i],
                ratio=ratio[d_i], delta=delta[d_i],
                rotation_deg=rotation[d_i], was_PD=bool(was_PD[d_i]),
            ))

        # Posterior summaries
        for qty_name, arr_vals in [
            ("lambda_1", lambda_1), ("lambda_2", lambda_2),
            ("ratio", ratio), ("delta", delta),
            ("rotation_deg", rotation),
        ]:
            s = posterior_summary(arr_vals)
            summary_rows.append(dict(
                environment=env, population=pop,
                n_measurement=N_MEASUREMENT[pop],
                n_draws_used=s["n_used"],
                quantity=qty_name,
                median=s["median"], lo95=s["lo95"], hi95=s["hi95"],
            ))

        # Console report for the ratio and rotation (readable at a glance)
        s_ratio = posterior_summary(ratio)
        s_rot = posterior_summary(rotation)
        print(f"        lambda_1/lambda_2 = {s_ratio['median']:.2f}  "
              f"[{s_ratio['lo95']:.2f}, {s_ratio['hi95']:.2f}]")
        print(f"        rotation (deg)    = {s_rot['median']:.1f}  "
              f"[{s_rot['lo95']:.1f}, {s_rot['hi95']:.1f}]")


# ----------------------------------------------------------------------
# WRITE OUTPUTS
# ----------------------------------------------------------------------
per_draw_df = pd.DataFrame(per_draw_rows)
per_draw_out = PROCESSED / "celegans_posterior_perdraw.csv"
per_draw_df.to_csv(per_draw_out, index=False)
print(f"\nWrote {per_draw_out}  ({len(per_draw_df):,} rows)")

summary_df = pd.DataFrame(summary_rows)
summary_out = TABLES / "celegans_posterior_summary.csv"
summary_df.to_csv(summary_out, index=False)
print(f"Wrote {summary_out}  ({len(summary_df)} rows)")

# Also write a pretty-printed summary for humans
print("\nFinal summary:")
display = summary_df.query("quantity in ['ratio','rotation_deg']").copy()
display["ci"] = display.apply(
    lambda r: f"[{r.lo95:.2f}, {r.hi95:.2f}]"
              if r.quantity == "ratio"
              else f"[{r.lo95:.1f}, {r.hi95:.1f}]", axis=1)
print(display[["environment", "population", "quantity", "median", "ci"]]
      .to_string(index=False))
