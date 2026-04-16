#!/usr/bin/env python3
"""
Stage 3: Summary figures and tables for the Senecio channel-operator pipeline.

Reads CSV outputs from Stages 1 and 2 and produces publication-quality
summary assets. This is the thin derived layer — all analytical decisions
happen in the R stages upstream. If a number looks wrong here, fix it in the
R stage, not here.

Outputs
-------
results/figures/senecio_channel_summary.pdf
    Three-panel figure: (a) G-matrix anisotropy per population, (b) pairwise
    gmax angles between populations, (c) S channel angles with and without
    the soil-proxy cross-validation.
results/tables/channel_operator_table.csv
    Flat table combining Stage 1 G-geometry and Stage 2 channel angles for
    all population pairs.

Usage
-----
python3 scripts/03_figures/make_figures.py            # run from pipeline root
python3 make_figures.py                               # run from this directory

Both invocations work — the script resolves the pipeline root from its own
location, not from the working directory.

Author: Ortiz-Barrientos, James, Burrage
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# Matplotlib configuration — Nature Genetics house style
# ──────────────────────────────────────────────────────────────────────────────
# Vector PDF output, Arial/Helvetica, Okabe-Ito colour palette (colour-blind
# safe). We embed fonts so the PDF is portable.
plt.rcParams.update({
    "pdf.fonttype": 42,           # embed TrueType fonts, editable in Illustrator
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})

# Okabe-Ito palette, indexed by ecotype.
ECOTYPE_COLOUR = {
    "Dune":     "#0072B2",  # blue
    "Head":     "#D55E00",  # vermillion
    "Headland": "#D55E00",
    "Table":    "#009E73",  # green
    "Tableland": "#009E73",
    "Wood":     "#CC79A7",  # purple
    "Woodland": "#CC79A7",
}
POPULATION_ECOTYPE = {
    "D1": "Dune", "D3": "Dune",
    "H1": "Head", "H2": "Head",
    "T1": "Table", "T9": "Table",
    "W2": "Wood",  "W4": "Wood",
}
POPULATION_ORDER = ["D1", "D3", "H1", "H2", "T1", "T9", "W2", "W4"]


# ──────────────────────────────────────────────────────────────────────────────
# Path resolution
# ──────────────────────────────────────────────────────────────────────────────
# Resolve the pipeline root from the script's own file location, so the script
# works regardless of the caller's working directory.
SCRIPT_DIR    = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent.parent          # two levels up
STAGE1_DIR    = PIPELINE_ROOT / "data" / "processed" / "population_G_D_results"
STAGE2_DIR    = PIPELINE_ROOT / "data" / "processed" / "channel_operator_results"
FIG_DIR       = PIPELINE_ROOT / "results" / "figures"
TABLE_DIR     = PIPELINE_ROOT / "results" / "tables"


def require_file(path: Path) -> Path:
    """Fail early with a clear message if an expected upstream file is missing."""
    if not path.is_file():
        sys.stderr.write(
            f"\nMISSING upstream file: {path}\n"
            f"Run Stages 1 and 2 first:\n"
            f"  bash run_pipeline.sh\n\n"
        )
        sys.exit(2)
    return path


# ──────────────────────────────────────────────────────────────────────────────
# Panel (a): G anisotropy
# ──────────────────────────────────────────────────────────────────────────────
def panel_anisotropy(ax, leaf_summary: pd.DataFrame, veg_summary: pd.DataFrame):
    """
    Paired bar chart: for each population, plot lambda1/lambda2 for the leaf
    block and the vegetative block side by side.

    High lambda1/lambda2 = strongly anisotropic G = one direction dominates
    heritable variation.
    """
    populations = [p for p in POPULATION_ORDER if p in leaf_summary["population"].values]

    leaf_anis = (
        leaf_summary.set_index("population")
        .reindex(populations)["lambda1_over_lambda2"]
        .values
    )
    veg_anis = (
        veg_summary.set_index("population")
        .reindex(populations)["lambda1_over_lambda2"]
        .values
    )

    x     = np.arange(len(populations))
    width = 0.38

    colours = [ECOTYPE_COLOUR[POPULATION_ECOTYPE[p]] for p in populations]

    # Leaf block (solid), vegetative block (hatched)
    ax.bar(x - width / 2, leaf_anis,  width, color=colours,
           edgecolor="black", linewidth=0.5, label="Leaf PCs")
    ax.bar(x + width / 2, veg_anis,   width, color=colours, alpha=0.55,
           edgecolor="black", linewidth=0.5, hatch="///", label="Vegetative")

    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":", zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(populations)
    ax.set_ylabel(r"$\lambda_1 / \lambda_2$")
    ax.set_xlabel("Population")
    ax.set_title("(a) G-matrix anisotropy")
    ax.legend(frameon=False, loc="upper right")


# ──────────────────────────────────────────────────────────────────────────────
# Panel (b): pairwise gmax angles (Schluter-limit proxy)
# ──────────────────────────────────────────────────────────────────────────────
def panel_gmax_angles(ax, geom_leaf: pd.DataFrame):
    """
    Strip plot of pairwise gmax angles grouped by within/between ecotype.
    The small angles within-ecotype and larger angles between-ecotype are the
    Schluter-limit signal — the baseline against which S-channel angles
    should be compared.
    """
    # Classify each pair as within-ecotype or between-ecotype
    def ecotype_pair_kind(row):
        e1 = POPULATION_ECOTYPE.get(row["pop1"])
        e2 = POPULATION_ECOTYPE.get(row["pop2"])
        return "Within ecotype" if e1 == e2 else "Between ecotype"

    geom_leaf = geom_leaf.copy()
    geom_leaf["pair_kind"] = geom_leaf.apply(ecotype_pair_kind, axis=1)

    categories = ["Within ecotype", "Between ecotype"]
    xpositions = {"Within ecotype": 0, "Between ecotype": 1}

    rng = np.random.default_rng(42)
    for cat in categories:
        sub = geom_leaf[geom_leaf["pair_kind"] == cat]
        jitter = rng.uniform(-0.12, 0.12, size=len(sub))
        ax.scatter(
            np.full(len(sub), xpositions[cat]) + jitter,
            sub["gmax_angle_deg"],
            s=26, c="#444444", alpha=0.75, edgecolor="none",
        )
        # Median marker
        if len(sub) > 0:
            ax.hlines(
                sub["gmax_angle_deg"].median(),
                xpositions[cat] - 0.22, xpositions[cat] + 0.22,
                colors="black", linewidth=1.2,
            )

    ax.set_xticks(list(xpositions.values()))
    ax.set_xticklabels(categories, fontsize=6.5)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel(r"$g_{\mathrm{max}}$ angle (°)")
    ax.set_title("(b) Pairwise $g_{\\mathrm{max}}$ angles (leaf PCs)")


# ──────────────────────────────────────────────────────────────────────────────
# Panel (c): S channel angles (Route A vs Route B)
# ──────────────────────────────────────────────────────────────────────────────
def panel_S_angles(ax, s_la: pd.DataFrame, s_soil: pd.DataFrame | None):
    """
    Compare Route A (Lande-Arnold quadratic gradients) and Route B
    (soil-based rank-1 proxy) channel angles. If both routes agree on the
    ranking of population pairs, the channel angle is a robust geometric
    signal rather than a regression artefact.

    Only leaf-block Dune x Headland comparisons are plotted, since that is
    the contrast the geometric-snowball framework predicts.
    """
    # Filter to leaf-block Dune x Headland pairs
    def is_DH_leaf(df, block_col, block_value):
        mask_block = df[block_col] == block_value
        mask_eco = (
            ((df["eco1"] == "Dune") & (df["eco2"] == "Head"))
            | ((df["eco1"] == "Head") & (df["eco2"] == "Dune"))
        )
        return df[mask_block & mask_eco].copy()

    la_dh = is_DH_leaf(s_la, "block", "leaf")

    if la_dh.empty:
        ax.text(0.5, 0.5,
                "No Dune × Headland pairs\nconverged in Stage 2",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="grey")
        ax.set_title("(c) S channel angles")
        ax.set_xticks([])
        ax.set_yticks([])
        return

    # Build a pair key and sort for consistent plotting
    la_dh["pair"] = la_dh.apply(
        lambda r: f"{r['pop1']}-{r['pop2']} ({r['env']})", axis=1
    )
    la_dh = la_dh.sort_values("S_channel_angle_deg").reset_index(drop=True)

    y = np.arange(len(la_dh))
    ax.barh(y, la_dh["S_channel_angle_deg"],
            height=0.55, color="#0072B2", alpha=0.85,
            edgecolor="black", linewidth=0.5,
            label="Route A (Lande–Arnold)")

    # Overlay soil proxy if available, matched by pair
    if s_soil is not None and not s_soil.empty:
        soil_dh = is_DH_leaf(s_soil, "block", "leaf_soil")
        if not soil_dh.empty:
            soil_dh["pair"] = soil_dh.apply(
                lambda r: f"{r['pop1']}-{r['pop2']} ({r['env']})", axis=1
            )
            # Match to the same ordering as la_dh
            merged = la_dh[["pair"]].merge(
                soil_dh[["pair", "S_channel_angle_deg"]],
                on="pair", how="left",
                suffixes=("", "_soil"),
            )
            ax.scatter(
                merged["S_channel_angle_deg"], y,
                s=34, color="#D55E00", zorder=5,
                label="Route B (soil proxy)",
                edgecolor="black", linewidth=0.4,
            )

    ax.set_yticks(y)
    ax.set_yticklabels(la_dh["pair"], fontsize=6)
    ax.set_xlabel("S channel angle (°)")
    ax.set_title("(c) Channel angles, Dune × Headland")
    ax.legend(frameon=False, loc="lower right", fontsize=6)


# ──────────────────────────────────────────────────────────────────────────────
# Summary table
# ──────────────────────────────────────────────────────────────────────────────
def build_summary_table(
    leaf_summary: pd.DataFrame,
    geom_leaf: pd.DataFrame,
    s_la: pd.DataFrame,
) -> pd.DataFrame:
    """
    Flat table combining per-population G anisotropy, pairwise gmax angles,
    and pairwise S channel angles for the leaf block.
    """
    g_small = leaf_summary[
        ["population", "ecotype", "lambda1", "lambda2",
         "lambda1_over_lambda2", "leading_plane_prop"]
    ].copy()
    g_small.columns = [
        "Population", "Ecotype",
        "λ1", "λ2", "λ1/λ2", "Leading plane prop."
    ]

    # Pairwise summary: one row per pair
    pair_rows = geom_leaf[["pop1", "pop2", "gmax_angle_deg"]].copy()
    pair_rows["ecotype_pair"] = pair_rows.apply(
        lambda r: f"{POPULATION_ECOTYPE.get(r['pop1'])} × "
                  f"{POPULATION_ECOTYPE.get(r['pop2'])}",
        axis=1,
    )

    # Attach S angle (mean across environments if multiple)
    s_la_leaf = s_la[s_la["block"] == "leaf"].copy()
    if not s_la_leaf.empty:
        s_per_pair = (
            s_la_leaf.groupby(["pop1", "pop2"])["S_channel_angle_deg"]
            .mean()
            .reset_index()
        )
        pair_rows = pair_rows.merge(s_per_pair, on=["pop1", "pop2"], how="left")
    else:
        pair_rows["S_channel_angle_deg"] = np.nan

    pair_rows = pair_rows[
        ["pop1", "pop2", "ecotype_pair", "gmax_angle_deg", "S_channel_angle_deg"]
    ]
    pair_rows.columns = [
        "Pop 1", "Pop 2", "Ecotype pair",
        "gmax angle (°)", "Mean S channel angle (°)"
    ]

    # Build a single stacked report: G summary first, then pairwise
    g_small.insert(0, "Block", "Per-population G (leaf)")
    pair_rows.insert(0, "Block", "Pairwise (leaf)")

    # Align columns with pd.concat → writes NaNs where absent
    combined = pd.concat([g_small, pair_rows], ignore_index=True, sort=False)
    return combined


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main() -> int:
    # Verify upstream stages ran
    leaf_summary_f = require_file(STAGE1_DIR / "leaf_PCs_G_summary.csv")
    veg_summary_f  = require_file(STAGE1_DIR / "vegetative_G_summary.csv")
    geom_leaf_f    = require_file(STAGE1_DIR / "pairwise_G_geometry_leafPC.csv")
    s_la_f         = require_file(STAGE2_DIR / "S_channel_angles_LA.csv")

    # Optional: soil-based angles
    s_soil_f = STAGE2_DIR / "S_channel_angles_soil.csv"
    s_soil_df: pd.DataFrame | None = (
        pd.read_csv(s_soil_f) if s_soil_f.is_file() else None
    )

    # Load
    leaf_summary = pd.read_csv(leaf_summary_f)
    veg_summary  = pd.read_csv(veg_summary_f)
    geom_leaf    = pd.read_csv(geom_leaf_f)
    s_la         = pd.read_csv(s_la_f)

    # Ensure output directories
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    # Build figure
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.8),
                             gridspec_kw={"width_ratios": [1.0, 0.65, 1.15]})
    panel_anisotropy(axes[0], leaf_summary, veg_summary)
    panel_gmax_angles(axes[1], geom_leaf)
    panel_S_angles(axes[2], s_la, s_soil_df)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig_out = FIG_DIR / "senecio_channel_summary.pdf"
    fig.savefig(fig_out, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote figure: {fig_out}")

    # Build table
    table = build_summary_table(leaf_summary, geom_leaf, s_la)
    table_out = TABLE_DIR / "channel_operator_table.csv"
    table.to_csv(table_out, index=False)
    print(f"Wrote table : {table_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
