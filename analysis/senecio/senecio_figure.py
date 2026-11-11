#!/usr/bin/env python3
"""Plot the descriptive Senecio pair summary with scientifically cautious labels."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SUMMARY = ROOT / "data" / "derived" / "senecio_pair_summary.csv"


def _column(frame: pd.DataFrame, corrected: str, legacy: str) -> pd.Series:
    if corrected in frame:
        return frame[corrected]
    if legacy in frame:
        return frame[legacy]
    raise ValueError(f"summary has neither {corrected!r} nor {legacy!r}")


def make_figure(summary_path: Path, output: Path) -> None:
    frame = pd.read_csv(summary_path).sort_values("Pair")
    morphology = _column(frame, "morphology_distance", "Pheno_dist")
    soil = _column(frame, "soil_distance", "Env_dist")
    colours = ["C0" if clade == "eastern" else "C1" for clade in frame["Clade"]]
    x = np.arange(len(frame))

    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.45))
    axes[0].bar(x, morphology, color=colours)
    axes[0].set_xticks(x, frame["Pair"])
    axes[0].set(xlabel="Dune–Headland pair", ylabel="standardised morphology distance")
    axes[0].set_title("a  Descriptive phenotype separation", fontsize=10, loc="left")

    axes[1].bar(x, soil, color=colours)
    axes[1].set_xticks(x, frame["Pair"])
    axes[1].set(xlabel="Dune–Headland pair", ylabel="standardised soil-chemistry distance")
    axes[1].set_title("b  Descriptive habitat separation", fontsize=10, loc="left")

    axes[2].scatter(soil, morphology, c=colours, s=42)
    for _, row in frame.iterrows():
        soil_value = row["soil_distance"] if "soil_distance" in row else row["Env_dist"]
        morphology_value = row["morphology_distance"] if "morphology_distance" in row else row["Pheno_dist"]
        axes[2].annotate(str(int(row["Pair"])), (soil_value, morphology_value), xytext=(4, 3), textcoords="offset points", fontsize=8)
    r, p = pearsonr(soil, morphology)
    axes[2].set(xlabel="soil-chemistry distance", ylabel="morphology distance")
    axes[2].set_title("c  Descriptive association only", fontsize=10, loc="left")
    axes[2].text(0.04, 0.96, f"Pearson r={r:.2f}; p={p:.3f}; n={len(frame)}", transform=axes[2].transAxes, va="top", fontsize=8)
    for ax in axes:
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(9)
        ax.yaxis.label.set_size(9)
    fig.text(0.73, 0.015, "Neither axis estimates H, a channel angle, or reproductive isolation.", ha="center", fontsize=7.5)

    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    if output.suffix.lower() == ".pdf":
        fig.savefig(output.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--outfile", type=Path, default=ROOT / "analysis" / "senecio" / "results" / "senecio_pair_descriptive.pdf")
    args = parser.parse_args()
    make_figure(args.summary, args.outfile)


if __name__ == "__main__":
    main()
