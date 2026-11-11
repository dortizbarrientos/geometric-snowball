#!/usr/bin/env python3
"""Descriptive Dune–Headland distances for motivating, not testing, the theory.

The outputs are deliberately not called estimates or proxies of H, channel angle,
or reproductive isolation. Direct tests require population-specific G matrices,
selection-curvature estimates, recombinant covariance, and hybrid fitness.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]

TRAITS = [
    "TotalHeight",
    "MSD",
    "MSL_W",
    "Area",
    "Circularity",
    "Nindents.Peri",
    "IndentWidth",
    "IndentDepth",
    "SB",
]
SOIL_VARIABLES = [
    "pH Value",
    "Electrical Conductivity @ 25_C",
    "Total Organic Carbon",
    "Cation Exchange Capacity",
    "Phosphorus",
    "Total Nitrogen as N",
    "Exchangeable Calcium",
    "Exchangeable Magnesium",
]
PAIRS = [
    (1, "eastern", "DN0", "HN0", "D00", "H00"),
    (2, "eastern", "DN1", "HN1", "D01", "H01"),
    (3, "eastern", "DN3", "HN2", "D03", "H02"),
    (4, "eastern", "DN4", "HN5", "D04", "H05"),
    (5, "eastern", "DN12", "HN14", "D12", "H14"),
    (6, "southeastern", "DN23", "HS21", "D23", "H21"),
    (7, "southeastern", "DN32", "HN12", "D32", "H12"),
    (8, "southeastern", "DN14", "HN15", "D14", "H15"),
]


def _read_morphology(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path)
    missing = sorted(set(["Ecotype", "Population", *TRAITS]) - set(frame.columns))
    if missing:
        raise ValueError(f"morphology file is missing columns: {missing}")
    return frame


def _read_soil(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name="All soil data")
    missing = sorted(set(["Ecotype", "Population", *SOIL_VARIABLES]) - set(frame.columns))
    if missing:
        raise ValueError(f"soil file is missing columns: {missing}")
    return frame


def analyse(morphology_path: Path, soil_path: Path, outdir: Path) -> None:
    morphology = _read_morphology(morphology_path)
    soil = _read_soil(soil_path)
    outdir.mkdir(parents=True, exist_ok=True)

    focal = morphology[morphology["Ecotype"].isin(["D.Dune.N", "G.Head_N", "H.Head_S"])].dropna(subset=TRAITS)
    scale = focal[TRAITS].std(ddof=1).replace(0.0, np.nan)
    population_means = morphology.dropna(subset=TRAITS).groupby(["Ecotype", "Population"])[TRAITS].mean()

    def morphology_vector(population: str) -> np.ndarray:
        mask = population_means.index.get_level_values("Population") == population
        if mask.sum() != 1:
            raise ValueError(f"expected exactly one morphology population named {population}")
        return population_means.loc[mask, TRAITS].to_numpy(dtype=float)[0] / scale.to_numpy(dtype=float)

    soil_focal = soil[soil["Ecotype"].isin(["Dune", "Headland"])].copy()
    soil_numeric = soil_focal[SOIL_VARIABLES].apply(pd.to_numeric, errors="coerce")
    soil_scale = soil_numeric.std(ddof=1).replace(0.0, np.nan)
    soil_lookup = soil_focal.set_index(["Ecotype", "Population"])

    rows: list[dict[str, float | str | int]] = []
    for pair, clade, dune_morph, head_morph, dune_soil, head_soil in PAIRS:
        morphology_difference = morphology_vector(dune_morph) - morphology_vector(head_morph)
        dune_values = pd.to_numeric(soil_lookup.loc[("Dune", dune_soil), SOIL_VARIABLES], errors="coerce")
        head_values = pd.to_numeric(soil_lookup.loc[("Headland", head_soil), SOIL_VARIABLES], errors="coerce")
        soil_difference = (dune_values - head_values) / soil_scale
        rows.append(
            {
                "Pair": pair,
                "Clade": clade,
                "Dune_pop": dune_morph,
                "Head_pop": head_morph,
                "morphology_distance": float(np.linalg.norm(morphology_difference)),
                "soil_distance": float(np.linalg.norm(soil_difference.dropna())),
                "delta_height_sd": float(morphology_difference[TRAITS.index("TotalHeight")]),
                "delta_area_sd": float(morphology_difference[TRAITS.index("Area")]),
                "delta_circularity_sd": float(morphology_difference[TRAITS.index("Circularity")]),
                "delta_pH": float(dune_values["pH Value"] - head_values["pH Value"]),
                "delta_organic_carbon": float(dune_values["Total Organic Carbon"] - head_values["Total Organic Carbon"]),
                "delta_total_nitrogen": float(dune_values["Total Nitrogen as N"] - head_values["Total Nitrogen as N"]),
            }
        )
    pair_summary = pd.DataFrame(rows).sort_values("Pair")
    pair_summary.to_csv(outdir / "senecio_pair_summary.csv", index=False)

    all_population_means = morphology.dropna(subset=TRAITS).groupby(["Ecotype", "Population"])[TRAITS].mean()
    pca = PCA(n_components=2)
    scores = pca.fit_transform(StandardScaler().fit_transform(all_population_means[TRAITS]))
    pca_scores = pd.DataFrame(scores, columns=["PC1", "PC2"], index=all_population_means.index).reset_index()
    pca_scores.to_csv(outdir / "senecio_morphology_pca_scores.csv", index=False)
    loadings = pd.DataFrame(pca.components_.T, index=TRAITS, columns=["PC1_loading", "PC2_loading"])
    loadings.to_csv(outdir / "senecio_morphology_pca_loadings.csv")
    variance = pd.DataFrame({"component": ["PC1", "PC2"], "variance_explained": pca.explained_variance_ratio_[:2]})
    variance.to_csv(outdir / "senecio_morphology_pca_variance.csv", index=False)
    print(pair_summary.to_string(index=False))
    print(f"\nWrote descriptive outputs to {outdir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--morphology", type=Path, required=True)
    parser.add_argument("--soil", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=ROOT / "analysis" / "senecio" / "results")
    args = parser.parse_args()
    analyse(args.morphology, args.soil, args.outdir)


if __name__ == "__main__":
    main()
