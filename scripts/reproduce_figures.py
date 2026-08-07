#!/usr/bin/env python3
"""Reproduce the canonical quantitative and conceptual theory figures."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geometric_snowball.channels import channel_basis, davis_kahan_bound  # noqa: E402
from geometric_snowball.dynamics import simulate_standing_variation  # noqa: E402
from geometric_snowball.linalg import principal_angles  # noqa: E402
from geometric_snowball.simulations import (  # noqa: E402
    dmi_handoff_experiment,
    local_quadratic_experiment,
    ridge_load_experiment,
    selection_vs_drift_experiment,
)

plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig: plt.Figure, outdir: Path, *stems: str) -> None:
    for stem in stems:
        fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
        fig.savefig(outdir / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def local_quadratic_figure(outdir: Path) -> None:
    frame = local_quadratic_experiment()
    fig, ax = plt.subplots(figsize=(4.8, 3.55))
    ax.loglog(frame["time"], frame["exact_chordal_distance_sq"], linewidth=2.1, label=r"exact $d_C^2(t)$")
    ax.loglog(frame["time"], frame["quadratic_prediction"], "--", linewidth=1.5, label=r"$t^2\|\mathcal{K}\|_F^2$")
    ax.set(xlabel="divergence time", ylabel=r"squared channel distance, $d_C^2$", title="Deterministic local geometric snowball")
    ax.legend(frameon=False)
    fig.tight_layout(); save(fig, outdir, "01_local_quadratic_law")


def hybrid_load_figure(outdir: Path) -> None:
    frame = ridge_load_experiment()
    fig, ax = plt.subplots(figsize=(4.8, 3.55))
    ax.plot(frame["ridge_theorem_prediction"], frame["exact_covariance_load"], "o", markersize=3.8, label="matrix trace calculation")
    maximum = float(frame["exact_covariance_load"].max())
    ax.plot([0, maximum], [0, maximum], "--", linewidth=1.1, label="exact equality")
    ax.set(xlabel=r"$(\Delta g\,\Delta h/4)\sum_j\sin^2\psi_j$", ylabel=r"direct covariance log-load, $L_{\rm cov}$", title="Exact ridge calibration of geometry to fitness")
    ax.legend(frameon=False)
    fig.tight_layout(); save(fig, outdir, "03_hybrid_load_calibration")


def combined_geometry_and_fitness(outdir: Path) -> None:
    local = local_quadratic_experiment()
    ridge = ridge_load_experiment()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.15))
    axes[0].loglog(local["time"], local["exact_chordal_distance_sq"], linewidth=2, label=r"exact $d_C^2(t)$")
    axes[0].loglog(local["time"], local["quadratic_prediction"], "--", linewidth=1.4, label=r"$t^2\|\mathcal{K}\|_F^2$")
    axes[0].set(xlabel="divergence time", ylabel=r"squared channel distance, $d_C^2$", title="a  Deterministic local snowball")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(ridge["ridge_theorem_prediction"], ridge["exact_covariance_load"], "o", markersize=3.3)
    maximum = float(ridge["exact_covariance_load"].max())
    axes[1].plot([0, maximum], [0, maximum], "--", linewidth=1)
    axes[1].set(xlabel="angle-calibrated ridge prediction", ylabel="direct covariance–curvature load", title="b  Geometry requires fitness calibration")
    fig.tight_layout(); save(fig, outdir, "figure1_local_geometry_and_fitness")


def selection_vs_drift_figure(outdir: Path) -> None:
    frame = selection_vs_drift_experiment()
    fig, ax = plt.subplots(figsize=(4.8, 3.55))
    for process, group in frame.groupby("process"):
        ax.loglog(group["time"], group["mean_squared_channel_distance"], linewidth=2, label=process)
    ax.set(xlabel="time", ylabel="expected squared channel displacement", title="Deterministic selection and diffusion predict different exponents")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, outdir, "02_selection_vs_drift_scaling", "figure2_selection_vs_drift")


def standing_variation_figure(outdir: Path) -> None:
    frame = simulate_standing_variation()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.15))
    axes[0].plot(frame["time"], frame["channel_angle_deg"], linewidth=2)
    axes[0].set(xlabel="generation", ylabel="common-metric channel angle (degrees)", title="a  Sorting standing variation rotates channels")
    axes[1].plot(frame["time"], frame["gap_1"], label="lineage 1")
    axes[1].plot(frame["time"], frame["gap_2"], label="lineage 2")
    axes[1].set(xlabel="generation", ylabel="leading spectral gap", title="b  Spectral susceptibility also evolves")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, outdir, "05_standing_variation_rotation", "figure3_standing_variation")


def dmi_handoff_figure(outdir: Path) -> None:
    frame = dmi_handoff_experiment()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.15))
    for scenario, group in frame.groupby("scenario"):
        label = scenario.replace("_", " ")
        axes[0].plot(group["time"], group["effective_migration"], linewidth=2, label=label)
        axes[1].plot(group["time"], group["intrinsic_isolation"], linewidth=2, label=label)
    axes[0].set(xlabel="time", ylabel=r"effective migration, $m_{\rm eff}$", title="a  Early load reduces exchange")
    axes[1].set(xlabel="time", ylabel="intrinsic isolation", title="b  Phenomenological DMI handoff")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, outdir, "04_geometric_to_dmi_handoff", "figure4_geometric_to_dmi_handoff")


def vector_angle(u: np.ndarray, v: np.ndarray) -> float:
    uu = u / np.linalg.norm(u); vv = v / np.linalg.norm(v)
    return float(np.arccos(np.clip(abs(uu @ vv), 0.0, 1.0)))


def channel_rotation_figure(outdir: Path) -> None:
    ratios = np.linspace(1.0, 5.0, 180)
    fig, ax = plt.subplots(figsize=(4.8, 3.55))
    for correlation in [0.1, 0.4, 0.7, 0.9]:
        g = np.array([[1.0, correlation], [correlation, 1.0]])
        values, vectors = np.linalg.eigh(g)
        gmax = vectors[:, np.argmax(values)]
        rotations = []
        for ratio in ratios:
            h = np.diag([1.0, ratio])
            channel = channel_basis(g, h, k=1).phenotype_basis[:, 0]
            rotations.append(np.degrees(vector_angle(gmax, channel)))
        ax.plot(ratios, rotations, linewidth=1.8, label=rf"$r={correlation}$")
    ax.axhline(0.0, color="0.4", linestyle="--", linewidth=0.9, label=r"$\mathbf{H}\propto\mathbf{I}$")
    ax.set(xlabel=r"curvature ratio, $h_{22}/h_{11}$", ylabel=r"rotation from $\mathbf{g}_{\max}$ (degrees)", title="Two-trait response-channel rotation")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, outdir, "figS1_channel_rotation_checked")


def dk_bound_figure(outdir: Path) -> None:
    rng = np.random.default_rng(1402)
    observed: list[float] = []
    bounds: list[float] = []
    commuting_observed: list[float] = []
    commuting_bounds: list[float] = []
    for _ in range(180):
        eigenvalues = np.sort(rng.uniform(0.3, 3.0, size=4))[::-1]
        q, _ = np.linalg.qr(rng.normal(size=(4, 4)))
        s1 = q @ np.diag(eigenvalues) @ q.T
        raw = rng.normal(size=(4, 4)); perturbation = 0.5 * (raw + raw.T)
        perturbation *= 0.09 / np.linalg.norm(perturbation, ord=2)
        s2 = s1 + perturbation
        if np.linalg.eigvalsh(s2).min() <= 0:
            continue
        q1 = np.linalg.eigh(s1)[1][:, -1:]
        q2 = np.linalg.eigh(s2)[1][:, -1:]
        observed.append(float(np.sin(principal_angles(q1, q2)[0])))
        bounds.append(davis_kahan_bound(s1, s2, k=1))
        commuting = q @ np.diag([0.06, -0.03, 0.02, -0.01]) @ q.T
        s3 = s1 + commuting
        q3 = np.linalg.eigh(s3)[1][:, -1:]
        commuting_observed.append(float(np.sin(principal_angles(q1, q3)[0])))
        commuting_bounds.append(davis_kahan_bound(s1, s3, k=1))
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.15))
    limit = 1.05 * max(max(bounds), max(observed))
    axes[0].scatter(bounds, observed, s=17, alpha=0.7)
    axes[0].plot([0, limit], [0, limit], "--", color="0.3", linewidth=1)
    axes[0].set(xlim=(0, limit), ylim=(0, limit), xlabel=r"bound $\|\Delta\mathbf{S}\|_2/\delta$", ylabel=r"observed $\sin\psi$", title="a  Davis–Kahan is an upper bound")
    axes[1].scatter(commuting_bounds, commuting_observed, s=17, alpha=0.7)
    axes[1].axhline(0.0, linestyle="--", color="0.3", linewidth=1)
    axes[1].set(xlabel=r"positive bound $\|\Delta\mathbf{S}\|_2/\delta$", ylabel=r"observed $\sin\psi$", title="b  Commuting perturbations do not rotate")
    fig.tight_layout(); save(fig, outdir, "figS2_DK_bound_checked")


def flow_box(ax: plt.Axes, x: float, text: str, *, width: float) -> None:
    box = FancyBboxPatch((x - width / 2, 0.45), width, 0.25, boxstyle="round,pad=0.015", facecolor="0.96", edgecolor="0.35", zorder=2)
    ax.add_patch(box); ax.text(x, 0.575, text, ha="center", va="center", fontsize=8, zorder=3)


def conceptual_handoff_figure(outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 2.7)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    nodes = [(0.07, "divergent\nselection"), (0.23, "channel\nrotation"), (0.42, "curvature-weighted\nhybrid load"), (0.61, "reduced effective\ngene flow"), (0.80, "lineage-specific\nsubstitutions"), (0.95, "DMIs")]
    widths = [0.12, 0.12, 0.16, 0.15, 0.16, 0.07]
    for index, ((x, label), width) in enumerate(zip(nodes, widths)):
        flow_box(ax, x, label, width=width)
        if index < len(nodes) - 1:
            next_x, _ = nodes[index + 1]
            ax.add_patch(FancyArrowPatch((x + width / 2 + 0.006, 0.575), (next_x - widths[index + 1] / 2 - 0.006, 0.575), arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color="0.25", zorder=4))
    ax.text(0.36, 0.21, "early: ecological and potentially reversible", ha="center", fontsize=9)
    ax.text(0.82, 0.21, "later: combinatorial and path dependent", ha="center", fontsize=9)
    fig.tight_layout(); save(fig, outdir, "figure_geometric_snowball_concept")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=ROOT / "figures")
    args = parser.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    local_quadratic_figure(args.outdir)
    hybrid_load_figure(args.outdir)
    combined_geometry_and_fitness(args.outdir)
    selection_vs_drift_figure(args.outdir)
    standing_variation_figure(args.outdir)
    dmi_handoff_figure(args.outdir)
    channel_rotation_figure(args.outdir)
    dk_bound_figure(args.outdir)
    conceptual_handoff_figure(args.outdir)
    print(f"Figures written to {args.outdir}")


if __name__ == "__main__":
    main()
