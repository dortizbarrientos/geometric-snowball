"""Command-line interface for reproducible numerical experiments."""
from __future__ import annotations

import argparse
from pathlib import Path

from .dynamics import simulate_standing_variation
from .simulations import (
    dmi_handoff_experiment,
    local_quadratic_experiment,
    ridge_load_experiment,
    selection_vs_drift_experiment,
)


def reproduce(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    experiments = {
        "local_quadratic.csv": local_quadratic_experiment(),
        "ridge_load.csv": ridge_load_experiment(),
        "selection_vs_drift.csv": selection_vs_drift_experiment(),
        "standing_variation.csv": simulate_standing_variation(),
        "dmi_handoff.csv": dmi_handoff_experiment(),
    }
    for filename, frame in experiments.items():
        frame.to_csv(outdir / filename, index=False)
    print(f"Wrote {len(experiments)} reproducible experiment tables to {outdir}")


def main() -> None:
    parser = argparse.ArgumentParser(prog="geometric-snowball")
    sub = parser.add_subparsers(dest="command", required=True)
    reproduce_parser = sub.add_parser("reproduce", help="run all canonical numerical experiments")
    reproduce_parser.add_argument("--outdir", type=Path, default=Path("build/simulations"))
    args = parser.parse_args()
    if args.command == "reproduce":
        reproduce(args.outdir)


if __name__ == "__main__":
    main()
