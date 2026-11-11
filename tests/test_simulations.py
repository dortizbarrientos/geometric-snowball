import numpy as np

from geometric_snowball.simulations import (
    local_quadratic_experiment,
    ridge_load_experiment,
    selection_vs_drift_experiment,
)


def test_local_experiment_converges_to_quadratic_prediction() -> None:
    frame = local_quadratic_experiment()
    ratio = frame["exact_chordal_distance_sq"].iloc[0] / frame["quadratic_prediction"].iloc[0]
    np.testing.assert_allclose(ratio, 1.0, rtol=2e-3)


def test_ridge_experiment_matches_exactly() -> None:
    frame = ridge_load_experiment()
    np.testing.assert_allclose(frame["exact_covariance_load"], frame["ridge_theorem_prediction"], atol=1e-12)


def test_selection_and_drift_exponents() -> None:
    frame = selection_vs_drift_experiment(replicates=10_000)
    slopes = {}
    for process, group in frame.groupby("process"):
        slopes[process] = np.polyfit(np.log(group["time"]), np.log(group["mean_squared_channel_distance"]), 1)[0]
    assert 1.98 < slopes["deterministic selection"] < 2.02
    assert 0.96 < slopes["unbiased diffusion"] < 1.04
