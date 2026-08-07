import numpy as np

from geometric_snowball.dynamics import (
    coupled_isolation,
    effective_migration,
    mean_response_identity,
    simulate_dmi_handoff,
    simulate_standing_variation,
)


def test_mean_response_including_migration() -> None:
    rng = np.random.default_rng(21)
    p = rng.uniform(0.1, 0.9, size=100)
    q = rng.uniform(0.1, 0.9, size=100)
    effects = rng.normal(scale=0.05, size=(100, 4))
    beta = rng.normal(size=4)
    direct, formula = mean_response_identity(p, q, effects, beta, migration=0.03)
    np.testing.assert_allclose(direct, formula, atol=1e-12)


def test_standing_variation_rotates_channels() -> None:
    frame = simulate_standing_variation(generations=10.0, dt=0.2, loci=80)
    assert frame["channel_angle_deg"].iloc[-1] > frame["channel_angle_deg"].iloc[0]
    assert np.isfinite(frame.to_numpy()).all()


def test_geometric_precursor_reduces_migration_and_increases_intrinsic_isolation() -> None:
    frame = simulate_dmi_handoff(duration=5.0, dt=0.02)
    final = frame.groupby("scenario").tail(1).set_index("scenario")
    assert final.loc["geometric_precursor", "effective_migration"] < final.loc["no_geometric_precursor", "effective_migration"]
    assert final.loc["geometric_precursor", "intrinsic_isolation"] > final.loc["no_geometric_precursor", "intrinsic_isolation"]


def test_combined_barrier_identity() -> None:
    le, d, eta = 0.4, 2.3, 0.2
    ie = 1.0 - np.exp(-le)
    idi = 1.0 - np.exp(-eta * d)
    np.testing.assert_allclose(coupled_isolation(le, d, eta=eta), 1.0 - (1.0 - ie) * (1.0 - idi))
    assert effective_migration(0.2, le) < 0.2
