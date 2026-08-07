"""Public API for the corrected geometric snowball framework."""
from .channels import (
    ChannelBasis,
    channel_basis,
    channel_operator,
    channel_operator_derivative,
    compare_channels,
    davis_kahan_bound,
    local_chordal_coefficient,
    local_rotation_matrix,
    sqrt_derivative,
)
from .hybrid_fitness import (
    bounded_isolation_from_log_load,
    covariance_cross_load,
    expected_negative_log_fitness,
    gaussian_mean_fitness,
    ridge_covariance_load,
)
from .linalg import (
    chordal_distance_sq,
    geodesic_distance,
    metric_principal_angles,
    principal_angles,
    volume_overlap_deficit,
)

__all__ = [
    "ChannelBasis",
    "bounded_isolation_from_log_load",
    "channel_basis",
    "channel_operator",
    "channel_operator_derivative",
    "chordal_distance_sq",
    "compare_channels",
    "covariance_cross_load",
    "davis_kahan_bound",
    "expected_negative_log_fitness",
    "gaussian_mean_fitness",
    "geodesic_distance",
    "local_chordal_coefficient",
    "local_rotation_matrix",
    "metric_principal_angles",
    "principal_angles",
    "ridge_covariance_load",
    "sqrt_derivative",
    "volume_overlap_deficit",
]

__version__ = "0.2.0"
