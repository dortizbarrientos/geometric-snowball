"""Minimal channel comparison under a shared G matrix."""
import numpy as np

from geometric_snowball import channel_basis, chordal_distance_sq, compare_channels

G = np.array([[1.0, 0.25], [0.25, 0.6]])
H1 = np.array([[0.5, 0.0], [0.0, 1.4]])
H2 = np.array([[1.3, 0.0], [0.0, 0.6]])

channel_1 = channel_basis(G, H1, k=1)
angles = compare_channels(G, H1, G, H2, k=1)

print("leading S eigenvalues:", channel_1.eigenvalues)
print("principal angle (degrees):", np.degrees(angles[0]))
print("squared chordal distance:", chordal_distance_sq(angles))
