"""Show that geometry and hybrid fitness are separate calculations."""
import numpy as np

from geometric_snowball import compare_channels, covariance_cross_load

G = np.array([[1.0, 0.25], [0.25, 0.6]])
H1 = np.array([[0.5, 0.0], [0.0, 1.4]])
H2 = np.array([[1.3, 0.0], [0.0, 0.6]])

angles = compare_channels(G, H1, G, H2, k=1)
load = covariance_cross_load(G, G, H1, H2)

print("channel angle (degrees):", np.degrees(angles[0]))
print("covariance-cross hybrid load:", load)
print("The angle is nonzero, but equal parental G matrices make this load exactly zero.")
