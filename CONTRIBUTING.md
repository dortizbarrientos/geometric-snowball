# Contributing

Contributions are welcome through focused pull requests.

Before submitting a change:

1. State the mathematical or biological claim being changed.
2. Identify the assumptions under which it is valid.
3. Add or update a numerical identity test.
4. Run `pytest`.
5. Run `python scripts/run_derivation_audit.py`.
6. Do not call a geometric distance reproductive isolation unless an explicit fitness mapping is supplied.

For unequal \(\mathbf G\) matrices, declare the shared metric. For perturbation claims, distinguish bounds from realised rotation and report cross-channel forcing as well as spectral gaps. For empirical analyses, distinguish measured quantities from model-dependent latent quantities.
