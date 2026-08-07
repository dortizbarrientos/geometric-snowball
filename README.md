# The Geometric Snowball of Speciation

A corrected mathematical and computational framework for asking how divergent natural selection can generate an early, accelerating component of reproductive isolation and thereby facilitate the later accumulation of intrinsic Dobzhansky–Muller incompatibilities.

**Daniel Ortiz-Barrientos, Maddie E. James, and Kevin Burrage**

> **Status:** active theory development. The mathematics and software are tested. The proposed biological mechanism is a set of explicit, falsifiable hypotheses, not yet an established empirical conclusion.

## Central argument

Local multivariate response around an optimum depends jointly on additive genetic covariance, \(\mathbf G\), and the curvature of minus log fitness, \(\mathbf H\). Although the local dynamics involve \(\mathbf G\mathbf H\), that product is generally non-symmetric. The channel geometry is therefore defined by the symmetric similar operator

\[
\mathbf S=\mathbf G^{1/2}\mathbf H\mathbf G^{1/2}.
\]

The leading invariant subspace of \(\mathbf S\) is the local additive-response channel in genetically whitened coordinates. If two lineages begin from the same channel and experience different smooth deterministic perturbations,

\[
d_C^2(t)=t^2\lVert\mathcal K\rVert_F^2+O(t^3),
\]

where \(\mathcal K\) contains the between-lineage cross-channel forcing divided by the relevant spectral gaps. Channel angles grow to first order, so squared subspace discrepancy grows to second order: the local geometric snowball.

A channel angle is not reproductive isolation by definition. For Gaussian phenotypes and quadratic selection, the covariance component of average hybrid log-load is

\[
L_{\mathrm{cov}}
=-\frac18\operatorname{tr}\!\left[(\mathbf H_2-\mathbf H_1)(\mathbf G_2-\mathbf G_1)\right]
+\frac14\operatorname{tr}\!\left[(\mathbf H_1+\mathbf H_2)\mathbf R\right].
\]

Under the balanced channel-ridge model, this fitness bridge becomes exact:

\[
L_{\mathrm{cov}}
=\frac{\Delta g\,\Delta h}{4}
\sum_{j=1}^k\sin^2\psi_j.
\]

The proposed sequence is therefore

\[
\text{divergent selection}
\rightarrow \text{channel rotation}
\rightarrow \text{curvature-weighted hybrid load}
\rightarrow \text{reduced effective gene flow}
\rightarrow \text{DMI accumulation}.
\]

The last arrow is implemented as an explicitly phenomenological extension. It is not a consequence of the channel operator alone.

## What this rebuild corrects

The former repository treated

\[
1-\prod_j\cos^2\psi_j
\]

as reproductive isolation or hybrid fitness cost. In v2 it is retained only as a **geometric volume-overlap deficit**. Fitness is calculated separately.

The rebuild also makes four structural corrections:

1. Channel geometry uses \(\mathbf S=\mathbf G^{1/2}\mathbf H\mathbf G^{1/2}\), not Euclidean eigenspace geometry of \(\mathbf G\mathbf H\).
2. When parental \(\mathbf G\) matrices differ, phenotype-space channels are compared under one declared common metric.
3. A small spectral gap is a susceptibility, not a cause. Rotation additionally requires cross-channel forcing.
4. Soil and morphology distances in the exploratory *Senecio* analysis are descriptive distances, not estimates of \(\Delta\mathbf H\), channel angle, or reproductive isolation.

The complete old-to-new map is in [`docs/MIGRATION_FROM_V1.md`](docs/MIGRATION_FROM_V1.md).

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[all]"

pytest
python scripts/run_derivation_audit.py --outdir build/audit
geometric-snowball reproduce --outdir build/simulations
python scripts/reproduce_figures.py --outdir build/figures
```

A minimal calculation is:

```python
import numpy as np
from geometric_snowball import (
    channel_basis,
    chordal_distance_sq,
    compare_channels,
    covariance_cross_load,
)

G = np.array([[1.0, 0.25], [0.25, 0.6]])
H1 = np.array([[0.5, 0.0], [0.0, 1.4]])
H2 = np.array([[1.3, 0.0], [0.0, 0.6]])

channel = channel_basis(G, H1, k=1)
angles = compare_channels(G, H1, G, H2, k=1)
distance = chordal_distance_sq(angles)

# Fitness is a separate calculation. This example has zero covariance-cross
# load even though the response channels differ.
load = covariance_cross_load(G, G, H1, H2)
```

## Repository map

```text
src/geometric_snowball/
  linalg.py              SPD geometry, shared metrics, principal angles
  channels.py            response operator, channels, perturbation theory
  hybrid_fitness.py      Gaussian fitness and covariance–curvature bridge
  dynamics.py            standing variation and phenomenological DMI handoff
  simulations.py         canonical numerical experiments

tests/                   unit and mathematical-identity tests
scripts/                 derivation audit and figure reproduction
examples/                compact worked examples

docs/
  derivations/           full teaching monograph and technical note sources
  THEORY.md              compact formal statement
  MATHEMATICAL_CONVENTIONS.md
  BIOLOGICAL_INTERPRETATION.md
  EMPIRICAL_PREDICTIONS.md
  MIGRATION_FROM_V1.md
  DATA_PROVENANCE.md

analysis/senecio/         cautious descriptive pair analysis
data/derived/             small reproducible summary only
figures/                  generation instructions, not binary outputs
archive/README.md         recovery of the former repository through history
```

## Validation

The automated suite tests:

- whitening and SPD square-root identities;
- spectral equivalence of \(\mathbf G\mathbf H\) and \(\mathbf S\);
- recovery of \(\mathbf g_{\max}\) when \(\mathbf H\propto\mathbf I\);
- coordinate invariance under a shared metric;
- the local \(t^2\) channel-distance coefficient;
- zero first-order rotation for commuting perturbations;
- derivatives of \(\mathbf G^{1/2}\) and \(\mathbf S\);
- expected quadratic load and exact Gaussian mean fitness;
- the general covariance–curvature identity;
- the exact ridge theorem;
- standing-variation response with migration; and
- the direction of the geometric-to-DMI feedback.

Continuous integration repeats the tests, numerical audit, package build, and figure smoke tests on every push and pull request.

## Empirical status

The *Senecio lautus* radiation motivates the theory but does not yet directly test it. Existing work establishes divergence in multivariate genetic variance, strong habitat-dependent fitness differences, F1 heterosis, severe F2 fertility loss, and later loss of fitness variance and habitat trade-offs. A decisive test must estimate population-specific \(\mathbf G\), local \(\mathbf H\), recombinant covariance, and hybrid fitness in the same replicated population pairs.

The committed pair summary contains only standardised distances in measured morphology and soil chemistry. It is retained as an exploratory design aid. It is not converted into a geometric isolation score.

## Documentation and citation

The most complete mathematical teaching treatment is [`docs/derivations/full_derivations.tex`](docs/derivations/full_derivations.tex). The compact theory is in [`docs/THEORY.md`](docs/THEORY.md), with empirical predictions in [`docs/EMPIRICAL_PREDICTIONS.md`](docs/EMPIRICAL_PREDICTIONS.md). Generated PDFs and figure files are deliberately excluded from the clean source tree. With a LaTeX distribution and `latexmk` installed, rebuild the full monograph, technical note, audits, and figures with:

```bash
make docs
```

Please cite the repository using [`CITATION.cff`](CITATION.cff), together with the relevant primary empirical papers listed in [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md).

## License

Software is released under the MIT License. Manuscript text, teaching documents, and original figures remain © the authors unless a file states otherwise. Published empirical data retain the terms and citation requirements of their original sources; see [`NOTICE.md`](NOTICE.md).
