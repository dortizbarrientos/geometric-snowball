# Empirical predictions and experimental design

## Primary predictions

### Prospective susceptibility

For an ancestral channel of dimension \(k\), estimate

\[
\chi_k^2=\|\mathcal K\|_F^2.
\]

Pairs with larger \(\chi_k^2\) should show faster early channel divergence per unit of measured perturbation. A small spectral gap alone is insufficient because the ecological perturbation may fail to couple the channel to its complement.

### Temporal exponent

Repeated smooth directional selection should produce a reproducible quadratic component in squared channel distance. Unbiased channel diffusion should produce an approximately linear expectation. A hierarchical time-series analysis can estimate both components:

\[
E[d_C^2(t)]\approx at+bt^2,
\]

while allowing replicate-specific paths and uncertainty in \(\mathbf G\) and \(\mathbf H\).

### Fitness calibration

Hybrid fitness should be predicted better by the exact Gaussian load or its covariance–curvature component than by:

- neutral genetic distance;
- Euclidean distance between \(\mathbf G\) matrices;
- environmental Euclidean distance;
- phenotype-mean distance; or
- principal angle alone.

The comparison should be prospective or cross-validated.

### Generation and environment specificity

The additive covariance mechanism should be clearer in F2 and advanced recombinant generations than in F1s. The early barrier should be strongest in parental habitats and may weaken in a benign common garden. A later intrinsic DMI component should remain in common conditions.

### Handoff

Conditional on divergence time and opportunity for contact, stronger early ecological load should predict lower effective gene flow and a greater later accumulation of intrinsic recombinant failure.

## A decisive replicated Senecio design

1. Estimate population-pair-specific \(\mathbf G_i\), with posterior draws rather than point estimates alone.
2. Estimate local quadratic selection surfaces from trait–fitness data in each parental habitat, yielding posterior distributions of \(\mathbf H_i\).
3. Select channel dimension from supported eigengaps and biological interpretation, not a fixed \(k=1\) convention.
4. Map channels to phenotype space and compare them under a prespecified shared metric; repeat under defensible alternatives.
5. Generate F1, F2, and advanced recombinant hybrids for each Dune–Headland pair.
6. Measure hybrid means, covariance, segregation excess, and fitness in both parental habitats and a common garden.
7. Partition mean mismatch, covariance load, segregation load, extrinsic isolation, and intrinsic isolation.
8. Compare preregistered predictive models using held-out population pairs or posterior predictive scoring.
9. Include within-ecotype crosses to distinguish repeated ecological forcing from pair-specific mutation-order processes.
10. Track genomic ancestry and recombinant viability to connect the quantitative-genetic precursor to later locus-specific incompatibilities.

## Estimating H

A local model can be written as

\[
\log W=\alpha+\boldsymbol\beta^\top\mathbf z-rac12\mathbf z^\top\mathbf H\mathbf z+\varepsilon.
\]

In practice, binary survival, seed set, or lifetime fitness may require appropriate generalized models. The Hessian must be transformed consistently with the trait scale and estimated with uncertainty. Environmental distance does not estimate \(\mathbf H\): it describes habitat difference, whereas \(\mathbf H\) describes the local curvature of fitness in phenotype space.

## Decision criteria

The theory gains support only if channel geometry adds reproducible, out-of-sample information about hybrid load and if the predicted temporal and generation-specific signatures are observed. A visually appealing correlation among distances is not sufficient.
