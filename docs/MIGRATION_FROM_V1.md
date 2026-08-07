# Migration from the former repository

The v2 tree is a deliberate scientific reset. The previous files remain recoverable from the parent commit and complete Git history, but they are not retained in the canonical working tree because several names and interpretations were no longer defensible.

| Former practice | v2 replacement | Reason |
|---|---|---|
| Euclidean channel geometry from \(\mathbf G\mathbf H\) | Eigenspaces of \(\mathbf S=\mathbf G^{1/2}\mathbf H\mathbf G^{1/2}\) | \(\mathbf G\mathbf H\) is generally non-symmetric; \(\mathbf S\) supplies the correct orthonormal geometry |
| Direct angle comparison between separately whitened populations | Phenotype-space channels under one declared metric | Separate whitenings do not define a shared Euclidean coordinate system |
| Davis–Kahan bound interpreted as rotation law | Exact local coefficient \(\mathcal K\), with Davis–Kahan retained as a bound | Spectral gap is susceptibility; realised rotation also requires cross-channel forcing |
| \(I_{\rm geom}=1-\prod\cos^2\psi_j\) labelled reproductive isolation or hybrid fitness cost | `volume_overlap_deficit` labelled geometric overlap only | Angle is not fitness without a fitness model |
| Phenotype distance labelled a proxy for channel angle | Descriptive morphology distance | Mean divergence and covariance-subspace divergence are different quantities |
| Soil distance labelled a proxy for \(\Delta\mathbf H\) | Descriptive soil-chemistry distance | Habitat descriptors do not estimate curvature of fitness in trait space |
| Quadratic curve imposed by an assumed linear angle path | Local perturbation theorem deriving \(d_C^2=t^2\|\mathcal K\|_F^2+O(t^3)\) | Separates assumption, coefficient, and scope of the early law |
| Hybrid cost inferred directly from overlap deficit | Exact Gaussian fitness and covariance–curvature trace identity | Identifies when geometry produces breakdown, no effect, or heterosis |
| DMI acceleration presented as a direct implication | Separately labelled phenomenological handoff | Intrinsic incompatibility requires an additional substitution-and-epistasis model |
| Hard-coded local filesystem paths | Command-line arguments and repository-relative defaults | Makes analyses portable and auditable |
| Mixed manuscripts, raw data, generated binaries, and exploratory scripts at root | Installable package, tests, docs, optional analysis, generated outputs ignored | Makes the repository reviewable and reproducible |

## Recovering the old tree

The commit immediately preceding the v2 replacement is

```text
ff9b0cf196af823cef372e55c9b91cf6e0c0aebf
```

To inspect it without changing the current branch:

```bash
git show ff9b0cf196af823cef372e55c9b91cf6e0c0aebf:README.md
git worktree add ../geometric-snowball-v1 ff9b0cf196af823cef372e55c9b91cf6e0c0aebf
```

Do not copy old figures or scripts into v2 without reviewing their mathematical labels against the table above.
