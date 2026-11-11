# Mathematical conventions

## Spaces and metrics

Three spaces must not be conflated.

**Phenotype coordinates** contain measured traits. Covariances transform as \(\mathbf G' = \mathbf T\mathbf G\mathbf T^\top\), while curvature transforms contragrediently as \(\mathbf H'=\mathbf T^{-\top}\mathbf H\mathbf T^{-1}\).

**Population-specific whitened coordinates** use \(\mathbf x_i=\mathbf G_i^{-1/2}\mathbf z\). They are natural for diagonalising \(\mathbf S_i\), but two different whitenings are not automatically one common Euclidean space.

**A shared comparison geometry** is supplied by a single positive-definite metric \(\mathbf M\) in phenotype space. Plausible choices include ancestral \(\mathbf G_0^{-1}\), inverse pooled covariance, or a measurement-error-adjusted biological metric. The choice must be declared and subjected to sensitivity analysis.

## Positive definiteness

The current implementation assumes positive-definite \(\mathbf G\) and \(\mathbf H\). Empirical estimates close to singularity require regularisation or restriction to an estimable subspace. Silently replacing a low-rank matrix with an arbitrary ridge can change the inferred channel dimension and should be avoided.

## Channel dimension

A one-dimensional leading vector is justified only when the first eigengap is biologically and statistically supported. If \(\lambda_1\approx\lambda_2\) but \(\lambda_2\gg\lambda_3\), the two-dimensional channel may be stable even though its leading vector is not. The relevant boundary is \(\delta_k=\lambda_k-\lambda_{k+1}\).

## Principal angles

For two \(k\)-dimensional subspaces, there are \(k\) principal angles. The largest angle is useful for worst-direction bounds; squared chordal distance uses all angles; determinant overlap is sensitive to the joint volume. Report the full angle spectrum whenever possible.

## Perturbation notation

The local theorem uses

\[
\mathbf S_i(t)=\mathbf S_0+t\mathbf A_i+O(t^2).
\]

The coefficient \(\mathcal K\) is relative: it depends on \(\mathbf A_2-\mathbf A_1\), not on either perturbation alone. Davis–Kahan gives an upper bound. The matrix \(\mathcal K\) gives the realised first-order rotational coefficient.

## Evolving G

If \(\mathbf R_G=\mathbf G^{1/2}\) and \(\mathbf X=D\sqrt{\mathbf G}[\dot{\mathbf G}]\), then

\[
\mathbf R_G\mathbf X+\mathbf X\mathbf R_G=\dot{\mathbf G},
\]

and

\[
\dot{\mathbf S}
=\mathbf X\mathbf H\mathbf R_G
+\mathbf R_G\dot{\mathbf H}\mathbf R_G
+\mathbf R_G\mathbf H\mathbf X.
\]

This separates contributions of evolving standing variation and evolving fitness curvature without pretending that matrix square roots differentiate element by element.

## Fitness conventions

\(\mathbf H\) is the curvature of **minus log fitness**, so it is positive definite near a local maximum of fitness. Log-loads add naturally. A bounded isolation index such as \(1-e^{-L}\) is optional and should be introduced only after the load is shown to be nonnegative for the comparison being made.
