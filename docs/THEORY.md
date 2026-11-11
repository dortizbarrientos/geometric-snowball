# Compact theory

## 1. Local response geometry

Around optimum \(\boldsymbol\theta_i\), let

\[
\log W_i(\mathbf z)=c_i-rac12(\mathbf z-\boldsymbol\theta_i)^\top\mathbf H_i(\mathbf z-\boldsymbol\theta_i),
\qquad \mathbf H_i\succ0.
\]

Lande's equation gives

\[
\dot{\bar{\mathbf z}}_i
=\mathbf G_i\mathbf H_i(\boldsymbol\theta_i-\bar{\mathbf z}_i).
\]

The generally non-symmetric matrix \(\mathbf G_i\mathbf H_i\) has the correct dynamical eigenvalues but not the Euclidean geometry needed for orthonormal channel bases. Its symmetric similar operator is

\[
\mathbf S_i=\mathbf G_i^{1/2}\mathbf H_i\mathbf G_i^{1/2}.
\]

If \(\mathbf S_i\mathbf v_{ij}=\lambda_{ij}\mathbf v_{ij}\), the corresponding phenotype-space direction is \(\mathbf u_{ij}=\mathbf G_i^{1/2}\mathbf v_{ij}\), and

\[
\mathbf G_i\mathbf H_i\mathbf u_{ij}=\lambda_{ij}\mathbf u_{ij}.
\]

When \(\mathbf H_i=h\mathbf I\), \(\mathbf S_i=h\mathbf G_i\), so the leading channel recovers \(\mathbf g_{\max}\).

## 2. Comparing channels

For equal \(\mathbf G\), compare the orthonormal leading eigenvectors of \(\mathbf S_1\) and \(\mathbf S_2\) directly in their shared whitened space.

For unequal \(\mathbf G_1\) and \(\mathbf G_2\), map the channels to phenotype space,

\[
\mathbf U_i=\mathbf G_i^{1/2}\mathbf V_i,
\]

choose one biologically defensible common metric \(\mathbf M\succ0\), and form

\[
\mathbf Q_i=\mathbf U_i(\mathbf U_i^\top\mathbf M\mathbf U_i)^{-1/2}.
\]

The singular values of \(\mathbf Q_1^\top\mathbf M\mathbf Q_2\) are \(\cos\psi_j\). This construction is invariant to every nonsingular linear transformation of the measured traits, provided the matrices and metric transform consistently.

Useful discrepancies are

\[
d_C^2=\sum_j\sin^2\psi_j
=\frac12\|\mathbf P_1-\mathbf P_2\|_F^2,
\]

\[
d_G=\left(\sum_j\psi_j^2\right)^{1/2},
\]

and

\[
I_{\rm overlap}=1-\prod_j\cos^2\psi_j.
\]

The last object is a volume-overlap deficit. None is a fitness measure without an explicit biological calibration.

## 3. Exact local geometric snowball

Suppose

\[
\mathbf S_i(t)=\mathbf S_0+t\mathbf A_i+O(t^2).
\]

Let \(\mathbf v_1,\ldots,\mathbf v_k\) span the ancestral channel and \(\mathbf w_1,\ldots,\mathbf w_{n-k}\) its complement. Define

\[
\mathcal K_{ab}
=\frac{\mathbf w_a^\top(\mathbf A_2-\mathbf A_1)\mathbf v_b}
{\lambda_b-\lambda_{k+a}}.
\]

Then

\[
d_C^2(t)=t^2\|\mathcal K\|_F^2+O(t^3).
\]

The numerator is cross-channel forcing. The denominator is spectral resistance. A small eigengap amplifies forcing but cannot create rotation when the perturbation commutes with \(\mathbf S_0\).

Smooth deterministic rotation therefore gives \(d_C^2=O(t^2)\). Unbiased diffusive channel motion gives \(E[d_C^2]=O(t)\) locally. A mixed process motivates, but does not force, the empirical form \(E[d_C^2]\approx at+bt^2\).

## 4. From geometry to hybrid fitness

For a Gaussian hybrid phenotype \(\mathbf Z_H\sim N(\boldsymbol\mu_H,\mathbf P_H)\), expected negative log fitness in environment \(e\) is

\[
L_{H,e}=\frac12(\boldsymbol\mu_H-\boldsymbol\theta_e)^\top\mathbf H_e(\boldsymbol\mu_H-\boldsymbol\theta_e)
+\frac12\operatorname{tr}(\mathbf H_e\mathbf P_H).
\]

Arithmetic mean fitness also has a closed form, implemented in `hybrid_fitness.gaussian_mean_fitness`.

Write

\[
\mathbf P_H=\frac12(\mathbf G_1+\mathbf G_2)+\mathbf R.
\]

Relative to the two parental covariance benchmarks, the mean covariance contribution is

\[
L_{\rm cov}
=-\frac18\operatorname{tr}[(\mathbf H_2-\mathbf H_1)(\mathbf G_2-\mathbf G_1)]
+\frac14\operatorname{tr}[(\mathbf H_1+\mathbf H_2)\mathbf R].
\]

Thus channel divergence may increase load, have no effect, or contribute to heterosis. The sign is biological, not guaranteed by subspace distance.

## 5. Exact ridge calibration

For rank-\(k\) channel projectors \(\mathbf P_i\), let

\[
\mathbf G_i=g_\perp\mathbf I+(g_\parallel-g_\perp)\mathbf P_i,
\]

\[
\mathbf H_i=h_\perp\mathbf I-(h_\perp-h_\parallel)\mathbf P_i,
\]

with \(g_\parallel>g_\perp\) and \(h_\perp>h_\parallel\). With \(\mathbf R=0\),

\[
L_{\rm cov}
=\frac{(g_\parallel-g_\perp)(h_\perp-h_\parallel)}{4}
\sum_{j=1}^k\sin^2\psi_j.
\]

This theorem is the cleanest angle-to-fitness bridge in the framework. Its assumptions are strong and testable.

## 6. Standing variation and later intrinsic isolation

For additive allelic effects \(\mathbf a_\ell\),

\[
\bar{\mathbf z}_i=2\sum_\ell(p_{i\ell}-p_{\ell0})\mathbf a_\ell,
\qquad
\mathbf G_i=2\sum_\ell p_{i\ell}(1-p_{i\ell})\mathbf a_\ell\mathbf a_\ell^\top.
\]

Under weak selection and migration,

\[
\dot p_{i\ell}=p_{i\ell}(1-p_{i\ell})\mathbf a_\ell^\top\boldsymbol\beta_i
+m(p_{j\ell}-p_{i\ell}),
\]

which implies

\[
\dot{\bar{\mathbf z}}_i=\mathbf G_i\boldsymbol\beta_i+m(\bar{\mathbf z}_j-\bar{\mathbf z}_i).
\]

The early ecological load may reduce effective migration through a phenomenological relation such as \(m_{\rm eff}=m_0e^{-L_E}\). A separate substitution-and-epistasis model is then required to generate DMIs. The repository keeps these two stages mathematically distinct.
