# Optional descriptive Senecio analysis

This directory retains a small exploratory analysis because the replicated Dune–Headland pairs are useful for designing a direct test. It does not claim to test the geometric snowball with the currently supplied variables.

The analysis calculates:

- standardised Euclidean distance between population mean morphologies;
- standardised Euclidean distance between measured soil-chemistry profiles; and
- a morphology PCA for visualisation.

These are descriptive quantities. Morphology-mean distance is not a channel angle, and soil distance is not \(\Delta\mathbf H\). A test of the theory requires population-specific \(\mathbf G\), selection curvature \(\mathbf H\), recombinant covariance, and hybrid fitness.

Run with local copies of the published project spreadsheets:

```bash
python analysis/senecio/senecio_pair_analysis.py \
  --morphology /path/to/data_AL_Arch\&Leaf_ARVersion.xlsx \
  --soil /path/to/All\ S.\ lautus\ soil\ data.xlsx \
  --outdir analysis/senecio/results

python analysis/senecio/senecio_figure.py \
  --summary analysis/senecio/results/senecio_pair_summary.csv \
  --outfile analysis/senecio/results/senecio_pair_descriptive.pdf
```

The plotting script can also use the committed small summary under `data/derived/`, which permits a CI smoke test without redistributing raw spreadsheets.
