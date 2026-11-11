# Data provenance and empirical scope

## Core theory

The mathematical package and canonical numerical experiments are synthetic and have no empirical-data dependency.

## Optional Senecio analysis

The descriptive pair analysis expects two previously used project spreadsheets, supplied locally as command-line inputs:

- `data_AL_Arch&Leaf_ARVersion.xlsx`: plant architecture and leaf traits;
- `All S. lautus soil data.xlsx`: population soil chemistry.

They are not duplicated in the v2 Git tip. The former repository history and original data archives remain the appropriate provenance trail. The committed `data/derived/senecio_pair_summary.csv` is a small descriptive output retained to exercise the plotting workflow.

Run the analysis with:

```bash
python analysis/senecio/senecio_pair_analysis.py \
  --morphology /path/to/data_AL_Arch\&Leaf_ARVersion.xlsx \
  --soil /path/to/All\ S.\ lautus\ soil\ data.xlsx \
  --outdir analysis/senecio/results
```

## Primary empirical papers

The current biological motivation draws on:

- Walter, G. M. et al. (2016), *Diversification across a heterogeneous landscape*, **Evolution** 70:1979–1992, doi:10.1111/evo.13009. Data archive: doi:10.5061/dryad.d9f36.
- Walter, G. M. et al. (2018), *Evolution of genetic variance during adaptive radiation*, **The American Naturalist** 191:E108–E128.
- Walter, G. M. et al. (2020), *Loss of ecologically important genetic variation in late generation hybrids reveals links between adaptation and speciation*, **Evolution Letters** 4:302–316, doi:10.1002/evl3.187. Data archive: doi:10.5061/dryad.sbcc2fr3t.

These studies motivate population-specific predictions but do not, by themselves, estimate the geometric snowball quantities. In particular, the pairwise descriptive distances in this repository must not be relabelled as \(\Delta\mathbf H\), channel angle, or reproductive isolation.
