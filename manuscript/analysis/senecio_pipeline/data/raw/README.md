# Raw data — provenance and column descriptions

These files are **immutable inputs** to the Senecio channel-operator pipeline.
Do not modify them; if corrections are needed, document them in the pipeline
scripts instead.

---

## Data_Exp1_Dmatrix.csv

**Source:** Common-garden experiment (Experiment 1) from the
Walter et al. breeding-design study.

**Rows:** 242 individuals across 8 populations.

**Used by:** Stage 1 — computing pairwise D among population means and
building the global leaf PCA basis.

| Column | Description |
|--------|-------------|
| Ecotype | Dune, Head, Table, or Wood |
| Population | D1, D3, H1, H2, T1, T9, W2, W4 |
| VegHeight | Plant height in mm (renamed to `Height` internally) |
| MSL_W | Mean specific leaf weight |
| SB | Stem branching count |
| MSD | Mean stem diameter |
| Area | Leaf area |
| P2A2 | Perimeter² / Area² (leaf complexity) |
| Circularity | Leaf circularity |
| Nindents.Peri | Number of indentations per unit perimeter |
| IndentWidth | Mean leaf indentation width |
| IndentDepth | Mean leaf indentation depth |

---

## Data_Exp2_Gvariance.csv

**Source:** Breeding-design offspring (Experiment 2) — full-sib families
from a sire-dam nested design within each population.

**Rows:** 934 offspring across 8 populations.

**Used by:** Stage 1 — fitting population-specific G matrices via MCMCglmm.

| Column | Description |
|--------|-------------|
| Type | Population identifier (= Block) |
| Block | Population: D1, D3, H1, H2, T1, T9, W2, W4 |
| ID | Individual identifier |
| Replicate | Replicate number |
| Family | Full-sibling family ID |
| dam | Dam identifier |
| sire | Sire identifier |
| Height | Plant height in mm |
| MSL_W | Mean specific leaf weight |
| SB | Stem branching count |
| MSD | Mean stem diameter |
| Area | Leaf area |
| P2A2 | Leaf complexity |
| Circularity | Leaf circularity |
| Nindents.Peri | Indentations per unit perimeter |
| IndentWidthMean | Mean indentation width (renamed to `IndentWidth` internally) |
| IndentDepthMean | Mean indentation depth (renamed to `IndentDepth` internally) |

---

## Walter_Data_Exp2.csv

**Source:** Field transplant experiment from Walter et al. (2016) *Evolution*
and Walter et al. (2020) *Evolution Letters*.

**Rows:** 26 246 seeds/seedlings transplanted into four environments.

**Used by:** Stage 2 — linking glasshouse morphology to field fitness via
shared family IDs.

| Column | Description |
|--------|-------------|
| ENV | Transplant environment: DUNE, HEAD, TABLE, WOOD |
| block | Experimental block |
| GRID | Grid reference in the field |
| GENO | Generation (parental, F1, F4) |
| Type | General cross-type identifier |
| TypeS | Specific identifier (R, Y, P for F4 crossing lines; population codes for parentals) |
| Family | Full-sibling family — shared with Data_Exp2_Gvariance.csv |
| Dam, Sire | F3 parental identifiers |
| GE | Germination/emergence (binary) |
| D100 | Alive at day 100 (binary) |
| L10 | Reached seedling establishment, 10 leaves (binary) — **default fitness measure** |
| Bud | Reached maturity, produced a bud (binary) |
| DayL10–DayL10Bud | Days to reach each developmental stage |

---

## Walter_Data_Exp1.csv

**Source:** Parental ecotype morphology from Walter et al. (2020) *Evolution
Letters*.

**Rows:** 1 139 individuals (parental ecotypes + F4 hybrid).

**Used by:** Not directly used by the pipeline, but archived for reference
and potential future analyses.

---

## All_S__lautus_soil_data.xlsx

**Source:** Soil chemistry sampled at each population site. Measurements
from commercial soil testing (methods described in published supplements).

**Sheet used:** "All soil data"

**Used by:** Stage 2 — building the soil-based rank-1 H proxy. Three
composite axes (Salinity, Nutrients, Acidity) are extracted from PCA on
subsets of the soil variables.

Key variable groups:

- **Salinity:** Electrical Conductivity, Chloride, Exchangeable Na/Mg/Ca
- **Nutrients:** Total N, Total Organic C, Colwell P, CEC
- **Acidity:** pH, Exchangeable Al, Aluminium

The `Population` column uses site codes (D01, D03, H01, H02, T1, T9, I5,
I24), which Stage 2 maps to the breeding-design codes (D1, D3, H1, H2, T1,
T9, W2, W4) internally.

---

## CrossesDH-ReciprocalsAveraged.txt

**Source:** James et al. (2025) *Molecular Ecology* — reproductive
isolation measurements between Dune and Headland population pairs,
reciprocals averaged.

**Rows:** 18 population pairs.

**Used by:** Downstream analysis (not directly consumed by the pipeline's
R scripts, but available for merging with channel-angle outputs to test
whether S-based angles predict RI better than gmax angles).

| Column | Description |
|--------|-------------|
| CrossType | Population pair (e.g., D01-H01) |
| RIseedSet-SobelChen | Reproductive isolation index — seed set |
| RIviability-SobelChen | Reproductive isolation index — viability |
| Ecotype | Same or Different |
| Clade | Within or Between |
| Divergence time | Estimated divergence time |
| EnvDist-MDS-5PCs | Environmental distance (MDS) |
| PhenoDist-MDS-7PCs | Phenotypic distance (MDS) |

---

## IntrinsicRIdata.txt

**Source:** James et al. (2025) *Molecular Ecology* — raw individual-cross
RI data (both directions, not averaged).

**Rows:** 33 directed crosses.

**Used by:** Same downstream purpose as CrossesDH-ReciprocalsAveraged.txt.
Provides per-direction RI estimates (mother → father) for crosses that
show reciprocal asymmetry.

---

## Data citations

- Walter GM, Aguirre JD, Blows MW, Ortiz-Barrientos D (2016)
  Diversification across a heterogeneous landscape. *Evolution* 70: 1979–1992.
- Walter GM, Wilkinson MJ, James ME, et al. (2020) Environmentally induced
  development costs underlie fitness tradeoffs. *Ecology* 99: 1–11.
- Walter GM, et al. (2020) Data from: *Evolution Letters*.
  Dryad. doi:10.5061/dryad.sbcc2fr3t
- James ME, Melo MC, Roda F, et al. (2025) Ecological and mutation-order
  speciation in *Senecio*. *Molecular Ecology* [doi].
