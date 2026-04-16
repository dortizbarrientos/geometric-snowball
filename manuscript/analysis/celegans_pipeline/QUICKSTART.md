# Quickstart

Three steps. Should take about 2 minutes on a laptop.

## 1. Install dependencies

```bash
pip install -r requirements.txt
```

You also need R (≥ 4.0) on your PATH. No R packages are required.

```bash
Rscript --version
```

should print something like `R scripting front-end version 4.x.x`.

## 2. Run the pipeline

From the repository root:

```bash
bash run_pipeline.sh
```

The script runs six stages in sequence and prints progress as it goes.
Expect ~30 seconds per environment for the R chain extraction and
another 30 seconds for the Python posterior propagation.

## 3. Look at the outputs

- `results/figures/celegans_posterior_violins.pdf` — the main figure
- `results/tables/table_celegans_summary.docx` — the summary table
- `results/tables/celegans_posterior_summary.csv` — raw numbers

Compare your numbers to `docs/expected_outputs.md`. They should match
exactly to the reported precision. If they don't, see that document's
troubleshooting section.

## Re-running selectively

Each pipeline stage is idempotent. You can comment out earlier stages
in `run_pipeline.sh` if their outputs already exist. For example, the
slow R extraction (stage 01c) only needs to run once per clone.
