# HFDR Pipeline

A small, CI-friendly pipeline step for mass-univariate IDP analysis + 2-level Hierarchical FDR (HFDR).

## Input format (long table)

Your CSV should look like your `subset_data_exp` table. Required columns:

- subject
- sex
- phenotype (must include "Control" plus one or more diseases)
- IDP
- interview_age (months)

Optional (but recommended):

- value_z (preferred outcome; otherwise value)
- modality (if missing, derived from IDP prefix before first underscore)
- measure

## Run locally

From the repo root:

```bash
Rscript inst/scripts/run_analysis.R --input data/example.csv --out_dir artifacts --alpha 0.05
```

Outputs are written under:

```
artifacts/<DATASET>/<DISEASE>/
  raw_outcomes.csv
  hfdr_results.csv
  hfdr_parent_table.csv
  hfdr_summary.png
artifacts/<DATASET>/
  raw_outcomes_ALL.csv
  hfdr_results_ALL.csv
```

## Reproducibility

Use `renv`:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Then CI will run `renv::restore()`.

## GitHub Actions

Workflow is in `.github/workflows/ci.yaml` and runs on every push/PR.
It restores dependencies, runs tests, runs the pipeline on `data/example.csv`,
and uploads `artifacts/` as a build artifact.
