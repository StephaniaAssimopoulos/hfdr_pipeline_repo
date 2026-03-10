# HFDR Pipeline

Reproducible data analysis pipeline in R with a CLI entrypoint, dependency locking (renv), and GitHub Actions CI that runs the pipeline and produces results (tables + plots).

Runs mass-univariate IDP analysis + 2-level Hierarchical FDR (HFDR).

## Input format (long table)

Your CSV should look like your `subset_data_exp` table. Required columns:

- subject
- sex
- phenotype (must include "Control" plus one or more diseases)
- IDP
- interview_age (consistent format across subjects)

Optional (but recommended):

- value_z (preferred outcome; otherwise value)
- modality (if missing, derived from IDP prefix before first underscore)
- measure

## Run locally

Clone the repository:

```bash
git clone https://github.com/StephaniaAssimopoulos/hfdr_pipeline_repo.git
cd hfdr_pipeline_repo
```

Install project dependencies:

```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")'
Rscript -e 'renv::restore(prompt = FALSE)'
```

Run the pipeline from the repo root:

```bash
Rscript inst/scripts/run_analysis.R --input data/example.csv --out_dir output --alpha 0.05
```

Outputs are written under:

```
output/<DATASET>/<DISEASE>/
  raw_outcomes.csv
  hfdr_results.csv
  hfdr_parent_table.csv
  hfdr_summary.png
output/<DATASET>/
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
It restores dependencies, runs tests, runs the pipeline on `data/example.csv` (as provided),
and uploads `output/` as a build output.

