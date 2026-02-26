#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Path to input long-format CSV"),
  make_option(c("-o","--out_dir"), type="character", default="artifacts", help="Output directory"),
  make_option(c("--alpha"), type="double", default=0.05, help="HFDR alpha"),
  make_option(c("--outcome_col"), type="character", default=NULL,
              help="Outcome column (default: value_z if present else value)")
)

options(error = function() { traceback(2); quit(status = 1) })

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Please provide --input path/to/file.csv")

# Expect to be run from repo root (works in CI and typical local usage)
source("R/hfdr.R")
source("R/plots.R")
source("R/pipeline.R")

run_hfdr_pipeline(
  input_csv = opt$input,
  out_dir = opt$out_dir,
  alpha = opt$alpha,
  outcome_col = opt$outcome_col
)
