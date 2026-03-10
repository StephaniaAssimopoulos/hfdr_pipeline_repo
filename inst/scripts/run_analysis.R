#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Path to input long-format CSV"),
  make_option(c("-o","--out_dir"), type="character", default="output", help="Results directory [default: %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="HFDR alpha [default: %default]"),
  make_option(c("--outcome_col"), type="character", default=NULL,
              help="Outcome column (default: value_z if present else value)")
)

options(error = function() { traceback(2); quit(status = 1) })

parser <- OptionParser(
  usage = "%prog --input data/example.csv --out_dir output --alpha 0.05",
  option_list = option_list,
  description = "Run mass-univariate IDP analysis followed by two-level hierarchical FDR (HFDR)."
)

opt <- parse_args(parser)

if (is.null(opt$input)) {
  stop("Please provide --input path/to/file.csv", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call. = FALSE)
}

# Increase buffer for large CSV lines
Sys.setenv(VROOM_CONNECTION_SIZE = 10485760)

# Ensure output directory exists
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

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
