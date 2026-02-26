# ------------------------------------------------------------------
# Hierarchical FDR (HFDR) utilities
# ------------------------------------------------------------------

# ---- Parent combination methods ----------------------------------

simes_p <- function(pvalues) {
  pvalues <- pvalues[is.finite(pvalues)]
  if (length(pvalues) == 0) return(NA_real_)

  m <- length(pvalues)
  ps <- sort(pvalues)
  min((m / seq_len(m)) * ps)
}

bonf_minp <- function(pvalues) {
  pvalues <- pvalues[is.finite(pvalues)]
  if (length(pvalues) == 0) return(NA_real_)

  min(1, length(pvalues) * min(pvalues))
}

# ---- Two-level HFDR ----------------------------------------------

hfdr_two_level <- function(outcome_table,
                           alpha = 0.05,
                           parent_col = "modality",
                           p_col = "pval",
                           parent_method = c("bonf_minp", "simes")) {
  
  parent_method <- match.arg(parent_method)
  
  if (!all(c(parent_col, p_col) %in% names(outcome_table))) {
    stop("outcome_table must contain columns: ", parent_col, " and ", p_col)
  }
  
  parent_fun <- switch(
    parent_method,
    bonf_minp = bonf_minp,
    simes    = simes_p
  )
  
  # Parent p per group
  parent_split <- split(outcome_table[[p_col]], outcome_table[[parent_col]])
  parent_p <- vapply(parent_split, parent_fun, numeric(1))
  
  parent_table <- data.frame(
    parent = names(parent_p),
    parent_p = as.numeric(parent_p),
    stringsAsFactors = FALSE
  )
  parent_table$parent_q <- p.adjust(parent_table$parent_p, method = "BH")
  parent_table$reject_parent <- parent_table$parent_q < alpha
  
  rejected_parents <- parent_table$parent[parent_table$reject_parent]
  
  # Child testing
  outcome_table$tested_child <- outcome_table[[parent_col]] %in% rejected_parents
  outcome_table$qval_child <- NA_real_
  outcome_table$reject_child <- FALSE
  
  if (length(rejected_parents) > 0) {
    for (grp in rejected_parents) {
      idx <- which(outcome_table[[parent_col]] == grp)
      child_q <- p.adjust(outcome_table[[p_col]][idx], method = "BH")
      outcome_table$qval_child[idx] <- child_q
      outcome_table$reject_child[idx] <- child_q < alpha
    }
  }
  
  # Build "results" in the shape your pipeline expects:
  # - a column named `parent` (the group label)
  # - parent_p/parent_q/reject_parent
  # - the original outcome rows + tested/qval/reject
  results <- outcome_table
  results$parent <- results[[parent_col]]
  
  # merge parent stats onto each row
  results <- merge(
    results,
    parent_table,
    by = "parent",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Return BOTH styles for compatibility
  list(
    results = results,
    outcome_table = outcome_table,
    parent_table = parent_table
  )
}