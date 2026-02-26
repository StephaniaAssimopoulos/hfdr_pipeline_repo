#' Run HFDR pipeline on one long-format CSV
#'
#' Expects columns like: subject, sex, phenotype, IDP, interview_age, modality,
#' plus an outcome column (`value_z` preferred, otherwise `value`).
#'
#' @param input_csv path to CSV
#' @param out_dir output directory
#' @param alpha HFDR/FDR level
#' @param outcome_col optional outcome column name; defaults to value_z if present else value
#' @param dataset optional dataset label; defaults to filename prefix before first underscore
#'
#' @return data.frame of HFDR-annotated results
#' @export
run_hfdr_pipeline <- function(input_csv, out_dir = "artifacts", alpha = 0.05,
                              outcome_col = NULL, dataset = NULL) {

  rlang::check_installed(c("readr","dplyr","ggplot2","patchwork","ggrepel"))

  if (!file.exists(input_csv)) rlang::abort(paste0("Input CSV not found: ", input_csv))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  df <- readr::read_csv(input_csv, show_col_types = FALSE)

  required <- c("subject","sex","phenotype","IDP","interview_age")
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) rlang::abort(paste0("Missing required columns: ", paste(miss, collapse=", ")))

  if (is.null(outcome_col)) {
    outcome_col <- if ("value_z" %in% names(df)) "value_z" else if ("value" %in% names(df)) "value" else NULL
  }
  if (is.null(outcome_col) || !(outcome_col %in% names(df))) {
    rlang::abort("Could not find outcome column. Provide outcome_col, or include value_z/value.")
  }

  if (is.null(dataset)) {
    nm <- basename(input_csv)
    dataset <- strsplit(nm, "_", fixed = TRUE)[[1]][1]
  }

  # modality if missing: derive from IDP prefix
  if (!("modality" %in% names(df))) {
    df <- df |> dplyr::mutate(modality = vapply(strsplit(.data$IDP, "_", fixed = TRUE), `[[`, character(1), 1))
  }

  # prep covariates
  df <- df |>
    dplyr::mutate(
      phenotype_fac = stats::relevel(factor(.data$phenotype), ref = "Control"),
      sex_fac = stats::relevel(factor(.data$sex), ref = "M"),
      age_yrs = .data$interview_age / 12.0
    )

  diseases <- df |>
    dplyr::filter(.data$phenotype != "Control") |>
    dplyr::pull(.data$phenotype) |>
    unique()

  all_hfdr <- list()
  all_raw  <- list()

  for (grp in diseases) {
    sub <- df |> dplyr::filter(.data$phenotype %in% c("Control", grp))

    idps <- unique(sub$IDP)

    p_values <- numeric(length(idps))
    coefs <- numeric(length(idps))
    tstats <- numeric(length(idps))
    measures <- if ("measure" %in% names(sub)) character(length(idps)) else rep(NA_character_, length(idps))

    for (k in seq_along(idps)) {
      idp <- idps[[k]]
      sidp <- sub |> dplyr::filter(.data$IDP == idp)

      if ("measure" %in% names(sidp)) measures[[k]] <- unique(sidp$measure)[1]

      # drop rows with missing outcome or covariates
      sidp2 <- sidp |> dplyr::filter(!is.na(.data[[outcome_col]]),
                                     !is.na(.data$age_yrs),
                                     !is.na(.data$sex_fac),
                                     !is.na(.data$phenotype_fac))

      # guard: too few observations
      if (nrow(sidp2) < 10) {
        p_values[[k]] <- NA_real_
        coefs[[k]] <- NA_real_
        tstats[[k]] <- NA_real_
        next
      }

      fml <- stats::as.formula(paste0(outcome_col, " ~ phenotype_fac + age_yrs + sex_fac"))
      fit <- stats::lm(fml, data = sidp2)
      sm <- summary(fit)

      # coefficient for phenotype_fac<grp> is row 2 when Control is ref
      coef_row <- 2
      p_values[[k]] <- sm$coefficients[coef_row, 4]
      coefs[[k]] <- sm$coefficients[coef_row, 1]
      tstats[[k]] <- sm$coefficients[coef_row, 3]
    }

    outcomes <- data.frame(
      Dataset = dataset,
      Disease = grp,
      IDP = idps,
      measure = measures,
      pval = p_values,
      Coefficient = coefs,
      T_stat = tstats,
      modality = vapply(strsplit(idps, "_", fixed = TRUE), `[[`, character(1), 1),
      stringsAsFactors = FALSE
    )

    all_raw[[grp]] <- outcomes

    # HFDR (parent=modality)
    hf <- hfdr_two_level(outcomes, alpha = alpha, parent_col = "modality", p_col = "pval", parent_method = "bonf_minp")
    out_hf <- hf$results |>
      dplyr::rename(parent = .data$parent,
                    parent_p = .data$parent_p,
                    parent_q = .data$parent_q,
                    reject_parent = .data$reject_parent)

    all_hfdr[[grp]] <- out_hf

    # write per-disease outputs
    out_subdir <- file.path(out_dir, dataset, grp)
    dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

    readr::write_csv(outcomes, file.path(out_subdir, "raw_outcomes.csv"))
    readr::write_csv(out_hf, file.path(out_subdir, "hfdr_results.csv"))
    readr::write_csv(hf$parent_table, file.path(out_subdir, "hfdr_parent_table.csv"))

    # plot
    plot_path <- file.path(out_subdir, "hfdr_summary.png")
    make_hfdr_summary_plot_labeled(out_hf, dataset, grp, plot_path)
  }

  # combined outputs
  raw_all <- dplyr::bind_rows(all_raw)
  hfdr_all <- dplyr::bind_rows(all_hfdr)
  readr::write_csv(raw_all, file.path(out_dir, dataset, "raw_outcomes_ALL.csv"))
  readr::write_csv(hfdr_all, file.path(out_dir, dataset, "hfdr_results_ALL.csv"))

  hfdr_all
}
