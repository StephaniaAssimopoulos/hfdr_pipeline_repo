#' @keywords internal
parent_p_bonf_minp <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & !is.na(pvals)]
  if (length(pvals) == 0) return(NA_real_)
  min(1, length(pvals) * min(pvals))
}

#' Make HFDR summary plot with labels for top discoveries
#' @keywords internal
make_hfdr_summary_plot_labeled <- function(df, dataset, disease, out_path,
                                          label_top_n = 10,
                                          status_colors = c(
                                            "Not tested" = "grey80",
                                            "Tested (NS)" = "#377eb8",
                                            "HFDR significant" = "#e41a1c"
                                          )) {

  plot_df <- df |>
    dplyr::filter(.data$Dataset == dataset, .data$Disease == disease) |>
    dplyr::mutate(
      status = dplyr::case_when(
        .data$tested_child == FALSE ~ "Not tested",
        .data$tested_child == TRUE & .data$reject_child == FALSE ~ "Tested (NS)",
        .data$reject_child == TRUE ~ "HFDR significant",
        TRUE ~ "Not tested"
      )
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df$modality <- factor(plot_df$modality, levels = c("GM", "WM"))

  parent_df <- plot_df |>
    dplyr::group_by(.data$modality) |>
    dplyr::summarise(
      m = dplyr::n(),
      parent_p = parent_p_bonf_minp(.data$pval),
      tested_any = any(.data$tested_child),
      rejected_any = any(.data$reject_child),
      .groups = "drop"
    )

  plot_df <- plot_df |>
    dplyr::group_by(.data$modality) |>
    dplyr::arrange(dplyr::desc(abs(.data$T_stat)), .by_group = TRUE) |>
    dplyr::mutate(IDP_rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # label: top n significant by abs(T)
  sig_lab <- plot_df |>
    dplyr::filter(.data$reject_child == TRUE) |>
    dplyr::slice_max(order_by = abs(.data$T_stat), n = label_top_n, with_ties = FALSE)

  p_parent <- ggplot2::ggplot(parent_df, ggplot2::aes(x = .data$modality, y = 1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$rejected_any), color = "black", width = 0.85, height = 0.75) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "grey95"), guide = "none") +
    ggplot2::geom_text(ggplot2::aes(label = paste0("parent p = ",
                                                   ifelse(is.na(.data$parent_p), "NA",
                                                          format(.data$parent_p, scientific = TRUE, digits = 2)))),
                       size = 4) +
    ggplot2::coord_cartesian(ylim = c(0.7, 1.3), clip = "off") +
    ggplot2::labs(
      title = paste0(dataset, " — Control vs ", disease, " — HFDR summary"),
      subtitle = "Top: parent gatekeeping; Bottom: IDP evidence (children) with HFDR outcome",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5),
      axis.text.x = ggplot2::element_text(face = "bold")
    )

  p_child <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$IDP_rank, y = .data$T_stat, color = .data$status)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::facet_grid(. ~ modality, scales = "free_x", space = "free_x") +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::labs(x = "IDPs (ranked by |T| within modality)", y = "T-statistic", color = "HFDR outcome") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )

  # add labels on child panel if any significant
  if (nrow(sig_lab) > 0) {
    p_child <- p_child +
      ggrepel::geom_text_repel(
        data = sig_lab,
        ggplot2::aes(label = .data$IDP),
        size = 3,
        max.overlaps = Inf
      )
  }

  fig <- p_parent / p_child + patchwork::plot_layout(heights = c(1, 3))
  ggplot2::ggsave(out_path, fig, width = 11, height = 7, dpi = 300)
  invisible(fig)
}
