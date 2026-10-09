# R/plotting.R
# ---------------------------------------------------------------------------
# Figure helpers for the paper.
#
#   plot_qdiabetes_density()       : density of QDiabetes risk with threshold
#   plot_clinical_features()       : 4-panel clinical-feature figure (Figure 1)
#   plot_love_balance()            : Love plot of covariate balance after matching
#   plot_cox_forest()              : forest plot of adjusted hazard ratios (Figure 3)
#   plot_treatment_km()            : supplementary KM plot with risk table
#   replace_risk_table_times()     : use custom times in a KM risk table
#   plot_km_single_outcome()       : KM plot and risk table for one outcome
#   plot_km_complication_panels()  : grouped KM panels for Figure 2
# ---------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(cowplot)
library(purrr)

plot_qdiabetes_density <- function(df,
                                   score_col = "qdiabetes_10y_risk_score",
                                   threshold = 5.6,
                                   out_file = NULL,
                                   width = 18,
                                   height = 10,
                                   dpi = 300) {
  x <- df[[score_col]]
  d <- density(x, na.rm = TRUE)
  y_top <- max(d$y)

  p <- ggplot(df, aes(x = .data[[score_col]])) +
    geom_area(
      stat = "density",
      aes(
        y = after_stat(ifelse(x <= threshold, density, 0)),
        fill = "Atypical type 2 diabetes"
      ),
      alpha = 0.5
    ) +
    geom_area(
      stat = "density",
      aes(
        y = after_stat(ifelse(x > threshold, density, 0)),
        fill = "Typical type 2 diabetes"
      ),
      alpha = 0.5
    ) +
    geom_density(colour = "black", linewidth = 1) +
    annotate(
      "segment",
      x = threshold, xend = threshold,
      y = 0, yend = y_top * 1.09,
      linetype = "dashed",
      linewidth = 1.2
    ) +
    annotate(
      "text",
      x = threshold,
      y = y_top * 1.10,
      label = paste0(
        "NHS Health Check QDiabetes\nscreening threshold (",
        threshold,
        "%)"
      ),
      hjust = 0.5,
      vjust = 0,
      size = 10,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "Atypical type 2 diabetes" = "#009E73",
        "Typical type 2 diabetes"  = "#D55E00"
      ),
      breaks = c(
        "Atypical type 2 diabetes",
        "Typical type 2 diabetes"
      )
    ) +
    labs(
      x = "QDiabetes 10-year risk (%)",
      y = "Density",
      fill = "QDiabetes risk score"
    ) +
    theme_minimal(base_size = 28) +
    scale_y_continuous(limits = c(0, y_top * 1.20), expand = c(0, 0)) +
    theme(
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    coord_cartesian(clip = "off")

  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = dpi)
  }

  p
}



plot_love_balance <- function(m,
                              threshold = 0.1,
                              abs = TRUE,
                              var_order = "unadjusted",
                              out_file = NULL,
                              width = 10,
                              height = 6,
                              dpi = 300) {
  p <- cobalt::love.plot(
    m,
    stats = "mean.diffs",
    abs = abs,
    thresholds = c(m = threshold),
    var.order = var_order,
    drop.distance = TRUE,
    size = 7
  ) +
    theme_minimal(base_size = 26) +
    labs(
      x = "Standardized mean difference",
      y = NULL
    )

  if (!is.null(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  p
}

# 4-panel clinical features violin plot (unmatched cohort)
plot_clinical_features <- function(df,
                                   out_file = "outputs/figures/clinical_features_unmatched.png",
                                   width = 20,
                                   height = 7,
                                   dpi = 300) {
  library(ggplot2)
  library(dplyr)
  library(cowplot)

  cols <- c(low = "#009E73", high = "#D55E00")

  theme_clean <- theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 25, face = "bold"),
      axis.text.y = element_text(size = 22),
      plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 8, r = 16, b = 8, l = 16),
      axis.line.y = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.x = element_blank()
    )

  # Create individual violin plots with low/high risk groups
  create_feature_plot <- function(var, title, ylim = NULL) {
    plot_df <- df %>%
      filter(!is.na(.data[[var]]), !is.na(qdiabetes_risk_cat)) %>%
      mutate(
        risk_group = factor(
          as.character(qdiabetes_risk_cat),
          levels = c("low", "high")
        ),
        # Tiny horizontal offset so the two x-axis labels do not touch.
        risk_group_x = dplyr::if_else(risk_group == "low", 1, 2.12)
      )

    p <- ggplot(
      plot_df,
      aes(x = risk_group_x, y = .data[[var]], fill = risk_group, group = risk_group)
    ) +
      geom_violin(
        colour = "black",
        linewidth = 0.6,
        trim = TRUE,
        alpha = 0.6
      ) +
      geom_boxplot(
        width = 0.12,
        colour = "black",
        linewidth = 0.7,
        outlier.shape = NA,
        median.linewidth = 0
      ) +
      stat_summary(
        fun.min = median,
        fun.max = median,
        geom = "errorbar",
        width = 0.12,
        colour = "black",
        linewidth = 0.8
      ) +
      scale_fill_manual(
        values = cols,
        name = "QDiabetes risk score",
        labels = c(
          low = "Low diabetes risk score",
          high = "High diabetes risk score"
        )
      ) +
      scale_x_continuous(
        limits = c(0.5, 2.62),
        breaks = c(1, 2.12),
        labels = NULL,
        expand = c(0, 0)
      ) +
      labs(title = title) +
      theme_clean +
      theme(
        axis.text.x = element_blank(),
        plot.margin = margin(t = 8, r = 16, b = 0, l = 16),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.direction = "horizontal"
      )

    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }

    label_strip <- ggplot() +
      annotate(
        "text",
        x = 1,
        y = 1,
        label = "Atypical\ntype 2\ndiabetes",
        colour = unname(cols["low"]),
        size = 22 / ggplot2::.pt,
        fontface = "bold",
        lineheight = 0.9
      ) +
      annotate(
        "text",
        x = 2.12,
        y = 1,
        label = "Typical\ntype 2\ndiabetes",
        colour = unname(cols["high"]),
        size = 22 / ggplot2::.pt,
        fontface = "bold",
        lineheight = 0.9
      ) +
      scale_x_continuous(limits = c(0.5, 2.62), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0.5, 1.5), clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t = 0, r = 16, b = 4, l = 16))

    cowplot::plot_grid(
      p + theme(legend.position = "none"),
      label_strip,
      ncol = 1,
      rel_heights = c(1, 0.24),
      align = "v",
      axis = "lr"
    )
  }

  # Create 4 plots
  p1 <- create_feature_plot("dm_diag_age", "Age (years)")
  p2 <- create_feature_plot("prebmi", "BMI (kg/m\u00B2)", ylim = c(20, 50))
  p3 <- create_feature_plot("prehba1c", "HbA1c (mmol/mol)", ylim = c(48, 130))
  p4 <- create_feature_plot("pretghdl_ratio", "TG/HDL ratio", ylim = c(0, 5))

  panel_grid <- cowplot::plot_grid(
    p1,
    p2,
    p3,
    p4,
    ncol = 4,
    align = "hv",
    axis = "tblr"
  )

  combined <- panel_grid

  if (!is.null(out_file)) {
    ggsave(out_file, combined, width = width, height = height, dpi = dpi)
  }

  combined
}

# ---- Forest plot for adjusted Cox hazard ratios -------------------------------
plot_cox_forest <- function(cox_results,
                            outcome_categories = list(
                              "Microvascular" = c(
                                "retinopathy_severe",
                                "neuropathy_severe",
                                "nephropathy_severe"
                              ),
                              "Macrovascular" = c(
                                "mi_fatal_nonfatal",
                                "stroke_fatal_nonfatal",
                                "hf_fatal_nonfatal"
                              ),
                              "Acute" = c("dka", "hypoglycaemia", "emerg_hosp")
                            ),
                            point_colour = "#1f77b4",
                            outcome_labels = list(
                              retinopathy_severe = "Severe retinopathy",
                              neuropathy_severe = "Severe neuropathy",
                              nephropathy_severe = "Severe nephropathy",
                              mi_fatal_nonfatal = "Myocardial infarction",
                              stroke_fatal_nonfatal = "Stroke",
                              hf_fatal_nonfatal = "Heart failure",
                              dka = "Diabetic ketoacidosis",
                              hypoglycaemia = "Hypoglycaemia",
                              emerg_hosp = "Emergency hospitalisation"
                            ),
                            out_file = NULL,
                            width = 22,
                            height = 15,
                            dpi = 300,
                            base_size = 24,
                            show_group_counts = TRUE,
                            axis_text_x_size = NULL,
                            axis_title_x_size = base_size,
                            axis_tick_length = 4,
                            axis_tick_linewidth = 0.4,
                            gap_after_block = 0.9,
                            gap_header_to_row = -0.3,
                            point_size = 6.5,
                            errorbar_linewidth = 2.8,
                            hr_header = "Adjusted hazard ratio\n(95% CI)",
                            hr_column_hjust = 0.5,
                            label_plot_gap = 0.22,
                            visible_categories = names(outcome_categories),
                            x_breaks = c(0.65, 0.8, 1.0, 1.2, 1.4, 1.6),
                            x_limits = NULL) {
  library(dplyr)
  library(ggplot2)

  plot_data <- cox_results %>%
    dplyr::filter(!is.na(HR), !is.na(LCL), !is.na(UCL))

  # Keep one row per outcome (in case multiple models are present)
  plot_data <- plot_data %>%
    dplyr::group_by(outcome) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  has_grp <- all(c("n_low", "ev_low", "n_high", "ev_high") %in% names(plot_data))

  fmt_int <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")
  fmt_grp <- function(ev, n) {
    if (is.na(ev) || is.na(n)) {
      return(NA_character_)
    }
    sprintf("%s / %s", fmt_int(ev), fmt_int(n))
  }

  # Vertical layout spacing. gap_header_to_row is negative so the first row of a
  # block sits closer to its category header than rows are to each other.
  # gap_header_to_row is supplied as a function argument (space between a
  # category header and its first row; negative pulls the first row closer)
  row_step <- 1.05 # spacing between successive complication rows
  # gap_after_block is supplied as a function argument (extra space after each block)

  forest_df <- data.frame()
  zebra_raw <- list()
  zebra_toggle <- FALSE
  # Seed so the first block header gets the same extra space as later ones.
  # (Default gap_after_block = 0.9 keeps the original layout unchanged.)
  y_counter <- gap_after_block - 0.9

  for (category in names(outcome_categories)) {
    outcomes <- outcome_categories[[category]]
    outcomes <- outcomes[outcomes %in% unique(plot_data$outcome)]
    if (length(outcomes) == 0) next

    # Category header
    y_counter <- y_counter + 1
    forest_df <- dplyr::bind_rows(forest_df, data.frame(
      outcome = NA_character_,
      outcome_category = category,
      HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_,
      y_pos = y_counter, outcome_label = category, is_category = TRUE,
      low_label = NA_character_, high_label = NA_character_, hr_label = NA_character_,
      stringsAsFactors = FALSE
    ))
    y_counter <- y_counter + gap_header_to_row

    for (outcome in outcomes) {
      row_data <- plot_data %>% dplyr::filter(outcome == !!outcome)
      if (nrow(row_data) == 0) next

      low_str <- if (has_grp) fmt_grp(row_data$ev_low[1], row_data$n_low[1]) else NA_character_
      high_str <- if (has_grp) fmt_grp(row_data$ev_high[1], row_data$n_high[1]) else NA_character_

      y_counter <- y_counter + row_step
      hr_str <- sprintf("%.2f (%.2f\u2013%.2f)", row_data$HR, row_data$LCL, row_data$UCL)

      forest_df <- dplyr::bind_rows(forest_df, data.frame(
        outcome = row_data$outcome,
        outcome_category = category,
        HR = row_data$HR,
        LCL = row_data$LCL,
        UCL = row_data$UCL,
        p = row_data$p,
        y_pos = y_counter,
        outcome_label = outcome_labels[[outcome]] %||% outcome,
        is_category = FALSE,
        low_label = low_str,
        high_label = high_str,
        hr_label = hr_str,
        stringsAsFactors = FALSE
      ))

      zebra_toggle <- !zebra_toggle
      zebra_raw[[length(zebra_raw) + 1]] <- list(
        ymin = y_counter - row_step / 2,
        ymax = y_counter + row_step / 2,
        stripe = zebra_toggle,
        outcome_category = category
      )
    }
    y_counter <- y_counter + gap_after_block
  }

  y_max_total <- max(y_counter) + 1
  forest_df <- forest_df %>%
    dplyr::mutate(y_pos = y_max_total - y_pos)

  zebra_df <- if (length(zebra_raw) > 0) {
    do.call(dplyr::bind_rows, lapply(zebra_raw, function(z) {
      data.frame(
        ymin = y_max_total - z$ymax,
        ymax = y_max_total - z$ymin,
        stripe = z$stripe,
        outcome_category = z$outcome_category,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    NULL
  }

  visible_categories <- intersect(visible_categories, names(outcome_categories))
  plot_data_only <- forest_df %>%
    dplyr::filter(
      !is_category,
      !is.na(HR),
      outcome_category %in% visible_categories
    )

  # X-axis range on the HR scale, derived from the requested breaks.
  breaks_exp <- x_breaks
  breaks_log <- log(breaks_exp)
  if (is.null(x_limits)) {
    xlim_l <- log(min(breaks_exp))
    xlim_r <- log(max(breaks_exp))
    max_bar_right <- suppressWarnings(max(log(plot_data_only$UCL), na.rm = TRUE))
    if (is.finite(max_bar_right)) xlim_r <- max(xlim_r, max_bar_right + 0.04)
  } else {
    stopifnot(
      length(x_limits) == 2,
      all(is.finite(x_limits)),
      all(x_limits > 0),
      x_limits[2] > x_limits[1]
    )
    xlim_l <- log(x_limits[1])
    xlim_r <- log(x_limits[2])
  }
  x_range_forest <- xlim_r - xlim_l

  char_w <- base_size / 2.835 / 25.4 * 0.55
  lbl_max <- max(nchar(c(unlist(outcome_labels), names(outcome_categories))))
  hr_chars <- max(
    nchar(unlist(strsplit(hr_header, "\n", fixed = TRUE))),
    nchar("1.00 (1.00\u20131.00)")
  )
  # Half-width (inches) of the widest "ev / n" string, for column placement
  val_chars <- max(
    nchar(stats::na.omit(c(forest_df$low_label, forest_df$high_label))),
    nchar("High diabetes risk score (matched)")
  )
  val_half <- (val_chars / 2) * char_w

  # Right-hand margin (HR column)
  right_outer_pad_in <- 0.4
  m_r_in <- hr_chars * char_w + right_outer_pad_in
  m_r_pt <- ceiling(m_r_in * 72)

  # Left-hand columns: minimal spacing to compress layout
  left_outer_pad_in <- 0.15
  if (isTRUE(show_group_counts)) {
    off_high_in <- -0.15
    off_low_in <- off_high_in + 1.0 * val_half + 0.02
    off_label_in <- off_low_in + 0.5 * val_half + 0.02
  } else {
    off_high_in <- NA_real_
    off_low_in <- NA_real_
    off_label_in <- -0.15
  }

  # Gap between high-risk column and left edge of forest plot (log units)
  high_plot_gap <- label_plot_gap

  # Base left margin (columns + label text + outer pad), excluding the gap
  base_m_l_in <- off_label_in + lbl_max * char_w + left_outer_pad_in

  # Solve for inches-per-log-unit accounting for the extra gap shift so the
  # left margin expands to fit the shifted label column (prevents clipping).
  ipl <- (width - base_m_l_in - m_r_in) / (x_range_forest + high_plot_gap)
  m_l_in <- base_m_l_in + high_plot_gap * ipl
  m_l_pt <- ceiling(m_l_in * 72)

  # Width available for the actual forest plot panel
  panel_in <- max(width - m_l_in - m_r_in, 5)

  col_label <- xlim_l - off_label_in / ipl - high_plot_gap
  if (isTRUE(show_group_counts)) {
    col_high <- xlim_l - off_high_in / ipl - high_plot_gap
    col_low <- xlim_l - off_low_in / ipl - high_plot_gap
  }

  # Right-hand HR column, placed just past the (extended) panel edge.
  hr_col_width <- (hr_chars * char_w) / ipl
  col_hr_left <- xlim_r + 0.02
  col_hr_anchor <- col_hr_left + hr_col_width * hr_column_hjust

  y_top <- max(forest_df$y_pos, na.rm = TRUE) + 1.0

  scale_x_min <- col_label - 0.02
  scale_x_max <- col_hr_left + hr_col_width * 1.08 + 0.01
  band_x_max <- col_hr_left + hr_col_width * 1.08 + 0.01
  band_x_min <- col_label - (lbl_max * char_w) / ipl - 0.005
  cat_df <- forest_df %>%
    dplyr::filter(is_category, outcome_category %in% visible_categories)
  out_df <- forest_df %>%
    dplyr::filter(!is_category, outcome_category %in% visible_categories)

  group_count_layers <- if (isTRUE(show_group_counts)) {
    list(
      annotate("text",
        x = col_low, y = y_top, label = "Atypical type 2\ndiabetes\nEvents / N",
        hjust = 0.5, vjust = 0.5, size = base_size / 2.835,
        fontface = "bold", colour = "grey15", lineheight = 0.9
      ),
      geom_text(
        data = out_df,
        aes(x = col_low, y = y_pos, label = low_label),
        hjust = 0.5, vjust = 0.5, size = (base_size - 2) / 2.835,
        colour = "grey25", inherit.aes = FALSE
      ),
      annotate("text",
        x = col_high, y = y_top, label = "Typical type 2\ndiabetes\nEvents / N",
        hjust = 0.5, vjust = 0.5, size = base_size / 2.835,
        fontface = "bold", colour = "grey15", lineheight = 0.9
      ),
      geom_text(
        data = out_df,
        aes(x = col_high, y = y_pos, label = high_label),
        hjust = 0.5, vjust = 0.5, size = (base_size - 2) / 2.835,
        colour = "grey25", inherit.aes = FALSE
      )
    )
  } else {
    list()
  }

  p <- ggplot(forest_df, aes(x = log(HR), y = y_pos)) +
    geom_blank(
      data = forest_df,
      aes(x = 0, y = y_pos),
      inherit.aes = FALSE
    )

  if (!is.null(zebra_df)) {
    stripe_rows <- zebra_df[
      zebra_df$stripe & zebra_df$outcome_category %in% visible_categories, ,
      drop = FALSE
    ]
    if (nrow(stripe_rows) > 0) {
      p <- p + geom_rect(
        data = stripe_rows,
        aes(xmin = band_x_min, xmax = band_x_max, ymin = ymin, ymax = ymax),
        fill = "#f6f6f6", colour = NA, inherit.aes = FALSE
      )
    }
  }

  p <- p +
    annotate(
      "segment",
      x = band_x_min, xend = band_x_max,
      y = y_top - 0.55, yend = y_top - 0.55,
      colour = "grey75", linewidth = 0.5
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = "grey55") +
    geom_errorbarh(
      data = plot_data_only,
      aes(xmin = log(LCL), xmax = log(UCL)),
      linewidth = errorbar_linewidth, height = 0.0, colour = point_colour
    ) +
    geom_point(data = plot_data_only, size = point_size, shape = 16, colour = point_colour) +

    # Category headers + outcome labels
    geom_text(
      data = cat_df,
      aes(x = col_label, y = y_pos, label = outcome_label),
      hjust = 1, vjust = 0.5, size = base_size / 2.835,
      fontface = "bold", colour = "grey15", inherit.aes = FALSE
    ) +
    geom_text(
      data = out_df,
      aes(x = col_label, y = y_pos, label = outcome_label),
      hjust = 1, vjust = 0.5, size = base_size / 2.835,
      fontface = "plain", colour = "grey15", inherit.aes = FALSE
    ) +

    # Low-risk column header + values
    # High-risk-matched column header + values
    group_count_layers +

    # Hazard ratio column (right of plot)
    # Hazard ratio column (right of plot)
    # Hazard ratio column (right of plot)
    annotate(
      "text",
      x = col_hr_anchor,
      y = y_top,
      label = hr_header,
      hjust = hr_column_hjust,
      vjust = 0.5,
      size = base_size / 2.835,
      fontface = "bold",
      colour = "grey15",
      lineheight = 0.9
    ) +
    geom_text(
      data = plot_data_only,
      aes(x = col_hr_anchor, y = y_pos, label = hr_label),
      hjust = hr_column_hjust,
      vjust = 0.5,
      size = (base_size - 2) / 2.835,
      colour = "grey25",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(band_x_min, band_x_max),
      breaks = breaks_log,
      labels = format(breaks_exp, nsmall = 1),
      name = "Adjusted hazard ratio (95% CI)",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = NULL,
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.title.x = element_text(
        margin = margin(t = 12),
        colour = "grey15",
        size = axis_title_x_size
      ),
      axis.title.y = element_blank(),
      axis.text.x = element_text(colour = "grey25", size = axis_text_x_size),
      axis.text.y = element_blank(),
      axis.ticks.x = element_line(colour = "grey60", linewidth = axis_tick_linewidth),
      axis.ticks.length.x = unit(axis_tick_length, "pt"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(l = max(m_l_pt - 25, 5), r = m_r_pt, t = 30, b = 20)
    ) +
    coord_cartesian(xlim = c(xlim_l, xlim_r), clip = "off")

  if (!is.null(out_file)) {
    ggsave(out_file,
      plot = p, width = width, height = height, dpi = dpi,
      limitsize = FALSE
    )
  }

  p
}


# Replace the default number-at-risk times in a survminer risk table.
replace_risk_table_times <- function(km_obj, fit, risk_times) {
  risk_summary <- summary(fit, times = risk_times, extend = TRUE)
  risk_counts <- data.frame(
    time = as.numeric(risk_summary$time),
    strata = as.character(risk_summary$strata),
    n.risk = round(as.numeric(risk_summary$n.risk)),
    stringsAsFactors = FALSE
  )

  table_data <- km_obj$table$data
  table_data_new <- purrr::map_dfr(
    seq_len(nrow(risk_counts)),
    function(i) {
      strata_i <- risk_counts$strata[i]
      table_strata <- as.character(table_data$strata)
      row_i <- match(strata_i, table_strata)

      if (is.na(row_i)) {
        group_i <- if (grepl("low", strata_i, ignore.case = TRUE)) {
          "low"
        } else {
          "high"
        }
        row_i <- which(grepl(group_i, table_strata, ignore.case = TRUE))[1]
      }

      if (is.na(row_i)) {
        stop("Could not match a survival stratum to the risk table.", call. = FALSE)
      }

      template_row <- table_data[row_i, , drop = FALSE]
      template_row$time <- risk_counts$time[i]
      template_row$n.risk <- risk_counts$n.risk[i]

      if ("llabels" %in% names(template_row)) {
        template_row$llabels <- as.character(risk_counts$n.risk[i])
      }
      if ("label" %in% names(template_row)) {
        template_row$label <- as.character(risk_counts$n.risk[i])
      }

      template_row
    }
  )

  km_obj$table$data <- table_data_new
  for (i in seq_along(km_obj$table$layers)) {
    if (
      inherits(km_obj$table$layers[[i]]$geom, "GeomText") &&
        !is.null(km_obj$table$layers[[i]]$data)
    ) {
      km_obj$table$layers[[i]]$data <- table_data_new
    }
  }

  km_obj
}


# Kaplan-Meier plot and number-at-risk table used for supplementary outcomes.
plot_treatment_km <- function(
    surv_df,
    title,
    out_file = NULL,
    xlim,
    x_breaks,
    y_limits,
    y_breaks,
    break_time_by = 2,
    risk_times = NULL,
    ahr_label = NULL,
    ahr_pos = "bottomright",
    width = 12,
    height = 9,
    dpi = 600,
    base_size = 25,
    risk_table_base_size = 22,
    risk_table_height = 0.26,
    show_y_axis_label = TRUE,
    x_pad_left = 0.10,
    x_pad_right = 0.05) {
  if (!ahr_pos %in% c("topleft", "bottomright")) {
    stop("`ahr_pos` must be either 'topleft' or 'bottomright'.", call. = FALSE)
  }

  plot_df <- surv_df %>%
    dplyr::mutate(group = factor(group, levels = c("low", "high")))

  x_limits_padded <- c(
    xlim[1] - x_pad_left * diff(xlim),
    xlim[2] + x_pad_right * diff(xlim)
  )

  fit <- survival::survfit(
    survival::Surv(time_years, event) ~ group,
    data = plot_df
  )

  km <- survminer::ggsurvplot(
    fit,
    data = plot_df,
    fun = "event",
    conf.int = TRUE,
    conf.int.alpha = 0.15,
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.y.text = FALSE,
    risk.table.y.text.col = TRUE,
    xlim = x_limits_padded,
    break.time.by = break_time_by,
    legend.title = "QDiabetes risk score",
    legend.labs = c(
      "Low diabetes risk score",
      "High diabetes risk score (matched)"
    ),
    palette = c("#009E73", "#D55E00"),
    xlab = "Years since diagnosis",
    ylab = if (show_y_axis_label) "Cumulative incidence (%)" else "",
    title = title,
    tables.theme = survminer::theme_cleantable(),
    risk.table.height = risk_table_height,
    size = 0.6
  )

  km$plot <- km$plot +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits_padded,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      labels = scales::percent,
      breaks = y_breaks
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0),
      legend.position = "bottom"
    )

  if (!show_y_axis_label) {
    km$plot <- km$plot +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }

  if (!is.null(ahr_label)) {
    annotation_position <- if (ahr_pos == "topleft") {
      list(
        x = -Inf,
        y = y_limits[2] - 0.02 * diff(y_limits),
        hjust = 0,
        vjust = 1
      )
    } else {
      list(
        x = xlim[2],
        y = y_limits[1] + 0.02 * diff(y_limits),
        hjust = 1,
        vjust = 0
      )
    }

    km$plot <- km$plot +
      ggplot2::annotate(
        "text",
        x = annotation_position$x,
        y = annotation_position$y,
        label = ahr_label,
        hjust = annotation_position$hjust,
        vjust = annotation_position$vjust,
        size = base_size / 3.4,
        fontface = "plain"
      )
  }

  km$table <- km$table +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits_padded,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    survminer::theme_cleantable() +
    ggplot2::theme(
      text = ggplot2::element_text(size = risk_table_base_size),
      plot.title = ggplot2::element_text(
        size = risk_table_base_size,
        hjust = 0
      ),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  risk_table_text_size <- risk_table_base_size / 2.845276
  for (i in seq_along(km$table$layers)) {
    if (inherits(km$table$layers[[i]]$geom, "GeomText")) {
      km$table$layers[[i]]$aes_params$size <- risk_table_text_size
    }
  }

  if (!is.null(risk_times)) {
    km <- replace_risk_table_times(km, fit, risk_times)
  }

  final_plot <- cowplot::plot_grid(
    km$plot + ggplot2::theme(legend.position = "none"),
    km$table,
    ncol = 1,
    rel_heights = c(1 - risk_table_height, risk_table_height),
    align = "v",
    axis = "lr"
  )

  if (!is.null(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot = final_plot,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  final_plot
}


# =============================================================================
# Kaplan-Meier cumulative-incidence plots for the matched cohort
# =============================================================================

km_outcome_settings <- list(
  # microvascular
  retinopathy_severe = list(label = "Severe retinopathy", y_max = 0.12),
  neuropathy_severe = list(label = "Severe neuropathy", y_max = 0.12),
  nephropathy_severe = list(label = "Severe nephropathy", y_max = 0.12),
  # macrovascular (fatal or non-fatal)
  mi_fatal_nonfatal = list(label = "Myocardial infarction", y_max = 0.12),
  stroke_fatal_nonfatal = list(label = "Stroke", y_max = 0.12),
  hf_fatal_nonfatal = list(label = "Heart failure", y_max = 0.12),

  # mortality
  cv_death_primary = list(label = "Cardiovascular mortality", y_max = 0.12),
  cancer_death_primary = list(label = "Cancer mortality", y_max = 0.12),

  # Secondary outcomes clinical outcomes
  dka = list(label = "Diabetic ketoacidosis", y_max = 0.12),
  hypoglycaemia = list(label = "Hypoglycaemia", y_max = 0.12),

  # treatment outcomes
  insulin_gt1y = list(
    label = "Insulin initiation",
    y_max = 0.16,
    x_min = 1,
    x_breaks = c(1, 3, 5, 7, 9)
  ),
  first_line_treatment = list(label = "First-line therapy initiation", y_max = 1.00),
  second_line_treatment = list(label = "Second-line therapy initiation", y_max = 0.70)
)



y_break_step <- function(y_max) {
  ymax_pct <- y_max * 100

  if (ymax_pct < 4) {
    step <- 1
  } else if (ymax_pct <= 10) {
    step <- 2
  } else if (ymax_pct > 10 && ymax_pct <= 28) {
    step <- 4
  } else if (ymax_pct > 28 && ymax_pct <= 55) {
    step <- 10
  } else if (ymax_pct > 55 && ymax_pct <= 100) {
    step <- 20
  } else {
    step <- 20
  }
  step / 100
}



# High-quality PNG save
save_hi_res <- function(out_file, plot, width, height, dpi) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(out_file,
      plot = plot, width = width, height = height,
      dpi = dpi, device = ragg::agg_png
    )
  } else {
    ggplot2::ggsave(out_file,
      plot = plot, width = width, height = height,
      dpi = dpi, type = "cairo"
    )
  }
}

extract_legend_safe <- function(p) {
  p <- p + theme(legend.position = "bottom")


  guides <- tryCatch(
    cowplot::get_plot_component(p, pattern = "guide-box", return_all = TRUE),
    error = function(e) list()
  )
  guides <- Filter(function(g) !is.null(g) && !inherits(g, "zeroGrob"), guides)
  if (length(guides) > 0) {
    return(guides[[1]])
  }

  lg <- tryCatch(cowplot::get_legend(p), error = function(e) NULL)
  if (!is.null(lg) && !inherits(lg, "zeroGrob")) {
    return(lg)
  }

  grid::nullGrob()
}


plot_km_single_outcome <- function(df,
                                   title = NULL,
                                   y_max = 0.15,
                                   cols = c(low = "#009E73", high = "#D55E00"),
                                   legend_labels = c(
                                     "Low diabetes risk score",
                                     "High diabetes risk score (matched)"
                                   ),
                                   base_size = 32,
                                   x_min = 0,
                                   x_max = 10,
                                   x_break = 2,
                                   x_breaks = NULL,
                                   curve_linewidth = 0.2,
                                   conf_int_alpha = 0.15,
                                   risk_table_size_offset = 2) {
  df <- df %>%
    mutate(group = factor(group, levels = c("low", "high")))

  fit <- survfit(
    Surv(time_years, event) ~ group,
    data = df
  )


  rt_base <- max(base_size - risk_table_size_offset, 1)

  p <- ggsurvplot(
    fit,
    data = df,
    fun = "event",
    conf.int = TRUE,
    conf.int.alpha = conf_int_alpha,
    size = curve_linewidth,
    censor = FALSE,
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.y.text = FALSE,
    risk.table.y.text.col = TRUE,
    fontsize = rt_base / 3.2,
    xlim = c(x_min, x_max),
    break.time.by = x_break,
    legend.title = "QDiabetes-defined type 2 diabetes phenotype",
    legend.labs = legend_labels,
    palette = unname(cols[c("low", "high")]),
    xlab = "Years since diagnosis",
    ylab = "Cumulative incidence (%)",
    title = title,
    tables.theme = theme_cleantable(),
    risk.table.height = 0.30
  )


  text_layer_size <- rt_base / 3.1

  y_step <- y_break_step(y_max)

  p$plot <- p$plot +
    theme_bw(base_size = base_size) +
    scale_y_continuous(
      limits = c(0, y_max),
      labels = scales::percent,
      breaks = seq(0, y_max, by = y_break_step(y_max))
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.title.x = element_text(size = base_size),
      axis.text.x = element_text(size = base_size),
      axis.title.y = element_text(size = base_size * 1.10, face = "plain"),
      axis.text.y = element_text(size = base_size * 1.05, face = "plain"),
      plot.title = element_text(size = base_size + 2, hjust = 0),
      legend.text = element_text(size = base_size * 1.30),
      legend.title = element_text(size = base_size * 1.05, face = "bold"),
      legend.key.size = unit(base_size * 1.4, "pt"),
      plot.margin = margin(
        t = base_size * 0.9, r = base_size * 0.4,
        b = base_size * 0.2, l = base_size * 0.4
      )
    )


  for (i in seq_along(p$plot$layers)) {
    lyr <- p$plot$layers[[i]]
    if (inherits(lyr$geom, "GeomStep") || inherits(lyr$geom, "GeomLine")) {
      lyr$aes_params$linewidth <- curve_linewidth
      lyr$aes_params$size <- curve_linewidth
      p$plot$layers[[i]] <- lyr
    }
  }

  if (!is.null(x_breaks)) {
    p$plot <- p$plot +
      scale_x_continuous(breaks = x_breaks) +
      coord_cartesian(xlim = c(x_min, x_max))
  }

  # Keep risk table minimal and uncluttered.
  p$table <- p$table +
    theme_cleantable() +
    theme(
      plot.title = element_text(
        size = text_layer_size * ggplot2::.pt,
        face = "plain",
        hjust = 0,
        margin = margin(0, 0, 2, 0)
      ),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 55)
    ) +
    coord_cartesian(
      xlim = c(x_min - 0.2, x_max),
      clip = "off"
    )

  if (!is.null(x_breaks)) {
    sfit <- summary(fit, times = x_breaks, extend = TRUE)
    sdf <- data.frame(
      strata = as.character(sfit$strata),
      time = as.numeric(sfit$time),
      nrisk = round(as.numeric(sfit$n.risk)),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(p$table$layers)) {
      if (!inherits(p$table$layers[[i]]$geom, "GeomText")) next
      ld <- p$table$layers[[i]]$data
      if (is.null(ld) || !"strata" %in% names(ld) || !"time" %in% names(ld)) next
      rows <- lapply(seq_len(nrow(sdf)), function(r) {
        tmpl <- ld[match(sdf$strata[r], as.character(ld$strata)), , drop = FALSE]
        tmpl$time <- sdf$time[r]
        if ("n.risk" %in% names(tmpl)) tmpl$n.risk <- sdf$nrisk[r]
        if ("llabels" %in% names(tmpl)) tmpl$llabels <- as.character(sdf$nrisk[r])
        if ("label" %in% names(tmpl)) tmpl$label <- as.character(sdf$nrisk[r])
        tmpl
      })
      p$table$layers[[i]]$data <- do.call(rbind, rows)
    }

    p$table <- p$table + scale_x_continuous(breaks = x_breaks)
  }

  risk_row_compact <- 0.72
  marker_points <- data.frame(
    x = numeric(),
    y = numeric(),
    strata = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(p$table$layers)) {
    if (inherits(p$table$layers[[i]]$geom, "GeomText")) {
      layer_data <- p$table$layers[[i]]$data
      if (!is.null(layer_data) && "y" %in% names(layer_data)) {
        y_vals <- suppressWarnings(as.numeric(layer_data$y))
        if (length(y_vals) > 0 && any(is.finite(y_vals))) {
          y_center <- mean(y_vals[is.finite(y_vals)], na.rm = TRUE)
          layer_data$y <- y_center + (y_vals - y_center) * risk_row_compact
        }
      }
      label_field <- if ("label" %in% names(layer_data)) {
        "label"
      } else if ("llabels" %in% names(layer_data)) {
        "llabels"
      } else {
        NULL
      }

      if (!is.null(label_field)) {
        lbl <- as.character(layer_data[[label_field]])
        marker_idx <- grepl("^\\s*-\\s*\\d", lbl)
        if (any(marker_idx) && all(c("x", "y") %in% names(layer_data))) {
          strata_vals <- if ("strata" %in% names(layer_data)) {
            as.character(layer_data$strata[marker_idx])
          } else {
            rep(NA_character_, sum(marker_idx))
          }

          marker_points <- rbind(
            marker_points,
            data.frame(
              x = as.numeric(layer_data$x[marker_idx]),
              y = as.numeric(layer_data$y[marker_idx]),
              strata = strata_vals,
              stringsAsFactors = FALSE
            )
          )

          # Remove leading marker dash from first timepoint values.
          lbl[marker_idx] <- sub("^\\s*-\\s*", "", lbl[marker_idx])
        }

        if ("label" %in% names(layer_data)) {
          layer_data$label <- lbl
        }
        if ("llabels" %in% names(layer_data)) {
          layer_data$llabels <- lbl
        }

        p$table$layers[[i]]$data <- layer_data
      }
      if ("size" %in% names(p$table$layers[[i]]$mapping)) {
        p$table$layers[[i]]$mapping$size <- NULL
      }
      p$table$layers[[i]]$aes_params$size <- text_layer_size
      p$table$layers[[i]]$aes_params$fontface <- "plain"
    }
  }

  if (nrow(marker_points) > 0) {
    marker_points <- marker_points[
      is.finite(marker_points$x) & is.finite(marker_points$y), ,
      drop = FALSE
    ]
    marker_points <- marker_points[!duplicated(marker_points$y), , drop = FALSE]
    marker_points <- marker_points[order(marker_points$y, decreasing = TRUE), , drop = FALSE]

    marker_cols <- ifelse(
      grepl("low", marker_points$strata, ignore.case = TRUE),
      unname(cols["low"]),
      ifelse(
        grepl("high", marker_points$strata, ignore.case = TRUE),
        unname(cols["high"]),
        rep(c(unname(cols["low"]), unname(cols["high"])), length.out = nrow(marker_points))
      )
    )


    dash_start <- if (!is.null(x_breaks)) -0.55 else -1.45
    dash_end <- if (!is.null(x_breaks)) -0.12 else -0.40

    for (k in seq_len(nrow(marker_points))) {
      x_mid <- marker_points$x[k]
      p$table <- p$table +
        ggplot2::annotate(
          "segment",
          x = x_mid + dash_start,
          xend = x_mid + dash_end,
          y = marker_points$y[k],
          yend = marker_points$y[k],
          colour = marker_cols[k],
          linewidth = 2.4,
          lineend = "round"
        )
    }

    p$table <- p$table +
      ggplot2::coord_cartesian(
        xlim = if (!is.null(x_breaks)) c(x_min - 0.7, x_max) else c(x_min - 1.6, x_max),
        clip = "off"
      )
  }

  p
}


km_complication_panels <- list(
  "A) Microvascular" = c(
    "retinopathy_severe",
    "neuropathy_severe",
    "nephropathy_severe"
  ),
  "B) Macrovascular" = c(
    "mi_fatal_nonfatal",
    "stroke_fatal_nonfatal",
    "hf_fatal_nonfatal"
  )
)


plot_km_complication_panels <- function(surv_list,
                                        groups,
                                        cols = c(low = "#009E73", high = "#D55E00"),
                                        legend_labels = c(
                                          "Low diabetes risk score",
                                          "High diabetes risk score (matched)"
                                        ),
                                        base_size = 18,
                                        ncol_each = 3,
                                        out_file = NULL,
                                        width = 18,
                                        height = 17,
                                        dpi = 600,
                                        curve_linewidth = 0.2,
                                        conf_int_alpha = 0.15,
                                        risk_table_size_offset = 2,
                                        y_max_overrides = NULL) {
  library(cowplot)
  library(purrr)
  library(ggplot2)

  blank_panel <- ggplot() +
    theme_void()

  group_panels <- list()
  row_weights <- c()

  for (group_name in names(groups)) {
    outcomes <- groups[[group_name]]
    outcomes <- outcomes[outcomes %in% names(surv_list)]

    plots <- purrr::imap(outcomes, function(nm, idx) {
      y_use <- km_outcome_settings[[nm]][["y_max"]]
      if (!is.null(y_max_overrides) && nm %in% names(y_max_overrides)) {
        y_use <- as.numeric(y_max_overrides[[nm]])
      }

      full_plot <- plot_km_single_outcome(
        df = surv_list[[nm]],
        title = km_outcome_settings[[nm]][["label"]],
        y_max = y_use,
        cols = cols,
        legend_labels = legend_labels,
        base_size = base_size,
        x_min = km_outcome_settings[[nm]][["x_min"]] %||% 0,
        x_max = km_outcome_settings[[nm]][["x_max"]] %||% 10,
        x_break = km_outcome_settings[[nm]][["x_break"]] %||% 2,
        x_breaks = km_outcome_settings[[nm]][["x_breaks"]] %||% NULL,
        curve_linewidth = curve_linewidth,
        conf_int_alpha = conf_int_alpha,
        risk_table_size_offset = risk_table_size_offset
      )

      plot_panel <- full_plot$plot + theme(legend.position = "none")
      if (((idx - 1) %% ncol_each) != 0) {
        plot_panel <- plot_panel + theme(axis.title.y = element_blank())
      }

      # Combine plot with risk table
      p <- cowplot::plot_grid(
        plot_panel,
        full_plot$table,
        ncol = 1,
        rel_heights = c(0.70, 0.30),
        align = "v",
        axis = "lr"
      )

      p
    })

    n_rows <- ceiling(length(plots) / ncol_each)
    n_needed <- n_rows * ncol_each
    if (length(plots) < n_needed) {
      plots <- c(plots, rep(list(blank_panel), n_needed - length(plots)))
    }

    grid <- cowplot::plot_grid(
      plotlist = plots,
      ncol = ncol_each,
      align = "hv",
      axis = "bt"
    )

    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        group_name,
        fontface = "bold",
        size = base_size + 8,
        x = 0,
        hjust = 0
      )

    title_spacer <- cowplot::ggdraw()

    panel <- cowplot::plot_grid(
      title,
      title_spacer,
      grid,
      ncol = 1,
      rel_heights = c(0.08, 0.035, 1)
    )

    group_panels[[group_name]] <- panel
    row_weights <- c(row_weights, n_rows + 0.22)
  }

  first_outcome <- unlist(groups)[unlist(groups) %in% names(surv_list)][1]

  legend_plot_full <- plot_km_single_outcome(
    df = surv_list[[first_outcome]],
    title = km_outcome_settings[[first_outcome]][["label"]],
    y_max = km_outcome_settings[[first_outcome]][["y_max"]],
    cols = cols,
    legend_labels = legend_labels,
    base_size = base_size,
    x_min = km_outcome_settings[[first_outcome]][["x_min"]] %||% 0,
    x_max = km_outcome_settings[[first_outcome]][["x_max"]] %||% 10,
    x_break = km_outcome_settings[[first_outcome]][["x_break"]] %||% 2,
    x_breaks = km_outcome_settings[[first_outcome]][["x_breaks"]] %||% NULL,
    curve_linewidth = curve_linewidth,
    conf_int_alpha = conf_int_alpha,
    risk_table_size_offset = risk_table_size_offset
  )

  legend <- extract_legend_safe(legend_plot_full$plot)

  spacer <- cowplot::ggdraw()

  panels_with_spacing <- unlist(
    lapply(group_panels, function(p) list(p, spacer)),
    recursive = FALSE
  )


  panels_with_spacing <- panels_with_spacing[-length(panels_with_spacing)]


  spacing_weight <- 0.15
  weights_with_spacing <- as.vector(rbind(row_weights, rep(spacing_weight, length(row_weights))))
  weights_with_spacing <- weights_with_spacing[-length(weights_with_spacing)]

  final <- cowplot::plot_grid(
    plotlist = c(panels_with_spacing, list(legend)),
    ncol = 1,
    rel_heights = c(weights_with_spacing, 0.15)
  )

  if (!is.null(out_file)) {
    save_hi_res(out_file, final, width = width, height = height, dpi = dpi)
  }

  final
}
