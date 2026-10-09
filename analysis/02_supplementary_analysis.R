# =============================================================================
# 02_supplementary_analysis.R
# -----------------------------------------------------------------------------
# Supplementary analyses for the diabetes discordance paper:
#
#   Figure 1. Acute (DKA, hypoglycaemia) and mortality (cardiovascular,
#      non-cardiovascular) KM curves with the minimally-adjusted HR.
#   Figure 2. Treatment KM curves (insulin, first-line, second-line
#      therapy initiation) with the minimally-adjusted HR.
#   Subgroup analyses for the main complication hazard ratios:
#        - Minimally adjusted HR
#        - Age at diagnosis <50 years
#        - Age at diagnosis >=50 years
#        - Male
#        - Female
#
# =============================================================================

if (!file.exists("R/setup.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}

# ---- Load packages and analysis functions ----------
source("R/setup.R")
source("R/cohort_definition.R")
source("R/calculators/qdiabetes_2018.R")
source("R/calculators/qrisk2.R")
source("R/matching.R")
source("R/tables.R")
source("R/survival_outcomes.R")
source("R/cox_models.R")
source("R/plotting.R")

# ---- Output folders ---------------------------------------------------------
dir.create("outputs/supplementary/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/supplementary/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load CPRD data
# =============================================================================
cprd <- connect_cprd()
at_diag <- load_at_diagnosis_cohort(cprd)
hba1c_clean <- load_hba1c_data(cprd, at_diag$patid)
drug_start_stop <- load_drug_start_stop(cprd)


# =============================================================================
# 2. Build the analysis cohort and compute risk scores
# =============================================================================
# Apply inclusion criteria, derive variables used in the tables/matching/models
cohort <- prepare_main_dataset(
  at_diag_raw = at_diag,
  hba1c_data = hba1c_clean,
  qdiabetes_risk_horizon = 10L,
  qdiabetes_risk_threshold = 5.6,
  drug_classes = drug_start_stop
)

# =============================================================================
# 3. Match and build treatment survival datasets
# =============================================================================
# Matching
matched_df <- match_low_high(
  cohort,
  treated_level = "low",
  control_level = "high",
  ratio = 4
)$matched


# =============================================================================
# 4. Supplementary Figure 1: acute and mortality KM curves (minimally adjusted HR)
# =============================================================================

plot_base_size <- 26
risk_table_base_size <- 20

# A) Acute: DKA, hypoglycaemia. B) Mortality: cardiovascular, non-cardiovascular.
acute_mortality_specs <- c(
  outcome_specs[c("dka", "hypoglycaemia", "cv_death_primary")],
  list(
    non_cv_death = list(
      event_col      = "non_cv_death_date",
      pre_col        = "pre_index_date_death",
      censor_profile = "mortality"
    )
  )
)

surv_list_acute_mortality <- make_surv_list(
  matched_df    = matched_df,
  horizon_years = 10,
  specs         = acute_mortality_specs
)

# Minimally-adjusted HR labels (matching variables only).
acute_mortality_hr_labels <- run_all_cox_matching_only(
  surv_list_acute_mortality
) %>%
  mutate(
    hr_label = sprintf(
      "Minimally adjusted HR (95%% CI):\n%.2f (%.2f-%.2f)",
      HR,
      LCL,
      UCL
    )
  ) %>%
  select(outcome, hr_label) %>%
  tibble::deframe()


km_dka <- plot_treatment_km(
  surv_df = surv_list_acute_mortality$dka,
  title = "Diabetic ketoacidosis",
  xlim = c(0, 10), x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 0.04), y_breaks = seq(0, 0.04, by = 0.01),
  ahr_label = acute_mortality_hr_labels[["dka"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size
)

km_hypo <- plot_treatment_km(
  surv_df = surv_list_acute_mortality$hypoglycaemia,
  title = "Hypoglycaemia",
  xlim = c(0, 10), x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 0.04), y_breaks = seq(0, 0.04, by = 0.01),
  ahr_label = acute_mortality_hr_labels[["hypoglycaemia"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size,
  show_y_axis_label = FALSE
)

km_cvd <- plot_treatment_km(
  surv_df = surv_list_acute_mortality$cv_death_primary,
  title = "Cardiovascular mortality",
  xlim = c(0, 10), x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 0.14), y_breaks = seq(0, 0.14, by = 0.02),
  ahr_label = acute_mortality_hr_labels[["cv_death_primary"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size
)

km_noncvd <- plot_treatment_km(
  surv_df = surv_list_acute_mortality$non_cv_death,
  title = "Non-cardiovascular mortality",
  xlim = c(0, 10), x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 0.14), y_breaks = seq(0, 0.14, by = 0.02),
  ahr_label = acute_mortality_hr_labels[["non_cv_death"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size,
  show_y_axis_label = FALSE
)

panel_label <- function(lbl) {
  cowplot::ggdraw() +
    cowplot::draw_label(lbl, fontface = "bold", size = 30, x = 0, hjust = 0)
}

supp_figure1_legend <- extract_legend_safe(
  ggplot2::ggplot(
    data.frame(
      x = rep(c(0, 1), 2),
      y = rep(c(0, 1), 2),
      phenotype = rep(c("Atypical type 2 diabetes", "Typical type 2 diabetes"), each = 2)
    ),
    ggplot2::aes(x, y, colour = phenotype)
  ) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::scale_colour_manual(
      name = "QDiabetes-defined type 2 diabetes phenotype",
      values = c(
        "Atypical type 2 diabetes" = "#009E73",
        "Typical type 2 diabetes" = "#D55E00"
      ),
      breaks = c("Atypical type 2 diabetes", "Typical type 2 diabetes")
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(title.position = "left", nrow = 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 26, face = "bold", hjust = 0.5),
      legend.text = ggplot2::element_text(size = 26),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.8, "cm")
    )
)

supp_figure1 <- cowplot::plot_grid(
  panel_label("A) Acute"),
  cowplot::plot_grid(km_dka, km_hypo, ncol = 2, align = "hv", axis = "tblr"),
  cowplot::ggdraw(),
  panel_label("B) Mortality"),
  cowplot::plot_grid(km_cvd, km_noncvd, ncol = 2, align = "hv", axis = "tblr"),
  supp_figure1_legend,
  ncol = 1,
  rel_heights = c(0.07, 1, 0.12, 0.07, 1, 0.14)
)

ggplot2::ggsave(
  filename = "outputs/supplementary/figures/supplementary_figure1.png",
  plot = supp_figure1,
  width = 17,
  height = 17,
  dpi = 600
)
supp_figure1


# =============================================================================
# 5. Supplementary Figure 2: treatment-outcome KM curves (minimally adjusted HR)
# =============================================================================
treatment_outcome_names <- c(
  "insulin_gt1y",
  "first_line_treatment",
  "second_line_treatment"
)

treatment_outcome_labels <- c(
  insulin_gt1y = "Insulin initiation",
  first_line_treatment = "First-line therapy initiation",
  second_line_treatment = "Second-line therapy initiation"
)

surv_list_treatment <- make_surv_list(
  matched_df = matched_df,
  horizon_years = 10,
  specs = outcome_specs[treatment_outcome_names]
)

# Minimally-adjusted hazard ratio (matching variables only: age, sex, year).
treatment_cox_results <- run_all_cox_matching_only(surv_list_treatment)

treatment_ahr_labels <- treatment_cox_results %>%
  mutate(
    ahr_label = sprintf(
      " Minimally adjusted HR (95%% CI):\n %.2f (%.2f-%.2f)",
      HR,
      LCL,
      UCL
    )
  ) %>%
  select(outcome, ahr_label) %>%
  tibble::deframe()



km_insulin_initiation_final <- plot_treatment_km(
  surv_df = surv_list_treatment$insulin_gt1y,
  title = "Insulin initiation",
  xlim = c(1, 10),
  x_breaks = seq(1, 9, by = 2),
  y_limits = c(0, 0.16),
  y_breaks = seq(0, 0.16, by = 0.04),
  break_time_by = 2,
  risk_times = seq(1, 9, by = 2),
  ahr_label = treatment_ahr_labels[["insulin_gt1y"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size,
  show_y_axis_label = FALSE
)

km_first_line_initiation_final <- plot_treatment_km(
  surv_df = surv_list_treatment$first_line_treatment,
  title = "First-line therapy initiation",
  xlim = c(0, 10),
  x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 1.00),
  y_breaks = seq(0, 1.00, by = 0.20),
  break_time_by = 2,
  ahr_label = treatment_ahr_labels[["first_line_treatment"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size
)

km_second_line_initiation_final <- plot_treatment_km(
  surv_df = surv_list_treatment$second_line_treatment,
  title = "Second-line therapy initiation",
  xlim = c(0, 10),
  x_breaks = seq(0, 10, by = 2),
  y_limits = c(0, 0.70),
  y_breaks = seq(0, 0.70, by = 0.10),
  break_time_by = 2,
  ahr_label = treatment_ahr_labels[["second_line_treatment"]],
  ahr_pos = "topleft",
  base_size = plot_base_size,
  risk_table_base_size = risk_table_base_size,
  show_y_axis_label = FALSE
)

# First-line, second-line and insulin initiation.
supp_figure2 <- cowplot::plot_grid(
  cowplot::plot_grid(
    km_first_line_initiation_final,
    km_second_line_initiation_final,
    km_insulin_initiation_final,
    ncol = 3,
    align = "hv",
    axis = "tblr"
  ),
  supp_figure1_legend,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

ggplot2::ggsave(
  filename = "outputs/supplementary/figures/supplementary_figure2.png",
  plot = supp_figure2,
  width = 22,
  height = 9,
  dpi = 600
)

supp_figure2


# =============================================================================
# 6. Subgroup analyses for the complication hazard ratios
# =============================================================================
# Analyses:
#   1. Age at diagnosis <50 years
#   2. Age at diagnosis >=50 years
#   3. Male
#   4. Female
#
# Each analysis re-matches atypical vs typical subgroups,
# builds survival datasets, imputes missing baseline HbA1c and fits
# the pooled Cox model.

# -----------------------------------------------------------------------------
# Outcomes included in this supplementary subgroup table
# -----------------------------------------------------------------------------

sensitivity_outcome_labels <- c(
  retinopathy_severe    = "Severe retinopathy",
  neuropathy_severe     = "Severe neuropathy",
  nephropathy_severe    = "Severe nephropathy",
  mi_fatal_nonfatal     = "Myocardial infarction",
  stroke_fatal_nonfatal = "Stroke",
  hf_fatal_nonfatal     = "Heart failure"
)

# Keep only these outcomes from the full outcome_specs list
sensitivity_specs <- outcome_specs[names(sensitivity_outcome_labels)]

# -----------------------------------------------------------------------------
# Run match -> survival datasets -> MICE -> pooled Cox for one analysis cohort
# -----------------------------------------------------------------------------

run_sensitivity_hr <- function(cohort_df,
                               label,
                               horizon_years = 10,
                               m = 10,
                               seed = 1,
                               ratio = 4) {
  matched <- match_low_high(
    cohort_df,
    treated_level = "low",
    control_level = "high",
    ratio = ratio,
    seed = seed
  )$matched

  survival_data <- make_surv_list(
    matched_df = matched,
    horizon_years = horizon_years,
    specs = sensitivity_specs
  )

  imputed_data <- run_mice_surv_list(
    surv_list = survival_data,
    m = m,
    seed = seed
  )

  # Fully adjusted (pooled, imputed HbA1c) and minimally adjusted (matching only).
  results <- dplyr::bind_rows(
    run_all_cox_pooled(
      surv_list = survival_data,
      imp_by_outcome = imputed_data,
      m = m,
      seed = seed
    ) %>%
      dplyr::mutate(scenario = label, model = "Fully adjusted"),
    run_all_cox_matching_only(survival_data) %>%
      dplyr::mutate(scenario = label, model = "Minimally adjusted")
  )

  attr(results, "overall_counts") <- tibble::tibble(
    scenario        = label,
    n_low_matched   = sum(matched$group == "low"),
    n_high_matched  = sum(matched$group == "high"),
    n_low_distinct  = dplyr::n_distinct(matched$patid[matched$group == "low"]),
    n_high_distinct = dplyr::n_distinct(matched$patid[matched$group == "high"])
  )

  results
}

# -----------------------------------------------------------------------------
# Collate HRs into wide supplementary table
# -----------------------------------------------------------------------------

collate_sensitivity_hrs <- function(..., labels = sensitivity_outcome_labels,
                                    model_label = "Fully adjusted") {
  dplyr::bind_rows(...) %>%
    dplyr::filter(outcome %in% names(labels), model == model_label) %>%
    dplyr::mutate(
      Complication = factor(
        labels[as.character(outcome)],
        levels = unname(labels)
      ),
      hr_cell = sprintf("%.2f (%.2f\u2013%.2f)", HR, LCL, UCL),
      scenario = factor(scenario, levels = unique(scenario))
    ) %>%
    dplyr::select(Complication, scenario, hr_cell) %>%
    tidyr::pivot_wider(
      names_from = scenario,
      values_from = hr_cell
    ) %>%
    dplyr::arrange(Complication)
}

# -----------------------------------------------------------------------------
# Per-analysis counts file: matched group sizes and events
# -----------------------------------------------------------------------------

write_sensitivity_counts <- function(cox_tbl,
                                     label,
                                     labels = sensitivity_outcome_labels,
                                     dir = "outputs/supplementary/tables") {
  counts <- cox_tbl %>%
    dplyr::filter(outcome %in% names(labels), model == "Fully adjusted") %>%
    dplyr::mutate(
      Complication = factor(
        labels[as.character(outcome)],
        levels = unname(labels)
      )
    ) %>%
    dplyr::arrange(Complication) %>%
    dplyr::transmute(
      Complication,
      `Lower QDiabetes score N`       = n_low,
      `Lower QDiabetes score events`  = ev_low,
      `Higher QDiabetes score N`      = n_high,
      `Higher QDiabetes score events` = ev_high
    )

  slug <- gsub("[^A-Za-z0-9]+", "_", tolower(label))
  slug <- gsub("^_|_$", "", slug)

  out_file <- file.path(dir, paste0("counts_", slug, ".csv"))

  write.csv(
    counts,
    out_file,
    row.names = FALSE
  )

  invisible(counts)
}

# -----------------------------------------------------------------------------
# Age-stratified analyses
# -----------------------------------------------------------------------------

cohort_under50 <- cohort %>%
  dplyr::filter(dm_diag_age < 50)

cohort_50plus <- cohort %>%
  dplyr::filter(dm_diag_age >= 50)

sa_age_under50 <- run_sensitivity_hr(
  cohort_df = cohort_under50,
  label = "Age at diagnosis <50 years"
)

sa_age_50plus <- run_sensitivity_hr(
  cohort_df = cohort_50plus,
  label = "Age at diagnosis \u226550 years"
)

# -----------------------------------------------------------------------------
# Sex-stratified analyses
# -----------------------------------------------------------------------------

cohort_male <- cohort %>%
  dplyr::filter(gender == 1)

cohort_female <- cohort %>%
  dplyr::filter(gender == 2)

sa_male <- run_sensitivity_hr(
  cohort_df = cohort_male,
  label = "Male"
)

sa_female <- run_sensitivity_hr(
  cohort_df = cohort_female,
  label = "Female"
)

# -----------------------------------------------------------------------------
# Collated HR table
# -----------------------------------------------------------------------------

sensitivity_hr_table <- collate_sensitivity_hrs(
  sa_age_under50,
  sa_age_50plus,
  sa_male,
  sa_female,
  model_label = "Fully adjusted"
)

write.csv(
  sensitivity_hr_table,
  "outputs/supplementary/tables/sensitivity_hazard_ratios.csv",
  row.names = FALSE
)

sensitivity_hr_table_min_adj <- collate_sensitivity_hrs(
  sa_age_under50,
  sa_age_50plus,
  sa_male,
  sa_female,
  model_label = "Minimally adjusted"
)

write.csv(
  sensitivity_hr_table_min_adj,
  "outputs/supplementary/tables/sensitivity_hazard_ratios_minimally_adjusted.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Per-analysis group sizes and N/events per complication
# -----------------------------------------------------------------------------

write_sensitivity_counts(sa_age_under50, "Age at diagnosis under 50 years")
write_sensitivity_counts(sa_age_50plus, "Age at diagnosis 50 years plus")
write_sensitivity_counts(sa_male, "Male")
write_sensitivity_counts(sa_female, "Female")

# -----------------------------------------------------------------------------
# Overall matched lower / higher score group sizes for every analysis
# -----------------------------------------------------------------------------

overall_counts_table <- dplyr::bind_rows(
  attr(sa_age_under50, "overall_counts"),
  attr(sa_age_50plus, "overall_counts"),
  attr(sa_male, "overall_counts"),
  attr(sa_female, "overall_counts")
)

write.csv(
  overall_counts_table,
  "outputs/supplementary/tables/sensitivity_overall_counts.csv",
  row.names = FALSE
)
