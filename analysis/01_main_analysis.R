# =============================================================================
# 01_main_analysis.R
# -----------------------------------------------------------------------------
#
# Compares prognosis of people who develop diabetes despite low predicted risk
# ("atypical T2D"; QDiabetes risk score <=5.6%) with age- and sex-matched
# individuals who develop diabetes with high predicted risk
# ("typical T2D"; QDiabetes risk score >5.6%).
#
# Pipeline:
#   1. Load the at-diagnosis cohort from CPRD.
#   2. Build the analysis cohort and compute the QDiabetes risk score.
#   3. Figure 1  - clinical features of atypical vs typical T2D group.
#   4. Match atypical- vs typical patients and build Table 1.
#   5. Build the primary-outcome survival datasets.
#   6. Figure 2  - Kaplan-Meier cumulative incidence (micro/macrovascular).
#   7. Figure 3  - adjusted Cox hazard ratios.
#   8. Table 2   - secondary outcomes (acute, mortality, treatment).
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
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/supplementary/figures", recursive = TRUE, showWarnings = FALSE)

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
# Apply the inclusion criteria and derive variables used for complication
# analyses by joining the diagnosis, HbA1c, and drug start/stop data.
cohort <- prepare_main_dataset(
  at_diag_raw = at_diag,
  hba1c_data = hba1c_clean,
  qdiabetes_risk_horizon = 10L,
  qdiabetes_risk_threshold = 5.6,
  drug_classes = drug_start_stop
)

# =============================================================================
# 3. Describe the unmatched cohort
# =============================================================================
# Figure 1: 4-panel clinical-feature comparison by risk group.
figure1 <- plot_clinical_features(
  df       = cohort,
  out_file = "outputs/figures/figure1.png"
)


# =============================================================================
# 4. Match atypical- vs typical patients
# =============================================================================
# Exact matching on sex, age at diagnosis and diagnosis year (4:1, with
# replacement).
matched_obj <- match_low_high(
  cohort,
  treated_level = "low",
  control_level = "high",
  ratio = 4
)

matched_df <- matched_obj$matched


# Covariate balance after matching (Love plot)
love_plot <- plot_love_balance(
  matched_obj$matchit,
  threshold = 0.1,
  out_file  = "outputs/supplementary/figures/love_plot.png"
)


# Table 1: baseline characteristics, unmatched and matched, combined into one
# table with standardised mean differences.
tab1_unmatched <- build_table1(
  df = cohort,
  strata_col = "qdiabetes_risk_cat",
  strata_levels = c("low", "high"),
  strata_labels = c(
    low = "Low diabetes risk score",
    high = "High diabetes risk score"
  ),
  test = TRUE,
  smd = TRUE
)

tab1_matched <- build_table1(
  df = matched_df,
  strata_col = "group",
  strata_levels = c("low", "high"),
  strata_labels = c(
    low = "Low diabetes risk score",
    high = "High diabetes risk score (matched)"
  ),
  test = FALSE,
  smd = TRUE
)
tab1_combined <- combine_table1(
  unmatched_obj = tab1_unmatched,
  matched_obj   = tab1_matched
)
write.csv(tab1_combined, "outputs/tables/table1.csv", row.names = FALSE)


# =============================================================================
# 5. Primary-outcome survival datasets
# =============================================================================
# Microvascular and macrovascular complications
complications <- c(
  "retinopathy_severe", "neuropathy_severe", "nephropathy_severe",
  "mi_fatal_nonfatal", "stroke_fatal_nonfatal", "hf_fatal_nonfatal"
)

surv_list <- make_surv_list(
  matched_df    = matched_df,
  horizon_years = 10,
  specs         = outcome_specs[complications]
)


# =============================================================================
# 6. Figure 2: Kaplan-Meier cumulative incidence
# =============================================================================
figure2 <- plot_km_complication_panels(
  surv_list       = surv_list,
  groups          = km_complication_panels,
  legend_labels   = c("Atypical type 2 diabetes", "Typical type 2 diabetes"),
  ncol_each       = 3,
  base_size       = 25,
  curve_linewidth = 2.0,
  conf_int_alpha  = 0.35,
  out_file        = "outputs/figures/figure2.png",
  width           = 26,
  height          = 16,
  dpi             = 600
)

# 10-year Kaplan-Meier cumulative incidence for the primary outcomes.
primary_outcomes_cumulative_incidence_10y <- summarise_km_at_time(
  surv_list,
  time_years     = 10,
  outcome_labels = purrr::map_chr(km_outcome_settings[complications], "label")
)
write.csv(
  primary_outcomes_cumulative_incidence_10y,
  "outputs/tables/primary_outcomes_cumulative_incidence_10y.csv",
  row.names = FALSE
)


# =============================================================================
# 7. Figure 3: adjusted Cox hazard ratios
# =============================================================================
# Forest plot of risk.
# Reference = high risk. Adjusted for the matching variables + baseline HbA1c +
# IMD quintile + ethnicity. Missing baseline HbA1c is imputed.
hba1c_imputations <- run_mice_surv_list(surv_list, m = 10, seed = 1)

cox_results <- run_all_cox_pooled(
  surv_list,
  imp_by_outcome = hba1c_imputations, m = 10, seed = 1
)

write.csv(cox_results, "outputs/tables/cox_results.csv", row.names = FALSE)

figure3 <- plot_cox_forest(
  cox_results = cox_results,
  outcome_categories = list(
    "A) Microvascular" = c("retinopathy_severe", "neuropathy_severe", "nephropathy_severe"),
    "B) Macrovascular" = c("mi_fatal_nonfatal", "stroke_fatal_nonfatal", "hf_fatal_nonfatal")
  ),
  out_file = "outputs/figures/figure3.png",
  width = 40,
  height = 18,
  dpi = 600,
  base_size = 36
)


# =============================================================================
# 8. Table 2: secondary outcomes
# =============================================================================
# Acute complications, mortality and treatment initiation. Outputs per-group
# events, person-years, 10-year cumulative incidence, and the adjusted HR
# (reference = typical T2D). Acute and mortality outcomes adjust for the
# matching variables + imputed HbA1c + IMD + ethnicity. Treatment outcomes use
# the same adjustment without baseline HbA1c.

secondary_specs <- c(
  outcome_specs[c("dka", "hypoglycaemia", "cv_death_primary")],
  list(
    non_cv_death = list(
      event_col      = "non_cv_death_date",
      pre_col        = "pre_index_date_death",
      censor_profile = "mortality"
    )
  ),
  outcome_specs[c("first_line_treatment", "second_line_treatment", "insulin_gt1y")]
)

secondary_labels <- c(
  dka                   = "Diabetic ketoacidosis",
  hypoglycaemia         = "Hypoglycaemia",
  cv_death_primary      = "Cardiovascular mortality",
  non_cv_death          = "Non-cardiovascular mortality",
  first_line_treatment  = "First-line therapy initiation",
  second_line_treatment = "Second-line therapy initiation",
  insulin_gt1y          = "Insulin initiation"
)

surv_list_secondary <- make_surv_list(
  matched_df    = matched_df,
  horizon_years = 10,
  specs         = secondary_specs
)

# HbA1c-adjusted HR for acute/mortality (HbA1c excluded for treatment outcomes)
hba1c_adjusted_outcomes <- c("dka", "hypoglycaemia", "cv_death_primary", "non_cv_death")
treatment_outcomes_tab2 <- c("first_line_treatment", "second_line_treatment", "insulin_gt1y")

secondary_hba1c_imp <- run_mice_surv_list(
  surv_list_secondary[hba1c_adjusted_outcomes],
  m = 10, seed = 1
)

secondary_ahr <- dplyr::bind_rows(
  run_all_cox_pooled(
    surv_list_secondary[hba1c_adjusted_outcomes],
    imp_by_outcome = secondary_hba1c_imp, m = 10, seed = 1
  ),
  run_all_cox_matching_imd_ethnicity(surv_list_secondary[treatment_outcomes_tab2])
) %>%
  dplyr::transmute(outcome, ahr = sprintf("%.2f (%.2f-%.2f)", HR, LCL, UCL))

# Per-group events and person-years.
secondary_counts <- purrr::imap_dfr(surv_list_secondary, function(df, nm) {
  df %>%
    dplyr::group_by(grp = as.character(group)) %>%
    dplyr::summarise(
      events       = sum(event, na.rm = TRUE),
      person_years = round(sum(time_years, na.rm = TRUE)),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(outcome = nm)
})

# Per-group 10-year cumulative incidence
secondary_risk <- summarise_km_at_time(surv_list_secondary, time_years = 10) %>%
  dplyr::mutate(grp = dplyr::if_else(grepl("Low", group), "low", "high")) %>%
  dplyr::select(outcome, grp, risk = km_95_ci_pct)

table2 <- secondary_counts %>%
  dplyr::left_join(secondary_risk, by = c("outcome", "grp")) %>%
  tidyr::pivot_wider(
    names_from  = grp,
    values_from = c(events, person_years, risk),
    names_glue  = "{grp}_{.value}"
  ) %>%
  dplyr::right_join(
    tibble::tibble(
      outcome = names(secondary_labels),
      Outcome = unname(secondary_labels)
    ),
    by = "outcome"
  ) %>%
  dplyr::left_join(secondary_ahr, by = "outcome") %>%
  dplyr::arrange(match(outcome, names(secondary_labels))) %>%
  dplyr::transmute(
    Outcome,
    `Atypical events`               = low_events,
    `Atypical person-years`         = low_person_years,
    `Atypical 10y risk % (95% CI)`  = low_risk,
    `Typical events`                = high_events,
    `Typical person-years`          = high_person_years,
    `Typical 10y risk % (95% CI)`   = high_risk,
    `aHR (95% CI)`                  = ahr
  )

write.csv(table2, "outputs/tables/table2.csv", row.names = FALSE)
