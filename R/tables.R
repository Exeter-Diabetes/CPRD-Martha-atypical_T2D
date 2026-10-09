# R/tables.R
# ---------------------------------------------------------------------------
#   build_table1()           : Table 1 for either unmatched or matched cohort
#   combine_table1()         : merge unmatched + matched into one wide table
#   summarise_outcomes()     : per-outcome events / follow-up by group
#   summarise_km_at_time()   : KM cumulative incidence at a fixed time point
# ---------------------------------------------------------------------------

library(dplyr)
library(tableone)
library(survey)
library(tibble)

# Define Table 1 variables
table1_spec <- function() {
  vars <- c(
    "dm_diag_age",
    "gender",
    "ethnicity_5cat",
    "agecat",
    "imd_quintile",
    "prebmi",
    "prehba1c",
    "bmi_cat",
    "hba1c_cat",
    "fh_diabetes",
    "preegfr",
    "prehdl",
    "preldl",
    "prenonhdl",
    "pretriglyceride",
    "pretghdl_ratio",
    "qrisk2_10y_risk_score",
    "qrisk2_risk_cat",
    "smoking_cat",
    "alcohol_cat",
    "pre_index_date_cvd",
    "pre_index_date_heartfailure",
    "htn",
    "pre_index_date_ckd",
    "pre_index_date_solid_cancer",
    "pre_index_date_haem_cancer",
    "pre_index_date_acutepancreatitis",
    "pre_index_date_chronicpancreatitis",
    "pre_index_date_pancreaticcancer",
    "pre_index_date_surgicalpancreaticresection",
    "pre_index_date_qresearch_polycystic_ovaries",
    "pre_index_date_qresearch_learning_disability",
    "pre_index_date_qresearch_bipolarschiz",
    "ster",
    "stat",
    "apsy",
    "pre_index_date_severe_retinopathy",
    "pre_index_date_severe_neuropathy",
    "pre_index_date_severe_nephropathy",
    "pre_index_date_incident_mi",
    "pre_index_date_incident_stroke",
    "pre_index_date_dka",
    "pre_index_date_hypoglycaemia",
    "pre_index_date_myocardialinfarction"
  )

  factor_vars <- c(
    "gender",
    "ethnicity_5cat",
    "agecat",
    "imd_quintile",
    "ins_in_1_year",
    "bmi_cat",
    "type1_any",
    "hba1c_cat",
    "fh_diabetes",
    "qrisk2_risk_cat",
    "smoking_cat",
    "alcohol_cat",
    "pre_index_date_cvd",
    "pre_index_date_heartfailure",
    "htn",
    "pre_index_date_ckd",
    "pre_index_date_solid_cancer",
    "pre_index_date_haem_cancer",
    "pre_index_date_acutepancreatitis",
    "pre_index_date_chronicpancreatitis",
    "pre_index_date_pancreaticcancer",
    "pre_index_date_surgicalpancreaticresection",
    "pre_index_date_qresearch_polycystic_ovaries",
    "pre_index_date_qresearch_learning_disability",
    "pre_index_date_qresearch_bipolarschiz",
    "ster",
    "stat",
    "apsy",
    "pre_index_date_severe_retinopathy",
    "pre_index_date_severe_neuropathy",
    "pre_index_date_severe_nephropathy",
    "pre_index_date_incident_mi",
    "pre_index_date_incident_stroke",
    "pre_index_date_dka",
    "pre_index_date_hypoglycaemia",
    "pre_index_date_myocardialinfarction"
  )

  list(vars = vars, factor_vars = factor_vars)
}

tableone_to_df <- function(tab,
                           show_all_levels = TRUE,
                           missing = TRUE,
                           smd = FALSE) {
  out <- print(
    tab,
    showAllLevels = show_all_levels,
    missing = missing,
    smd = smd,
    printToggle = FALSE,
    noSpaces = TRUE
  )

  as.data.frame(out, check.names = FALSE) %>%
    tibble::rownames_to_column("variable")
}


# -------------------------------------------------------------------------
# Baseline table: generic Table 1 builder for unmatched or matched cohorts
# -------------------------------------------------------------------------
build_table1 <- function(df,
                         strata_col = "qdiabetes_risk_cat",
                         strata_levels = c("low", "high"),
                         strata_labels = c(
                           low = "Low diabetes risk score",
                           high = "High diabetes risk score"
                         ),
                         test = FALSE,
                         smd = FALSE) {
  spec <- table1_spec()
  needed <- unique(c(strata_col, spec$vars, spec$factor_vars))

  dat_tab <- df %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    dplyr::mutate(
      strata = factor(
        .data[[strata_col]],
        levels = strata_levels,
        labels = c(strata_labels[["low"]], strata_labels[["high"]])
      )
    )

  tab <- tableone::CreateTableOne(
    vars = spec$vars,
    strata = "strata",
    data = dat_tab,
    factorVars = spec$factor_vars,
    test = test
  )

  list(
    tab = tab,
    df = tableone_to_df(tab, smd = smd)
  )
}

# -------------------------------------------------------------------------
# Merge unmatched + matched Table 1 into one wide table
# -------------------------------------------------------------------------
combine_table1 <- function(unmatched_obj, matched_obj) {
  unmatched_df <- unmatched_obj$df
  matched_df <- matched_obj$df

  if (!"variable" %in% names(unmatched_df)) {
    unmatched_df <- tibble::rownames_to_column(as.data.frame(unmatched_df), "variable")
  }
  if (!"variable" %in% names(matched_df)) {
    matched_df <- tibble::rownames_to_column(as.data.frame(matched_df), "variable")
  }

  unmatched_cols <- setdiff(names(unmatched_df), "variable")
  matched_cols <- setdiff(names(matched_df), "variable")

  unmatched_df <- unmatched_df %>%
    dplyr::rename_with(~ paste0("Unmatched: ", .x), dplyr::all_of(unmatched_cols))

  matched_df <- matched_df %>%
    dplyr::rename_with(~ paste0("Matched: ", .x), dplyr::all_of(matched_cols))

  dplyr::full_join(unmatched_df, matched_df, by = "variable")
}

summarise_outcomes <- function(surv_list) {
  results <- purrr::imap_dfr(
    surv_list,
    function(df, outcome_name) {
      df %>%
        dplyr::group_by(treated) %>%
        dplyr::summarise(
          n_events = sum(event, na.rm = TRUE),
          n_subjects = n(),
          event_rate = mean(event, na.rm = TRUE),
          mean_followup = mean(time_years, na.rm = TRUE),
          person_years = sum(time_years, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          group = dplyr::if_else(
            treated == 1,
            "Low diabetes risk score",
            "High diabetes risk score (matched)"
          ),
          complication = outcome_name,
          .before = treated
        ) %>%
        dplyr::select(
          complication,
          group,
          n_subjects,
          n_events,
          event_rate,
          mean_followup,
          person_years
        )
    }
  )

  results %>%
    dplyr::mutate(
      event_rate = round(event_rate * 100, 2),
      across(c(mean_followup, person_years), ~ round(., 2)),
      complication = stringr::str_replace_all(complication, "_", " ") %>%
        stringr::str_to_sentence()
    ) %>%
    dplyr::arrange(complication, group)
}

# Kaplan-Meier cumulative incidence at a fixed time point.
summarise_km_at_time <- function(surv_list,
                                 time_years = 10,
                                 outcome_labels = NULL,
                                 group_labels = c(
                                   low = "Low diabetes risk score",
                                   high = "High diabetes risk score (matched)"
                                 )) {
  purrr::imap_dfr(
    surv_list,
    function(df, outcome_name) {
      complication_label <- stringr::str_to_sentence(
        stringr::str_replace_all(outcome_name, "_", " ")
      )

      if (!is.null(outcome_labels) && outcome_name %in% names(outcome_labels)) {
        complication_label <- unname(outcome_labels[outcome_name])
      }

      dat <- df %>%
        dplyr::mutate(group = factor(group, levels = c("low", "high")))

      fit <- survival::survfit(
        survival::Surv(time_years, event) ~ group,
        data = dat
      )

      km <- summary(fit, times = time_years, extend = TRUE)

      km_tbl <- tibble::tibble(
        outcome = outcome_name,
        complication = complication_label,
        group_raw = sub("^group=", "", as.character(km$strata)),
        time_years = as.numeric(km$time),
        n_at_risk = as.integer(round(km$n.risk)),
        km_cumulative_incidence = 1 - as.numeric(km$surv),
        lower_95_ci = 1 - as.numeric(km$upper),
        upper_95_ci = 1 - as.numeric(km$lower)
      )

      event_counts <- dat %>%
        dplyr::group_by(group_raw = as.character(group)) %>%
        dplyr::summarise(
          n_subjects = dplyr::n(),
          n_events = sum(event, na.rm = TRUE),
          .groups = "drop"
        )

      km_tbl %>%
        dplyr::left_join(event_counts, by = "group_raw") %>%
        dplyr::mutate(
          group = dplyr::coalesce(
            unname(group_labels[group_raw]),
            group_raw
          ),
          km_cumulative_incidence_pct =
            round(100 * km_cumulative_incidence, 1),
          lower_95_ci_pct = round(100 * lower_95_ci, 1),
          upper_95_ci_pct = round(100 * upper_95_ci, 1),
          km_95_ci_pct = sprintf(
            "%.1f (%.1f-%.1f)",
            km_cumulative_incidence_pct,
            lower_95_ci_pct,
            upper_95_ci_pct
          )
        ) %>%
        dplyr::select(
          outcome,
          complication,
          group,
          time_years,
          n_subjects,
          n_events,
          n_at_risk,
          km_cumulative_incidence,
          lower_95_ci,
          upper_95_ci,
          km_cumulative_incidence_pct,
          lower_95_ci_pct,
          upper_95_ci_pct,
          km_95_ci_pct
        )
    }
  )
}
