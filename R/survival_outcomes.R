# R/survival_outcomes.R
# Build matched survival datasets for outcome variables.

library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)

# Censoring profiles
# If the outcome is death-related, don't censor at death (death is the event)
censor_profiles <- list(
  clinical  = c("gp_end_date", "hes_end_date", "death_date"),
  mortality = c("gp_end_date", "hes_end_date")
)

# Outcome specifications:
# - event_col: incident event date
# - pre_col: prevalent disease flag (used for baseline exclusion)
# - censor_profile: censoring rule key
outcome_specs <- list(

  # a) microvascular outcomes
  retinopathy_severe = list(
    event_col = "post_index_date_first_severe_retinopathy",
    pre_col = "pre_index_date_severe_retinopathy",
    censor_profile = "clinical"
  ),
  neuropathy_severe = list(
    event_col = "post_index_date_first_severe_neuropathy",
    pre_col = "pre_index_date_severe_neuropathy",
    censor_profile = "clinical"
  ),
  nephropathy_severe = list(
    event_col = "post_index_date_first_severe_nephropathy",
    pre_col = "pre_index_date_severe_nephropathy",
    censor_profile = "clinical"
  ),

  # macrovascular outcomes (fatal or non-fatal)
  mi_fatal_nonfatal = list(
    event_col = "post_index_date_first_mi_fatal_nonfatal",
    pre_col = "pre_index_date_mi_fatal_nonfatal",
    censor_profile = "mortality"
  ),
  stroke_fatal_nonfatal = list(
    event_col = "post_index_date_first_stroke_fatal_nonfatal",
    pre_col = "pre_index_date_stroke_fatal_nonfatal",
    censor_profile = "mortality"
  ),
  hf_fatal_nonfatal = list(
    event_col = "post_index_date_first_hf_fatal_nonfatal",
    pre_col = "pre_index_date_hf_fatal_nonfatal",
    censor_profile = "mortality"
  ),

  # --- Secondary outcomes ---

  dka = list(
    event_col = "post_index_date_first_dka",
    pre_col = "pre_index_date_dka",
    censor_profile = "clinical"
  ),
  hypoglycaemia = list(
    event_col = "post_index_date_first_hypoglycaemia",
    pre_col = "pre_index_date_hypoglycaemia",
    censor_profile = "clinical"
  ),
  cv_death_primary = list(
    event_col = "cv_death_primary_cause_date",
    pre_col = "pre_index_date_death",
    censor_profile = "mortality"
  ),


  # --- Treatment outcomes ---
  insulin_gt1y = list(
    event_col = "dm_diag_insdate_gt1y",
    pre_col = NULL,
    censor_profile = "clinical",
    landmark_years = 1
  ),
  first_line_treatment = list(
    event_col = "treatment_first_line_earliest",
    pre_col = "pre_index_date_diabetes_meds",
    censor_profile = "clinical",
    allow_index_date_events = TRUE
  ),
  second_line_treatment = list(
    event_col = "treatment_second_line_earliest",
    pre_col = "pre_index_date_diabetes_meds",
    censor_profile = "clinical"
  )
)



make_surv_df <- function(matched_df,
                         event_col,
                         pre_exclude_col,
                         index_col = "dm_diag_date",
                         censor_cols,
                         horizon_years = 10,
                         landmark_years = 0,
                         allow_index_date_events = FALSE) {
  dat0 <- matched_df %>%
    mutate(treated = as.integer(treated))

  if (is.null(pre_exclude_col)) {
    dat0 <- dat0 %>%
      mutate(pre_exclude = FALSE)
  } else {
    dat0 <- dat0 %>%
      mutate(pre_exclude = coalesce(as.integer(.data[[pre_exclude_col]]), 0L) == 1L)
  }

  # If treated has prevalent disease, drop the entire matched set
  bad_sets <- dat0 %>%
    filter(treated == 1L, pre_exclude) %>%
    distinct(match_id)

  dat1 <- dat0 %>%
    anti_join(bad_sets, by = "match_id") %>%
    # If a control has prevalent disease, drop that control only
    filter(!(treated == 0L & pre_exclude))

  # Keep sets with >=1 treated and >=1 control
  keep_sets <- dat1 %>%
    group_by(match_id) %>%
    summarise(
      n_treated = sum(treated == 1L),
      n_control = sum(treated == 0L),
      .groups = "drop"
    ) %>%
    filter(n_treated >= 1, n_control >= 1) %>%
    pull(match_id)

  dat2 <- dat1 %>%
    filter(match_id %in% keep_sets)

  dat3 <- dat2 %>%
    mutate(
      index = as.Date(.data[[index_col]]),
      event_date = as.Date(.data[[event_col]]),
      # Use %m+% so 29-Feb diagnoses roll to 28-Feb instead of returning NA
      # (plain + years() yields NA for non-existent dates like 2026-02-29).
      horizon_cap = index %m+% years(horizon_years),
      censor_date = do.call(
        pmin,
        c(as.list(across(all_of(censor_cols))), list(horizon_cap), list(na.rm = TRUE))
      ),
      fu_end_date = pmin(event_date, censor_date, na.rm = TRUE),
      reached_horizon = fu_end_date >= horizon_cap,
      index_date_event_with_follow_up =
        isTRUE(allow_index_date_events) &
        !is.na(event_date) &
        event_date == index &
        event_date <= censor_date &
        censor_date > index,
      time_years = if_else(
        index_date_event_with_follow_up,
        0.5 / 365.25,
        if_else(
          reached_horizon,
          as.numeric(horizon_years),
          as.numeric(fu_end_date - index) / 365.25
        )
      ),
      event = as.integer(!is.na(event_date) & event_date <= censor_date)
    ) %>%
    filter(
      !is.na(time_years),
      time_years >= 0
    ) %>%
    mutate(
      time_years = if_else(time_years == 0, 0.5 / 365.25, time_years)
    )

  # Keep only individuals who are still event-free
  # and under follow-up at the landmark (i.e. their follow-up extends beyond it);
  if (landmark_years > 0) {
    dat3 <- dat3 %>%
      filter(time_years > landmark_years)
  }

  valid_sets <- dat3 %>%
    count(match_id, treated) %>%
    pivot_wider(names_from = treated, values_from = n, values_fill = 0) %>%
    filter(`1` >= 1, `0` >= 1) %>%
    pull(match_id)

  dat3 %>%
    filter(match_id %in% valid_sets)
}

make_surv_list <- function(matched_df,
                           horizon_years = 10,
                           specs = outcome_specs,
                           profiles = censor_profiles) {
  purrr::imap(
    specs,
    function(spec, .name) {
      make_surv_df(
        matched_df = matched_df,
        event_col = spec$event_col,
        pre_exclude_col = spec$pre_col,
        censor_cols = profiles[[spec$censor_profile]],
        horizon_years = horizon_years,
        landmark_years = spec$landmark_years %||% 0,
        allow_index_date_events = spec$allow_index_date_events %||% FALSE
      )
    }
  )
}
