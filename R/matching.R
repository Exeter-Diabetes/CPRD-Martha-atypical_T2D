# R/matching.R
# ---------------------------------------------------------------------------
# Exact matching of individuals with atypical vs typical T2D
#
# Atypical T2D (treated) are exact-matched to Typical T2D (control)
# on sex, age at diagnosis and year of diagnosis, 4:1 with replacement.
# Returns the MatchIt object and the matched data (with weights + match id).
# ---------------------------------------------------------------------------

library(dplyr)
library(MatchIt)

# Assumes the input data has:
# - qdiabetes_risk_cat ("low"/"high")
# - gender (1 = male / 2 = female), dm_diag_age, dx_year
match_low_high <- function(df,
                           treated_level = "low",
                           control_level = "high",
                           ratio = 4,
                           seed = 1,
                           replace = TRUE) {
  dat <- df %>%
    mutate(
      group = case_when(
        qdiabetes_risk_cat == treated_level ~ treated_level,
        qdiabetes_risk_cat == control_level ~ control_level,
        TRUE ~ NA_character_
      ),
      treated = as.integer(group == treated_level),
      sex_match = factor(if_else(gender == 1L, "Male", "Female")),
      age_match = as.integer(dm_diag_age),
      dx_year = as.integer(dx_year)
    ) %>%
    filter(
      !is.na(group),
      !is.na(sex_match),
      !is.na(age_match),
      !is.na(dx_year)
    )

  set.seed(seed)

  m <- MatchIt::matchit(
    treated ~ 1,
    data    = dat,
    method  = "nearest",
    exact   = ~ sex_match + age_match + dx_year,
    ratio   = ratio,
    replace = replace,
    m.order = "random"
  )

  matched <- MatchIt::get_matches(m) %>%
    mutate(
      weight   = weights,
      match_id = subclass,
      group    = if_else(treated == 1L, treated_level, control_level)
    ) %>%
    select(-weights) %>%
    relocate(patid, group, treated, match_id, weight)

  list(
    matchit = m,
    matched = matched
  )
}
