# R/cox_models.R
# ---------------------------------------------------------------------------
# Cox proportional-hazards modelling for the matched survival datasets.
#
# The primary analysis fits a single adjusted Cox model per outcome on the
# matched data, with missing baseline HbA1c handled by
# multiple imputation (MICE)
#   impute_prehba1c_mice() : impute missing prehba1c (once per patient)
#   run_mice_surv_list()   : run the imputation for every outcome
#   fit_cox_pooled()       : fit the adjusted Cox model + pool over imputations
#   run_all_cox_pooled()   : apply fit_cox_pooled() across all outcomes
#
# Adjusted model:
#   Surv(time, event) ~ group + age_match + sex_match + dx_year +
#                       prehba1c + imd_quintile + ethnicity_5cat
# (reference group = typical T2D).
# ---------------------------------------------------------------------------

library(dplyr)
library(survival)
library(broom)
library(purrr)

# ------------------------------------------------------------------------------
# MICE: multiple imputation of missing baseline HbA1c (prehba1c)
#
# Imputes prehba1c for a single survival dataset. Imputation is done once per
# patient (distinct patid) so that rows duplicated by matching-with-replacement
# do not bias the imputation model
#
# ------------------------------------------------------------------------------

impute_prehba1c_mice <- function(surv_df, m = 10, seed = 1, return_mids = FALSE) {
  imp_base <- surv_df %>%
    distinct(patid, .keep_all = TRUE) %>%
    transmute(
      patid,
      prehba1c       = as.numeric(prehba1c),
      age_match      = as.numeric(age_match),
      sex_match      = factor(sex_match),
      ethnicity_5cat = factor(ethnicity_5cat),
      imd_quintile   = factor(imd_quintile),
      prebmi         = as.numeric(prebmi),
      dx_year        = as.numeric(dx_year),
      group          = factor(group),
      time_years     = as.numeric(time_years),
      event          = as.integer(event)
    )

  imp_base$na_cumhaz <- mice::nelsonaalen(imp_base, time_years, event)

  imp_vars <- imp_base %>%
    select(
      prehba1c, age_match, sex_match, ethnicity_5cat,
      imd_quintile, prebmi, dx_year, group, event, na_cumhaz
    )

  meth <- mice::make.method(imp_vars)
  meth[setdiff(names(meth), "prehba1c")] <- "" # impute prehba1c only
  pred <- mice::make.predictorMatrix(imp_vars)
  pred["prehba1c", "prehba1c"] <- 0

  mids <- mice::mice(imp_vars,
    m = m, method = meth, predictorMatrix = pred,
    seed = seed, printFlag = FALSE
  )

  if (isTRUE(return_mids)) {
    return(mids)
  }

  lapply(seq_len(m), function(i) {
    tibble::tibble(
      patid        = imp_base$patid,
      prehba1c_imp = mice::complete(mids, i)$prehba1c
    )
  })
}

# Run MICE for every outcome in a surv_list
run_mice_surv_list <- function(surv_list, m = 10, seed = 1) {
  purrr::imap(
    surv_list,
    function(surv_df, .name) {
      if (sum(is.na(surv_df$prehba1c)) == 0) {
        return(NULL)
      }
      impute_prehba1c_mice(surv_df, m = m, seed = seed)
    }
  )
}

# ------------------------------------------------------------------------------
# Adjusted Cox model fitted on MICE-imputed data
#   group + matching variables + baseline HbA1c (imputed) +
#   IMD quintile (imd_quintile) + ethnicity (ethnicity_5cat)
#
# ------------------------------------------------------------------------------
fit_cox_pooled <- function(surv_df, outcome_name,
                           imp_list = NULL,
                           m = 10, seed = 1,
                           extra_covariates = NULL,
                           model_label = NULL) {
  dat <- surv_df %>%
    mutate(
      group          = factor(group, levels = c("high", "low")), # reference = high
      ethnicity_5cat = factor(ethnicity_5cat),
      imd_quintile   = factor(imd_quintile)
    ) %>%
    droplevels()

  extra_covariates <- unique(extra_covariates)
  extra_covariates <- extra_covariates[!is.na(extra_covariates)]

  missing_extra_covariates <- setdiff(extra_covariates, names(dat))
  if (length(missing_extra_covariates) > 0) {
    warning(
      "Dropping extra covariates not found in survival data: ",
      paste(missing_extra_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  extra_covariates <- intersect(extra_covariates, names(dat))

  if (length(extra_covariates) > 0) {
    dat <- dat %>%
      mutate(
        across(
          all_of(extra_covariates),
          ~ if (is.character(.x) || is.factor(.x) || is.logical(.x)) factor(.x) else .x
        )
      ) %>%
      filter(if_all(all_of(extra_covariates), ~ !is.na(.x))) %>%
      droplevels()
  }

  if (is.null(model_label)) {
    model_label <- "Matching + HbA1c (imputed) + IMD + ethnicity"
    if (length(extra_covariates) > 0) {
      model_label <- paste0(
        model_label,
        " + ",
        paste(extra_covariates, collapse = " + ")
      )
    }
  }

  cox_terms <- c(
    "group",
    "age_match",
    "dx_year",
    "prehba1c",
    "imd_quintile",
    "ethnicity_5cat"
  )

  if (
    "sex_match" %in% names(dat) &&
      dplyr::n_distinct(dat$sex_match[!is.na(dat$sex_match)]) > 1
  ) {
    cox_terms <- append(cox_terms, "sex_match", after = 2)
  }

  cox_terms <- c(cox_terms, extra_covariates)

  cox_formula <- stats::as.formula(
    paste("Surv(time_years, event) ~", paste(cox_terms, collapse = " + "))
  )


  grp_counts <- dat %>%
    group_by(group) %>%
    summarise(n = dplyr::n(), ev = sum(event), .groups = "drop")
  n_low <- sum(grp_counts$n[grp_counts$group == "low"])
  ev_low <- sum(grp_counts$ev[grp_counts$group == "low"])
  n_high <- sum(grp_counts$n[grp_counts$group == "high"])
  ev_high <- sum(grp_counts$ev[grp_counts$group == "high"])

  n_miss <- sum(is.na(dat$prehba1c))

  # ---- No missingness -> ordinary clustered Cox fit (no imputation needed) ----
  if (n_miss == 0) {
    fit <- coxph(cox_formula,
      data = dat,
      ties = "efron", robust = TRUE, cluster = patid
    )
    est <- tidy(fit, conf.int = TRUE) %>% filter(term == "grouplow")
    return(tibble::tibble(
      outcome = outcome_name, model = model_label,
      HR = exp(est$estimate), LCL = exp(est$conf.low), UCL = exp(est$conf.high),
      p = est$p.value, BIC = BIC(fit),
      n_total = fit$n, n_events = fit$nevent,
      n_low = n_low, ev_low = ev_low, n_high = n_high, ev_high = ev_high,
      m_imputations = 0L
    ))
  }

  # ---- Impute prehba1c (m completed datasets)----
  if (is.null(imp_list)) {
    imp_list <- impute_prehba1c_mice(dat, m = m, seed = seed)
  }
  m <- length(imp_list)

  # ---- Fit clustered Cox models on each completed dataset ----
  ests <- vector("list", m)
  for (i in seq_len(m)) {
    dat_i <- dat %>%
      select(-prehba1c) %>%
      left_join(imp_list[[i]], by = "patid") %>%
      rename(prehba1c = prehba1c_imp)

    fit_i <- coxph(cox_formula,
      data = dat_i,
      ties = "efron", robust = TRUE, cluster = patid
    )
    co <- summary(fit_i)$coefficients
    ests[[i]] <- list(
      est = co["grouplow", "coef"],
      se  = co["grouplow", "robust se"],
      n   = fit_i$n,
      nev = fit_i$nevent
    )
  }

  # ---- Pool log-HR with Rubin's rules ----
  q_values <- vapply(ests, function(x) x$est, numeric(1))
  u_values <- vapply(ests, function(x) x$se^2, numeric(1))
  q_bar <- mean(q_values)
  u_bar <- mean(u_values)
  between_variance <- stats::var(q_values)
  total_variance <- u_bar + (1 + 1 / m) * between_variance
  se_pool <- sqrt(total_variance)

  r <- (1 + 1 / m) * between_variance / u_bar
  df_old <- (m - 1) * (1 + 1 / r)^2
  crit <- stats::qt(0.975, df = df_old)

  tibble::tibble(
    outcome = outcome_name, model = model_label,
    HR = exp(q_bar),
    LCL = exp(q_bar - crit * se_pool),
    UCL = exp(q_bar + crit * se_pool),
    p = 2 * stats::pt(-abs(q_bar / se_pool), df = df_old),
    BIC = NA_real_,
    n_total = ests[[1]]$n,
    n_events = ests[[1]]$nev,
    n_low = n_low, ev_low = ev_low, n_high = n_high, ev_high = ev_high,
    m_imputations = as.integer(m)
  )
}

# Run the pooled Cox model across all outcomes.
run_all_cox_pooled <- function(surv_list, imp_by_outcome = NULL,
                               m = 10, seed = 1,
                               extra_covariates = NULL,
                               model_label = NULL) {
  purrr::imap_dfr(
    surv_list,
    ~ fit_cox_pooled(.x,
      outcome_name = .y,
      imp_list = imp_by_outcome[[.y]],
      m = m, seed = seed,
      extra_covariates = extra_covariates,
      model_label = model_label
    )
  )
}

# Cox model without baseline HbA1c. This supports both the minimally adjusted
# model and the treatment-outcome model, which additionally includes IMD and
# ethnicity.
fit_cox_without_hba1c <- function(
    surv_df,
    outcome_name,
    extra_covariates = character(),
    model_label = "Age + sex + diagnosis year") {
  dat <- surv_df %>%
    dplyr::mutate(group = factor(group, levels = c("high", "low"))) %>%
    droplevels()

  extra_covariates <- unique(extra_covariates)
  missing_covariates <- setdiff(extra_covariates, names(dat))
  if (length(missing_covariates) > 0L) {
    stop(
      "Missing Cox-model covariates: ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(extra_covariates) > 0L) {
    dat <- dat %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(extra_covariates),
          ~ if (is.character(.x) || is.logical(.x)) factor(.x) else .x
        )
      )
  }

  cox_terms <- c("group", "age_match", "dx_year")
  if (
    "sex_match" %in% names(dat) &&
      dplyr::n_distinct(dat$sex_match[!is.na(dat$sex_match)]) > 1L
  ) {
    cox_terms <- append(cox_terms, "sex_match", after = 2L)
  }
  cox_terms <- c(cox_terms, extra_covariates)

  model_vars <- c("time_years", "event", "patid", cox_terms)
  dat <- dat %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(model_vars), ~ !is.na(.x))
    ) %>%
    droplevels()

  cox_formula <- stats::as.formula(
    paste(
      "survival::Surv(time_years, event) ~",
      paste(cox_terms, collapse = " + ")
    )
  )

  group_counts <- dat %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n = dplyr::n(), events = sum(event), .groups = "drop")

  fit <- survival::coxph(
    cox_formula,
    data = dat,
    ties = "efron",
    robust = TRUE,
    cluster = patid
  )

  coefficients <- summary(fit)$coefficients
  robust_se_col <- if ("robust se" %in% colnames(coefficients)) {
    "robust se"
  } else {
    "se(coef)"
  }
  beta <- coefficients["grouplow", "coef"]
  standard_error <- coefficients["grouplow", robust_se_col]

  tibble::tibble(
    outcome = outcome_name,
    model = model_label,
    HR = exp(beta),
    LCL = exp(beta - 1.96 * standard_error),
    UCL = exp(beta + 1.96 * standard_error),
    p = coefficients["grouplow", "Pr(>|z|)"],
    BIC = stats::BIC(fit),
    n_total = fit$n,
    n_events = fit$nevent,
    n_low = sum(group_counts$n[group_counts$group == "low"]),
    ev_low = sum(group_counts$events[group_counts$group == "low"]),
    n_high = sum(group_counts$n[group_counts$group == "high"]),
    ev_high = sum(group_counts$events[group_counts$group == "high"])
  )
}

fit_cox_matching_only <- function(
    surv_df,
    outcome_name,
    model_label = "Age + sex + diagnosis year") {
  fit_cox_without_hba1c(
    surv_df = surv_df,
    outcome_name = outcome_name,
    model_label = model_label
  )
}

fit_cox_matching_imd_ethnicity <- function(
    surv_df,
    outcome_name,
    model_label = "Matching + IMD + ethnicity") {
  fit_cox_without_hba1c(
    surv_df = surv_df,
    outcome_name = outcome_name,
    extra_covariates = c("imd_quintile", "ethnicity_5cat"),
    model_label = model_label
  )
}

run_all_cox_matching_only <- function(surv_list) {
  purrr::imap_dfr(
    surv_list,
    ~ fit_cox_matching_only(.x, outcome_name = .y)
  )
}

run_all_cox_matching_imd_ethnicity <- function(surv_list) {
  purrr::imap_dfr(
    surv_list,
    ~ fit_cox_matching_imd_ethnicity(.x, outcome_name = .y)
  )
}
