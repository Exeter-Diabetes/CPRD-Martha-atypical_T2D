# R/cohort_definition.R


is_integer64 <- function(x) identical(class(x), "integer64")

# Safe min across Date columns (returns NA if all inputs are NA)
pmin_date <- function(...) {
  x <- list(...)
  x <- lapply(x, function(z) if (inherits(z, "Date")) z else as.Date(z))
  out <- do.call(pmin, c(x, list(na.rm = TRUE)))
  out[is.infinite(as.numeric(out))] <- NA
  out
}

# ------------------------------------------------------------------------------
# Cohort definition
# ------------------------------------------------------------------------------
define_cohort_at_diag <- function(at_diag_raw,
                                  hba1c_data = NULL,
                                  study_start_date = as.Date("2013-01-01"),
                                  study_end_date = as.Date("2023-03-31"),
                                  age_min = 25,
                                  age_max = 85, # QDiabetes valid age < 85
                                  hba1c_min = 48,
                                  require_hes = TRUE,
                                  verbose = TRUE) {
  .n <- function(df) format(nrow(df), big.mark = ",")
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  # ---------------------------------------------------------------------------
  # Step 1: Baseline type 2 diabetes cohort (2013-2023)
  # ---------------------------------------------------------------------------
  cohort <- at_diag_raw %>%
    dplyr::filter(gender != 3) %>%
    dplyr::mutate(
      patid = as.character(patid),
      dplyr::across(dplyr::where(is_integer64), as.integer)
    ) %>%
    dplyr::filter(
      diabetes_type == "type 2",
      # diagnosed after GP registration
      dm_diag_date > regstartdate,
      # Diagnosed between study start and end dates
      dm_diag_date >= study_start_date,
      dm_diag_date <= study_end_date
    )

  # HbA1c evidence of diabetes.
  #   - Keep patients with an at-diagnosis HbA1c (prehba1c) >= hba1c_min.
  #   - Keep patients with MISSING prehba1c only if they have at least one
  #     HbA1c record > 47 mmol/mol at some point.
  #   - Exclude patients with MISSING prehba1c who either never had any HbA1c
  #     record, or had HbA1c record(s) but never > 47.
  if (!is.null(hba1c_data)) {
    qualifying_hba1c_patids <- hba1c_data %>%
      dplyr::mutate(
        patid     = as.character(patid),
        testvalue = as.numeric(testvalue)
      ) %>%
      dplyr::filter(!is.na(testvalue)) %>%
      dplyr::group_by(patid) %>%
      dplyr::summarise(any_hba1c_gt47 = any(testvalue > 47), .groups = "drop") %>%
      dplyr::filter(any_hba1c_gt47) %>%
      dplyr::pull(patid)

    n_pre <- nrow(cohort)
    cohort <- cohort %>%
      dplyr::filter(
        !is.na(prehba1c) | # HbA1c present at diagnosis: keep
          (patid %in% qualifying_hba1c_patids) # missing at dx but HbA1c > 47 elsewhere: keep
      ) %>%
      dplyr::filter(is.na(prehba1c) | prehba1c >= hba1c_min)
    .msg(
      "  Excluded missing at-dx HbA1c with no HbA1c > 47:      N = ",
      format(n_pre - nrow(cohort), big.mark = ",")
    )
  } else {
    cohort <- cohort %>% dplyr::filter(is.na(prehba1c) | prehba1c >= hba1c_min)
  }

  .msg(
    "Newly diagnosed type 2 diabetes ",
    format(study_start_date, "%Y"), "\u2013", format(study_end_date, "%Y"),
    " with HbA1c \u226548 mmol/mol:   N = ", .n(cohort)
  )

  # ---------------------------------------------------------------------------
  # Step 2a: HES linkage
  # ---------------------------------------------------------------------------
  if (isTRUE(require_hes)) {
    n_pre <- nrow(cohort)
    cohort <- cohort %>% dplyr::filter(with_hes == 1)
    .msg(
      "  Excluded without HES linkage:                         N = ",
      format(n_pre - nrow(cohort), big.mark = ",")
    )
  }

  # ---------------------------------------------------------------------------
  # Step 2b: Possible non-type 2 / secondary diabetes
  #   - insulin within 1 year of diagnosis
  #   - any type 1 diabetes code
  #   - acute/chronic pancreatitis, pancreatic resection, pancreatic cancer,
  #     haemochromatosis, or cystic fibrosis
  # ---------------------------------------------------------------------------
  # Type 1 reasons: insulin within 1 year of diagnosis or any type 1 code.
  type1_cols <- c("ins_in_1_year")
  # Type 3c reasons: pancreatic / secondary diabetes conditions.
  type3c_cols <- c(
    "pre_index_date_acutepancreatitis",
    "pre_index_date_chronicpancreatitis",
    "pre_index_date_surgicalpancreaticresection",
    "pre_index_date_pancreaticcancer",
    "pre_index_date_haemochromatosis",
    "pre_index_date_cysticfibrosis"
  )
  n_pre <- nrow(cohort)
  cohort <- cohort %>%
    dplyr::mutate(
      .type1_reason =
        (rowSums(dplyr::across(
          dplyr::any_of(type1_cols),
          ~ dplyr::coalesce(as.integer(.x), 0L)
        )) > 0) |
        dplyr::coalesce(type1_code_count > 0, FALSE),
      .type3c_reason =
        rowSums(dplyr::across(
          dplyr::any_of(type3c_cols),
          ~ dplyr::coalesce(as.integer(.x), 0L)
        )) > 0,
      .possible_non_t2d = .type1_reason | .type3c_reason
    )
  .msg(
    "  Excluded for type 1 reasons (type 1 code or insulin \u22641y): N = ",
    format(sum(cohort$.type1_reason), big.mark = ",")
  )
  .msg(
    "  Excluded for type 3c reasons (pancreatic/secondary):  N = ",
    format(sum(cohort$.type3c_reason), big.mark = ",")
  )
  cohort <- cohort %>%
    dplyr::filter(!.possible_non_t2d) %>%
    dplyr::select(-.type1_reason, -.type3c_reason, -.possible_non_t2d)
  .msg(
    "  Excluded possible non-type 2 / secondary diabetes:    N = ",
    format(n_pre - nrow(cohort), big.mark = ",")
  )

  # ---------------------------------------------------------------------------
  # Step 2c: Age at diagnosis (QDiabetes valid range: 25 to <85)
  # ---------------------------------------------------------------------------
  n_pre <- nrow(cohort)
  cohort <- cohort %>%
    dplyr::filter(dm_diag_age >= age_min, dm_diag_age < age_max)
  .msg(
    "  Excluded diagnosed <", age_min, " or \u2265", age_max, " years:              N = ",
    format(n_pre - nrow(cohort), big.mark = ",")
  )

  cohort <- cohort %>%
    dplyr::mutate(dx_year = lubridate::year(dm_diag_date))

  .msg("Assessed for QDiabetes eligibility:                      N = ", .n(cohort))

  cohort
}

# ------------------------------------------------------------------------------
# Derived variables
# ------------------------------------------------------------------------------

derive_features_at_diag <- function(cohort) {
  out <- cohort %>%
    mutate(
      pretghdl_ratio = pretriglyceride / prehdl,
      # Patients with gestational diabetes are removed before final merge
      pre_index_date_gdm = 0,
      pre_index_date_diabetes_meds = 0,
      pre_index_date_qresearch_polycystic_ovaries =
        if_else(gender == 1,
          0L,
          pre_index_date_qresearch_polycystic_ovaries
        ),
      # ethnicity (5-cat)
      ethnicity_5cat = dplyr::recode(
        as.character(ethnicity_5cat),
        "0" = "White",
        "1" = "South Asian",
        "2" = "Black",
        "3" = "Other",
        "4" = "Mixed",
        .default = NA_character_
      ),
      ethnicity_5cat = factor(ethnicity_5cat,
        levels = c("White", "South Asian", "Black", "Other", "Mixed")
      ),

      # IMD quintile (from deciles)
      imd_quintile = case_when(
        imd_decile %in% 1:2 ~ "1-2",
        imd_decile %in% 3:4 ~ "3-4",
        imd_decile %in% 5:6 ~ "5-6",
        imd_decile %in% 7:8 ~ "7-8",
        imd_decile %in% 9:10 ~ "9-10",
        TRUE ~ NA_character_
      ),
      imd_quintile = factor(imd_quintile,
        levels = c("1-2", "3-4", "5-6", "7-8", "9-10"),
        ordered = TRUE
      ),

      # age bands (used in Table 1)
      agecat = cut(
        dm_diag_age,
        breaks = c(25, 40, 50, 60, 70, 80, 85),
        labels = c(
          "25\u201339", "40\u201349", "50\u201359",
          "60\u201369", "70\u201379", "80\u201384"
        ),
        right = FALSE
      ),

      # HbA1c category at diagnosis
      hba1c_cat = case_when(
        prehba1c < 58 ~ "<58",
        prehba1c >= 58 & prehba1c < 75 ~ "58\u201374",
        prehba1c >= 75 & prehba1c < 86 ~ "75\u201385",
        prehba1c >= 86 ~ "\u226586",
        TRUE ~ NA_character_
      ),
      hba1c_cat = factor(
        hba1c_cat,
        levels = c("<58", "58\u201374", "75\u201385", "\u226586")
      ),

      # BMI category
      bmi_cat = case_when(
        prebmi < 25 ~ "<25",
        prebmi >= 25 & prebmi < 30 ~ "25-<30",
        prebmi >= 30 & prebmi < 40 ~ "30-<40",
        prebmi >= 40 ~ ">=40",
        TRUE ~ NA_character_
      ),
      bmi_cat = factor(bmi_cat, levels = c("<25", "25-<30", "30-<40", ">=40")),

      # type 1 code present (used for sensitivity)
      type1_any = as.integer(type1_code_count >= 1),

      # non-HDL cholesterol
      prenonhdl = pretotalcholesterol - prehdl
    ) %>%
    # frailty category (age + eFI category)
    mutate(
      frail_elderly_cat = case_when(
        dm_diag_age <= 70 ~ "Age \u2264 70",
        dm_diag_age > 70 & pre_index_date_efi_cat %in% c("fit", "mild") ~ "Non-frail > 70",
        dm_diag_age > 70 & pre_index_date_efi_cat %in% c("moderate", "severe") ~ "Frail > 70",
        TRUE ~ NA_character_
      ),
      frail_elderly_cat = factor(frail_elderly_cat,
        levels = c("Age \u2264 70", "Non-frail > 70", "Frail > 70")
      )
    ) %>%
    # composite CVD (baseline)
    {
      cvd_vars <- paste0(
        "pre_index_date_",
        c(
          "myocardialinfarction", "stroke", "angina", "ihd",
          "pad", "revasc", "tia", "primary_incident_mi"
        )
      )
      mutate(
        .,
        pre_index_date_cvd = if_else(
          rowSums(select(., all_of(cvd_vars)), na.rm = TRUE) > 0,
          1L,
          0L
        )
      )
    } %>%
    # CKD baseline (stage 3a+)
    mutate(
      pre_index_date_ckd = as.integer(
        preckdstage %in% c("stage_3a", "stage_3b", "stage_4", "stage_5")
      )
    ) %>%
    mutate(
      pre_index_date_mace3 = as.integer(
        rowSums(
          dplyr::select(
            .,
            pre_index_date_myocardialinfarction,
            pre_index_date_stroke,
            pre_index_date_incident_mi,
            pre_index_date_incident_stroke
          ),
          na.rm = TRUE
        ) > 0
      )
    ) %>%
    # cause-specific death dates (primary cause)
    mutate(
      pre_index_date_death = 0,
      cv_death_primary_cause_date = if_else(
        cv_death_primary_cause == 1L & !is.na(death_date),
        as.Date(death_date),
        as.Date(NA)
      ),
      cancer_death_primary_cause_date = if_else(
        cancer_death_primary_cause == 1L & !is.na(death_date),
        as.Date(death_date),
        as.Date(NA)
      ),


      # non-CVD primary death
      non_cv_death = if_else(
        !is.na(death_date) &
          (is.na(cv_death_primary_cause) | cv_death_primary_cause != 1L),
        1L, 0L
      ),
      non_cv_death_date = if_else(
        non_cv_death == 1L,
        as.Date(death_date),
        as.Date(NA)
      ),

      # other primary death = not cardiovascular and not cancer
      other_death = if_else(
        !is.na(death_date) &
          (is.na(cv_death_primary_cause) | cv_death_primary_cause != 1L) &
          (is.na(cancer_death_primary_cause) | cancer_death_primary_cause != 1L),
        1L, 0L
      ),
      other_death_date = if_else(
        other_death == 1L,
        as.Date(death_date),
        as.Date(NA)
      ),
      hf_death_primary_cause_date = if_else(
        hf_death_primary_cause == 1L & !is.na(death_date),
        as.Date(death_date),
        as.Date(NA)
      ),

      # MI death date
      mi_death_primary_cause_date = if_else(
        mi_death_primary_cause == 1L & !is.na(death_date),
        as.Date(death_date),
        as.Date(NA)
      ),

      # Stroke death date
      stroke_death_primary_cause_date = if_else(
        stroke_death_primary_cause == 1L & !is.na(death_date),
        as.Date(death_date),
        as.Date(NA)
      )
    ) %>%
    # Flag a primary incident MI before the index date.
    mutate(
      pre_index_date_index_date_primary_incident_mi = as.integer(
        !is.na(pre_index_date_latest_primary_incident_mi)
      )
    ) %>%
    # Pre-index: fatal or non-fatal MI
    mutate(
      pre_index_date_mi_fatal_nonfatal = as.integer(
        rowSums(
          dplyr::select(
            .,
            # pre_index_date_myocardialinfarction,
            pre_index_date_index_date_primary_incident_mi
          ),
          na.rm = TRUE
        ) > 0
      )
    ) %>%
    # Pre-index: fatal or non-fatal stroke
    mutate(
      pre_index_date_stroke_fatal_nonfatal = as.integer(
        rowSums(
          dplyr::select(
            .,
            # pre_index_date_stroke,
            pre_index_date_incident_stroke
          ),
          na.rm = TRUE
        ) > 0
      )
    ) %>%
    # Post-index: first fatal or non-fatal MI date
    mutate(
      post_index_date_first_mi_fatal_nonfatal = pmin_date(
        #  post_index_date_first_myocardialinfarction,
        post_index_date_first_primary_incident_mi,
        mi_death_primary_cause_date
      ),
      post_index_date_mi_fatal_nonfatal = as.integer(
        !is.na(post_index_date_first_mi_fatal_nonfatal)
      )
    ) %>%
    # Post-index: first fatal or non-fatal stroke date
    mutate(
      post_index_date_first_stroke_fatal_nonfatal = pmin_date(
        post_index_date_first_primary_incident_stroke,
        stroke_death_primary_cause_date
      ),
      post_index_date_stroke_fatal_nonfatal = as.integer(
        !is.na(post_index_date_first_stroke_fatal_nonfatal)
      )
    ) %>%
    mutate(
      post_index_date_first_mace3 = pmin_date(
        post_index_date_first_myocardialinfarction,
        post_index_date_first_stroke,
        cv_death_primary_cause_date,
        post_index_date_first_primary_incident_mi,
        post_index_date_first_primary_incident_stroke
      ),
      post_index_date_mace3 = as.integer(!is.na(post_index_date_first_mace3))
    ) %>%
    # Pre-index: fatal or non-fatal heart failure
    mutate(
      pre_index_date_hf_fatal_nonfatal = as.integer(
        rowSums(
          dplyr::select(
            .,
            pre_index_date_primary_hhf,
            pre_index_date_heartfailure
          ),
          na.rm = TRUE
        ) > 0
      )
    ) %>%
    # Post-index: first fatal or non-fatal heart failure date
    mutate(
      post_index_date_first_hf_fatal_nonfatal = pmin_date(
        post_index_date_first_primary_hhf,
        post_index_date_first_heartfailure,
        hf_death_primary_cause_date
      ),
      post_index_date_hf_fatal_nonfatal = as.integer(
        !is.na(post_index_date_first_hf_fatal_nonfatal)
      )
    ) %>%
    # heart failure (baseline: hospital HF or GP HF)
    mutate(
      pre_index_date_hf = as.integer(
        rowSums(
          dplyr::select(
            .,
            pre_index_date_primary_hhf,
            pre_index_date_heartfailure
          ),
          na.rm = TRUE
        ) > 0
      )
    ) %>%
    # heart failure (post-index: earliest hospital HF or GP HF)
    mutate(
      post_index_date_first_hf = pmin_date(
        post_index_date_first_primary_hhf,
        post_index_date_first_heartfailure
      ),
      post_index_date_hf = as.integer(!is.na(post_index_date_first_hf))
    ) %>%
    # insulin timing flags
    mutate(
      pre_index_date_insulin_le1y = as.integer(
        !is.na(dm_diag_insdate) &
          dm_diag_insdate <= dm_diag_date + years(1)
      ),
      dm_diag_insdate_gt1y = if_else(
        !is.na(dm_diag_insdate) & dm_diag_insdate > dm_diag_date + years(1),
        as.Date(dm_diag_insdate),
        as.Date(NA)
      )
    ) %>%
    # ESKD composite date (CKD5 code OR kidney-failure primary death)
    mutate(
      kf_death_event = !is.na(death_date) & (kf_death_primary_cause == 1L),
      post_index_date_first_eskd = pmin_date(
        post_index_date_first_ckd5_code,
        if_else(kf_death_event, as.Date(death_date), as.Date(NA))
      )
    ) %>%
    mutate(post_index_date_eskd = as.integer(!is.na(post_index_date_first_eskd)))


  out
}


build_tx_lines_oha <- function(drug_classes, cohort_df) {
  dx_dates <- cohort_df %>%
    mutate(patid = as.character(patid)) %>%
    select(patid, dm_diag_date)

  tx <- drug_classes %>%
    mutate(
      patid = as.character(patid)
    ) %>%
    left_join(dx_dates, by = "patid")

  pre_meds <- tx %>%
    filter(!is.na(dstartdate), !is.na(dm_diag_date), dstartdate < dm_diag_date) %>%
    distinct(patid) %>%
    mutate(pre_index_date_diabetes_meds_drug = 1L)

  tx_post_diag <- tx %>%
    filter(
      !is.na(dm_diag_date),
      !is.na(dstartdate),
      dstartdate >= dm_diag_date
    )

  first_line <- tx_post_diag %>%
    filter(drugline_all == 1L, drug_instance == 1L) %>%
    group_by(patid) %>%
    summarise(
      treatment_first_line_earliest = min(dstartdate, na.rm = TRUE),
      .groups = "drop"
    )

  second_line <- tx_post_diag %>%
    filter(drugline_all == 2L, drug_instance == 1L) %>%
    group_by(patid) %>%
    summarise(
      treatment_second_line_earliest = min(dstartdate, na.rm = TRUE),
      .groups = "drop"
    )

  first_line %>%
    full_join(second_line, by = "patid") %>%
    full_join(pre_meds, by = "patid")
}

add_tx_lines_to_cohort <- function(cohort_df, tx_lines_df) {
  tx_lines_df <- tx_lines_df %>%
    mutate(patid = as.character(patid))

  cohort_df %>%
    mutate(patid = as.character(patid)) %>%
    left_join(tx_lines_df, by = "patid") %>%
    mutate(
      pre_index_date_diabetes_meds = coalesce(
        pre_index_date_diabetes_meds_drug,
        as.integer(pre_index_date_diabetes_meds),
        0L
      )
    ) %>%
    select(-any_of("pre_index_date_diabetes_meds_drug"))
}


# ------------- wrapper for main pipeline --------------------------

prepare_main_dataset <- function(at_diag_raw,
                                 hba1c_data = NULL,
                                 qdiabetes_risk_horizon = 10L,
                                 qdiabetes_risk_threshold = 5.6,
                                 drug_classes = NULL,
                                 verbose = TRUE) {
  .n <- function(df) format(nrow(df), big.mark = ",")
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  cohort <- define_cohort_at_diag(at_diag_raw,
    hba1c_data = hba1c_data,
    verbose = verbose
  )
  cohort <- derive_features_at_diag(cohort)

  # ---------------------------------------------------------------------------
  # QDiabetes 2018 (10-year risk score + risk category)
  # Individuals with missing or out-of-range input variables (BMI, ethnicity,
  # smoking status, Townsend score) receive NA scores and are excluded.
  # ---------------------------------------------------------------------------
  n_pre <- nrow(cohort)
  cohort <- calculate_qdiabetes_2018(cohort,
    surv_years = qdiabetes_risk_horizon,
    threshold  = qdiabetes_risk_threshold
  ) %>%
    dplyr::filter(!is.na(qdiabetes_10y_risk_score))
  .msg(
    "  Excluded missing/out-of-range QDiabetes variables:    N = ",
    format(n_pre - nrow(cohort), big.mark = ",")
  )
  .msg("With valid QDiabetes risk score:                         N = ", .n(cohort))

  if (!is.null(qdiabetes_risk_threshold)) {
    n_low_risk <- sum(cohort$qdiabetes_risk_cat == "low", na.rm = TRUE)
    n_high_risk <- sum(cohort$qdiabetes_risk_cat == "high", na.rm = TRUE)
    .msg(
      "  Low QDiabetes risk score (\u2264", qdiabetes_risk_threshold, "%):              N = ",
      format(n_low_risk, big.mark = ",")
    )
    .msg(
      "  High QDiabetes risk score (>", qdiabetes_risk_threshold, "%):               N = ",
      format(n_high_risk, big.mark = ",")
    )
  }

  # QRISK2 supplementary risk measure
  cohort <- add_qrisk2(cohort, surv_years = qdiabetes_risk_horizon)

  if (!is.null(drug_classes)) {
    tx_lines <- build_tx_lines_oha(drug_classes, cohort_df = cohort)
    cohort <- add_tx_lines_to_cohort(cohort, tx_lines)

    .msg(
      "  With first-line treatment date:                       N = ",
      format(sum(!is.na(cohort$treatment_first_line_earliest)), big.mark = ",")
    )
    .msg(
      "  With second-line treatment after first-line:          N = ",
      format(sum(!is.na(cohort$treatment_second_line_earliest)), big.mark = ",")
    )
  }

  cohort
}
