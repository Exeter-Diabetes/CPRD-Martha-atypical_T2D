# R/setup.R
# ---------------------------------------------------------------------------
# Packages, CPRD database connection, and data loaders.
# ---------------------------------------------------------------------------

required_packages <- c(
  "aurum", "dplyr", "tidyr", "purrr", "EHRBiomarkr", "lubridate",
  "tableone", "survey", "MatchIt", "cobalt", "survival", "broom",
  "ggplot2", "survminer", "cowplot", "scales", "stringr", "tibble",
  "QDiabetes", "mice"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running the analysis.",
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)


# ---- CPRD connection --------------------------------------------------------

connect_cprd <- function(
    cprd_env = Sys.getenv("CPRD_ENV", unset = "diabetes-2024"),
    yaml_conf = Sys.getenv("AURUM_CONFIG", unset = "~/.aurum.yaml")) {
  aurum::CPRDData$new(cprdEnv = cprd_env, cprdConf = yaml_conf)
}


# ---- Data loaders -----------------------------------------------------------

# At-diagnosis cohort (one row per patient at type 2 diabetes diagnosis).
load_at_diagnosis_cohort <- function(cprd,
                                     analysis_name = "at_diag",
                                     table_name = "final_20260813") {
  analysis <- cprd$analysis(analysis_name)
  analysis$cached(name = table_name) %>%
    dplyr::collect()
}

# All cleaned HbA1c records for a set of patients (long: patid, date, testvalue).
load_hba1c_data <- function(cprd,
                            patient_ids,
                            analysis_name = "all_patid",
                            table_name = "all_patid_hba1c_clean_medcodes") {
  analysis <- cprd$analysis(analysis_name)

  cat(paste("Loading HbA1c data for", length(patient_ids), "patients...\n"))

  hba1c_data <- analysis$cached(name = table_name) %>%
    dplyr::filter(patid %in% patient_ids) %>%
    dplyr::select(patid, date, testvalue) %>%
    dplyr::collect() %>%
    dplyr::mutate(patid = as.character(patid))

  cat(sprintf("\u2713 Loaded %d HbA1c records\n\n", nrow(hba1c_data)))

  hba1c_data
}

# Glucose-lowering drug start/stop episodes (one row per drug episode).
load_drug_start_stop <- function(cprd,
                                 analysis_name = "mm",
                                 table_name = "drug_start_stop") {
  analysis <- cprd$analysis(analysis_name)
  analysis$cached(name = table_name) %>%
    dplyr::collect()
}
