# A risk model-based approach to identifying atypical type 2 diabetes at diagnosis using routine clinical features: a matched UK cohort study

**Analysis code repository**

## Study overview
This repostiry contains the R code for "A risk model-based approach to identifying atypical type 2 diabetes at diagnosis using routine clinical features: a matched UK cohort study". This study uses QDiabetes to define a type 2 diabete phenotype at diagnosis. Adults newly diagnossed with T2D between January 2013 to March 2023 are classified as:

- **Atypical T2D:** QDiabetes 10-year risk score ≤5.6%.
- **Typical T2D:** QDiabetes 10-year risk score >5.6%.

Eligibility requires a recorded T2D diagnosis and HbA1c evidence of diabetes (at-diagnosis HbA1c ≥48 mmol/mol, or another HbA1c record ≥48 mmol/mol when the at-diagnosis value is missing). People with evidence of type 1 diabetes, type 3c or specified secondary diabetes are excluded. A valid QDiabetes score is required. Cohort construction details are implemented in `R/cohort_definition.R` and `R/calculators/qdiabetes_2018.R`.

## Outcomes

### Primary outcomes

Severe complications after T2D diagnosis, identified using primary care, hospital and mortality records:

**Microvascular**
- Severe retinopathy: vitreous or pre-retinal haemorrhage, photocoagulation, proliferative retinopathy, blindness or visual impairment.
- Severe neuropathy: foot ulcer, gastroparesis, Charcot foot, painful peripheral neuropathy, neuropathic pain or lower-limb amputation.
- Severe nephropathy: sustained ≥40% decline in eGFR, end-stage kidney disease (renal replacement therapy or sustained eGFR <15 mL/min/1.73 m²), or death from kidney failure.

**Macrovascular**
- Fatal or non-fatal myocardial infarction, stroke and heart failure.

### Secondary outcomes

- **Acute:** diabetic ketoacidosis and severe hypoglycaemia.
- **Mortality:** cardiovascular and non-cardiovascular mortality.
- **Treatment initiation:** first-line therapy, second-line therapy and insulin.

Glucose-lowering therapies were identified from primary care prescribing records and included metformin, sulphonylureas, DPP-4 inhibitors, SGLT2 inhibitors, GLP-1 receptor agonists, thiazolidinediones, glinides and acarbose.

## Statistical analysis

Each person with atypical T2D was matched to up to four people with typical T2D, with replacement, exactly on age at diagnosis, sex and calendar year. Follow-up began at diagnosis and continued until the outcome, practice deregistration, death, or the end of available follow-up (March 2023). Ten-year cumulative incidence was estimated using Kaplan–Meier methods. Cox models adjusted for matching variables, ethnicity, deprivation and baseline HbA1c; missing baseline HbA1c was multiply imputed. Treatment-outcome models omitted baseline HbA1c.

## Repository layout

| Path | Purpose |
|---|---|
| `analysis/01_main_analysis.R` | Cohort construction, matching, primary figures, Table 1 and Table 2. |
| `analysis/02_supplementary_analysis.R` | Supplementary acute/mortality and treatment figures, and subgroup analyses. |
| `R/setup.R` | Package checks, CPRD connection and data loaders. |
| `R/cohort_definition.R` | Cohort eligibility and derived variables. |
| `R/calculators/` | QDiabetes and QRISK2 calculations. |
| `R/matching.R` | Exact matching. |
| `R/survival_outcomes.R` | Outcome specifications, censoring and survival-data construction. |
| `R/cox_models.R` | Cox models and HbA1c multiple imputation. |
| `R/tables.R`, `R/plotting.R` | Tables and figures. |
| `outputs/figures/` | Main-analysis figures. |
| `outputs/tables/` | Main-analysis tables |
| `outputs/supplementary/figures/` | Supplementary figures |
| `outputs/supplementary/tables/` | Supplementary subgroup results and counts. |

## Software requirements

The analyses were conducted using R 4.4.0, with `mice` 3.18.0 and `tidyverse` 2.0.0. Required package names are listed in `R/setup.R`; running an analysis script reports any missing packages.
