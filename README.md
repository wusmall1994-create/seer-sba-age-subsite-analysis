# SEER Small Bowel Adenocarcinoma Analysis Code

This repository contains the statistical analysis code for an observational study of age- and anatomic-subsite-specific incidence trends in small bowel adenocarcinoma using SEER data. It includes cohort construction, incidence and Joinpoint preparation, descriptive analyses, survival analyses, sensitivity analyses, and figure generation.

## Data access

No manuscript files, patient-level records, SEER exports, derived datasets, or analysis results are included. SEER Research Data are available to registered users through the US National Cancer Institute under the SEER data-use agreement. Users must obtain their own authorized data and reproduce the required SEER*Stat exports.

The analysis used SEER*Stat 9.0.43.0. The scripts expect the export filenames documented below and do not download or redistribute SEER data.

## Repository structure

All executable analysis files are in `R/`:

1. `01_clean_define_cohorts.R` cleans the SEER 17 case listing and defines the descriptive, stage-harmonized, and survival cohorts.
2. `02_baseline_table1.R` produces cohort characteristics and missingness summaries.
3. `03_survival_competing_risk.R` fits overall- and cancer-specific survival models.
4. `04_prepare_joinpoint_inputs.R` validates and prepares incidence files for the NCI Joinpoint Regression Program.
5. `05_survival_nonPH_robust.R` fits interval-specific Cox models and initial RMST analyses.
6. `06_build_manuscript_outputs.R` creates analysis tables and figures from preceding results.
7. `07_enhancement_analyses.R` runs age-stratified, time-varying, weighted RMST, and imputation analyses.
8. `08_cebp_sensitivity_analyses.R` runs COVID-era, multiplicity, stage, treatment, and other prespecified sensitivity analyses.
9. `09_cohort_flow_figure.R` generates the analytic-cohort flow diagram.
10. `10_surveillance_revision_analyses.R` generates the age-by-subsite surveillance figure, descriptive absolute incidence-rate changes and confidence intervals, follow-up support summaries, restricted-year RMST analysis, and strict-histology sensitivity analyses.
11. `11_anatomically_specified_incidence_sensitivity.R` combines mutually exclusive, identically standardized subsite rates to compare full-period age-specific trends for all SBA with trends after excluding small-intestine NOS tumors.
12. `12_terminal_year_and_reporting_sensitivities.R` recalculates fixed-structure trends and descriptive absolute rate changes after truncation at 2021 and 2019, and extracts official segment-specific APC and sex-stratified AAPC results for reporting.
13. `13_30_39_definition_bridge.R` compares matched SEER 8 annual series at ages 30–39 years for all malignant histologies, adenocarcinoma including NOS, and adenocarcinoma excluding NOS in the 2010–2019 and 2001–2021 windows. It retains zero-count years with a common Poisson population-offset model, inflates variance for Pearson overdispersion, and calculates descriptive absolute changes from exported age-adjusted rates.

## Required local inputs

Place the following authorized SEER*Stat exports in a local project directory. These files are deliberately excluded by `.gitignore`:

- `01_SEER8_annual_age_subsite_rates.txt`
- `02_SEER8_annual_sex_subsite_rates.txt`
- `05_SEER17_annual_race_subsite_rates.txt`
- `07_SEER17_SBA_master_case_listing.txt`
- `00_SEER8_histology_code_audit_1975_2023.txt`
- `09A_SEER8_annual_30_39_allhist_including_NOS.txt`
- `09B_SEER8_annual_30_39_adeno_including_NOS.txt`
- `09C_SEER8_annual_30_39_adeno_excluding_NOS.txt`

Joinpoint analyses additionally require the exported APC and AAPC files referenced by scripts 07, 08, and 10. The NCI Joinpoint Regression Program is separate software and is not bundled here.

## Running the analysis

Scripts accept paths through command-line arguments or the environment variables `SBA_INPUT_DIR`, `SBA_CLEAN_DIR`, `SBA_RESULTS_DIR`, `SBA_PROJECT_ROOT`, and `SBA_FIGURE_DIR`. A typical sequence is:

```powershell
Rscript R/01_clean_define_cohorts.R "D:/authorized/seer_exports" "D:/sba_project/cleaned"
Rscript R/02_baseline_table1.R "D:/sba_project/cleaned" "D:/sba_project/results/table1"
Rscript R/03_survival_competing_risk.R "D:/sba_project/cleaned" "D:/sba_project/results/survival"
Rscript R/04_prepare_joinpoint_inputs.R "D:/authorized/seer_exports" "D:/sba_project/results/joinpoint_input"
Rscript R/05_survival_nonPH_robust.R "D:/sba_project/cleaned" "D:/sba_project/results/survival_robust"
Rscript R/06_build_manuscript_outputs.R "D:/sba_project/results"
Rscript R/07_enhancement_analyses.R "D:/sba_project"
Rscript R/08_cebp_sensitivity_analyses.R "D:/sba_project"
Rscript R/09_cohort_flow_figure.R "D:/sba_project/results/manuscript"
Rscript R/10_surveillance_revision_analyses.R "D:/sba_project" "D:/sba_project/results/cebp_surveillance_revision"
Rscript R/11_anatomically_specified_incidence_sensitivity.R "D:/sba_project" "D:/sba_project/results/ccc_sensitivity"
Rscript R/12_terminal_year_and_reporting_sensitivities.R "D:/sba_project" "D:/sba_project/results/ccc_terminal_year"
Rscript R/13_30_39_definition_bridge.R "D:/authorized/seer_exports" "D:/sba_project/results/ccc_definition_bridge"
```

## Software dependencies

The code uses R packages `broom`, `cmprsk`, `dplyr`, `ggplot2`, `mice`, `patchwork`, `ragg`, `readr`, `scales`, `stringr`, `survival`, `svglite`, and `tidyr`, plus base R packages including `grid`, `grDevices`, `splines`, and `tools`.

## Reproducibility boundary

The repository supports transparent inspection of the statistical workflow. Exact reproduction requires authorized access to the same SEER submissions, the documented SEER*Stat exports, and the corresponding Joinpoint output files. The strict ICD-O-3 8140/3 trend analysis is a count-based population-offset sensitivity analysis and is not numerically interchangeable with the age-adjusted Joinpoint estimates. The terminal-year sensitivity preserves the full-period selected model structure; it does not reselect joinpoints after truncation. The 30–39-year definition bridge uses separate SEER*Stat exports with identical population denominators; it is a SEER 8 classification sensitivity analysis, not a reproduction of nationwide USCS estimates.

## License

The code is released under the MIT License. SEER data remain subject to the NCI SEER data-use agreement and are not covered by this software license.
