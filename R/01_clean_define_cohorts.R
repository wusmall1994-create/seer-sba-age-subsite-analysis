#!/usr/bin/env Rscript

# Small bowel adenocarcinoma (SBA) SEER cleaning and cohort-definition pipeline
# Input:  SEER*Stat Case Listing export 07_SEER17_SBA_master_case_listing.txt
# Output: cleaned master data, prespecified cohorts, flow/QC tables, codebook

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_INPUT_DIR")
if (!nzchar(input_dir)) stop("Supply the authorized SEER export directory as the first argument or SBA_INPUT_DIR.")
output_dir <- if (length(args) >= 2) args[[2]] else Sys.getenv("SBA_CLEAN_DIR", file.path(input_dir, "cleaned"))
input_file <- file.path(input_dir, "07_SEER17_SBA_master_case_listing.txt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing input file: ", input_file)

# Locked protocol definitions -------------------------------------------------
SBA_HISTOLOGY <- c(8140L,8143L,8144L,8145L,8210L,8211L,8220L,8221L,8255L,
                   8260L,8261L,8262L,8263L,8310L,8480L,8481L,8490L)
SBA_SITES <- c("C17.0-Duodenum","C17.1-Jejunum","C17.2-Ileum",
               "C17.3-Meckels diverticulum","C17.8-Overlapping lesion of small intestine",
               "C17.9-Small intestine, NOS")
RACE_LEVELS <- c("Non-Hispanic White","Non-Hispanic Black",
                 "Non-Hispanic American Indian/Alaska Native",
                 "Non-Hispanic Asian or Pacific Islander","Hispanic (All Races)",
                 "Non-Hispanic Unknown Race")

clean_missing <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "Blank(s)")] <- NA_character_
  x
}
parse_age <- function(x) {
  x <- as.character(x)
  ifelse(str_detect(x, "^90\\+"), 90,
         suppressWarnings(as.integer(str_extract(x, "^[0-9]+"))))
}
site4 <- function(x) case_when(
  x == "C17.0-Duodenum" ~ "Duodenum",
  x %in% c("C17.1-Jejunum","C17.2-Ileum") ~ "Jejunum/Ileum",
  x %in% c("C17.3-Meckels diverticulum","C17.8-Overlapping lesion of small intestine") ~ "Other specified",
  x == "C17.9-Small intestine, NOS" ~ "NOS",
  TRUE ~ NA_character_)
stage4 <- function(x) case_when(
  x == "Localized only" ~ "Localized",
  str_detect(x, "^Regional") ~ "Regional",
  x == "Distant site(s)/node(s) involved" ~ "Distant",
  x == "Unknown/unstaged/unspecified/DCO" ~ "Unknown/unstaged",
  TRUE ~ NA_character_)
grade4 <- function(old, clinical, pathological) {
  z <- coalesce(clean_missing(pathological), clean_missing(clinical), clean_missing(old))
  case_when(
    str_detect(z, "Well differentiated|Grade I(?!I)") ~ "Grade I",
    str_detect(z, "Moderately differentiated|Grade II(?!I)") ~ "Grade II",
    str_detect(z, "Poorly differentiated|Grade III") ~ "Grade III",
    str_detect(z, "Undifferentiated|anaplastic|Grade IV") ~ "Grade IV",
    TRUE ~ "Unknown")
}

# Import and retain every source field ---------------------------------------
raw <- read_tsv(input_file, show_col_types = FALSE, na = "NA", quote = '"',
                progress = FALSE, name_repair = "minimal")
required <- c("Patient ID","Year of diagnosis","Survival months","Sex",
              "Age recode with single ages and 90+","Primary Site - labeled",
              "Histologic Type ICD-O-3","Behavior code ICD-O-3","Diagnostic Confirmation",
              "Combined Summary Stage with Expanded Regional Codes (2004+)",
              "RX Summ--Surg Prim Site (1998-2022)","Reason no cancer-directed surgery",
              "Chemotherapy recode (yes, no/unk)","Vital status recode (study cutoff used)",
              "SEER cause-specific death classification","Sequence number","Type of Reporting Source")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse=", "))

dat <- raw %>% mutate(
  record_id = row_number(),
  patient_id = as.character(`Patient ID`),
  diagnosis_year = as.integer(`Year of diagnosis`),
  age = parse_age(`Age recode with single ages and 90+`),
  age_topcoded_90 = str_detect(`Age recode with single ages and 90+`, "^90\\+"),
  age_group3 = factor(case_when(age <= 49 ~ "20-49", age <= 69 ~ "50-69", TRUE ~ "70+"),
                      levels=c("20-49","50-69","70+")),
  early_onset = age >= 20 & age <= 49,
  sex = factor(Sex, levels=c("Female","Male")),
  race_ethnicity = factor(`Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`, levels=RACE_LEVELS),
  race_known = !is.na(race_ethnicity) & race_ethnicity != "Non-Hispanic Unknown Race",
  site_code = str_extract(`Primary Site - labeled`, "^C17\\.[0-9]"),
  site_main = factor(site4(`Primary Site - labeled`), levels=c("Duodenum","Jejunum/Ileum","Other specified","NOS")),
  histology_code = as.integer(`Histologic Type ICD-O-3`),
  histology_group = factor(case_when(
    histology_code %in% c(8480L,8481L) ~ "Mucinous adenocarcinoma",
    histology_code == 8490L ~ "Signet-ring cell carcinoma",
    histology_code %in% c(8260L:8263L) ~ "Papillary adenocarcinoma",
    TRUE ~ "Other adenocarcinoma")),
  malignant = `Behavior code ICD-O-3` == "Malignant",
  microscopically_confirmed = str_detect(`Diagnostic Confirmation`, "Positive histology|Positive exfoliative cytology|Positive microscopic confirm"),
  stage4 = factor(stage4(`Combined Summary Stage with Expanded Regional Codes (2004+)`),
                  levels=c("Localized","Regional","Distant","Unknown/unstaged")),
  grade4 = factor(grade4(`Grade Recode (thru 2017)`,`Grade Clinical (2018+)`,`Grade Pathological (2018+)`),
                  levels=c("Grade I","Grade II","Grade III","Grade IV","Unknown")),
  tumor_size_raw = clean_missing(`Tumor Size Summary (2016+)`),
  nodes_examined_raw = clean_missing(`Regional nodes examined (1988+)`),
  nodes_positive_raw = clean_missing(`Regional nodes positive (1988+)`),
  surgery = factor(case_when(
    str_detect(`Reason no cancer-directed surgery`, "^Surgery performed$") ~ "Yes",
    str_detect(`Reason no cancer-directed surgery`, "Not recommended|Recommended but not performed|Not performed") ~ "No",
    TRUE ~ "Unknown"), levels=c("No","Yes","Unknown")),
  surgery_available = diagnosis_year <= 2022 & !is.na(clean_missing(`RX Summ--Surg Prim Site (1998-2022)`)),
  chemotherapy = factor(if_else(`Chemotherapy recode (yes, no/unk)` == "Yes", "Yes", "No/unknown"),
                        levels=c("No/unknown","Yes")),
  radiation = factor(if_else(`Radiation recode` %in% c("None/Unknown",NA_character_), "None/unknown", "Yes"),
                     levels=c("None/unknown","Yes")),
  survival_months = suppressWarnings(as.numeric(`Survival months`)),
  survival_known = !is.na(survival_months),
  analysis_time_months = if_else(survival_known, pmax(survival_months,0.5), NA_real_),
  os_event = as.integer(`Vital status recode (study cutoff used)` == "Dead"),
  css_event = as.integer(`SEER cause-specific death classification` == "Dead (attributable to this cancer dx)"),
  dco_or_autopsy = str_detect(`Type of Reporting Source`, regex("death certificate|autopsy",ignore_case=TRUE)),
  first_primary = `Sequence number` %in% c("One primary only","1st of 2 or more primaries"),
  multiple_primary_record = !(`Sequence number` == "One primary only")
)

# Membership flags: all are explicit and non-destructive ---------------------
dat <- dat %>% mutate(
  eligible_master = diagnosis_year %in% 2000:2023 & age >= 20 &
                    `Primary Site - labeled` %in% SBA_SITES & histology_code %in% SBA_HISTOLOGY & malignant,
  cohort_descriptive = eligible_master,
  cohort_incidence_race = eligible_master & race_known,
  cohort_stage = eligible_master & diagnosis_year %in% 2004:2023 & !is.na(stage4),
  cohort_survival_primary = eligible_master & diagnosis_year %in% 2004:2022 &
                            !dco_or_autopsy & survival_known & !is.na(stage4),
  cohort_survival_first_primary = cohort_survival_primary & first_primary,
  cohort_survival_micro = cohort_survival_primary & microscopically_confirmed,
  cohort_treatment = cohort_survival_primary & surgery_available,
  cohort_early_onset = eligible_master & early_onset,
  cohort_early_onset_survival = cohort_survival_primary & early_onset
)

# Hard QC: fail rather than silently change the study population -------------
stopifnot(nrow(dat) == 14158L)
stopifnot(all(dat$eligible_master))
stopifnot(sum(dat$cohort_incidence_race) == 14135L)
stopifnot(sum(dat$cohort_stage) == 12391L)
stopifnot(sum(dat$cohort_early_onset) == 1436L)
stopifnot(sum(dat$dco_or_autopsy) == 35L)
stopifnot(all(dat$diagnosis_year %in% 2000:2023))
stopifnot(all(dat$histology_code %in% SBA_HISTOLOGY))
stopifnot(all(dat$`Primary Site - labeled` %in% SBA_SITES))

# Flow, QC, codebook ----------------------------------------------------------
flow <- tibble(
  cohort=c("Master descriptive","Race/ethnicity incidence","Stage analysis",
           "Primary survival","First-primary survival sensitivity",
           "Microscopically confirmed survival sensitivity","Treatment analysis",
           "Early-onset descriptive","Early-onset survival"),
  definition=c(
    "Adult malignant SBA, 2000-2023",
    "Master plus known race/ethnicity",
    "Master, 2004-2023, valid stage4",
    "Master, 2004-2022, valid stage4, known survival, excluding DCO/autopsy",
    "Primary survival plus one/first primary",
    "Primary survival plus microscopic confirmation",
    "Primary survival plus surgery field available",
    "Master, age 20-49",
    "Primary survival, age 20-49"),
  n=c(sum(dat$cohort_descriptive),sum(dat$cohort_incidence_race),sum(dat$cohort_stage),
      sum(dat$cohort_survival_primary),sum(dat$cohort_survival_first_primary),
      sum(dat$cohort_survival_micro),sum(dat$cohort_treatment),
      sum(dat$cohort_early_onset),sum(dat$cohort_early_onset_survival)))

qc <- bind_rows(
  tibble(check="Input MD5", value=unname(tools::md5sum(input_file)), expected="Recorded"),
  tibble(check="Input records",value=as.character(nrow(dat)),expected="14158"),
  tibble(check="Unique patient IDs",value=as.character(n_distinct(dat$patient_id)),expected="14063"),
  tibble(check="Unknown race",value=as.character(sum(!dat$race_known)),expected="23"),
  tibble(check="Stage cohort",value=as.character(sum(dat$cohort_stage)),expected="12391"),
  tibble(check="Early-onset records",value=as.character(sum(dat$cohort_early_onset)),expected="1436"),
  tibble(check="DCO/autopsy records",value=as.character(sum(dat$dco_or_autopsy)),expected="35"),
  tibble(check="2023 surgery raw missing",value=as.character(sum(dat$diagnosis_year==2023 & is.na(clean_missing(dat$`RX Summ--Surg Prim Site (1998-2022)`)))),expected="799")
) %>% mutate(status=if_else(value==expected | expected=="Recorded","PASS","FAIL"))

codebook <- tribble(
  ~variable,~definition,
  "record_id","Sequential tumor-record identifier generated at import",
  "patient_id","SEER patient identifier; may repeat for multiple primary tumors",
  "age","Age in years; 90 denotes top-coded 90+",
  "early_onset","Age 20-49 years",
  "site_main","Duodenum; Jejunum/Ileum; Other specified; NOS",
  "histology_group","Mucinous; signet-ring; papillary; other adenocarcinoma",
  "stage4","Localized; regional; distant; unknown/unstaged, valid 2004+",
  "grade4","Pathological 2018+ preferred, then clinical 2018+, then grade through 2017",
  "surgery","Derived primarily from Reason no cancer-directed surgery",
  "surgery_available","Surgery field nonblank and diagnosis year <=2022",
  "analysis_time_months","Survival months with zero months set to 0.5",
  "os_event","All-cause death indicator",
  "css_event","SEER cause-specific death indicator",
  "first_primary","One primary only or first of multiple primaries",
  "cohort_survival_primary","Diagnosis 2004-2022, valid stage and survival, excludes DCO/autopsy"
)

# Write R-native and analysis-friendly flat files -----------------------------
saveRDS(dat, file.path(output_dir,"sba_master_clean.rds"), compress="xz")
write_tsv(dat, file.path(output_dir,"sba_master_clean.tsv"), na="NA")
write_tsv(dat %>% filter(cohort_descriptive), file.path(output_dir,"cohort_descriptive.tsv"), na="NA")
write_tsv(dat %>% filter(cohort_stage), file.path(output_dir,"cohort_stage_2004_2023.tsv"), na="NA")
write_tsv(dat %>% filter(cohort_survival_primary), file.path(output_dir,"cohort_survival_primary_2004_2022.tsv"), na="NA")
write_tsv(dat %>% filter(cohort_treatment), file.path(output_dir,"cohort_treatment_2004_2022.tsv"), na="NA")
write_csv(flow, file.path(output_dir,"cohort_flow.csv"), na="NA")
write_csv(qc, file.path(output_dir,"quality_checks.csv"), na="NA")
write_csv(codebook, file.path(output_dir,"derived_variable_codebook.csv"), na="NA")

manifest <- c(
  paste0("Generated: ", Sys.time()),
  paste0("Input: ", normalizePath(input_file, winslash="/")),
  paste0("Input MD5: ", unname(tools::md5sum(input_file))),
  paste0("R version: ", R.version.string),
  "Analysis unit: tumor record",
  "Original source columns are preserved; derived variables and cohort flags are appended.",
  "Never interpret surgery blanks in 2023 as no surgery.",
  "Treatment variables represent observed associations, not causal effects."
)
writeLines(manifest,file.path(output_dir,"MANIFEST.txt"),useBytes=TRUE)

if (any(qc$status != "PASS")) stop("One or more quality checks failed. See quality_checks.csv")
cat("Cleaning complete. Output:", normalizePath(output_dir,winslash="/"), "\n")
print(flow)
print(qc)
