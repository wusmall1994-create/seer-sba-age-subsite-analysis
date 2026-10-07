#!/usr/bin/env Rscript

# Terminal-year, sex-stratified, and segment-reporting sensitivities.
#
# This script uses only aggregate SEER*Stat rate exports and official Joinpoint
# output files. It does not contain or redistribute patient-level SEER data.
# The 2019/2021 terminal-year estimates preserve the selected full-period model
# structure (including the 1995 duodenal age 50-69 joinpoint) and are therefore
# sensitivity estimates, not newly selected Joinpoint models.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_PROJECT_ROOT", unset = ".")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(project_root, "results", "ccc_terminal_year")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_dir <- Sys.getenv("SBA_INPUT_DIR", unset = project_root)
jp_dir <- file.path(project_root, "results", "joinpoint_input")

age_file <- file.path(input_dir, "01_SEER8_annual_age_subsite_rates.txt")
sex_file <- file.path(input_dir, "02_SEER8_annual_sex_subsite_rates.txt")
apc_file <- file.path(jp_dir, "JP01_SEER8_age_subsite_1975_2023.Export.APC.txt")
sex_aapc_file <- file.path(jp_dir, "JP02_SEER8_sex_subsite_1975_2023.Export.AAPC.txt")

needed <- c(age_file, sex_file, apc_file, sex_aapc_file)
missing_files <- needed[!file.exists(needed)]
if (length(missing_files)) stop("Missing required files: ", paste(missing_files, collapse = ", "))

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
age <- read_tab(age_file)
sex <- read_tab(sex_file)

clean_rates <- function(d) {
  names(d)[names(d) == "Age-Adjusted Rate"] <- "rate"
  names(d)[names(d) == "Standard Error"] <- "se"
  names(d)[names(d) == "YEAR_SINGLE"] <- "year"
  d$year <- as.integer(d$year)
  d$rate <- as.numeric(d$rate)
  d$se <- as.numeric(d$se)
  d[is.finite(d$rate) & d$rate > 0 & is.finite(d$se) & d$se > 0, , drop = FALSE]
}

age <- clean_rates(age)
sex <- clean_rates(sex)

fit_annual_change <- function(d, hinge_year = NA_integer_) {
  d$w <- (d$rate / d$se)^2
  if (is.na(hinge_year)) {
    fit <- lm(log(rate) ~ year, data = d, weights = w)
    b <- coef(fit)[["year"]]
    v <- vcov(fit)["year", "year"]
  } else {
    d$after_hinge <- pmax(d$year - hinge_year, 0)
    fit <- lm(log(rate) ~ year + after_hinge, data = d, weights = w)
    start <- min(d$year)
    finish <- max(d$year)
    post_fraction <- max(finish - hinge_year, 0) / (finish - start)
    contrast <- c(0, 1, post_fraction)
    b <- sum(contrast * coef(fit))
    v <- as.numeric(t(contrast) %*% vcov(fit) %*% contrast)
  }
  ci <- b + c(-1, 1) * 1.96 * sqrt(v)
  c(estimate = 100 * (exp(b) - 1), lower = 100 * (exp(ci[[1]]) - 1), upper = 100 * (exp(ci[[2]]) - 1))
}

absolute_change <- function(d) {
  d <- d[order(d$year), , drop = FALSE]
  first <- head(d, 5)
  last <- tail(d, 5)
  delta <- mean(last$rate) - mean(first$rate)
  delta_se <- sqrt(sum(last$se^2) / nrow(last)^2 + sum(first$se^2) / nrow(first)^2)
  c(
    estimate = delta,
    lower = delta - 1.96 * delta_se,
    upper = delta + 1.96 * delta_se,
    first_start = min(first$year),
    first_end = max(first$year),
    last_start = min(last$year),
    last_end = max(last$year)
  )
}

one_series <- function(d, cutoff, level, subsite, stratum) {
  d <- d[d$year <= cutoff, , drop = FALSE]
  hinge <- if (level == "Age by subsite" && subsite == "Duodenum" && stratum == "Age 50-69") 1995L else NA_integer_
  trend <- fit_annual_change(d, hinge)
  abs_change <- absolute_change(d)
  data.frame(
    Analysis = level,
    Subsite = subsite,
    Stratum = stratum,
    TerminalYear = cutoff,
    FirstObservedYear = min(d$year),
    LastObservedYear = max(d$year),
    AnnualChange = trend[["estimate"]],
    AnnualChangeLCL = trend[["lower"]],
    AnnualChangeUCL = trend[["upper"]],
    AbsoluteRateChangePerMillion = 10 * abs_change[["estimate"]],
    AbsoluteRateChangePerMillionLCL = 10 * abs_change[["lower"]],
    AbsoluteRateChangePerMillionUCL = 10 * abs_change[["upper"]],
    FirstWindow = paste0(abs_change[["first_start"]], "-", abs_change[["first_end"]]),
    LastWindow = paste0(abs_change[["last_start"]], "-", abs_change[["last_end"]]),
    ModelStructure = if (is.na(hinge)) "log-linear, no joinpoint" else "continuous log-linear, fixed joinpoint at 1995",
    stringsAsFactors = FALSE
  )
}

cutoffs <- c(2023L, 2021L, 2019L)
rows <- list()
k <- 0L
for (cutoff in cutoffs) {
  for (site in unique(age$SBA_SITE_MAIN)) {
    for (grp in unique(age$AGE_3GROUPS_RATE)) {
      k <- k + 1L
      rows[[k]] <- one_series(
        age[age$SBA_SITE_MAIN == site & age$AGE_3GROUPS_RATE == grp, , drop = FALSE],
        cutoff, "Age by subsite", site, grp
      )
    }
  }
  overall <- sex[sex$Sex == "Male and female", , drop = FALSE]
  for (site in unique(overall$SBA_SITE_MAIN)) {
    k <- k + 1L
    rows[[k]] <- one_series(
      overall[overall$SBA_SITE_MAIN == site, , drop = FALSE],
      cutoff, "Overall subsite", site, "All ages"
    )
  }
}
terminal <- do.call(rbind, rows)
terminal[] <- lapply(terminal, function(x) if (is.numeric(x)) round(x, 4) else x)
write.csv(terminal, file.path(output_dir, "terminal_year_trend_and_absolute_change.csv"), row.names = FALSE, na = "")

apc <- read_tab(apc_file)
segment <- data.frame(
  Subsite = apc$Subsite,
  AgeGroup = apc$AgeGroup,
  Segment = apc$Segment + 1L,
  StartYear = apc[["Segment Start"]],
  EndYear = apc[["Segment End"]],
  APC = apc$APC,
  LCL = apc[["APC 95% LCL"]],
  UCL = apc[["APC 95% UCL"]],
  PValue = apc[["P-Value"]],
  Significant = ifelse(apc[["APC Significant"]] == 1, "Yes", "No"),
  stringsAsFactors = FALSE
)
write.csv(segment, file.path(output_dir, "official_segment_specific_apc.csv"), row.names = FALSE, na = "")

sex_aapc <- read_tab(sex_aapc_file)
sex_full <- sex_aapc[sex_aapc[["AAPC Index"]] == "Full Range" & sex_aapc$Sex %in% c("Female", "Male"), , drop = FALSE]
sex_table <- data.frame(
  Subsite = sex_full$Subsite,
  Sex = sex_full$Sex,
  StartYear = sex_full[["Start Obs"]],
  EndYear = sex_full[["End Obs"]],
  AAPC = sex_full$AAPC,
  LCL = sex_full[["AAPC C.I. Low"]],
  UCL = sex_full[["AAPC C.I. High"]],
  PValue = sex_full[["P-Value"]],
  stringsAsFactors = FALSE
)
write.csv(sex_table, file.path(output_dir, "official_sex_stratified_full_period_aapc.csv"), row.names = FALSE, na = "")

cat("Wrote terminal-year sensitivity, segment APC, and sex-stratified AAPC outputs to", normalizePath(output_dir), "\n")
