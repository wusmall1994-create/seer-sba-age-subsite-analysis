#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("SBA_PROJECT_ROOT")
if (!nzchar(project_root)) {
  stop("Supply the authorized SEER analysis project root as the first argument or SBA_PROJECT_ROOT.")
}
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(project_root, "results", "ccc_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(project_root, "01_SEER8_annual_age_subsite_rates.txt")
annual <- read_tsv(source_file, show_col_types = FALSE, na = c("NA", ""))

required <- c(
  "SBA_SITE_MAIN", "YEAR_SINGLE", "AGE_3GROUPS_RATE", "Age-Adjusted Rate",
  "Standard Error", "Count"
)
stopifnot(all(required %in% names(annual)))

# Direct standardization is linear. The component rates are mutually exclusive
# and use the same standard population, so their annual age-adjusted rates can
# be summed. Under independent Poisson counts, the component variances add.
# SEER leaves the standard error blank for zero-count cells; these cells
# contribute zero to the summed point estimate and variance.
annual <- annual %>%
  mutate(
    rate = case_when(
      !is.na(`Age-Adjusted Rate`) ~ `Age-Adjusted Rate`,
      is.na(`Age-Adjusted Rate`) & Count == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    se = case_when(
      !is.na(`Standard Error`) ~ `Standard Error`,
      is.na(`Standard Error`) & Count == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

if (any(is.na(annual$rate)) || any(is.na(annual$se))) {
  stop("Non-zero cells with missing rates or standard errors were found.")
}

aggregate_series <- function(data, definition, exclude_nos) {
  selected <- if (exclude_nos) {
    filter(data, SBA_SITE_MAIN != "Small intestine, NOS")
  } else {
    data
  }

  selected %>%
    group_by(AGE_3GROUPS_RATE, YEAR_SINGLE) %>%
    summarise(
      rate = sum(rate),
      se = sqrt(sum(se^2)),
      cases = sum(Count),
      .groups = "drop"
    ) %>%
    mutate(definition = definition)
}

series <- bind_rows(
  aggregate_series(annual, "All SBA sites, including NOS", FALSE),
  aggregate_series(annual, "Anatomically specified SBA, excluding NOS", TRUE)
)

fit_full_period <- function(z) {
  stopifnot(all(z$rate > 0), all(z$se > 0))
  fit <- lm(log(rate) ~ YEAR_SINGLE, data = z, weights = (rate / se)^2)
  co <- summary(fit)$coefficients["YEAR_SINGLE", ]
  beta <- unname(co["Estimate"])
  beta_se <- unname(co["Std. Error"])
  tibble(
    start_year = min(z$YEAR_SINGLE),
    end_year = max(z$YEAR_SINGLE),
    total_cases = sum(z$cases),
    annual_change_percent = 100 * (exp(beta) - 1),
    ci_low = 100 * (exp(beta - 1.96 * beta_se) - 1),
    ci_high = 100 * (exp(beta + 1.96 * beta_se) - 1),
    p_value = unname(co["Pr(>|t|)"]),
    early_5y_rate_per_million = mean(head(z$rate, 5)) * 10,
    late_5y_rate_per_million = mean(tail(z$rate, 5)) * 10,
    absolute_change_per_million = (mean(tail(z$rate, 5)) - mean(head(z$rate, 5))) * 10
  )
}

results <- series %>%
  arrange(definition, AGE_3GROUPS_RATE, YEAR_SINGLE) %>%
  group_by(definition, AGE_3GROUPS_RATE) %>%
  group_modify(~fit_full_period(.x)) %>%
  ungroup() %>%
  mutate(
    age_group = sub("Age ", "", AGE_3GROUPS_RATE),
    model = "Weighted log-linear sensitivity model"
  ) %>%
  select(
    definition, age_group, start_year, end_year, total_cases,
    annual_change_percent, ci_low, ci_high, p_value,
    early_5y_rate_per_million, late_5y_rate_per_million,
    absolute_change_per_million, model
  )

write_csv(series, file.path(out_dir, "anatomically_specified_annual_rates.csv"))
write_csv(results, file.path(out_dir, "anatomically_specified_full_period_trends.csv"))

cat("Anatomically specified incidence sensitivity outputs:", normalizePath(out_dir, winslash = "/"), "\n")

