#!/usr/bin/env Rscript

# Definition-bridging analysis for small-intestine cancer at ages 30-39 years.
#
# Three independently exported SEER 8 annual series are compared within the
# published-study calendar windows:
#   A. all malignant histologies, including C17.9;
#   B. adenocarcinoma, including C17.9;
#   C. adenocarcinoma, excluding C17.9.
#
# Sparse adenocarcinoma series contain zero-count years. A common
# Poisson log-linear count model with a population offset is therefore used
# for all definitions, avoiding deletion or continuity correction of
# zero-count years. Pearson overdispersion inflates the model variance when it
# exceeds 1; underdispersion is not used to narrow confidence intervals.
# Descriptive absolute changes use the exported
# age-adjusted rates and their standard errors.

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("SBA_INPUT_DIR")
output_dir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("SBA_RESULTS_DIR")
if (!nzchar(input_dir)) stop("Supply the authorized SEER*Stat export directory.")
if (!nzchar(output_dir)) output_dir <- file.path(input_dir, "results", "ccc_definition_bridge")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spec <- data.frame(
  definition_id = c("09A", "09B", "09C"),
  definition = c(
    "All malignant histologies, including NOS",
    "Adenocarcinoma, including NOS",
    "Adenocarcinoma, excluding NOS"
  ),
  filename = c(
    "09A_SEER8_annual_30_39_allhist_including_NOS.txt",
    "09B_SEER8_annual_30_39_adeno_including_NOS.txt",
    "09C_SEER8_annual_30_39_adeno_excluding_NOS.txt"
  ),
  stringsAsFactors = FALSE
)

read_one <- function(path, definition_id, definition) {
  if (!file.exists(path)) stop("Missing SEER*Stat export: ", path)
  d <- read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
  required <- c(
    "YEAR_SINGLE", "Age-Adjusted Rate", "Standard Error",
    "Lower Confidence Interval", "Upper Confidence Interval", "Count", "Population"
  )
  if (!all(required %in% names(d))) stop("Unexpected columns in ", path)
  d <- data.frame(
    definition_id = definition_id,
    definition = definition,
    year = as.integer(d[["YEAR_SINGLE"]]),
    rate = as.numeric(d[["Age-Adjusted Rate"]]),
    se = as.numeric(d[["Standard Error"]]),
    count = as.integer(d[["Count"]]),
    population = as.numeric(d[["Population"]]),
    stringsAsFactors = FALSE
  )
  if (!identical(d$year, 1975:2023)) stop("Annual series is not complete for ", definition_id)
  if (any(!is.finite(d$rate)) || any(!is.finite(d$population)) || any(d$population <= 0)) {
    stop("Invalid rate or population value in ", definition_id)
  }
  # SEER leaves the SE undefined for a zero-count annual cell. For a mean of
  # independent annual directly standardized rates, its variance contribution
  # is treated as zero, matching the zero point estimate.
  d$se[is.na(d$se) & d$count == 0L] <- 0
  if (any(is.na(d$se))) stop("Non-zero cell has a missing standard error in ", definition_id)
  d
}

series <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  read_one(
    file.path(input_dir, spec$filename[[i]]),
    spec$definition_id[[i]],
    spec$definition[[i]]
  )
}))

# All definitions must use the same population denominator in every year.
pop_wide <- reshape(
  series[c("definition_id", "year", "population")],
  idvar = "year", timevar = "definition_id", direction = "wide"
)
if (!all(pop_wide$population.09A == pop_wide$population.09B) ||
    !all(pop_wide$population.09A == pop_wide$population.09C)) {
  stop("Population denominators differ across definitions.")
}

fit_window <- function(d, start_year, end_year, comparator) {
  z <- d[d$year >= start_year & d$year <= end_year, , drop = FALSE]
  if (!identical(z$year, start_year:end_year)) stop("Incomplete analysis window.")
  z$year_centered <- z$year - start_year
  fit <- glm(
    count ~ year_centered,
    offset = log(population),
    family = poisson(link = "log"),
    data = z
  )
  beta <- unname(coef(fit)[["year_centered"]])
  pearson_dispersion <- sum(residuals(fit, type = "pearson")^2) / fit$df.residual
  variance_multiplier <- max(1, pearson_dispersion)
  beta_se <- sqrt(vcov(fit)["year_centered", "year_centered"] * variance_multiplier)
  annual <- 100 * (exp(beta) - 1)
  annual_ci <- 100 * (exp(beta + c(-1, 1) * 1.96 * beta_se) - 1)

  first <- z[seq_len(5), , drop = FALSE]
  last <- z[(nrow(z) - 4L):nrow(z), , drop = FALSE]
  first_mean <- mean(first$rate)
  last_mean <- mean(last$rate)
  delta <- last_mean - first_mean
  delta_se <- sqrt(sum(first$se^2) / 25 + sum(last$se^2) / 25)
  delta_ci <- delta + c(-1, 1) * 1.96 * delta_se

  data.frame(
    definition_id = z$definition_id[[1]],
    definition = z$definition[[1]],
    comparator = comparator,
    start_year = start_year,
    end_year = end_year,
    years = nrow(z),
    cases = sum(z$count),
    annual_change_percent = annual,
    annual_change_lcl = annual_ci[[1]],
    annual_change_ucl = annual_ci[[2]],
    p_value = 2 * pnorm(abs(beta / beta_se), lower.tail = FALSE),
    pearson_dispersion = pearson_dispersion,
    variance_multiplier = variance_multiplier,
    first_5y_window = paste0(min(first$year), "-", max(first$year)),
    last_5y_window = paste0(min(last$year), "-", max(last$year)),
    first_5y_rate_per_100k = first_mean,
    last_5y_rate_per_100k = last_mean,
    absolute_change_per_million = 10 * delta,
    absolute_change_per_million_lcl = 10 * delta_ci[[1]],
    absolute_change_per_million_ucl = 10 * delta_ci[[2]],
    model = "Poisson log-linear count model with population offset; variance inflated for Pearson dispersion >1",
    stringsAsFactors = FALSE
  )
}

windows <- list(
  list(start = 2010L, end = 2019L, comparator = "Koh et al. 2023"),
  list(start = 2001L, end = 2021L, comparator = "Yasinzai et al. 2026")
)

results <- do.call(rbind, lapply(windows, function(w) {
  do.call(rbind, lapply(split(series, series$definition_id), function(d) {
    fit_window(d, w$start, w$end, w$comparator)
  }))
}))
results <- results[order(results$start_year, results$definition_id, decreasing = TRUE), ]
row.names(results) <- NULL

write.csv(
  series,
  file.path(output_dir, "age30_39_definition_bridge_annual_series.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  results,
  file.path(output_dir, "age30_39_definition_bridge_results.csv"),
  row.names = FALSE, na = ""
)

qa <- data.frame(
  check = c(
    "Rows per definition",
    "Population denominators identical",
    "Total count 09A",
    "Total count 09B",
    "Total count 09C",
    "Nested annual counts 09A >= 09B >= 09C",
    "Total NOS adenocarcinoma count (09B - 09C)"
  ),
  value = c(
    paste(as.integer(table(series$definition_id)), collapse = ","),
    "TRUE",
    sum(series$count[series$definition_id == "09A"]),
    sum(series$count[series$definition_id == "09B"]),
    sum(series$count[series$definition_id == "09C"]),
    all(
      series$count[series$definition_id == "09A"] >= series$count[series$definition_id == "09B"] &
      series$count[series$definition_id == "09B"] >= series$count[series$definition_id == "09C"]
    ),
    sum(series$count[series$definition_id == "09B"] - series$count[series$definition_id == "09C"])
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(output_dir, "age30_39_definition_bridge_QA.csv"), row.names = FALSE)

cat("Wrote 30-39-year definition-bridging outputs to", normalizePath(output_dir), "\n")
