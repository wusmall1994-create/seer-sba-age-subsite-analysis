#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(splines)
  library(broom)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("SBA_PROJECT_ROOT")
if (!nzchar(root)) {
  stop("Supply the analysis project root as the first command-line argument or SBA_PROJECT_ROOT.")
}
results <- file.path(root, "results")
sens <- file.path(results, "cebp_sensitivity")
out <- if (length(args) >= 2L) args[[2L]] else file.path(results, "cebp_surveillance_revision")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Figure 2: the primary surveillance result on relative and absolute scales.
aapc <- read_csv(file.path(sens, "age_subsite_AAPC_with_BH_FDR.csv"), show_col_types = FALSE)
annual_rates <- read_tsv(file.path(root, "01_SEER8_annual_age_subsite_rates.txt"), show_col_types = FALSE)
absolute <- annual_rates %>%
  transmute(
    subsite = SBA_SITE_MAIN,
    age_group = AGE_3GROUPS_RATE,
    year = YEAR_SINGLE,
    rate_per_million = `Age-Adjusted Rate` * 10,
    se_per_million = `Standard Error` * 10
  ) %>%
  filter(!is.na(rate_per_million), !is.na(se_per_million)) %>%
  group_by(subsite, age_group) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    early_rate_per_million = mean(head(rate_per_million, 5)),
    late_rate_per_million = mean(tail(rate_per_million, 5)),
    absolute_change_per_million = late_rate_per_million - early_rate_per_million,
    absolute_change_se = sqrt(sum(head(se_per_million, 5)^2) / 25 + sum(tail(se_per_million, 5)^2) / 25),
    absolute_change_ci_low = absolute_change_per_million - 1.96 * absolute_change_se,
    absolute_change_ci_high = absolute_change_per_million + 1.96 * absolute_change_se,
    .groups = "drop"
  )
figdat <- aapc %>%
  transmute(
    subsite = Subsite,
    age_group = sub("Age ", "", AgeGroup),
    start_year = `Start Obs`,
    end_year = `End Obs`,
    AAPC,
    ci_low = `AAPC C.I. Low`,
    ci_high = `AAPC C.I. High`,
    p_fdr = p_FDR,
    fdr_significant = significant_FDR
  ) %>%
  left_join(
    absolute %>% mutate(age_group = sub("Age ", "", age_group)),
    by = c("subsite", "age_group")
  ) %>%
  mutate(
    subsite = recode(subsite, `Small intestine, NOS` = "NOS"),
    subsite = factor(subsite, levels = c("Duodenum", "Jejunum/Ileum", "Other specified", "NOS")),
    age_group = factor(age_group, levels = c("20-49", "50-69", "70+")),
    row_id = interaction(subsite, age_group, sep = " | ", lex.order = TRUE),
    row_id = factor(row_id, levels = rev(unique(row_id)))
  )
write_csv(figdat, file.path(out, "Figure2_source_data.csv"))

palette_age <- c("20-49" = "#3B78A8", "50-69" = "#D28A32", "70+" = "#8B3F55")
row_labels <- setNames(
  paste0(
    as.character(figdat$subsite), "   ",
    ifelse(as.character(figdat$age_group) == "70+", "\u226570", as.character(figdat$age_group)),
    " y"
  ),
  as.character(figdat$row_id)
)

theme_pub <- function() {
  theme_classic(base_size = 8, base_family = "DejaVu Sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7.2, colour = "black"),
      legend.position = "none",
      plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 5)),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

p_a <- ggplot(figdat, aes(x = AAPC, y = row_id, colour = age_group)) +
  geom_vline(xintercept = 0, colour = "#8F8F8F", linewidth = 0.35, linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.55) +
  geom_point(aes(shape = fdr_significant), size = 2.1, stroke = 0.5) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  scale_colour_manual(values = palette_age) +
  scale_y_discrete(labels = row_labels) +
  scale_x_continuous(breaks = seq(-6, 4, 2), limits = c(-6, 4)) +
  labs(
    title = "Relative annual change",
    x = "Full-period AAPC, % (95% CI)",
    y = NULL
  ) +
  theme_pub()

p_b <- ggplot(figdat, aes(x = absolute_change_per_million, y = row_id, colour = age_group)) +
  geom_vline(xintercept = 0, colour = "#8F8F8F", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = absolute_change_ci_low, xmax = absolute_change_ci_high), height = 0, linewidth = 0.55) +
  geom_point(size = 2.1) +
  geom_text(
    aes(
      x = absolute_change_per_million,
      label = sprintf("%+.2f", absolute_change_per_million)
    ),
    hjust = 0.5, vjust = -0.85,
    colour = "black", size = 2.1, family = "DejaVu Sans"
  ) +
  scale_colour_manual(values = palette_age) +
  scale_y_discrete(labels = NULL) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(-3.0, 23.0), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Absolute incidence-rate change",
    x = "Absolute rate difference per million (95% CI)",
    y = NULL
  ) +
  theme_pub() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())

fig1 <- p_a + p_b + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

width_in <- 183 / 25.4
height_in <- 118 / 25.4
svglite::svglite(file.path(out, "Figure2_age_subsite_incidence_change.svg"), width = width_in, height = height_in)
print(fig1)
dev.off()
grDevices::cairo_pdf(file.path(out, "Figure2_age_subsite_incidence_change.pdf"), width = width_in, height = height_in, family = "Arial")
print(fig1)
dev.off()
ragg::agg_tiff(file.path(out, "Figure2_age_subsite_incidence_change.tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(fig1)
dev.off()
ragg::agg_png(file.path(out, "Figure2_age_subsite_incidence_change.png"), width = width_in, height = height_in, units = "in", res = 300)
print(fig1)
dev.off()

# Follow-up support for RMST, including a cohort with potential for 10-year follow-up.
dat <- readRDS(file.path(root, "cleaned", "sba_master_clean.rds"))
marital2 <- function(x) case_when(
  x == "Married (including common law)" ~ "Married",
  is.na(x) | x == "Unknown" ~ "Unknown",
  TRUE ~ "Not married"
)
period4 <- function(y) factor(case_when(
  y <= 2008 ~ "2004-2008",
  y <= 2013 ~ "2009-2013",
  y <= 2018 ~ "2014-2018",
  TRUE ~ "2019-2022"
), levels = c("2004-2008", "2009-2013", "2014-2018", "2019-2022"))

d <- dat %>% filter(cohort_survival_primary) %>% transmute(
  record_id,
  diagnosis_year,
  histology_code,
  age,
  time = analysis_time_months,
  event = os_event,
  early = factor(if_else(early_onset, "Age 20-49", "Age 50+"), levels = c("Age 50+", "Age 20-49")),
  site = factor(site_main, levels = c("Duodenum", "Jejunum/Ileum", "Other specified", "NOS")),
  sex = relevel(factor(sex), ref = "Female"),
  race = relevel(factor(case_when(
    race_ethnicity %in% c("Non-Hispanic American Indian/Alaska Native", "Non-Hispanic Unknown Race") ~ "Other/unknown",
    TRUE ~ as.character(race_ethnicity)
  )), ref = "Non-Hispanic White"),
  marital = factor(marital2(`Marital status at diagnosis`), levels = c("Married", "Not married", "Unknown")),
  period = period4(diagnosis_year),
  histology = relevel(factor(histology_group), ref = "Other adenocarcinoma"),
  stage = factor(stage4, levels = c("Localized", "Regional", "Distant", "Unknown/unstaged")),
  grade = factor(grade4)
) %>% filter(complete.cases(.))
stopifnot(nrow(d) == 11563L)

risk_support <- bind_rows(lapply(c(60, 120), function(h) {
  d %>% group_by(early) %>% summarise(
    horizon_months = h,
    N = n(),
    at_risk = sum(time >= h),
    events_before_horizon = sum(event == 1 & time < h),
    censored_before_horizon = sum(event == 0 & time < h),
    .groups = "drop"
  )
}))
write_csv(risk_support, file.path(out, "RMST_followup_support.csv"))

rmst_one <- function(time, event, w, tau) {
  sf <- survfit(Surv(time, event) ~ 1, weights = w, robust = TRUE)
  tt <- c(0, sf$time[sf$time < tau], tau)
  ss <- c(1, sf$surv[sf$time < tau])
  sum(diff(tt) * ss)
}

weighted_rmst <- function(data, tau = 120, B = 1000, seed = 20261003) {
  dd <- droplevels(data)
  dd$z <- as.integer(dd$early == "Age 20-49")
  calc <- function(x) {
    m <- glm(z ~ sex + race + marital + period + histology + stage, data = x, family = binomial())
    ps <- pmin(pmax(predict(m, type = "response"), 0.01), 0.99)
    pz <- mean(x$z)
    w <- ifelse(x$z == 1, pz / ps, (1 - pz) / (1 - ps))
    q <- quantile(w, c(0.01, 0.99))
    w <- pmin(pmax(w, q[1]), q[2])
    c(
      later = rmst_one(x$time[x$z == 0], x$event[x$z == 0], w[x$z == 0], tau),
      early = rmst_one(x$time[x$z == 1], x$event[x$z == 1], w[x$z == 1], tau)
    )
  }
  point <- calc(dd)
  set.seed(seed)
  boots <- replicate(B, {
    i0 <- sample(which(dd$z == 0), sum(dd$z == 0), replace = TRUE)
    i1 <- sample(which(dd$z == 1), sum(dd$z == 1), replace = TRUE)
    tryCatch(calc(dd[c(i0, i1), ]), error = function(e) c(later = NA_real_, early = NA_real_))
  })
  dif <- boots["early", ] - boots["later", ]
  tibble(
    cohort = "Diagnosed 2004-2012",
    N = nrow(dd),
    events = sum(dd$event),
    at_risk_120 = sum(dd$time >= 120),
    tau_months = tau,
    RMST_age_50_plus = point["later"],
    RMST_age_20_49 = point["early"],
    adjusted_difference = point["early"] - point["later"],
    CI_low = quantile(dif, 0.025, na.rm = TRUE),
    CI_high = quantile(dif, 0.975, na.rm = TRUE),
    bootstrap_replicates = B
  )
}

d_10y <- filter(d, diagnosis_year <= 2012)
rmst_10y <- weighted_rmst(d_10y)
write_csv(rmst_10y, file.path(out, "RMST_120_month_restricted_diagnosis_years.csv"))

# Strict conventional histology sensitivity (ICD-O-3 8140/3).
strict <- filter(d, histology_code == 8140L) %>% droplevels()
strict_fit <- coxph(
  Surv(time, event) ~ early + site + sex + race + marital + period + stage + grade,
  data = strict,
  ties = "efron",
  robust = TRUE
)
strict_cox <- tidy(strict_fit, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term %in% c("earlyAge 20-49", "siteJejunum/Ileum", "siteOther specified", "siteNOS")) %>%
  transmute(
    definition = "ICD-O-3 8140/3 only",
    N = strict_fit$n,
    events = strict_fit$nevent,
    term,
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p_value = p.value
  )
write_csv(strict_cox, file.path(out, "strict_histology_survival_sensitivity.csv"))

# Count-based long-term trend check for the strict 8140/3 definition. This is
# deliberately supportive: it uses population offsets but is not a replacement
# for the age-adjusted Joinpoint estimates.
broad_codes <- c(8140L, 8143L, 8144L, 8145L, 8210L, 8211L, 8220L, 8221L,
                 8255L, 8260L, 8261L, 8262L, 8263L, 8310L, 8480L, 8481L, 8490L)
audit <- read_tsv(file.path(root, "00_SEER8_histology_code_audit_1975_2023.txt"), show_col_types = FALSE) %>%
  transmute(
    year = as.integer(PERIOD_5Y),
    subsite = case_when(
      SBA_SITE_6 == "Duodenum" ~ "Duodenum",
      SBA_SITE_6 %in% c("Jejunum", "Ileum") ~ "Jejunum/Ileum",
      SBA_SITE_6 %in% c("Meckel diverticulum", "Overlapping lesion") ~ "Other specified",
      SBA_SITE_6 == "Small intestine, NOS" ~ "NOS",
      TRUE ~ NA_character_
    ),
    histology_code = as.integer(`Histologic Type ICD-O-3`)
  ) %>% filter(!is.na(subsite), histology_code %in% broad_codes)
population <- read_tsv(file.path(root, "01_SEER8_annual_age_subsite_rates.txt"), show_col_types = FALSE) %>%
  filter(SBA_SITE_MAIN == "Duodenum") %>%
  distinct(YEAR_SINGLE, AGE_3GROUPS_RATE, Population) %>%
  group_by(year = YEAR_SINGLE) %>%
  summarise(population = sum(Population), .groups = "drop")
trend_counts <- bind_rows(
  audit %>% count(year, subsite, name = "cases") %>% mutate(definition = "Protocol adenocarcinoma codes"),
  audit %>% filter(histology_code == 8140L) %>% count(year, subsite, name = "cases") %>% mutate(definition = "ICD-O-3 8140/3 only")
) %>%
  complete(year = 1975:2023, subsite = c("Duodenum", "Jejunum/Ileum", "Other specified", "NOS"), definition, fill = list(cases = 0)) %>%
  left_join(population, by = "year")
strict_trend <- trend_counts %>% group_by(definition, subsite) %>% group_modify(function(z, key) {
  fit <- glm(cases ~ year + offset(log(population)), data = z, family = quasipoisson())
  co <- summary(fit)$coef["year", ]
  tibble(
    years = "1975-2023",
    total_cases = sum(z$cases),
    APC = 100 * (exp(co[1]) - 1),
    CI_low = 100 * (exp(co[1] - 1.96 * co[2]) - 1),
    CI_high = 100 * (exp(co[1] + 1.96 * co[2]) - 1),
    p_value = co[4]
  )
}) %>% ungroup()
write_csv(strict_trend, file.path(out, "strict_histology_count_based_trend_sensitivity.csv"))

cat("Surveillance revision outputs:", normalizePath(out, winslash = "/"), "\n")
