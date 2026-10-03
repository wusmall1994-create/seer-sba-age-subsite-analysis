#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival); library(broom)
})

args <- commandArgs(trailingOnly=TRUE)
clean_dir <- if (length(args)>=1) args[[1]] else Sys.getenv("SBA_CLEAN_DIR")
if (!nzchar(clean_dir)) stop("Supply the cleaned-data directory as the first argument or SBA_CLEAN_DIR.")
out_dir <- if (length(args)>=2) args[[2]] else file.path(dirname(clean_dir), "results", "survival_robust")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dat <- readRDS(file.path(clean_dir,"sba_master_clean.rds"))

marital2 <- function(x) case_when(
  x=="Married (including common law)" ~ "Married",
  is.na(x) | x=="Unknown" ~ "Unknown", TRUE ~ "Not married")
period4 <- function(y) factor(case_when(
  y<=2008~"2004-2008", y<=2013~"2009-2013",
  y<=2018~"2014-2018", TRUE~"2019-2022"),
  levels=c("2004-2008","2009-2013","2014-2018","2019-2022"))

d <- dat %>% filter(cohort_survival_primary) %>% transmute(
  record_id, time=analysis_time_months, os_event,
  early_onset=factor(if_else(early_onset,"20-49","50+"),levels=c("50+","20-49")),
  site_main=factor(site_main,levels=c("Duodenum","Jejunum/Ileum","Other specified","NOS")),
  sex=relevel(factor(sex),ref="Female"),
  # Pool the two extremely sparse categories for numerical stability. This
  # affects adjustment only; descriptive race/ethnicity estimates retain the
  # original five-category variable in the incidence analyses.
  race_ethnicity=relevel(factor(case_when(
    race_ethnicity %in% c("Non-Hispanic American Indian/Alaska Native",
                          "Non-Hispanic Unknown Race") ~ "Other/unknown",
    TRUE ~ race_ethnicity)),ref="Non-Hispanic White"),
  marital=factor(marital2(`Marital status at diagnosis`),levels=c("Married","Not married","Unknown")),
  diagnosis_period=period4(diagnosis_year),
  histology_group=relevel(factor(histology_group),ref="Other adenocarcinoma"),
  stage4=relevel(factor(stage4),ref="Localized"),
  grade4=relevel(factor(grade4),ref="Grade I"),
  surgery=relevel(factor(surgery),ref="No"),
  chemotherapy=relevel(factor(chemotherapy),ref="No/unknown"),
  radiation=relevel(factor(radiation),ref="None/unknown")) %>%
  filter(complete.cases(.))

stopifnot(nrow(d)==11563L, all(d$time>0))

# Conditional interval-specific Cox models. At each landmark, only patients
# alive at the interval start enter the risk set; follow-up is administratively
# censored at the interval end. This yields directly interpretable conditional
# HRs without imposing one coefficient over the full follow-up period.
base_formula <- Surv(interval_time,interval_event) ~
  early_onset + site_main + sex + race_ethnicity + marital + diagnosis_period +
  histology_group + stage4 + grade4 + surgery + chemotherapy + radiation
intervals <- tibble(start=c(0,12,36,60),end=c(12,36,60,Inf),
                    label=c("0-12","13-36","37-60",">60"))
fit_interval <- function(start,end,label) {
  di <- d %>% filter(time>start) %>% mutate(
    interval_time=pmin(time,end)-start,
    interval_event=as.integer(os_event==1L & time<=end))
  fit <- coxph(base_formula,data=di,ties="efron",robust=TRUE)
  tidy(fit,conf.int=TRUE,exponentiate=TRUE) %>%
    filter(grepl("^(early_onset|site_main)",term)) %>%
    mutate(exposure=case_when(
      term=="early_onset20-49"~"Age 20-49 vs 50+",
      term=="site_mainJejunum/Ileum"~"Jejunum/Ileum vs duodenum",
      term=="site_mainOther specified"~"Other specified vs duodenum",
      term=="site_mainNOS"~"NOS vs duodenum",
      TRUE~term),followup_interval=label,
      risk_set_n=nrow(di),events=sum(di$interval_event)) %>%
    transmute(exposure,followup_interval,risk_set_n,events,
              HR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value)
}
pw <- bind_rows(Map(fit_interval,intervals$start,intervals$end,intervals$label))
write_csv(pw,file.path(out_dir,"piecewise_cox_OS.csv"))

# Restricted mean survival time calculated as the area under the KM curve.
rmst_one <- function(time,event,tau=60) {
  fit <- survfit(Surv(time,event)~1)
  tt <- c(0, fit$time[fit$time<tau], tau)
  ss <- c(1, fit$surv[fit$time<tau])
  sum(diff(tt)*ss)
}

rmst_boot <- function(data,var,tau=60,B=500,seed=20260922) {
  g <- droplevels(data[[var]]); lev <- levels(g)
  point <- vapply(lev,function(z) rmst_one(data$time[g==z],data$os_event[g==z],tau),numeric(1))
  set.seed(seed)
  boots <- replicate(B, {
    unlist(lapply(lev,function(z) {
      ii <- which(g==z); jj <- sample(ii,length(ii),replace=TRUE)
      rmst_one(data$time[jj],data$os_event[jj],tau)
    }))
  })
  if (is.null(dim(boots))) boots <- matrix(boots,nrow=length(lev))
  rows <- lapply(seq_along(lev),function(i) tibble(
    variable=var,level=lev[i],reference=lev[1],tau_months=tau,
    RMST_months=point[i],RMST_difference=point[i]-point[1],
    CI_low=if(i==1) NA_real_ else quantile(boots[i,]-boots[1,],.025,na.rm=TRUE),
    CI_high=if(i==1) NA_real_ else quantile(boots[i,]-boots[1,],.975,na.rm=TRUE)))
  bind_rows(rows)
}

rmst <- bind_rows(
  rmst_boot(d,"early_onset",tau=60,B=500,seed=20260922),
  rmst_boot(d,"site_main",tau=60,B=500,seed=20260923),
  rmst_boot(d,"early_onset",tau=120,B=500,seed=20260924),
  rmst_boot(d,"site_main",tau=120,B=500,seed=20260925))
write_csv(rmst,file.path(out_dir,"RMST_OS_unadjusted.csv"))

# Complete-case and event-count audit for every piecewise interval.
qc <- intervals %>% rowwise() %>% transmute(
  interval=label,records=sum(d$time>start),
  events=sum(d$os_event==1L & d$time>start & d$time<=end)) %>% ungroup()
write_csv(qc,file.path(out_dir,"piecewise_cox_QC.csv"))
cat("Robust non-PH survival analysis complete:",normalizePath(out_dir,winslash="/"),"\n")
print(pw); print(rmst); print(qc)
