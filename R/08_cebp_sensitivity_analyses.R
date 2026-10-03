#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival); library(splines); library(broom); library(cmprsk)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_PROJECT_ROOT")
if (!nzchar(project_root)) stop("Supply the project root as the first argument or SBA_PROJECT_ROOT.")
clean_dir <- file.path(project_root, "cleaned")
root <- file.path(project_root, "results")
out <- file.path(root,"cebp_sensitivity")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dat <- readRDS(file.path(clean_dir,"sba_master_clean.rds"))

marital2 <- function(x) case_when(
  x=="Married (including common law)" ~ "Married",
  is.na(x) | x=="Unknown" ~ "Unknown", TRUE ~ "Not married")
period4 <- function(y) factor(case_when(
  y<=2008~"2004-2008", y<=2013~"2009-2013", y<=2018~"2014-2018", TRUE~"2019-2022"),
  levels=c("2004-2008","2009-2013","2014-2018","2019-2022"))

d <- dat %>% filter(cohort_survival_primary) %>% transmute(
  record_id, diagnosis_year, age, time=analysis_time_months, event=os_event,
  early=factor(if_else(early_onset,"Age 20-49","Age 50+"),levels=c("Age 50+","Age 20-49")),
  site=factor(site_main,levels=c("Duodenum","Jejunum/Ileum","Other specified","NOS")),
  sex=relevel(factor(sex),ref="Female"),
  race=relevel(factor(case_when(
    race_ethnicity %in% c("Non-Hispanic American Indian/Alaska Native","Non-Hispanic Unknown Race") ~ "Other/unknown",
    TRUE ~ race_ethnicity)),ref="Non-Hispanic White"),
  marital=factor(marital2(`Marital status at diagnosis`),levels=c("Married","Not married","Unknown")),
  period=period4(diagnosis_year), histology=relevel(factor(histology_group),ref="Other adenocarcinoma"),
  stage=factor(stage4,levels=c("Localized","Regional","Distant","Unknown/unstaged")),
  grade=factor(grade4), surgery=relevel(factor(surgery),ref="No"),
  chemotherapy=relevel(factor(chemotherapy),ref="No/unknown"),
  radiation=relevel(factor(radiation),ref="None/unknown")) %>%
  filter(complete.cases(.))
stopifnot(nrow(d)==11563L, any(d$stage=="Unknown/unstaged"))

base_no_tx <- Surv(time,event) ~ early + site + sex + race + marital + period +
  histology + stage + grade
base_with_tx <- update(base_no_tx,. ~ . + surgery + chemotherapy + radiation)

extract_target <- function(fit,analysis) tidy(fit,conf.int=TRUE,exponentiate=TRUE) %>%
  filter(term %in% c("earlyAge 20-49","siteJejunum/Ileum","siteOther specified","siteNOS")) %>%
  transmute(analysis,term,HR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value,N=fit$n,events=fit$nevent)

# Treatment is post-diagnosis and may be a mediator/collider. Compare models
# with and without treatment adjustment.
f_no_tx <- coxph(base_no_tx,data=d,ties="efron",robust=TRUE)
f_with_tx <- coxph(base_with_tx,data=d,ties="efron",robust=TRUE)
tx <- bind_rows(extract_target(f_no_tx,"No treatment variables"),
                extract_target(f_with_tx,"Surgery, chemotherapy, and radiation adjusted"))
write_csv(tx,file.path(out,"cox_treatment_adjustment_sensitivity.csv"))

# Explicit unknown/unstaged cases are retained in the primary model. Two
# transparent scenarios test whether their handling drives the age contrast.
d_known <- filter(d,stage!="Unknown/unstaged") %>% droplevels()
d_worst <- d %>% mutate(stage=factor(if_else(stage=="Unknown/unstaged","Distant",as.character(stage)),
                                    levels=c("Localized","Regional","Distant")))
stage_sens <- bind_rows(
  extract_target(f_with_tx,"Primary: unknown stage as category"),
  extract_target(coxph(base_with_tx,data=d_known,ties="efron",robust=TRUE),"Exclude unknown stage"),
  extract_target(coxph(base_with_tx,data=d_worst,ties="efron",robust=TRUE),"Worst-stage recode"))
write_csv(stage_sens,file.path(out,"stage_unknown_sensitivity.csv"))

# Stage-stratified estimates answer whether the early-onset association is
# confined to a particular disease extent. Treatment is omitted here.
stage_strat <- bind_rows(lapply(c("Localized","Regional","Distant"),function(s){
  ds <- filter(d,stage==s) %>% droplevels()
  fit <- coxph(Surv(time,event)~early+site+sex+race+marital+period+histology+grade,
               data=ds,ties="efron",robust=TRUE)
  tidy(fit,conf.int=TRUE,exponentiate=TRUE) %>% filter(term=="earlyAge 20-49") %>%
    transmute(stage=s,N=fit$n,events=fit$nevent,HR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value)
}))
write_csv(stage_strat,file.path(out,"cox_early_onset_by_stage.csv"))

# Excluding the anatomically nonspecific category.
d_no_nos <- filter(d,site!="NOS") %>% droplevels()
nos_sens <- bind_rows(extract_target(f_no_tx,"Primary cohort"),
                      extract_target(coxph(base_no_tx,data=d_no_nos,ties="efron",robust=TRUE),"Exclude NOS"))
write_csv(nos_sens,file.path(out,"cox_exclude_NOS_sensitivity.csv"))

# Continuous age specification with restricted cubic basis (natural spline),
# including sex interaction. Predictions are standardized at reference
# covariates and expressed relative to age 50 within sex.
f_age_linear <- coxph(Surv(time,event)~age*sex+site+race+marital+period+histology+stage+grade,
                      data=d,ties="efron")
f_age_spline <- coxph(Surv(time,event)~ns(age,df=4)*sex+site+race+marital+period+histology+stage+grade,
                      data=d,ties="efron")
lrt <- anova(f_age_linear,f_age_spline,test="LRT")
write_csv(tibble(linear_df=attr(logLik(f_age_linear),"df"),spline_df=attr(logLik(f_age_spline),"df"),
                 LR_chisq=lrt$Chisq[2],df=lrt$Df[2],p_value=lrt$`Pr(>|Chi|)`[2]),
          file.path(out,"continuous_age_spline_LRT.csv"))
refrow <- d[1,]
refrow$site <- factor("Duodenum",levels=levels(d$site)); refrow$race <- factor("Non-Hispanic White",levels=levels(d$race))
refrow$marital <- factor("Married",levels=levels(d$marital)); refrow$period <- factor("2014-2018",levels=levels(d$period))
refrow$histology <- factor("Other adenocarcinoma",levels=levels(d$histology)); refrow$stage <- factor("Localized",levels=levels(d$stage))
refrow$grade <- factor("Grade I",levels=levels(d$grade))
pred <- bind_rows(lapply(levels(d$sex),function(sx){
  nd <- refrow[rep(1,length(20:90)),]; nd$age <- 20:90; nd$sex <- factor(sx,levels=levels(d$sex))
  X <- model.matrix(delete.response(terms(f_age_spline)),nd)
  b <- coef(f_age_spline); V <- vcov(f_age_spline); eta <- drop(X[,names(b),drop=FALSE]%*%b)
  i50 <- which(nd$age==50); C <- X[,names(b),drop=FALSE]-matrix(X[i50,names(b)],nrow=nrow(X),ncol=length(b),byrow=TRUE)
  se <- sqrt(rowSums((C%*%V)*C)); est <- eta-eta[i50]
  tibble(sex=sx,age=nd$age,HR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se))
}))
write_csv(pred,file.path(out,"continuous_age_spline_predictions.csv"))

# Official Joinpoint AAPCs and BH-FDR across the 12 prespecified age-by-site
# full-range tests. Recent 5- and 10-year estimates remain descriptive.
aapc_file <- file.path(root,"joinpoint_input","JP01_SEER8_age_subsite_1975_2023.Export.AAPC.txt")
aapc <- read_tsv(aapc_file,show_col_types=FALSE)
full <- aapc %>% filter(tolower(as.character(`AAPC Index`)) %in% c("full range","entire range") |
                        (`Start Obs`==1975 & `End Obs`==2023))
if(nrow(full)==0) full <- aapc %>% group_by(Subsite,AgeGroup) %>% slice(1) %>% ungroup()
full <- full %>% mutate(p_FDR=p.adjust(`P-Value`,method="BH"),significant_FDR=p_FDR<.05)
write_csv(full,file.path(out,"age_subsite_AAPC_with_BH_FDR.csv"))

# Report all prespecified AAPC windows exactly as exported by Joinpoint.
aapc_all <- aapc %>% filter(`AAPC Index` %in% c("Full Range","Last 10","Last 5")) %>%
  transmute(Subsite,AgeGroup,window=`AAPC Index`,start_year=`Start Obs`,end_year=`End Obs`,
            AAPC,AAPC_CI_low=`AAPC C.I. Low`,AAPC_CI_high=`AAPC C.I. High`,p_value=`P-Value`)
write_csv(aapc_all,file.path(out,"age_subsite_AAPC_prespecified_windows.csv"))

# Formal weighted log-linear tests of whether temporal slopes differ across the
# 12 prespecified age-by-subsite strata. These complement rather than replace
# the official Joinpoint estimates.
age_rate <- read_tsv(file.path(root,"joinpoint_input","JP01_SEER8_age_subsite_1975_2023.txt"),show_col_types=FALSE) %>%
  filter(is.finite(Rate),Rate>0,is.finite(StandardError),StandardError>0) %>%
  mutate(stratum=interaction(Subsite,AgeGroup,drop=TRUE),w=1/(StandardError/Rate)^2)
fit_common <- lm(log(Rate)~stratum+Year,data=age_rate,weights=w)
fit_hetero <- lm(log(Rate)~stratum+stratum:Year,data=age_rate,weights=w)
global_an <- anova(fit_common,fit_hetero)
hetero <- tibble(test="Age-by-subsite slope heterogeneity",
                 F_statistic=global_an$F[2],df_num=global_an$Df[2],df_den=df.residual(fit_hetero),
                 p_value=global_an$`Pr(>F)`[2])
by_site <- bind_rows(lapply(unique(age_rate$Subsite),function(s){
  z <- filter(age_rate,Subsite==s)
  f0 <- lm(log(Rate)~AgeGroup+Year,data=z,weights=w)
  f1 <- lm(log(Rate)~AgeGroup+AgeGroup:Year,data=z,weights=w)
  a <- anova(f0,f1)
  tibble(test=paste0("Age-group slope heterogeneity within ",s),F_statistic=a$F[2],
         df_num=a$Df[2],df_den=df.residual(f1),p_value=a$`Pr(>F)`[2])
}))
write_csv(bind_rows(hetero,by_site),file.path(out,"age_subsite_trend_heterogeneity_tests.csv"))

# Absolute rate changes compare the mean of the first and last five observed
# annual rates in each stratum; estimates are descriptive and expressed per
# million person-years for easier interpretation of this rare cancer.
abs_change <- age_rate %>% arrange(Subsite,AgeGroup,Year) %>% group_by(Subsite,AgeGroup) %>%
  summarise(first_year=min(Year),last_year=max(Year),
            early_5y_rate_per_million=10*mean(head(Rate,5)),
            late_5y_rate_per_million=10*mean(tail(Rate,5)),
            absolute_change_per_million=late_5y_rate_per_million-early_5y_rate_per_million,
            .groups="drop")
write_csv(abs_change,file.path(out,"age_subsite_absolute_rate_change.csv"))

# COVID-era sensitivity: weighted log-linear trend excluding 2020-2021. This
# is not a replacement Joinpoint model; it checks whether terminal disruption
# materially changes the overall log-linear slope.
jpdat <- read_tsv(file.path(root,"joinpoint_input","JP00_SEER8_overall_subsite_1975_2023.txt"),show_col_types=FALSE)
covid_fit <- function(x,exclude=FALSE){
  z <- x %>% filter(is.finite(Rate),Rate>0,is.finite(StandardError),StandardError>0)
  if(exclude) z <- filter(z,!Year %in% 2020:2021)
  fit <- lm(log(Rate)~Year,data=z,weights=1/(StandardError/Rate)^2)
  cf <- summary(fit)$coef["Year",]; apc <- 100*(exp(cf[1])-1); lo <- 100*(exp(cf[1]-1.96*cf[2])-1); hi <- 100*(exp(cf[1]+1.96*cf[2])-1)
  tibble(model=if_else(exclude,"Exclude 2020-2021","All years"),APC=apc,CI_low=lo,CI_high=hi,p_value=cf[4],N=nrow(z))
}
covid <- jpdat %>% group_by(Subsite) %>% group_modify(~bind_rows(covid_fit(.x,FALSE),covid_fit(.x,TRUE))) %>% ungroup()
write_csv(covid,file.path(out,"covid_year_exclusion_loglinear_sensitivity.csv"))

# Five-year absolute cancer-specific mortality in the presence of other-cause
# death. Cause 1 is SEER-attributed SBA death; cause 2 is other-cause death.
cs <- dat %>% filter(cohort_survival_primary) %>% transmute(
  time=analysis_time_months,
  status=case_when(css_event==1L~1L,os_event==1L~2L,TRUE~0L),
  early=if_else(early_onset,"Age 20-49","Age 50+"),site=as.character(site_main))
cif_at <- function(time,status,group,horizon=60){
  ci <- cuminc(time,status,group=group,cencode=0)
  lev <- unique(group)
  bind_rows(lapply(lev,function(g){
    obj <- ci[[paste(g,1)]]; j <- max(which(obj$time<=horizon))
    est <- obj$est[j]; se <- sqrt(obj$var[j])
    tibble(group=g,months=horizon,risk=est,CI_low=max(0,est-1.96*se),CI_high=min(1,est+1.96*se),variance=obj$var[j])
  }))
}
age_cif <- cif_at(cs$time,cs$status,cs$early) %>% mutate(contrast="Age at diagnosis")
site_cs <- filter(cs,site %in% c("Duodenum","Jejunum/Ileum"))
site_cif <- cif_at(site_cs$time,site_cs$status,site_cs$site) %>% mutate(contrast="Primary subsite")
css_abs <- bind_rows(age_cif,site_cif)
write_csv(css_abs,file.path(out,"five_year_absolute_CSS_risk.csv"))
make_diff <- function(tab,exposed,reference,label){
  a <- filter(tab,group==exposed); b <- filter(tab,group==reference)
  dif <- a$risk-b$risk; se <- sqrt(a$variance+b$variance)
  tibble(contrast=label,exposed,reference,risk_difference=dif,
         CI_low=dif-1.96*se,CI_high=dif+1.96*se)
}
css_diff <- bind_rows(
  make_diff(age_cif,"Age 20-49","Age 50+","Age at diagnosis"),
  make_diff(site_cif,"Jejunum/Ileum","Duodenum","Primary subsite"))
write_csv(css_diff,file.path(out,"five_year_absolute_CSS_risk_differences.csv"))

# Audit that primary stage handling is reported correctly.
write_csv(d %>% count(stage,name="N") %>% mutate(percent=100*N/sum(N)),file.path(out,"primary_cohort_stage_distribution.csv"))
cat("CEBP sensitivity analyses complete:",normalizePath(out,winslash="/"),"\n")
