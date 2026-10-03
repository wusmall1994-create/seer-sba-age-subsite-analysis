#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(survival); library(cmprsk); library(ggplot2); library(broom)
})

args <- commandArgs(trailingOnly=TRUE)
clean_dir <- if (length(args)>=1) args[[1]] else Sys.getenv("SBA_CLEAN_DIR")
if (!nzchar(clean_dir)) stop("Supply the cleaned-data directory as the first argument or SBA_CLEAN_DIR.")
out_dir <- if (length(args)>=2) args[[2]] else file.path(dirname(clean_dir), "results", "survival")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dat <- readRDS(file.path(clean_dir,"sba_master_clean.rds"))

fmt_p <- function(x) ifelse(is.na(x),NA_character_,ifelse(x<0.001,"<0.001",sprintf("%.3f",x)))
marital2 <- function(x) case_when(x=="Married (including common law)" ~ "Married",
                                  is.na(x) | x=="Unknown" ~ "Unknown", TRUE ~ "Not married")
period4 <- function(y) factor(case_when(y<=2008~"2004-2008",y<=2013~"2009-2013",
                                        y<=2018~"2014-2018",TRUE~"2019-2022"),
                              levels=c("2004-2008","2009-2013","2014-2018","2019-2022"))

make_analysis <- function(flag) dat %>% filter(.data[[flag]]) %>% transmute(
  record_id, patient_id, time=analysis_time_months, os_event,
  status_css=case_when(css_event==1L~1L,os_event==1L~2L,TRUE~0L),
  early_onset=factor(if_else(early_onset,"20-49","50+"),levels=c("50+","20-49")),
  site_main=factor(site_main,levels=c("Duodenum","Jejunum/Ileum","Other specified","NOS")),
  sex=relevel(factor(sex),ref="Female"), race_ethnicity=relevel(factor(race_ethnicity),ref="Non-Hispanic White"),
  marital=factor(marital2(`Marital status at diagnosis`),levels=c("Married","Not married","Unknown")),
  diagnosis_period=period4(diagnosis_year), histology_group=relevel(factor(histology_group),ref="Other adenocarcinoma"),
  stage4=relevel(factor(stage4),ref="Localized"), grade4=relevel(factor(grade4),ref="Grade I"),
  surgery=relevel(factor(surgery),ref="No"), chemotherapy=relevel(factor(chemotherapy),ref="No/unknown"),
  radiation=relevel(factor(radiation),ref="None/unknown"))

survdat <- make_analysis("cohort_survival_primary")
stopifnot(nrow(survdat)==11563L, all(survdat$time>0), !anyNA(survdat$os_event), !anyNA(survdat$status_css))

# Kaplan-Meier summary and log-rank tests ------------------------------------
km_summary <- function(data,var) {
  f <- as.formula(paste0("Surv(time,os_event)~",var)); fit <- survfit(f,data=data)
  s <- summary(fit,times=c(12,36,60),extend=TRUE)
  tibble(stratum=sub(paste0("^",var,"="),"",s$strata),months=s$time,
         survival=s$surv,lower=s$lower,upper=s$upper,n_risk=s$n.risk,n_event=s$n.event)
}
km_early <- km_summary(survdat,"early_onset") %>% mutate(comparison="Age at diagnosis")
km_site <- km_summary(survdat,"site_main") %>% mutate(comparison="Primary subsite")
write_csv(bind_rows(km_early,km_site),file.path(out_dir,"KM_landmark_OS.csv"))

logrank <- function(var) {
  z <- survdiff(as.formula(paste0("Surv(time,os_event)~",var)),data=survdat)
  tibble(comparison=var,chisq=unname(z$chisq),df=length(z$n)-1,p_value=pchisq(z$chisq,length(z$n)-1,lower.tail=FALSE))
}
write_csv(bind_rows(logrank("early_onset"),logrank("site_main")),file.path(out_dir,"logrank_tests.csv"))

# Cumulative incidence of SBA death, treating other death as competing -------
cif_summary <- function(data,var) {
  g <- data[[var]]; ci <- cuminc(data$time,data$status_css,group=g,cencode=0)
  nm <- names(ci)[!names(ci)%in%"Tests"]
  rows <- lapply(nm,function(k) {
    event <- sub(".* ","",k); if (event!="1") return(NULL)
    stratum <- sub(" 1$","",k); obj <- ci[[k]]
    bind_rows(lapply(c(12,36,60),function(tt) {
      idx <- max(which(obj$time<=tt)); tibble(stratum=stratum,months=tt,cif=obj$est[idx],variance=obj$var[idx])
    }))
  })
  out <- bind_rows(rows) %>% mutate(lower=pmax(0,cif-1.96*sqrt(variance)),upper=pmin(1,cif+1.96*sqrt(variance)))
  test <- if (!is.null(ci$Tests)) tibble(comparison=var,statistic=ci$Tests[1,"stat"],df=ci$Tests[1,"df"],p_value=ci$Tests[1,"pv"]) else tibble()
  list(summary=out,test=test,object=ci)
}
cif_e <- cif_summary(survdat,"early_onset"); cif_s <- cif_summary(survdat,"site_main")
write_csv(bind_rows(cif_e$summary%>%mutate(comparison="Age at diagnosis"),cif_s$summary%>%mutate(comparison="Primary subsite")),file.path(out_dir,"CIF_landmark_CSS.csv"))
write_csv(bind_rows(cif_e$test,cif_s$test),file.path(out_dir,"gray_tests.csv"))

# Multivariable models --------------------------------------------------------
clinical_formula <- Surv(time,os_event) ~ early_onset + site_main + sex + race_ethnicity + marital + diagnosis_period + histology_group + stage4 + grade4
treatment_formula <- update(clinical_formula,. ~ . + surgery + chemotherapy + radiation)
fit_cox <- coxph(clinical_formula,data=survdat,x=TRUE,model=TRUE)
fit_cox_tx <- coxph(treatment_formula,data=survdat,x=TRUE,model=TRUE)
tidy_cox <- function(fit,name) tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% transmute(model=name,term,HR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value)
write_csv(bind_rows(tidy_cox(fit_cox,"Clinical"),tidy_cox(fit_cox_tx,"Clinical + treatment")),file.path(out_dir,"cox_OS_models.csv"))

ph <- cox.zph(fit_cox_tx,transform="km")$table
write_csv(tibble(term=rownames(ph),chisq=ph[,"chisq"],df=ph[,"df"],p_value=ph[,"p"]),file.path(out_dir,"cox_PH_assumption.csv"))

fg_data <- survdat %>% drop_na(early_onset,site_main,sex,race_ethnicity,marital,diagnosis_period,histology_group,stage4,grade4,surgery,chemotherapy,radiation)
Xclin <- model.matrix(~ early_onset + site_main + sex + race_ethnicity + marital + diagnosis_period + histology_group + stage4 + grade4,fg_data)[,-1,drop=FALSE]
Xtx <- model.matrix(~ early_onset + site_main + sex + race_ethnicity + marital + diagnosis_period + histology_group + stage4 + grade4 + surgery + chemotherapy + radiation,fg_data)[,-1,drop=FALSE]
fg_tidy <- function(X,name) {
  f <- crr(fg_data$time,fg_data$status_css,cov1=X,failcode=1,cencode=0)
  se <- sqrt(diag(f$var)); tibble(model=name,term=colnames(X),sHR=exp(f$coef),CI_low=exp(f$coef-1.96*se),CI_high=exp(f$coef+1.96*se),p_value=2*pnorm(-abs(f$coef/se)))
}
write_csv(bind_rows(fg_tidy(Xclin,"Clinical"),fg_tidy(Xtx,"Clinical + treatment")),file.path(out_dir,"finegray_CSS_models.csv"))

# Sensitivity analyses with the same clinical specification ------------------
sensitivity <- function(flag,label) {
  d <- make_analysis(flag)
  fit <- coxph(clinical_formula,data=d)
  tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(term=="early_onset20-49") %>%
    transmute(cohort=label,n=nrow(d),events=sum(d$os_event),HR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value)
}
sens <- bind_rows(sensitivity("cohort_survival_primary","Primary survival cohort"),
                  sensitivity("cohort_survival_first_primary","First-primary only"),
                  sensitivity("cohort_survival_micro","Microscopically confirmed"))
write_csv(sens,file.path(out_dir,"sensitivity_early_onset_OS.csv"))

# Publication-ready KM and CIF curves ---------------------------------------
km_plot_data <- function(var) {
  fit <- survfit(as.formula(paste0("Surv(time,os_event)~",var)),data=survdat); s<-summary(fit)
  tibble(time=s$time,survival=s$surv,stratum=sub(paste0("^",var,"="),"",s$strata))
}
p1 <- ggplot(km_plot_data("early_onset"),aes(time,survival,color=stratum))+geom_step(linewidth=.8)+
  coord_cartesian(xlim=c(0,120),ylim=c(0,1))+scale_y_continuous(labels=scales::percent)+
  labs(x="Months since diagnosis",y="Overall survival",color="Age group")+theme_classic(base_size=11)
ggsave(file.path(out_dir,"Fig_KM_OS_age.png"),p1,width=6.5,height=4.8,dpi=300)
p2 <- ggplot(km_plot_data("site_main"),aes(time,survival,color=stratum))+geom_step(linewidth=.8)+
  coord_cartesian(xlim=c(0,120),ylim=c(0,1))+scale_y_continuous(labels=scales::percent)+
  labs(x="Months since diagnosis",y="Overall survival",color="Primary subsite")+theme_classic(base_size=11)
ggsave(file.path(out_dir,"Fig_KM_OS_subsite.png"),p2,width=6.5,height=4.8,dpi=300)

qc <- tibble(check=c("Primary survival cohort","Positive follow-up time","OS event observed","CSS status defined","Fine-Gray complete cases"),
             value=c(nrow(survdat),sum(survdat$time>0),sum(!is.na(survdat$os_event)),sum(!is.na(survdat$status_css)),nrow(fg_data)),
             expected=c(11563,11563,11563,11563,nrow(fg_data))) %>% mutate(status=if_else(value==expected,"PASS","FAIL"))
write_csv(qc,file.path(out_dir,"survival_QC.csv"))
if(any(qc$status!="PASS")) stop("Survival QC failed")
cat("Survival analysis complete:",normalizePath(out_dir,winslash="/"),"\n"); print(qc); print(sens)
