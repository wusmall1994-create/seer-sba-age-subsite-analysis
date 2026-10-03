#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(stringr)})

args <- commandArgs(trailingOnly=TRUE)
clean_dir <- if(length(args)>=1) args[[1]] else Sys.getenv("SBA_CLEAN_DIR")
if (!nzchar(clean_dir)) stop("Supply the cleaned-data directory as the first argument or SBA_CLEAN_DIR.")
out_dir <- if(length(args)>=2) args[[2]] else file.path(dirname(clean_dir), "results", "table1")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(file.path(clean_dir,"sba_master_clean.rds")) %>% filter(cohort_descriptive)

fmt_n_pct <- function(n,den) sprintf("%s (%.1f)",format(n,big.mark=","),100*n/den)
cramers_v <- function(tab){
  if(any(dim(tab)<2)||sum(tab)==0)return(NA_real_)
  ch<-suppressWarnings(chisq.test(tab,correct=FALSE)); sqrt(unname(ch$statistic)/(sum(tab)*min(nrow(tab)-1,ncol(tab)-1)))
}
p_cat <- function(tab){
  if(any(dim(tab)<2)) return(NA_real_)
  suppressWarnings(chisq.test(tab,correct=FALSE)$p.value)
}

make_cat <- function(data,group,var,label,levels=NULL,missing_label="Missing"){
  gv <- data[[group]]
  glev <- if(is.factor(gv)) levels(gv) else unique(as.character(gv))
  g <- factor(as.character(gv),levels=glev); x <- as.character(data[[var]])
  x[is.na(x)|x==""] <- missing_label
  if(!is.null(levels)) {
    xlev <- unique(c(levels,if(any(x==missing_label)) missing_label else NULL))
    x <- factor(x,levels=xlev)
  }
  tab <- table(x,g,useNA="no"); groups<-colnames(tab); den<-colSums(tab)
  testtab <- tab[rowSums(tab)>0,colSums(tab)>0,drop=FALSE]
  p<-p_cat(testtab); smd<-cramers_v(testtab)
  bind_rows(lapply(seq_len(nrow(tab)),function(i){
    tibble(Variable=ifelse(i==1,label,""),Level=rownames(tab)[i],
           Overall=fmt_n_pct(sum(tab[i,]),sum(tab)),
           !!groups[1]:=fmt_n_pct(tab[i,1],den[1]),
           !!groups[2]:=fmt_n_pct(tab[i,2],den[2]),
           P_value=ifelse(i==1,p,NA_real_),SMD=ifelse(i==1,smd,NA_real_))
  }))
}
make_cont <- function(data,group,var,label){
  groups<-if(is.factor(data[[group]])) levels(data[[group]]) else unique(as.character(data[[group]])); x<-data[[var]]
  f<-function(z)sprintf("%.0f [%.0f-%.0f]",median(z,na.rm=TRUE),quantile(z,.25,na.rm=TRUE),quantile(z,.75,na.rm=TRUE))
  p<-wilcox.test(data[[var]]~data[[group]])$p.value
  s<-abs(mean(x[data[[group]]==groups[1]],na.rm=T)-mean(x[data[[group]]==groups[2]],na.rm=T))/sqrt((var(x[data[[group]]==groups[1]],na.rm=T)+var(x[data[[group]]==groups[2]],na.rm=T))/2)
  tibble(Variable=label,Level="Median [IQR]",Overall=f(x),!!groups[1]:=f(x[data[[group]]==groups[1]]),!!groups[2]:=f(x[data[[group]]==groups[2]]),P_value=p,SMD=s)
}

period <- function(y) cut(y,c(1999,2005,2011,2017,2023),labels=c("2000-2005","2006-2011","2012-2017","2018-2023"))
d <- d %>% mutate(
  onset_group=factor(if_else(early_onset,"Early-onset (20–49)","Age 50+"),levels=c("Early-onset (20–49)","Age 50+")),
  period4=period(diagnosis_year),
  marital=case_when(str_detect(`Marital status at diagnosis`,"Married")~"Married",str_detect(`Marital status at diagnosis`,"Unknown")~"Unknown",TRUE~"Not married"),
  stage_display=as.character(stage4),grade_display=as.character(grade4),
  histology_display=as.character(histology_group),surgery_display=as.character(surgery),
  chemo_display=as.character(chemotherapy),radiation_display=as.character(radiation))

vars <- list(
  c("sex","Sex"),c("race_ethnicity","Race/ethnicity"),c("marital","Marital status"),
  c("period4","Diagnosis period"),c("site_main","Primary subsite"),c("histology_display","Histologic subtype"),
  c("stage_display","Summary stage (2004+)") ,c("grade_display","Grade"),
  c("microscopically_confirmed","Microscopic confirmation"),c("surgery_display","Surgery"),
  c("chemo_display","Chemotherapy"),c("radiation_display","Radiation"))
levs <- list(
  sex=c("Female","Male"),race_ethnicity=levels(d$race_ethnicity),marital=c("Married","Not married","Unknown"),
  period4=levels(d$period4),site_main=levels(d$site_main),histology_display=levels(d$histology_group),
  stage_display=levels(d$stage4),grade_display=levels(d$grade4),
  microscopically_confirmed=c("FALSE","TRUE"),surgery_display=levels(d$surgery),
  chemo_display=levels(d$chemotherapy),radiation_display=levels(d$radiation))

eo <- bind_rows(make_cont(d,"onset_group","diagnosis_year","Year of diagnosis, median [IQR]"),
                lapply(vars,\(z)make_cat(d,"onset_group",z[1],z[2],levs[[z[1]]])))
names(eo)[4:5] <- c("Early_onset_20_49","Age_50_plus")

site_d <- d %>% filter(site_main %in% c("Duodenum","Jejunum/Ileum")) %>%
  mutate(site_compare=factor(as.character(site_main),levels=c("Duodenum","Jejunum/Ileum")))
site_vars <- vars[!vapply(vars,\(x)x[1]=="site_main",logical(1))]
site_tab <- bind_rows(make_cont(site_d,"site_compare","age","Age, years, median [IQR]"),
                      lapply(site_vars,\(z)make_cat(site_d,"site_compare",z[1],z[2],levs[[z[1]]])))

missing <- tibble(variable=c("Stage (2000-2003 not applicable)","Grade","Tumor size (2016+)","Nodes examined: raw blank only","Nodes positive: raw blank only","Surgery","Survival months"),
                  missing_n=c(sum(is.na(d$stage4)),sum(d$grade4=="Unknown"),sum(is.na(d$tumor_size_raw)),sum(is.na(d$nodes_examined_raw)),sum(is.na(d$nodes_positive_raw)),sum(d$surgery=="Unknown"),sum(!d$survival_known))) %>%
  mutate(total=nrow(d),missing_percent=100*missing_n/total)

write_csv(eo,file.path(out_dir,"table1A_early_vs_older.csv"),na="")
write_csv(site_tab,file.path(out_dir,"table1B_duodenum_vs_jejunum_ileum.csv"),na="")
write_csv(missing,file.path(out_dir,"missingness_audit.csv"),na="")

qc <- tibble(check=c("Descriptive cohort","Early-onset","Age 50+","Duodenum vs jejunum/ileum subset"),
             value=c(nrow(d),sum(d$early_onset),sum(!d$early_onset),nrow(site_d)),
             expected=c(14158,1436,12722,12077)) %>% mutate(status=if_else(value==expected,"PASS","FAIL"))
write_csv(qc,file.path(out_dir,"table1_QC.csv"))
stopifnot(all(qc$status=="PASS"))
cat("Table 1 outputs written to",out_dir,"\n"); print(qc)
