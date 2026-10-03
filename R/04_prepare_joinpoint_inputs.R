#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE,scipen=999)
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(stringr)})

args<-commandArgs(trailingOnly=TRUE)
input_dir<-if(length(args)>=1) args[[1]] else Sys.getenv("SBA_INPUT_DIR")
if(!nzchar(input_dir)) stop("Supply the authorized SEER export directory as the first argument or SBA_INPUT_DIR.")
out_dir<-if(length(args)>=2) args[[2]] else file.path(input_dir,"results","joinpoint_input")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

clean_names<-function(d){names(d)<-names(d)%>%str_replace_all('"','')%>%str_trim();d}
valid_rate<-function(d) d%>%filter(!is.na(Rate),!is.na(StandardError),Rate>0,StandardError>0)

age<-read_tsv(file.path(input_dir,"01_SEER8_annual_age_subsite_rates.txt"),show_col_types=FALSE,na=c("NA",""))%>%clean_names()%>%
  transmute(Subsite=SBA_SITE_MAIN,AgeGroup=AGE_3GROUPS_RATE,Year=as.integer(YEAR_SINGLE),
            Rate=as.numeric(`Age-Adjusted Rate`),StandardError=as.numeric(`Standard Error`),Count=as.integer(Count),Population=as.numeric(Population))%>%
  arrange(Subsite,AgeGroup,Year)
sex<-read_tsv(file.path(input_dir,"02_SEER8_annual_sex_subsite_rates.txt"),show_col_types=FALSE,na=c("NA",""))%>%clean_names()%>%
  transmute(Subsite=SBA_SITE_MAIN,Sex,Year=as.integer(YEAR_SINGLE),Rate=as.numeric(`Age-Adjusted Rate`),
            StandardError=as.numeric(`Standard Error`),Count=as.integer(Count),Population=as.numeric(Population))%>%
  arrange(Subsite,Sex,Year)
race<-read_tsv(file.path(input_dir,"05_SEER17_annual_race_subsite_rates.txt"),show_col_types=FALSE,na=c("NA",""))%>%clean_names()%>%
  transmute(Subsite=SBA_SITE_MAIN,RaceEthnicity=Race_ETH5,Year=as.integer(YEAR_SINGLE_2000_2023),
            Rate=as.numeric(`Age-Adjusted Rate`),StandardError=as.numeric(`Standard Error`),Count=as.integer(Count),Population=as.numeric(Population))%>%
  arrange(Subsite,RaceEthnicity,Year)

# User-supplied rate files. Rows with zero rate or unavailable SE cannot be fit
# using the log-linear rate + SE model and are retained separately for audit.
write_tsv(valid_rate(age),file.path(out_dir,"JP01_SEER8_age_subsite_1975_2023.txt"),na="NA")
write_tsv(valid_rate(sex),file.path(out_dir,"JP02_SEER8_sex_subsite_1975_2023.txt"),na="NA")
write_tsv(valid_rate(race),file.path(out_dir,"JP03_SEER17_race_subsite_2000_2023.txt"),na="NA")
write_tsv(valid_rate(sex%>%filter(Sex=="Male and female"))%>%select(-Sex),file.path(out_dir,"JP00_SEER8_overall_subsite_1975_2023.txt"),na="NA")
excluded<-bind_rows(
  age%>%filter(is.na(StandardError)|Rate<=0|StandardError<=0)%>%mutate(source="Age"),
  sex%>%filter(is.na(StandardError)|Rate<=0|StandardError<=0)%>%mutate(source="Sex"),
  race%>%filter(is.na(StandardError)|Rate<=0|StandardError<=0)%>%mutate(source="Race"))
write_csv(excluded,file.path(out_dir,"excluded_zero_or_missing_SE.csv"),na="NA")

spec<-tribble(
  ~file,~by_variables,~independent,~dependent,~standard_error,~period,~primary_use,
  "JP00_SEER8_overall_subsite_1975_2023.txt","Subsite","Year","Rate","StandardError","1975-2023","Primary overall subsite trends",
  "JP01_SEER8_age_subsite_1975_2023.txt","Subsite, AgeGroup","Year","Rate","StandardError","1975-2023","Age-stratified trends",
  "JP02_SEER8_sex_subsite_1975_2023.txt","Subsite, Sex","Year","Rate","StandardError","1975-2023","Sex-stratified trends",
  "JP03_SEER17_race_subsite_2000_2023.txt","Subsite, RaceEthnicity","Year","Rate","StandardError","2000-2023","Race/ethnicity-stratified trends")
write_csv(spec,file.path(out_dir,"Joinpoint_analysis_specifications.csv"))

qc<-tibble(dataset=c("Overall","Age","Sex","Race"),
           rows=c(nrow(valid_rate(sex%>%filter(Sex=="Male and female"))),nrow(valid_rate(age)),nrow(valid_rate(sex)),nrow(valid_rate(race))),
           series=c(n_distinct(valid_rate(sex%>%filter(Sex=="Male and female"))$Subsite),
                    n_distinct(interaction(valid_rate(age)$Subsite,valid_rate(age)$AgeGroup)),
                    n_distinct(interaction(valid_rate(sex)$Subsite,valid_rate(sex)$Sex)),
                    n_distinct(interaction(valid_rate(race)$Subsite,valid_rate(race)$RaceEthnicity))),
           min_year=c(min(valid_rate(sex)$Year),min(valid_rate(age)$Year),min(valid_rate(sex)$Year),min(valid_rate(race)$Year)),
           max_year=c(max(valid_rate(sex)$Year),max(valid_rate(age)$Year),max(valid_rate(sex)$Year),max(valid_rate(race)$Year)))
write_csv(qc,file.path(out_dir,"Joinpoint_input_QC.csv"))
cat("Joinpoint inputs prepared:",normalizePath(out_dir,winslash="/"),"\n");print(qc)
