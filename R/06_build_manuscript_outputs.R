#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_RESULTS_DIR")
if (!nzchar(root)) stop("Supply the results directory as the first argument or SBA_RESULTS_DIR.")
out <- file.path(root,"manuscript")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# Tables ------------------------------------------------------------------
t1 <- read_csv(file.path(root,"table1/table1A_early_vs_older.csv"),show_col_types=FALSE)
write_csv(t1,file.path(out,"Table1_early_vs_later_onset.csv"))

jp_dir <- file.path(root,"joinpoint_input")
jp00 <- read_tsv(file.path(jp_dir,"JP00_SEER8_overall_subsite_1975_2023.Export.APC.txt"),show_col_types=FALSE)
table2 <- jp00 %>% transmute(
  Subsite,Period=paste(`Segment Start`,`Segment End`,sep="-"),
  APC=round(APC,2),CI=paste0(round(`APC 95% LCL`,2)," to ",round(`APC 95% UCL`,2)),
  P_value=`P-Value`,Significant=if_else(`APC Significant`==1,"Yes","No"))
write_csv(table2,file.path(out,"Table2_joinpoint_overall_subsite.csv"))

cox <- read_csv(file.path(root,"survival/cox_OS_models.csv"),show_col_types=FALSE) %>%
  filter(model=="Clinical + treatment",term %in% c("early_onset20-49","site_mainJejunum/Ileum")) %>%
  transmute(Analysis="Global Cox OS (secondary)",Exposure=case_when(
    term=="early_onset20-49"~"Age 20-49 vs 50+",TRUE~"Jejunum/Ileum vs duodenum"),
    Interval="Entire follow-up",Estimate=HR,CI_low,CI_high,p_value)
fg <- read_csv(file.path(root,"survival/finegray_CSS_models.csv"),show_col_types=FALSE) %>%
  filter(model=="Clinical + treatment",term %in% c("early_onset20-49","site_mainJejunum/Ileum")) %>%
  transmute(Analysis="Fine-Gray CSS (supportive)",Exposure=case_when(
    term=="early_onset20-49"~"Age 20-49 vs 50+",TRUE~"Jejunum/Ileum vs duodenum"),
    Interval="Entire follow-up",Estimate=sHR,CI_low,CI_high,p_value)
pw <- read_csv(file.path(root,"survival_robust/piecewise_cox_OS.csv"),show_col_types=FALSE) %>%
  filter(exposure %in% c("Age 20-49 vs 50+","Jejunum/Ileum vs duodenum")) %>%
  transmute(Analysis="Interval-specific Cox OS (primary)",
            Exposure=exposure,Interval=followup_interval,Estimate=HR,CI_low,CI_high,p_value)
table3 <- bind_rows(cox,fg,pw) %>% mutate(across(c(Estimate,CI_low,CI_high),~round(.x,3)))
write_csv(table3,file.path(out,"Table3_survival_relative_effects.csv"))

rmst <- read_csv(file.path(root,"survival_robust/RMST_OS_unadjusted.csv"),show_col_types=FALSE)
write_csv(rmst,file.path(out,"TableS1_RMST_OS.csv"))

# Figure 2: annual incidence trends --------------------------------------
inc <- read_tsv(file.path(jp_dir,"JP00_SEER8_overall_subsite_1975_2023.txt"),show_col_types=FALSE)
names(inc) <- make.names(names(inc))
sub_col <- grep("Subsite",names(inc),value=TRUE)[1]
year_col <- grep("Year",names(inc),value=TRUE)[1]
rate_col <- grep("Rate",names(inc),value=TRUE)[1]
plot_inc <- tibble(Subsite=inc[[sub_col]],Year=as.numeric(inc[[year_col]]),Rate=as.numeric(inc[[rate_col]]))
p2 <- ggplot(plot_inc,aes(Year,Rate,color=Subsite)) +
  geom_line(linewidth=.8,alpha=.9) + geom_point(size=.8,alpha=.6) +
  scale_color_brewer(palette="Dark2") +
  labs(x="Year of diagnosis",y="Age-adjusted incidence rate per 100,000",color=NULL) +
  theme_classic(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(out,"Figure2_incidence_by_subsite.png"),p2,width=7.2,height=5.2,dpi=600)
ggsave(file.path(out,"Figure2_incidence_by_subsite.pdf"),p2,width=7.2,height=5.2)

# Figure 4: interval-specific adjusted HRs -------------------------------
forest <- read_csv(file.path(root,"survival_robust/piecewise_cox_OS.csv"),show_col_types=FALSE) %>%
  filter(exposure %in% c("Age 20-49 vs 50+","Jejunum/Ileum vs duodenum")) %>%
  mutate(followup_interval=factor(followup_interval,levels=c("0-12","13-36","37-60",">60")))
p4 <- ggplot(forest,aes(HR,followup_interval,color=exposure)) +
  geom_vline(xintercept=1,linetype=2,color="grey45") +
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.16,position=position_dodge(width=.42)) +
  geom_point(size=2.2,position=position_dodge(width=.42)) +
  scale_x_log10(breaks=c(.25,.5,1,2),labels=label_number(accuracy=.01)) +
  scale_color_manual(values=c("#0072B2","#D55E00")) +
  labs(x="Adjusted hazard ratio (log scale)",y="Follow-up interval, months",color=NULL) +
  theme_classic(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(out,"Figure4_interval_specific_HR.png"),p4,width=7.2,height=4.8,dpi=600)
ggsave(file.path(out,"Figure4_interval_specific_HR.pdf"),p4,width=7.2,height=4.8)

cat("Manuscript tables and figures written to",normalizePath(out,winslash="/"),"\n")
