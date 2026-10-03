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

# Figure 3: observed annual incidence trends -----------------------------
inc <- read_tsv(file.path(jp_dir,"JP00_SEER8_overall_subsite_1975_2023.txt"),show_col_types=FALSE)
names(inc) <- make.names(names(inc))
sub_col <- grep("Subsite",names(inc),value=TRUE)[1]
year_col <- grep("Year",names(inc),value=TRUE)[1]
rate_col <- grep("Rate",names(inc),value=TRUE)[1]
plot_inc <- tibble(Subsite=inc[[sub_col]],Year=as.numeric(inc[[year_col]]),Rate=as.numeric(inc[[rate_col]]))
end_labels <- plot_inc %>%
  filter(!is.na(Rate)) %>%
  group_by(Subsite) %>%
  slice_max(Year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = recode(
      Subsite,
      "Small intestine, NOS" = "NOS",
      "Jejunum/Ileum" = "Jejunum/Ileum",
      .default = as.character(Subsite)
    )
  )
p2 <- ggplot(plot_inc,aes(Year,Rate,color=Subsite)) +
  geom_line(aes(linetype=Subsite),linewidth=.8,alpha=.9) +
  geom_text(
    data=end_labels,
    aes(x=Year+1.0,y=Rate,label=label,color=Subsite),
    hjust=0,size=3.2,fontface="bold",show.legend=FALSE
  ) +
  scale_color_brewer(palette="Dark2") +
  scale_linetype_manual(values=c("solid","22","42","13")) +
  scale_x_continuous(limits=c(1975,2035),breaks=seq(1980,2020,10),expand=expansion(mult=c(0.01,0))) +
  scale_y_continuous(breaks=seq(0,0.8,0.2),labels=label_number(accuracy=0.1),expand=expansion(mult=c(0.02,0.08))) +
  labs(x="Year of diagnosis",y="Observed age-adjusted incidence rate per 100,000",color=NULL,linetype=NULL) +
  theme_classic(base_size=11,base_family="Arial") +
  theme(legend.position="none",plot.margin=margin(5.5,34,5.5,5.5))
ragg::agg_tiff(file.path(out,"Figure3_incidence_by_subsite.tiff"),width=7.2,height=5.2,units="in",res=600,compression="lzw")
print(p2)
dev.off()
if (identical(Sys.getenv("SBA_QA_PREVIEW"), "1")) {
  qa_path <- Sys.getenv("SBA_QA_PREVIEW_PATH", file.path(tempdir(), "Figure3_QA_preview.png"))
  ragg::agg_png(qa_path,width=7.2,height=5.2,units="in",res=150)
  print(p2)
  dev.off()
}

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
