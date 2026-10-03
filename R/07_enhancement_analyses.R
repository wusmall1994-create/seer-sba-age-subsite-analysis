#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival); library(splines)
  library(broom); library(mice); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_PROJECT_ROOT")
if (!nzchar(project_root)) stop("Supply the project root as the first argument or SBA_PROJECT_ROOT.")
clean_dir <- file.path(project_root, "cleaned")
root <- file.path(project_root, "results")
out <- file.path(root,"enhancement")
fig_dir <- file.path(root,"manuscript")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
dat <- readRDS(file.path(clean_dir,"sba_master_clean.rds"))

marital2 <- function(x) case_when(
  x=="Married (including common law)" ~ "Married",
  is.na(x) | x=="Unknown" ~ "Unknown", TRUE ~ "Not married")
period4 <- function(y) factor(case_when(
  y<=2008~"2004-2008", y<=2013~"2009-2013",
  y<=2018~"2014-2018", TRUE~"2019-2022"),
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
  period=period4(diagnosis_year),
  histology=relevel(factor(histology_group),ref="Other adenocarcinoma"),
  stage=relevel(factor(stage4),ref="Localized"),
  grade=relevel(factor(grade4),ref="Grade I"),
  surgery=relevel(factor(surgery),ref="No"),
  chemotherapy=relevel(factor(chemotherapy),ref="No/unknown"),
  radiation=relevel(factor(radiation),ref="None/unknown")) %>%
  filter(complete.cases(across(-grade)))
stopifnot(nrow(d)==11563L)

# 1. Age-stratified Joinpoint results already estimated in the official program.
jp <- read_tsv(file.path(root,"joinpoint_input","JP01_SEER8_age_subsite_1975_2023.Export.APC.txt"),show_col_types=FALSE)
jp_tab <- jp %>% transmute(
  Subsite,AgeGroup,
  Period=paste(`Segment Start`,`Segment End`,sep="-"),
  APC=APC,CI_low=`APC 95% LCL`,CI_high=`APC 95% UCL`,P_value=`P-Value`,
  Significant=if_else(`APC Significant`==1,"Yes","No"))
write_csv(jp_tab,file.path(out,"Table_age_stratified_joinpoint.csv"))

# 2. Continuous time-varying coefficient on log follow-up time. This model
# avoids arbitrary interval cut points and uses the same adjustment set as the
# primary conditional interval models.
d$early_num <- as.integer(d$early=="Age 20-49")
tv_data <- survSplit(Surv(time,event)~.,data=d,cut=seq(3,117,by=3),
                     start="tstart",end="tstop",event="event",episode="episode")
tv_data$log_mid <- log1p((tv_data$tstart+tv_data$tstop)/2)
tv_fit <- coxph(Surv(tstart,tstop,event) ~ early_num + early_num:log_mid +
                  site + sex + race + marital + period + histology + stage +
                  grade + surgery + chemotherapy + radiation,
                data=tv_data,ties="efron",robust=TRUE,cluster=record_id)
times <- 1:120
bn <- names(coef(tv_fit)); b0 <- coef(tv_fit)["early_num"]; b1 <- coef(tv_fit)["early_num:log_mid"]
V <- vcov(tv_fit)[c("early_num","early_num:log_mid"),c("early_num","early_num:log_mid")]
tv <- bind_rows(lapply(times,function(tm){
  L <- c(1,log1p(tm)); est <- b0+b1*log1p(tm); se <- sqrt(drop(t(L)%*%V%*%L))
  tibble(month=tm,HR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se))
}))
write_csv(tv,file.path(out,"continuous_time_varying_HR_early_vs_later.csv"))

# Restricted-cubic-spline sensitivity analysis for the time-varying early-onset
# coefficient. Three degrees of freedom allow smooth nonlinearity without
# overfitting the smaller late risk sets.
sb <- ns(tv_data$log_mid,df=3)
tv_data$s1 <- sb[,1]; tv_data$s2 <- sb[,2]; tv_data$s3 <- sb[,3]
spline_fit <- coxph(Surv(tstart,tstop,event) ~ early_num +
                      early_num:s1 + early_num:s2 + early_num:s3 +
                      site + sex + race + marital + period + histology + stage +
                      grade + surgery + chemotherapy + radiation,
                    data=tv_data,ties="efron",robust=TRUE,cluster=record_id)
pred_times <- c(6,12,36,60,120)
new_basis <- predict(sb,newx=log1p(pred_times))
spline_terms <- c("early_num","early_num:s1","early_num:s2","early_num:s3")
spline_beta <- coef(spline_fit)[spline_terms]
spline_V <- vcov(spline_fit)[spline_terms,spline_terms]
spline_hr <- bind_rows(lapply(seq_along(pred_times),function(i){
  L <- c(1,new_basis[i,]); est <- sum(L*spline_beta)
  se <- sqrt(drop(t(L)%*%spline_V%*%L))
  tibble(month=pred_times[i],HR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se))
}))
write_csv(spline_hr,file.path(out,"spline_time_varying_HR_sensitivity.csv"))

# 3. Covariate-standardized RMST using stabilized inverse-probability weights.
rmst_one <- function(time,event,w,tau) {
  sf <- survfit(Surv(time,event)~1,weights=w,robust=TRUE)
  tt <- c(0,sf$time[sf$time<tau],tau); ss <- c(1,sf$surv[sf$time<tau])
  sum(diff(tt)*ss)
}
weighted_rmst <- function(data,var,tau=60,B=300,seed=1) {
  dd <- data
  dd$z <- as.integer(dd[[var]]==levels(dd[[var]])[2])
  psf <- if(var=="early")
    z ~ sex + race + marital + period + histology + stage
  else
    z ~ ns(age,df=4) + sex + race + marital + period + histology + stage
  calc <- function(x) {
    m <- glm(psf,data=x,family=binomial())
    ps <- pmin(pmax(predict(m,type="response"),.01),.99)
    pz <- mean(x$z); w <- ifelse(x$z==1,pz/ps,(1-pz)/(1-ps))
    q <- quantile(w,c(.01,.99)); w <- pmin(pmax(w,q[1]),q[2])
    c(ref=rmst_one(x$time[x$z==0],x$event[x$z==0],w[x$z==0],tau),
      exposed=rmst_one(x$time[x$z==1],x$event[x$z==1],w[x$z==1],tau))
  }
  point <- calc(dd); set.seed(seed)
  boots <- replicate(B,{
    ii0 <- sample(which(dd$z==0),sum(dd$z==0),replace=TRUE)
    ii1 <- sample(which(dd$z==1),sum(dd$z==1),replace=TRUE)
    tryCatch(calc(dd[c(ii0,ii1),]),error=function(e)c(ref=NA,exposed=NA))
  })
  dif <- boots["exposed",]-boots["ref",]
  tibble(variable=var,contrast=paste(levels(dd[[var]])[2],"vs",levels(dd[[var]])[1]),
         tau_months=tau,RMST_reference=point["ref"],RMST_exposed=point["exposed"],
         adjusted_difference=point["exposed"]-point["ref"],
         CI_low=quantile(dif,.025,na.rm=TRUE),CI_high=quantile(dif,.975,na.rm=TRUE),
         bootstrap_replicates=B)
}
rmst_adj <- bind_rows(
  weighted_rmst(d,"early",60,1000,20260922),weighted_rmst(d,"early",120,1000,20260923),
  weighted_rmst(filter(d,site %in% c("Duodenum","Jejunum/Ileum")),"site",60,1000,20260924),
  weighted_rmst(filter(d,site %in% c("Duodenum","Jejunum/Ileum")),"site",120,1000,20260925))
write_csv(rmst_adj,file.path(out,"adjusted_RMST_IPW.csv"))

make_weights <- function(data,var,type=c("iptw","overlap")) {
  type <- match.arg(type); x <- data
  x$z <- as.integer(x[[var]]==levels(x[[var]])[2])
  psf <- if(var=="early")
    z ~ sex + race + marital + period + histology + stage
  else
    z ~ ns(age,df=4) + sex + race + marital + period + histology + stage
  m <- glm(psf,data=x,family=binomial())
  ps <- pmin(pmax(predict(m,type="response"),.01),.99); pz <- mean(x$z)
  if(type=="iptw") {
    w <- ifelse(x$z==1,pz/ps,(1-pz)/(1-ps)); q <- quantile(w,c(.01,.99))
    w <- pmin(pmax(w,q[1]),q[2])
  } else w <- ifelse(x$z==1,1-ps,ps)
  x$weight <- w; x
}
smd_binary <- function(x,z,w=rep(1,length(z))) {
  p1 <- weighted.mean(x[z==1],w[z==1]); p0 <- weighted.mean(x[z==0],w[z==0])
  (p1-p0)/sqrt((p1*(1-p1)+p0*(1-p0))/2)
}
balance_table <- function(x,label) {
  vars <- c("sex","race","marital","period","histology","stage")
  bind_rows(lapply(vars,function(v) bind_rows(lapply(levels(x[[v]]),function(lev){
    xx <- as.numeric(x[[v]]==lev)
    tibble(analysis=label,variable=v,level=lev,
           unweighted_SMD=smd_binary(xx,x$z),weighted_SMD=smd_binary(xx,x$z,x$weight))
  }))))
}
diag_one <- function(x,label,method) tibble(
  analysis=label,method,N=nrow(x),ESS=(sum(x$weight)^2)/sum(x$weight^2),
  minimum=min(x$weight),p01=quantile(x$weight,.01),median=median(x$weight),
  mean=mean(x$weight),p99=quantile(x$weight,.99),maximum=max(x$weight))
age_w <- make_weights(d,"early","iptw")
site_w <- make_weights(filter(d,site %in% c("Duodenum","Jejunum/Ileum")),"site","iptw")
age_ow <- make_weights(d,"early","overlap")
site_ow <- make_weights(filter(d,site %in% c("Duodenum","Jejunum/Ileum")),"site","overlap")
write_csv(bind_rows(diag_one(age_w,"Age 20-49 vs Age 50+","Stabilized IPTW, 1st/99th truncation"),
                    diag_one(site_w,"Jejunum/Ileum vs Duodenum","Stabilized IPTW, 1st/99th truncation"),
                    diag_one(age_ow,"Age 20-49 vs Age 50+","Overlap weights"),
                    diag_one(site_ow,"Jejunum/Ileum vs Duodenum","Overlap weights")),
          file.path(out,"ipw_weight_diagnostics.csv"))
write_csv(bind_rows(balance_table(age_w,"Age 20-49 vs Age 50+ | IPTW"),
                    balance_table(site_w,"Jejunum/Ileum vs Duodenum | IPTW"),
                    balance_table(age_ow,"Age 20-49 vs Age 50+ | overlap"),
                    balance_table(site_ow,"Jejunum/Ileum vs Duodenum | overlap")),
          file.path(out,"ipw_balance_smd.csv"))

ow_rmst <- bind_rows(lapply(list(
  list(age_ow,"early",60),list(age_ow,"early",120),
  list(site_ow,"site",60),list(site_ow,"site",120)),function(a){
    x <- a[[1]]; var <- a[[2]]; tau <- a[[3]]
    r0 <- rmst_one(x$time[x$z==0],x$event[x$z==0],x$weight[x$z==0],tau)
    r1 <- rmst_one(x$time[x$z==1],x$event[x$z==1],x$weight[x$z==1],tau)
    tibble(variable=var,tau_months=tau,RMST_reference=r0,RMST_exposed=r1,difference=r1-r0)
  }))
write_csv(ow_rmst,file.path(out,"overlap_weighted_RMST_sensitivity.csv"))

# Weighted survival curves for the age contrast.
sf <- survfit(Surv(time,event)~early,data=age_w,weights=age_w$weight,robust=TRUE)
ss <- summary(sf,times=seq(0,120,by=3),extend=TRUE)
curve <- tibble(time=ss$time,surv=ss$surv,lower=ss$lower,upper=ss$upper,
                group=sub("early=","",ss$strata))
write_csv(curve,file.path(out,"adjusted_survival_curve_age.csv"))

# 4. Grade missingness and multiple-imputation sensitivity analysis.
grade_audit <- dat %>% filter(cohort_survival_primary) %>% mutate(
  grade_missing=grade4=="Unknown",period=period4(diagnosis_year)) %>%
  count(period,site_main,grade_missing,name="N") %>% group_by(period,site_main) %>%
  mutate(percent=100*N/sum(N)) %>% ungroup()
write_csv(grade_audit,file.path(out,"grade_missingness_by_period_site.csv"))

# Grade fields changed structurally after 2018; imputation is therefore
# restricted to 2004-2018, where the legacy grade construct is observable.
mi <- d %>% filter(diagnosis_year<=2018) %>%
  transmute(time,event,diagnosis_year,early,site,sex,race,marital,period,histology,stage,
            grade=if_else(as.character(grade)=="Unknown",NA_character_,as.character(grade)),
                      surgery,chemotherapy,radiation)
mi$grade <- factor(mi$grade)
bh <- basehaz(coxph(Surv(time,event)~1,data=mi),centered=FALSE)
mi$nelson_aalen <- approx(c(0,bh$time),c(0,bh$hazard),xout=mi$time,
                          method="constant",f=0,rule=2)$y
ini <- mice(mi,maxit=0,printFlag=FALSE)
method <- ini$method; method[] <- ""; method["grade"] <- "polyreg"
pred <- ini$predictorMatrix; pred[,] <- 0
pred["grade",c("event","nelson_aalen","diagnosis_year","early","site","sex","race","marital","period","histology","stage","surgery","chemotherapy","radiation")] <- 1
imp <- mice(mi,m=10,maxit=10,method=method,predictorMatrix=pred,seed=20260922,printFlag=FALSE)
mi_fit <- with(imp,coxph(Surv(time,event)~early+site+sex+race+marital+period+histology+stage+grade+surgery+chemotherapy+radiation))
mi_pool <- summary(pool(mi_fit),conf.int=TRUE,exponentiate=TRUE) %>%
  filter(term %in% c("earlyAge 20-49","siteJejunum/Ileum"))
write_csv(mi_pool,file.path(out,"multiple_imputation_grade_sensitivity.csv"))

# Publication figures (R-only workflow): source CSVs plus vector/raster outputs.
theme_pub <- theme_classic(base_size=9)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
p1 <- ggplot(tibble(x=c(1,1,1),y=c(3,2,1),
                    label=c("14,158 adults with SBA\nSEER 17, 2000-2023",
                            "12,391 stage-harmonized cases\nDiagnosed 2004-2023",
                            "11,563 primary survival cohort\nDiagnosed 2004-2022")),aes(x,y))+
  geom_label(aes(label=label),size=3.4,label.size=.45,fill=c("#E8F1F8","#EAF5EA","#FFF3DE"),label.padding=unit(.45,"lines"))+
  annotate("segment",x=1,xend=1,y=2.72,yend=2.28,arrow=arrow(length=unit(.16,"inches")),linewidth=.6)+
  annotate("segment",x=1,xend=1,y=1.72,yend=1.28,arrow=arrow(length=unit(.16,"inches")),linewidth=.6)+
  annotate("text",x=1.65,y=2.5,label="Exclude pre-2004 stage-incompatible records",hjust=0,size=2.7)+
  annotate("text",x=1.65,y=1.5,label="Exclude 2023 survival-truncation records\nand invalid follow-up",hjust=0,size=2.7)+
  coord_cartesian(xlim=c(.2,3.5),ylim=c(.65,3.35),clip="off")+theme_void()

p3 <- ggplot(curve,aes(time,surv,color=group,fill=group))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,color=NA)+geom_step(linewidth=.75)+
  scale_color_manual(values=c("Age 20-49"="#0072B2","Age 50+"="#D55E00"))+
  scale_fill_manual(values=c("Age 20-49"="#0072B2","Age 50+"="#D55E00"))+
  scale_y_continuous(labels=percent_format(accuracy=1),limits=c(0,1))+
  labs(x="Months after diagnosis",y="IPW-adjusted overall survival",color=NULL,fill=NULL)+theme_pub

p5 <- ggplot(tv,aes(month,HR))+geom_ribbon(aes(ymin=CI_low,ymax=CI_high),fill="#0072B2",alpha=.18)+
  geom_line(color="#0072B2",linewidth=.8)+geom_hline(yintercept=1,linetype=2,color="grey40")+
  scale_y_log10(breaks=c(.25,.5,1,2),labels=label_number(accuracy=.01))+
  labs(x="Months after diagnosis",y="Adjusted HR: age 20-49 vs age 50+ (log scale)")+theme_pub

save_fig <- function(p,stem,w=183/25.4,h=115/25.4){
  ggsave(file.path(fig_dir,paste0(stem,".png")),p,width=w,height=h,dpi=600,bg="white")
  ggsave(file.path(fig_dir,paste0(stem,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")
  ggsave(file.path(fig_dir,paste0(stem,".svg")),p,width=w,height=h,device=svglite::svglite,bg="white")
  ggsave(file.path(fig_dir,paste0(stem,".tiff")),p,width=w,height=h,dpi=600,compression="lzw",bg="white")
}
save_fig(p1,"Figure1_STROBE_flow",183/25.4,120/25.4)
save_fig(p3,"Figure3_adjusted_survival",183/25.4,115/25.4)
save_fig(p5,"Figure5_continuous_time_varying_HR",183/25.4,110/25.4)

write_csv(tibble(cohort=c("Descriptive master","Stage-harmonized","Primary survival"),N=c(14158,12391,11563)),
          file.path(fig_dir,"Figure1_source_data.csv"))
write_csv(curve,file.path(fig_dir,"Figure3_source_data.csv"))
write_csv(tv,file.path(fig_dir,"Figure5_source_data.csv"))

cat("Enhancement analyses complete:",normalizePath(out,winslash="/"),"\n")
