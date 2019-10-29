rm(list=ls())
options(stringsAsFactors = F)
# install.packages("survminer")
library(survminer)
a=read.table('PLEKHA5-BRCA.tsv',header = T,sep = '\t')
head(a)
dat=a[a$sample_type.samples=='Primary Tumor',4:6]
head(dat)
dat$vital_status.demographic=ifelse(dat$vital_status.demographic=='Alive',0,1)
surv_rnaseq.cut <- surv_cutpoint(
  dat,
  time = "OS.time",
  event = "vital_status.demographic",
  variables = c("ENSG00000052126.13")
)
summary(surv_rnaseq.cut)
plot(surv_rnaseq.cut, "ENSG00000052126.13", palette = "npg")

surv_rnaseq.cat <- surv_categorize(surv_rnaseq.cut) 

library(survival)
fit <- survfit(Surv(OS.time, vital_status.demographic) ~ ENSG00000052126.13,
               data = surv_rnaseq.cat)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,3000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1000,    # break X axis in time intervals by 500. 
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)


