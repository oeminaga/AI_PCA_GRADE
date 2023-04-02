#### LIBRARIES ####
library(rms)
library(classifierplots)
library(readr)
library(survival)
library(coxphw)
library(survivalROC)
library(survminer)
####  LOAD FILE ####
ModelName = "PlexusNET_BCR_10x_BEST_APPROACH_HUE"
filename_test_set = paste0("../",ModelName,"_test_set.csv")

data_test <- read_csv(filename_test_set)
data_test$BCR_status=0
data_test$BCR_status[data_test$X1st.BCR.Type != "-"]=1

table(data_test$PID)->count
match(data_test$PID, names(count))->m
match(names(count), data_test$PID)->fir
length(fir)
#### FOR 10x PATCHES ####
data_test$Patch_0= data_test$`0...47`
data_test$Patch_1= data_test$`1...48`
data_test$Patch_2= data_test$`2...49`
data_test$Patch_3= data_test$`3...50`
data_test$Patch_4= data_test$`4`
data_test$Patch_5= data_test$`5`
data_test$Patch_6= data_test$`6`
data_test$Patch_7= data_test$`7`
data_test$Patch_8= data_test$`8`

1/count[m]->data_test$W
#### TMA CORE LEVEL PREDICTION ####
patches_data=data_test[,c("Patch_0","Patch_1","Patch_2","Patch_3","Patch_4","Patch_5","Patch_6", "Patch_7", "Patch_8")]
data_test$PredictionScore=rowMeans(patches_data)
data_test=data_test[data_test$Interval.RP.to.BCR.or.last.contact.death>0,]

#### CASE-LEVEL PREDICTION
vl_case_pred = c()
for (f in names(count)){
  vl=data_test$PredictionScore[data_test$PID==f]
  vl_case_pred =c(vl_case_pred ,mean(vl))
}

data_test_cases=data_test[fir,]
data_test_cases$PredictionCaseLevel=vl_case_pred

#### TMA CORE LEVEL EVALUATION ####
cox_score_spots <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionScore,data=data_test,c(data_test$W),x=TRUE,y=TRUE, )
summary(cox_score_spots)

BIC(cox_score_spots)
AIC(cox_score_spots)

coxphw_score_spots<-coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~PredictionScore, data=data_test, caseweights = data_test$W,  robust=TRUE)
summary(coxphw_score_spots)

#### CASE LEVEL EVALUATION ####
cox_score_cases<- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionCaseLevel,data=data_test_cases,x=TRUE,y=TRUE, )
summary(cox_score_cases)

BIC(cox_score_cases)
AIC(cox_score_cases)

coxphw_score_cases=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionCaseLevel,data=data_test_cases,robust=TRUE)
summary(coxphw_score_cases)

survivalROC(data_test_cases$Interval.RP.to.BCR.or.last.contact.death, data_test_cases$BCR_status, marker = haz.aug,predict.time = 60, method = "KM")->t1
t1$AUC

survivalROC(data_test_cases$Interval.RP.to.BCR.or.last.contact.death, data_test_cases$BCR_status, marker = haz.aug,predict.time = 120, method = "KM")->t1
t1$AUC

###### CALIBRATION EVALUATION ######
f <- cph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status) ~ PredictionCaseLevel, x=TRUE, y=TRUE, surv=TRUE, time.inc=120, data=data_test_cases)
cal <- calibrate(f, u=120, cmethod='KM', m=50, B=2000)
plot(cal)

#### CATEGORIZATION #####
### USE the CHAID algorithms provided by SPSS and use the cutoffs to categorize the patients
### Categories are: [0,0.06),[0.06, 0.43), [0.43,0.75), [0.75,1]
groups_risk_=cut(data_test_cases$PredictionCaseLevel, breaks = c(0,0.06, 0.43,0.75,1),
                labels = c("0-5","6-42","43-74", "75-100"),
                include.lowest = TRUE)
data_test_cases$Risk_Group=groups_risk_

coxphw_score_cases_risk_group=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group,data=data_test_cases,robust=TRUE)
summary(coxphw_score_cases_risk_group)

cox_score_cases_risk_group<- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group,data=data_test_cases,x=TRUE,y=TRUE, )
summary(cox_score_cases_risk_group)

BIC(cox_score_cases_risk_group)
AIC(cox_score_cases_risk_group)
#### COMPARISON BETWEEE AI AND GG ####
########## COXPHW ##########
##---------- BCR --------###
#Univariate
#CPCBN VARIABLES#
UnivariateAnalysis.Risk_Group=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.AgeAtDiagnosis=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~AgeAtDiagnosis,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.IDC_status=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~IDC_status,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.PSM_status=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PSM_status,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.pT_stage=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~pT_stage,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.pN_stage=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~pN_stage,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.PSA_level=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PSA_level,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.GG_RPE=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~GG_RPE,data=data_test_cases,robust=TRUE)
UnivariateAnalysis.PredictionCaseLevel=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionCaseLevel,data=data_test_cases,robust=TRUE)
# NOTE: FOR PROCURE, we considerd PSA_Level, pN_stage, pT_stage, PSM_status in addition to GG_RPE, Risk_Group or PredictionCaseLevel#
summary(UnivariateAnalysis.AgeAtDiagnosis)
summary(UnivariateAnalysis.Risk_Group)
summary(UnivariateAnalysis.GG_RPE)
summary(UnivariateAnalysis.PredictionCaseLevel)
summary(UnivariateAnalysis.IDC_status)
summary(UnivariateAnalysis.PSM_status)
summary(UnivariateAnalysis.PSA_level)
summary(UnivariateAnalysis.pT_stage)
summary(UnivariateAnalysis.pN_stage)
#Multivariate ##
#Multivariate - CPCBN
MultivariateAnalyses_BCR_RISK_GROUP_CPCBN_w_pN_stage=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group+PSA_Level+IDC_status+pT_stage+pN_stage+RPE_GG,data=data_test_cases,robust=TRUE)
MultivariateAnalyses_BCR_RISK_GROUP_CPCBN_wo_pN_stage=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group+PSA_Level+IDC_status+pT_stage+RPE_GG,data=data_test_cases,robust=TRUE)
MultivariateAnalyses_BCR_PredictionCaseLevel_CPCBN=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionCaseLevel+PSA_Level+IDC_status+pT_stage+RPE_GG,data=data_test_cases,robust=TRUE)
summary(MultivariateAnalyses_BCR_RISK_GROUP_CPCBN_w_pN_stage)
summary(MultivariateAnalyses_BCR_RISK_GROUP_CPCBN_wo_pN_stage)
summary(MultivariateAnalyses_BCR_PredictionCaseLevel_CPCBN)

#Multivariate - PROCURE
# data_test_cases_cpcbn=data_test_cases
# data_test_cases=data_test_cases_procure
##PROIVIDED HERE FOR CODE COMPLETENESS - WE CHANGED THE DATA SOURCE TO THE PROCURE DATA SET WHEN RUNNING THE FOLLOWING TESTS
MultivariateAnalyses_BCR_RISK_GROUP_PROCURE=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group+PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,robust=TRUE)
MultivariateAnalyses_BCR_RISK_GROUP_PROCURE_W_RPE_GG=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~Risk_Group+PSA_Level+pT_stage+pN_stage+PSM_status+RPE_GG,data=data_test_cases,robust=TRUE)
MultivariateAnalyses_BCR_PredictionCaseLevel_PROCURE=coxphw(Surv(Interval.RP.to.BCR.or.last.contact.death,BCR_status)~PredictionCaseLevel+PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,robust=TRUE)
summary(MultivariateAnalyses_BCR_RISK_GROUP_PROCURE)
summary(MultivariateAnalyses_BCR_RISK_GROUP_PROCURE_W_RPE_GG)
summary(MultivariateAnalyses_BCR_PredictionCaseLevel_PROCURE)

########## NESTED PARTIAL LIKELIHOOD RATIO TEST -BCR- ###############
library(nonnestcox)
# --- CPCBN --- #

pltest_full_model <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~Risk_Group+pT_stage+GG_RPE+PSA_Level,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_full_model)
AIC(pltest_full_model)
BIC(pltest_full_model)

pltest_model_RPE_GG_excluded <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~Risk_Group+pT_stage+PSA_Level,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_model_RPE_GG_excluded)
AIC(pltest_model_RPE_GG_excluded)
BIC(pltest_model_RPE_GG_excluded)

pltest_model_RISK_GROUP_excluded <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~GG_RPE+pT_stage+PSA_Level,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_model_RISK_GROUP_excluded)
AIC(pltest_model_RISK_GROUP_excluded)
BIC(pltest_model_RISK_GROUP_excluded)

pltest_basemodel <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~pT_stage+PSA_Level,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_basemodel)
AIC(pltest_basemodel)
BIC(pltest_basemodel)

plrtest(pltest_full_model,pltest_model_RPE_GG_excluded, nested = T)
plrtest(pltest_full_model,pltest_model_RISK_GROUP_excluded, nested = T)
plrtest(pltest_full_model,pltest_basemodel, nested = T)

# --- PROCURE --- #
#PROIVIDED HERE FOR CODE COMPLETENESS - WE CHANGED THE DATA SOURCE TO THE PROCURE DATA SET WHEN RUNNING THE FOLLOWING TESTS
# data_test_cases_cpcbn=data_test_cases
# data_test_cases=data_test_cases_procure
pltest_full_model <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~Risk_Group+GG_RPE+PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_full_model)
AIC(pltest_full_model)
BIC(pltest_full_model)

pltest_model_RPE_GG_excluded <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~Risk_Group+PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_model_RPE_GG_excluded)
AIC(pltest_model_RPE_GG_excluded)
BIC(pltest_model_RPE_GG_excluded)

pltest_model_RISK_GROUP_excluded <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~GG_RPE+PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_model_RISK_GROUP_excluded)
AIC(pltest_model_RISK_GROUP_excluded)
BIC(pltest_model_RISK_GROUP_excluded)

pltest_basemodel <- coxph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status)~PSA_Level+pT_stage+pN_stage+PSM_status,data=data_test_cases,x=TRUE,y=TRUE)
summary(pltest_basemodel)
AIC(pltest_basemodel)
BIC(pltest_basemodel)

plrtest(pltest_full_model,pltest_model_RPE_GG_excluded, nested = T)
plrtest(pltest_full_model,pltest_model_RISK_GROUP_excluded, nested = T)
plrtest(pltest_full_model,pltest_basemodel, nested = T)


########## KM Curves ##########
### APPLIED TO CPCBN and PROCURE ###
# -------- Risk groups ------- #
km_test_set_Risk_Group <- survfit(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status) ~Risk_Group  , data=data_test_cases)
summary(km_test_set_Risk_Group)

ggsurv_km_test_set_Risk_Group <- ggsurvplot(
  km_test_set_Risk_Group,                     # survfit object with calculated statistics.
  data = data_test_cases,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = "jco",
  xlab = "Time from curative treatment in months",   # customize X axis label.
  ylab = "BCR-free Survival probability",   # customize Y axis label.
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.2, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.2,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
)
ggsurv_km_test_set_Risk_Group
# --------     GG      ------- #
km_test_set_GG_RPE <- survfit(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status) ~GG_RPE  , data=data_test_cases)
summary(km_test_set_GG_RPE)

ggsurv_km_test_set_GG <- ggsurvplot(
  km_test_set_GG_RPE,                     # survfit object with calculated statistics.
  data = data_test_cases,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = "jco",
  xlab = "Time from curative treatment in months",   # customize X axis label.
  ylab = "BCR-free Survival probability",   # customize Y axis label.
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.2, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.2,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
)
ggsurv_km_test_set_GG

##---------- CSS --------##
# -------- Risk groups ------- #
km_trt_fit_test_set_all_Risk_Group_CSS <- survfit(Surv(Follow.up.Interval.last.contact.or.death.and.RP, DEATH) ~Risk_Group  , data=data_test_cases)
summary(km_trt_fit_test_set_all_Risk_Group_CSS)

ggsurv_km_test_set_Risk_Group_CSS <- ggsurvplot(
  km_trt_fit_test_set_all_Risk_Group_CSS,                     # survfit object with calculated statistics.
  data = data_test_cases,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = "jco",
  xlab = "Time from curative treatment in months",   # customize X axis label.
  ylab = "Cancer-specific Survival probability",   # customize Y axis label.
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.2, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.2,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
)
ggsurv_km_test_set_Risk_Group_CSS
# --------     GG      ------- #
km_trt_fit_test_set_all_GG_CSS <- survfit(Surv(Follow.up.Interval.last.contact.or.death.and.RP, DEATH) ~GG_RPE  , data=data_test_cases)
summary(km_trt_fit_test_set_all_GG_CSS)

ggsurv_km_test_set_GG_CSS <- ggsurvplot(
  km_trt_fit_test_set_all_GG_CSS,                     # survfit object with calculated statistics.
  data = data_test_cases,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = "jco",
  xlab = "Time from curative treatment in months",   # customize X axis label.
  ylab = "Cancer-specific Survival probability",   # customize Y axis label.
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.2, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.2,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
)
ggsurv_km_test_set_GG_CSS