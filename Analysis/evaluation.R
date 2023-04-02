#### LIBRARIES ####
library(rms)
library(classifierplots)
library(readr)
library(survival)
library(coxphw)
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

###### CALIBRATION EVALUATION ######
f <- cph(Surv(Interval.RP.to.BCR.or.last.contact.death, BCR_status) ~ PredictionCaseLevel, x=TRUE, y=TRUE, surv=TRUE, time.inc=120, data=data_test_cases)
cal <- calibrate(f, u=120, cmethod='KM', m=50, B=2000)
plot(cal)

#### CATEGORIZATION #####
### USE the CHAID algorithms provided by SPSS and use the cutoffs to categorize the patients
### Cutoffs are: 0.06,0.43, 0.75

#### COMPARISON BETWEEE AI AND GG ####
