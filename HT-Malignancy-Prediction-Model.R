## Project: Prediction of 1-Year Malignancy After Heart Transplant
## Script purpose: Code used in Paper
##################################################

#Required Libraries
library(sas7bdat)                             # Package to read in SAS files (sas7bdat)
library(glmnet) 
library(mice)                                 # Package for evaluating and imputing missing data
library(lattice)
library(dplyr)
library(matrixStats)
library(rms) #
library(naniar)#
library(ggplot2)
library(DescTools)
library(haven)
library(pmsampsize)
library(DescTools)



####################################################################################### 
# Data pre-process 
####################################################################################### 


data_base <- read.csv("directory to data_base") #edit me 
data_base<-data.frame(data_base)
#any blanks now are NA 
data_base[data_base == ""] <- NA 
# Assign NA to weight less than 35 and heights less than 90
omit_weight=c(which(data_base$REC_WGT_KG<35))
data_base$REC_WGT_KG[omit_weight]=NA
omit_hgt=c(which(data_base$REC_HGT_CM<90))
data_base$REC_HGT_CM[omit_hgt]=NA

#recalculate REC_CLCR: 
#before REC_CLCR was calculated with NA's with unfiltered(based) on criteria for REC_WGT_KG
#REC_CLCR = (((140 - REC_AGE_AT_TX) * REC_WGT_KG)/(72 * REC_CREAT))*REC_CLCR_GEN
data_base$REC_CLCR_GEN=ifelse(data_base$CAN_GENDER=="M",1,0.85)
data_base$REC_CLCR=(((140 - data_base$REC_AGE_AT_TX) * data_base$REC_WGT_KG)/(72 * data_base$REC_CREAT))*
  data_base$REC_CLCR_GEN

vars_miss = c("OUT_MALIG", "Out_Death_1", "DON_AGE", "DON_HGT_CM", "DON_WGT_KG", "REC_HTN",
              "DON_CREAT", "CAN_GENDER", "REC_BMI", "REC_MALIG", "REC_ECMO", "REC_IABP", "REC_INOTROP",
              "REC_VENTILATOR", "REC_CARDIAC_OUTPUT", "REC_CREAT", "REC_TOT_BILI", "REC_HR_ISCH", 
              "REC_AGE_AT_TX", "CAN_TOT_ALBUMIN", "TX_YEAR", "REC_DRUG_INDUCTION", "REC_ABO", 
              "REC_ANGINA", "HF_ETIOL", "REC_DM", "REC_COPD", "REC_PVD", "CAN_SCD", 
              "DON_ABO1", "DON_ANTI_CMV_POS", "REC_CMV_IGG_POS", "REC_EBV_STAT_POS",
              "REC_HBV_ANTIBODY_POS", "REC_DR_MM_EQUIV_TX", "REC_VAD", "REC_CLCR", "REC_HIST_CIGARETTE",
              "REC_HGT_CM", "REC_WGT_KG", "REC_RACE_WHITE")

## Select desired variables: for outcome
vars_model = c("OUT_MALIG", "REC_DR_MM_EQUIV_TX", "CAN_GENDER", "REC_HGT_CM", "REC_CREAT", "REC_AGE_AT_TX", 
               "TX_YEAR", "REC_DRUG_INDUCTION", "REC_HTN", "REC_HIST_CIGARETTE", "REC_RACE_WHITE", "REC_WGT_KG",
               "DON_ANTI_CMV_POS", "REC_CMV_IGG_POS", "HF_ETIOL", "REC_EBV_STAT_POS", "REC_HBV_ANTIBODY_POS",
               "REC_MALIG", "REC_DM")

data_miss = data_base[, vars_miss]
data_model = data_base[, vars_model]



####################################################################################### 
# Sample size calculation using the 'pmsampsize' package
####################################################################################### 

#pmsampsize - Calculates the minimum sample size required for developing a multivariable prediction model
pmsampsize(type = "b",               # b = binary outcome
           rsquared = 0.15,          # expected value of (Cox-Snell) R-squared in new model
           parameters = 30,          # Number of candidate predictor parameters in model
           prevalence = 0.022)       # Overall outcome proportion/prevalence

## Limit data to desired variables and look at number of patients & variables
CountCompCases(data_model) 
## Look at 'data_model' file (n=38598, 19 variables)
dim(data_model)
# Use the 'transform' function to change the type of each feature
data_miss <- transform(OUT_MALIG=as.factor(OUT_MALIG), Out_Death_1=as.factor(Out_Death_1),
                       CAN_GENDER=as.factor(CAN_GENDER), REC_MALIG=as.factor(REC_MALIG), 
                       REC_ECMO=as.factor(REC_ECMO), REC_IABP=as.factor(REC_IABP), 
                       REC_INOTROP=as.factor(REC_INOTROP), REC_VENTILATOR=as.factor(REC_VENTILATOR),
                       REC_DRUG_INDUCTION=as.factor(REC_DRUG_INDUCTION), REC_ABO=as.factor(REC_ABO),
                       REC_ANGINA=as.factor(REC_ANGINA), HF_ETIOL=as.factor(HF_ETIOL),
                       REC_DM=as.factor(REC_DM), REC_COPD=as.factor(REC_COPD), REC_PVD=as.factor(REC_PVD),
                       CAN_SCD=as.factor(CAN_SCD), DON_ABO1=as.factor(DON_ABO1), 
                       DON_ANTI_CMV_POS=as.factor(DON_ANTI_CMV_POS), REC_CMV_IGG_POS=as.factor(REC_CMV_IGG_POS),
                       REC_DR_MM_EQUIV_TX=as.factor(REC_DR_MM_EQUIV_TX), REC_VAD=as.factor(REC_VAD),
                       REC_HIST_CIGARETTE=as.factor(REC_HIST_CIGARETTE),REC_HTN=as.factor(REC_HTN),
                       data_miss)

data_model <- transform(OUT_MALIG=as.factor(OUT_MALIG), 
                        REC_DRUG_INDUCTION=as.factor(REC_DRUG_INDUCTION), REC_HTN=as.factor(REC_HTN),
                        REC_RACE_WHITE=as.factor(REC_RACE_WHITE), DON_ANTI_CMV_POS=as.factor(DON_ANTI_CMV_POS),
                        REC_CMV_IGG_POS=as.factor(REC_CMV_IGG_POS), REC_EBV_STAT_POS=as.factor(REC_EBV_STAT_POS),
                        REC_HBV_ANTIBODY_POS=as.factor(REC_HBV_ANTIBODY_POS), REC_MALIG=as.factor(REC_MALIG),
                        REC_DM=as.factor(REC_DM),
                        data_model)
##Missing data_summary
miss_var_summary(data_model)
miss_var_summary(data_miss)



####################################################################################### 
# Identify Auxiliary Variables for Missing Relevant Covariates
####################################################################################### 

#Missing covariates: REC_HTN,CAN_HIST_CIGARETTE,OUT_MALIG,
#REC_DR_MM_EQUIV_TX,REC_DRUG_INDUCTION,REC_HGT_CM, REC_WGT_KG,REC_CREAT,REC_DM
#Run regression to identify factors associated with missingness
#lrm: Fit binary and proportional odds ordinal logistic regression models 
#using maximum likelihood estimation or penalized maximum likelihood estimation
library(rms)
set.seed(0621)
dd <- datadist(data_miss)
options(datadist="dd")

HTN_MISS_REG = lrm(is.na(REC_HTN) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                     rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                     CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                     REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                     rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                     rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                     REC_DRUG_INDUCTION + REC_ABO + REC_ANGINA + HF_ETIOL + 
                     REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                     DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                     REC_HBV_ANTIBODY_POS + REC_DR_MM_EQUIV_TX + REC_VAD + 
                     rcs(REC_CLCR,3) + REC_HIST_CIGARETTE + rcs(REC_HGT_CM,3) + 
                     rcs(REC_WGT_KG,3),
                   data=data_miss, x=TRUE, y=TRUE)

###Function Description: Puco ###
# Purpose: Calculate proportion of useable cases and outbound statistics 
#proportion of usable cases=puc
#o=outbound stat
##################################
puco=function(dat)
{
  p<-md.pairs(dat)
  #Missing variables: variables to include based on proportion of usable cases
  #proportion of usable cases
  #(Van Buuren, Boshuizen, and Knook 1999)
  puc=p$mr/(p$mr+p$mm) #can use to include imputed variables
  o=p$rm/(p$rm+p$rr)#outbound statistic   
  res=list(as.data.frame(puc),as.data.frame(o));
  return(res)
}

hmr_imp=puco(data_miss)
# View(hmr_imp[[1]])
# hmr_imp[[1]]["REC_HTN",]



#OUT_MALIG
MALIG_MISS_REG = lrm(is.na(OUT_MALIG) ~ REC_HIST_CIGARETTE + Out_Death_1 + rcs(DON_AGE,3) + 
                       rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                       CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                       REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                       rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                       REC_DRUG_INDUCTION + REC_ABO + REC_ANGINA + HF_ETIOL + 
                       REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_DR_MM_EQUIV_TX + REC_VAD + 
                       rcs(REC_CLCR,3) + REC_HTN + rcs(REC_HGT_CM,3) + 
                       rcs(REC_WGT_KG,3)+REC_RACE_WHITE,
                     data=data_miss, x=TRUE, y=TRUE)

#hmr_imp[[1]]["OUT_MALIG",]

#REC_HIST_CIGGARETTE
CIG_MISS_REG = lrm(is.na(REC_HIST_CIGARETTE) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                     rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                     CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                     REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                     rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                     rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + 
                     REC_DRUG_INDUCTION + REC_ABO + REC_ANGINA + HF_ETIOL + 
                     REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                     DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                     REC_HBV_ANTIBODY_POS + REC_DR_MM_EQUIV_TX + REC_VAD + 
                     rcs(REC_CLCR,3) + REC_HTN + rcs(REC_HGT_CM,3) + 
                     rcs(REC_WGT_KG,3)+REC_RACE_WHITE,
                   data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#hmr_imp[[1]]["REC_HIST_CIGARETTE",]

#REC_DR_MM_EQUIV_TX
DR_MISS_REG = lrm(is.na(REC_DR_MM_EQUIV_TX) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                    rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                    CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                    REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                    rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                    rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                    REC_DRUG_INDUCTION + REC_ABO + REC_ANGINA + HF_ETIOL + 
                    REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                    DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                    REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                    rcs(REC_CLCR,3) + REC_HTN + rcs(REC_HGT_CM,3) + 
                    rcs(REC_WGT_KG,3)+REC_RACE_WHITE,
                  data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#hmr_imp[[1]]["REC_DR_MM_EQUIV_TX",]

#REC_DRUG_INDUCTION
INDUC_MISS_REG = lrm(is.na(REC_DRUG_INDUCTION) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                       rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                       CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                       REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                       rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                       REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                       REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                       rcs(REC_CLCR,3) + REC_HTN + rcs(REC_HGT_CM,3) + 
                       rcs(REC_WGT_KG,3)+REC_RACE_WHITE,
                     data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#hmr_imp[[1]]["REC_DRUG_INDUCTION",]


#REC_HGT_CM
HGT_MISS_REG = lrm(is.na(REC_HGT_CM) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                     rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                     CAN_GENDER  + REC_MALIG + REC_ECMO + REC_IABP + 
                     REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                     rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                     rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                     REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                     REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                     DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                     REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                     rcs(REC_CLCR,3) + REC_HTN + REC_DRUG_INDUCTION + 
                     rcs(REC_WGT_KG,3),
                   data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#removed bmi in hgt_miss_reg1
HGT_MISS_REG1 = lrm(is.na(REC_HGT_CM) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                      rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                      CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                      REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                      rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                      rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                      REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                      REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                      DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                      REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                      rcs(REC_CLCR,3) + REC_HTN + REC_DRUG_INDUCTION + 
                      rcs(REC_WGT_KG,3),
                    data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#hmr_imp[[1]]["REC_HGT_CM",]

WGT_MISS_REG = lrm(is.na(REC_WGT_KG) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                     rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                     CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                     REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                     rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                     rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                     REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                     REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                     DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                     REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD  
                   + REC_HTN + REC_DRUG_INDUCTION + 
                     rcs(REC_HGT_CM,3)+REC_RACE_WHITE,
                   data=data_miss, x=TRUE, y=TRUE, maxit=1000)



#REC_CREAT
CREAT_MISS_REG = lrm(is.na(REC_CREAT) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                       rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                       CAN_GENDER + REC_MALIG + REC_ECMO + REC_IABP + 
                       REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                       rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                       REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                       rcs(REC_WGT_KG,3) + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                       REC_HTN + REC_DRUG_INDUCTION + 
                       rcs(REC_HGT_CM,3)+REC_RACE_WHITE,
                     data=data_miss, x=TRUE, y=TRUE, maxit=1000)

DM_MISS_REG = lrm(is.na(REC_DM) ~ OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                    rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                    CAN_GENDER  + REC_MALIG + REC_ECMO + REC_IABP + 
                    REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                    rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                    rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + TX_YEAR + 
                    REC_DR_MM_EQUIV_TX + REC_ABO + REC_ANGINA + HF_ETIOL + 
                    rcs(REC_WGT_KG,3) + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                    DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                    REC_HBV_ANTIBODY_POS + REC_HIST_CIGARETTE + REC_VAD + 
                    rcs(REC_CLCR,3) + REC_HTN + REC_DRUG_INDUCTION + 
                    rcs(REC_HGT_CM,3)+REC_RACE_WHITE,
                  data=data_miss, x=TRUE, y=TRUE, maxit=1000)



temp_vars=c("OUT_MALIG", "DON_AGE", "DON_HGT_CM", "DON_WGT_KG", "REC_HTN",
            "DON_CREAT", "CAN_GENDER", "REC_BMI", "REC_MALIG", "REC_IABP", "REC_INOTROP",
            "REC_VENTILATOR", "REC_CARDIAC_OUTPUT", "REC_TOT_BILI", "REC_HR_ISCH", 
            "REC_AGE_AT_TX", "CAN_TOT_ALBUMIN", "TX_YEAR", "REC_DRUG_INDUCTION", "REC_ABO", 
            "REC_ANGINA", "HF_ETIOL", "REC_COPD", "REC_PVD", "CAN_SCD", 
            "DON_ABO1", "DON_ANTI_CMV_POS", "REC_CMV_IGG_POS", "REC_EBV_STAT_POS",
            "REC_HBV_ANTIBODY_POS","REC_VAD", "REC_CLCR", "REC_HIST_CIGARETTE",
            "REC_RACE_WHITE","REC_DR_MM_EQUIV_TX","REC_HGT_CM","REC_WGT_KG","REC_CREAT","REC_DM","Out_Death_1")
temp_dat=data_miss[,temp_vars]
temp_imp=puco(temp_dat)
# View(temp_imp[[1]])

###Out_Death_1
#REC_HIST_CIGGARETTE
OUT_DEATH_miss = lrm(is.na(Out_Death_1) ~ REC_HIST_CIGARETTE+OUT_MALIG + Out_Death_1 + rcs(DON_AGE,3) + 
                       rcs(DON_HGT_CM,3) + rcs(DON_WGT_KG,3) + rcs(DON_CREAT,3) + 
                       CAN_GENDER + rcs(REC_BMI,3) + REC_MALIG + REC_ECMO + REC_IABP + 
                       REC_INOTROP + REC_VENTILATOR + rcs(REC_CARDIAC_OUTPUT,3) + 
                       rcs(REC_CREAT,3) + rcs(REC_TOT_BILI,3) + rcs(REC_HR_ISCH,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(CAN_TOT_ALBUMIN,3) + 
                       REC_DRUG_INDUCTION + REC_ABO + REC_ANGINA + HF_ETIOL + 
                       REC_DM + REC_COPD + REC_PVD + CAN_SCD + DON_ABO1 + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_DR_MM_EQUIV_TX + REC_VAD + 
                       rcs(REC_CLCR,3) + REC_HTN + rcs(REC_HGT_CM,3) + 
                       rcs(REC_WGT_KG,3)+REC_RACE_WHITE,
                     data=data_miss, x=TRUE, y=TRUE, maxit=1000)
#Error in fitter(X, Y, offset = offs, penalty.matrix = penalty.matrix,  : 
#NA/NaN/Inf in foreign function call (arg 1)
#In this case , chose auxiliary variables by including analysis model covariates 
#also used PUC. 
OUT_DEATH_miss = lrm(is.na(Out_Death_1) ~  OUT_MALIG+REC_DR_MM_EQUIV_TX+CAN_GENDER+rcs(REC_HGT_CM,3)
                     +rcs(REC_CREAT,3)+rcs(REC_AGE_AT_TX,3)+TX_YEAR+REC_DRUG_INDUCTION+REC_HTN
                     +REC_HIST_CIGARETTE+REC_RACE_WHITE+rcs(REC_WGT_KG,3)+DON_ANTI_CMV_POS+
                       REC_CMV_IGG_POS+HF_ETIOL+REC_EBV_STAT_POS+REC_HBV_ANTIBODY_POS+REC_MALIG+REC_DM,
                     data=data_miss, x=TRUE, y=TRUE, maxit=1000)


#Calculate the proportion of cases (rows) that contain missing or complete values.
naniar::prop_miss_case(data_model)

#using vars_miss to impute (contains all analysis model's covariates and imputation model's
#covariates)

#note vars_miss now includes REC_CLCR_GEN
vars_miss<-c("OUT_MALIG", "Out_Death_1", "DON_AGE", "DON_HGT_CM", "DON_WGT_KG", "REC_HTN",
             "DON_CREAT", "CAN_GENDER", "REC_BMI", "REC_MALIG", "REC_ECMO", "REC_IABP", "REC_INOTROP",
             "REC_VENTILATOR", "REC_CARDIAC_OUTPUT", "REC_CREAT", "REC_TOT_BILI", "REC_HR_ISCH", 
             "REC_AGE_AT_TX", "CAN_TOT_ALBUMIN", "TX_YEAR", "REC_DRUG_INDUCTION", "REC_ABO", 
             "REC_ANGINA", "HF_ETIOL", "REC_DM", "REC_COPD", "REC_PVD", "CAN_SCD", 
             "DON_ABO1", "DON_ANTI_CMV_POS", "REC_CMV_IGG_POS", "REC_EBV_STAT_POS",
             "REC_HBV_ANTIBODY_POS", "REC_DR_MM_EQUIV_TX", "REC_VAD", "REC_CLCR", "REC_HIST_CIGARETTE",
             "REC_HGT_CM", "REC_WGT_KG", "REC_RACE_WHITE","REC_CLCR_GEN")
data_mi= data_base[,vars_miss]
data_mi <- transform(OUT_MALIG=as.factor(OUT_MALIG), Out_Death_1=as.factor(Out_Death_1),
                     CAN_GENDER=as.factor(CAN_GENDER), REC_MALIG=as.factor(REC_MALIG), 
                     REC_ECMO=as.factor(REC_ECMO), REC_IABP=as.factor(REC_IABP), 
                     REC_INOTROP=as.factor(REC_INOTROP), REC_VENTILATOR=as.factor(REC_VENTILATOR),
                     REC_DRUG_INDUCTION=as.factor(REC_DRUG_INDUCTION), REC_ABO=as.factor(REC_ABO),
                     REC_ANGINA=as.factor(REC_ANGINA), HF_ETIOL=as.factor(HF_ETIOL),
                     REC_DM=as.factor(REC_DM), REC_COPD=as.factor(REC_COPD), REC_PVD=as.factor(REC_PVD),
                     CAN_SCD=as.factor(CAN_SCD), DON_ABO1=as.factor(DON_ABO1), 
                     DON_ANTI_CMV_POS=as.factor(DON_ANTI_CMV_POS), REC_CMV_IGG_POS=as.factor(REC_CMV_IGG_POS),
                     REC_DR_MM_EQUIV_TX=as.factor(REC_DR_MM_EQUIV_TX), REC_VAD=as.factor(REC_VAD),
                     REC_HIST_CIGARETTE=as.factor(REC_HIST_CIGARETTE),REC_HTN=as.factor(REC_HTN),
                     REC_CLCR_GEN=as.factor(REC_CLCR_GEN),
                     data_mi)


####################################################################################### 
# Run Multiple Imputation using 'MICE'package
####################################################################################### 

#Need to run imp once to then modify predictor matrix for imputation
imp<-mice(data_mi,m=2,meth='pmm',print=F) 
meth=imp$method
#meth["REC_BMI"]<- "~I(REC_WGT_KG/((REC_HGT_CM*0.01)^2))"
meth["REC_CLCR"]<-"~I((((140 - REC_AGE_AT_TX)*REC_WGT_KG)/(72*REC_CREAT)) *REC_CLCR_GEN)"
#write.csv(imp$predictorMatrix,"pred1.csv")
#manually adjust prediction matrix 
preds <- as.matrix(read.csv("", row.names=1)) #edit me : predictormatrix.csv 
pas.imp<-mice(data_mi,m=20,maxit=20,meth=meth,pred=preds,print=F,seed=2) 
#Rhat.mice: to investigate Rhat
miceadds::Rhat.mice(pas.imp)



#mid_filter: which observations had missing observations 
mid_filter=which(is.na(data_base$OUT_MALIG))
get_mid_dats<-function(imp_dats,imps)
{
  mid_dats=list();
  for(i in 1:imps)
  {
    complete(imp_dats,"long",include=T,i)
    dat1=complete(imp_dats,i)
    #For each complete data set, remove outcomes that are imputed 
    mid_dats[[i]]=dat1[-c(mid_filter),]
  }
  return(mid_dats)
}

####################################################################################### 
# Create model using imputed missing data
####################################################################################### 

#Model Fit 
f <- fit.mult.impute(OUT_MALIG ~ CAN_GENDER + rcs(REC_HGT_CM,3) + rcs(REC_CREAT,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(TX_YEAR,3) + REC_DRUG_INDUCTION + 
                       REC_HTN + REC_HIST_CIGARETTE + REC_RACE_WHITE + rcs(REC_WGT_KG,3) + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + HF_ETIOL + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_MALIG + REC_DM + REC_DR_MM_EQUIV_TX, 
                     fitter=lrm, xtrans=pas.imp, data=data_model, pr=FALSE, x = TRUE, y = TRUE)

#shows which subscripts are being tested

print(anova(f),'subscripts')
print(anova(f), table.env=TRUE) 
plot(anova(f))
# Metrics
print(f, coef = FALSE)

#Receive complete data sets from imputation
complete_dats=mice::complete(pas.imp,"long",include=T)
#omits observations that had missing outcome
complete_mid_dats= complete_dats %>% dplyr::filter(!.id %in% mid_filter)
#as.mids: convert to MIDS object to use in fmid for fit.mult.impute 
pas.imp.mid = as.mids(complete_mid_dats)

#Model fit with 
fmid=fit.mult.impute(OUT_MALIG ~ CAN_GENDER + rcs(REC_HGT_CM,3) + rcs(REC_CREAT,3) + 
                       rcs(REC_AGE_AT_TX,3) + rcs(TX_YEAR,3) + REC_DRUG_INDUCTION + 
                       REC_HTN + REC_HIST_CIGARETTE + REC_RACE_WHITE + rcs(REC_WGT_KG,3) + 
                       DON_ANTI_CMV_POS + REC_CMV_IGG_POS + HF_ETIOL + REC_EBV_STAT_POS + 
                       REC_HBV_ANTIBODY_POS + REC_MALIG + REC_DM + REC_DR_MM_EQUIV_TX, 
                     fitter=lrm, xtrans=pas.imp.mid, data=data_model, pr=FALSE, x = TRUE, y = TRUE)

print(anova(fmid),'subscripts')
print(anova(fmid), table.env=TRUE) 
plot(anova(fmid))
# Metrics
print(fmid, coef = FALSE)

#### Function: get_c_CI ############
# reps=number of replications
# B=number of bootstrap samples 
# n=number of observations in data
# example: c_stat_f=get_c_CI(reps=1,B=40,n=nrow(data_base))
# Warning: This is computationally intensive 
####################################

get_c_CI<-function(reps,B,n,mod_f)
{
  C = c()
  for (i in 1:reps)
  {
    g <- update(mod_f, subset = sample(1:n, n, replace = TRUE))
    v <- validate(g, B = B)
    C[i] <- v['Dxy', 'index.corrected'] / 2 + .5
  }
  return(quantile(C, c(.025, .975)))
}
#CI_f= get_c_CI(reps=500,B=500,n=nrow(data_base),mod_f=f)
#CI_mid=get_c_CI(reps=500,B=500,n=nrow(data_base),mod_f=fmid)

# Compute the C-statistic (a.k.a. area under the ROC curve)
(c_opt_corr <- 0.5 * (val_1[1, 5] + 1))
(c_opt_corr <- 0.5 * (valmid[1, 5] + 1))

####################################################################################### 
# Create calibration plots using 500 bootstrapped samples
####################################################################################### 

# Calibration Plots
cal_1 <- calibrate(f, B = 500)
plot(cal_1)
cal_1


# jpeg("calibration_plot.jpeg",res=300,height=6, width=7 , unit="in")
plot(cal_1,subtitles = T,xlab="Predicted Probability",legend=F)
legend(.12,0.1,legend=c("Apparent","Bias-Corrected","Ideal"),lty=c(3,1,2),box.lty=0)
# dev.off()


calmid=calibrate(fmid, B = 500)
plot(calmid,subtitles = T, xlab="Predicted Probability")


####################################################################################### 
# Run model validation using 500 bootstrapped samples
####################################################################################### 

# Validation
val_1 <- validate(f, B = 500)
val_1
c_index_optimisim_corrected= val_1['Dxy', 'index.corrected'] / 2 + .5
valmid<-validate(fmid, B = 500)
valmid



####################################################################################### 
# Create nomogram of model
####################################################################################### 

#### 
#png(filename="nomogram.png",height=6,width=8,units="in",res=300)
labels_f<-rms::Newlabels(f,c(CAN_GENDER="Male Sex",
                        REC_HGT_CM="Recipient Height (cm)",
                        REC_CREAT="Recipient SCr (mg/dL)",
                        REC_AGE_AT_TX="Recipient age (years)",
                        TX_YEAR="Transplant year",
                        REC_DRUG_INDUCTION="Induction therapy",
                        REC_HTN="Recipient hypertension",
                        REC_HIST_CIGARETTE="Recipient smoking history",
                        REC_RACE_WHITE="Recipient White race",
                        REC_WGT_KG="Recipient weight (kg)",
                        DON_ANTI_CMV_POS="Donor CMV positive",
                        REC_CMV_IGG_POS="Recipient CMV Positive",
                        HF_ETIOL="HF Etiology",
                        REC_EBV_STAT_POS="Recipient EBV positive",
                        REC_HBV_ANTIBODY_POS="Recipient HBV positive",
                        REC_MALIG="Recipient malignancy history",
                        REC_DM="Recipient diabetes",
                        REC_DR_MM_EQUIV_TX="Num of DR mismatches"))
# jpeg("nomogram.jpeg",res=300,height=8, width=7 , unit="in")
plot(nomogram(labels_f, fun=function(x)1/(1+exp(-x)),  # or fun=plogis
              fun.at=c(.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999),
              funlabel="Predicted Probability"),col.conf=c('red'),
     label.every=1,col.grid = gray(c(0.8, 0.95)),cex.axis=.5,lmgp=.1,cex.var=.65)
#dev.off()

####### End of script 












