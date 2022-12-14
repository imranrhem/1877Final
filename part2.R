library(tidyverse)
library(janitor)
library(lubridate)
library(corrplot)
library(psych)
library(dplyr)
library(MASS)

# Import data
transfusion_data <- read.csv("transfusion_data.csv", stringsAsFactors = T, na.strings = c("", "?", "#VALUE!"))

# Clean up column names
names(transfusion_data) <- tolower(names(transfusion_data))
transfusion_data <- clean_names(transfusion_data)

# Convert to date format
transfusion_data$or_date <- ymd(transfusion_data$or_date)
transfusion_data$death_date <- dmy(transfusion_data$death_date)

# Clean up data set
selected_data <- transfusion_data %>%
  # select variables of interest
  select(c(or_date:redo_lung_transplant, exvivo_lung_perfusion:ecls_cpb, intra_fresh_frozen_plasma:intra_cryoprecipitate,duration_of_ventilation:hospital_los, rbc_0_24hrs:cryo_72hr_total, total_24hr_rbc:massive_transfusion)) %>%
  select(-pre_fibrinogen) %>%
  # replace NA with 0s where appropriate (transfusion variables)
  mutate_at(vars(rbc_0_24hrs:cryo_72hr_total), ~replace_na(., 0))

## Create new variables that we need and remove old ones
new_vars <- selected_data %>%
  # Create peri variables
  mutate(peri_RBC = intra_packed_cells + rbc_0_24hrs) %>%
  mutate(peri_plasma = intra_fresh_frozen_plasma + ffp_0_24hrs) %>%
  mutate(peri_platelets = intra_platelets + plt_0_24hrs) %>%
  mutate(peri_cryoprecipitate = intra_cryoprecipitate + cryo_0_24hrs) %>%
  mutate(any_transfusion = ifelse(peri_RBC > 0 | peri_plasma > 0 | peri_platelets > 0 | peri_cryoprecipitate > 0, "Yes", "No")) %>%
  mutate(las_status = case_when(las_score < 30 ~ "20-29",
                                las_score >= 30 & las_score < 40 ~ "30-39",
                                las_score >= 40 ~ "40+")) %>%
  mutate(gender = case_when(gender_male == T ~ "male",
                            gender_male == F ~ "female")) %>%
  mutate(transplant_type = case_when(type == "Single Left Lung" ~ "Single",
                                     type == "Single Right Lung" ~ "Single",
                                     type == "Bilateral" ~ "Double")) %>%
  mutate(repeat_status = case_when(first_lung_transplant == T & redo_lung_transplant == F ~ "First",
                                   first_lung_transplant == F & redo_lung_transplant == F ~ "First",
                                   first_lung_transplant == F & redo_lung_transplant == T ~ "Redo",
                                   first_lung_transplant == T & redo_lung_transplant == T ~ "Redo")) %>%
  mutate(intra_ecls_type = case_when(ecls_ecmo == T & ecls_cpb == F ~ "ECMO",
                                     ecls_ecmo == F & ecls_cpb == T ~ "CPB",
                                     ecls_ecmo == F & ecls_cpb == F ~ "None")) %>% 
  mutate(death_days = death_date - or_date) %>% 
  mutate(transplant_reason = case_when(copd == T & alpha1_antitrypsin_deficiency == F & cystic_fibrosis == F & idiopathic_pulmonary_hypertension == F & interstitial_lung_disease == F & pulm_other == F ~ "COPD",
                                       copd == T & (alpha1_antitrypsin_deficiency == T | cystic_fibrosis == T | idiopathic_pulmonary_hypertension == T | interstitial_lung_disease == T | pulm_other == T) ~ "Multiple",
                                       copd == F & alpha1_antitrypsin_deficiency == T & cystic_fibrosis == F & idiopathic_pulmonary_hypertension == F & interstitial_lung_disease == F & pulm_other == F ~ "A1AD",
                                       alpha1_antitrypsin_deficiency == T & (copd == T | cystic_fibrosis == T | idiopathic_pulmonary_hypertension == T | interstitial_lung_disease == T | pulm_other == T) ~ "Multiple",
                                       copd == F & alpha1_antitrypsin_deficiency == F & cystic_fibrosis == T & idiopathic_pulmonary_hypertension == F & interstitial_lung_disease == F & pulm_other == F ~ "CF",
                                       cystic_fibrosis == T & (copd == T | alpha1_antitrypsin_deficiency == T | idiopathic_pulmonary_hypertension == T | interstitial_lung_disease == T | pulm_other == T) ~ "Multiple",
                                       copd == F & alpha1_antitrypsin_deficiency == F & cystic_fibrosis == F & idiopathic_pulmonary_hypertension == T & interstitial_lung_disease == F & pulm_other == F ~ "Pulmonary Hypertension",
                                       idiopathic_pulmonary_hypertension == T & (copd == T | alpha1_antitrypsin_deficiency == T | cystic_fibrosis == T | interstitial_lung_disease == T | pulm_other == T) ~ "Multiple",
                                       copd == F & alpha1_antitrypsin_deficiency == F & cystic_fibrosis == F & idiopathic_pulmonary_hypertension == F & interstitial_lung_disease == T & pulm_other == F ~ "Interstitial Lung Disease",
                                       interstitial_lung_disease == T & (copd == T | alpha1_antitrypsin_deficiency == T | cystic_fibrosis == T | idiopathic_pulmonary_hypertension == T | pulm_other == T) ~ "Multiple",
                                       copd == F & alpha1_antitrypsin_deficiency == F & cystic_fibrosis == F & idiopathic_pulmonary_hypertension == F & interstitial_lung_disease == F & pulm_other == T ~ "Other",
                                       pulm_other == T & (copd == T | alpha1_antitrypsin_deficiency == T | cystic_fibrosis == T | idiopathic_pulmonary_hypertension == T | interstitial_lung_disease == T) ~ "Multiple")) %>%
  mutate(comorbidity_score = coronary_artery_disease + hypertension + diabetes_insulin + diabetes_diet_ohgs + gerd_pud + renal_failure + stroke_cva + liver_disease + thyroid_disease)

new_vars$preoperative_ecls[new_vars$preoperative_ecls == T] <- "Yes"
new_vars$preoperative_ecls[new_vars$preoperative_ecls == F] <- "No"

new_vars$preoperative_ecls[new_vars$preoperative_ecls == T] <- "Yes"
new_vars$preoperative_ecls[new_vars$preoperative_ecls == F] <- "No"

# Create datasets
q1_data <- new_vars %>%
  select(peri_RBC, peri_plasma, peri_platelets, peri_cryoprecipitate, massive_transfusion, any_transfusion, gender, 
         height, weight, age, bmi, transplant_reason, comorbidity_score, pre_hb,
         pre_hct, pre_platelets, pre_pt, pre_inr, pre_ptt, las_status, transplant_type, 
         repeat_status, exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type)

q1_data <- q1_data[-(193:200),]

# removing NA values using complete case 
q1_dataclean <- q1_data[complete.cases(q1_data),]


q2_data <- new_vars %>%
  select(peri_RBC, peri_plasma, peri_platelets, peri_cryoprecipitate, massive_transfusion, any_transfusion, gender, 
         height, weight, age, bmi, transplant_reason, comorbidity_score, pre_hb,
         pre_hct, pre_platelets, pre_pt, pre_inr, pre_ptt, las_status, transplant_type, 
         repeat_status, exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type, 
         death_days, icu_los, hospital_los, alive_30days_yn, alive_90days_yn, alive_12mths_yn)

q2_data <- q2_data[-(193:200),]

test <- q1_data %>%
  filter(is.na(las_status))

##### Descriptive Statistics #####
summary(q2_data)
  #1 missing in pre_ptt -> impute, cannot classify as MNAR, treat as MAR?
  #12 missing in las_score -> impute, cannot classify as MNAR, treat as MAR?
  #160 missing in death_date -> haven't died yet

# Numerical variables
q2_data$death_days <- as.numeric(q2_data$death_days)
num_vars <- names(dplyr::select_if(q2_data, is.numeric))
describe(q2_data[num_vars])

# Categorical variables
cat_vars <- names(dplyr::select_if(q2_data, is.character))
lapply(q2_data[cat_vars],table)

factored_vars <- names(dplyr::select_if(q2_data, is.factor))
lapply(q2_data[cat_vars],table)

# Categorical variables

## Save as RDS

saveRDS(q1_data, "q1_data.rds")
saveRDS(q2_data, "q2_data.rds")

######## TEST CODE ######## 

### Pre-imputation EDA

# Check missingness and distributions
summary(selected_data)
  #1 missing in pre_ptt -> impute
  #12 missing in las_score -> impute
  #187 missing in pre_fibrinogen -> impute?
  #160 missing in death_date -> haven't died yet

# Check correlations
num_only <- q2_data %>% 
  dplyr::select(where(is.numeric))

cors <- cor(na.omit(num_only))

corrplot(cors)
  # height-weight are correlated
  # weight and BMI are correlated
  # pre_hb-rbc_0_24
  # pre_hb-pre_inr
  # all intras are corerlated
  # all intras and total_24rbc
  # all ttransfusion varibales are corelated for the most part
  # no correlations with LAS score that re of issue, assume no colinearity then

### Imputation

# Check pattern of missingness
md.pattern(selected_data, rotate.names = T)
  # Variables to include in imputation: variables of the intended model + variables related to missingness
  #12 observations where LAS is missing, 10 missing with death_date, no meaningful relationship
    # LAS is not MAR

# Imputation of LAS
  # factors used to calculate LAS: age, ht, wt, diangosis, functional status, ventilation, ECLS


### Multivariate Linear Regression
model1 <- lm(peri_RBC ~ 
               gender+ height+ weight+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_hct+ pre_platelets+ pre_pt+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model2 <- lm(peri_plasma ~ 
               gender+ height+ weight+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_hct+ pre_platelets+ pre_pt+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model3 <- lm(peri_platelets ~ 
               gender+ height+ weight+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_hct+ pre_platelets+ pre_pt+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model4 <- lm(peri_cryoprecipitate ~ 
               gender+ height+ weight+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_hct+ pre_platelets+ pre_pt+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

# check variance inflation factor (VIF) to check for collinearity, if vif > 5 model will be adjusted
# note: treating GVIF values as VIF for variables with df = 1, if df > 1, using the square of the scaled GVIF value
library(car)
vif(model1)
vif(model2) 
vif(model3)  
vif(model4) 

# adjusted models based on VIF results - excluding height, weight, pre_hct, and pre_pt
model1 <- lm(peri_RBC ~ 
               gender+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_platelets+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model2 <- lm(peri_plasma ~ 
               gender+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_platelets+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model3 <- lm(peri_platelets ~ 
               gender+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_platelets+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

model4 <- lm(peri_cryoprecipitate ~ 
               gender+ age+ bmi+ transplant_reason+ comorbidity_score+ pre_hb+
               pre_platelets+ pre_inr+ pre_ptt+ las_status+ transplant_type+ 
               repeat_status+ exvivo_lung_perfusion+ preoperative_ecls+ intra_ecls_type, 
             data = q1_dataclean)

# backwards stepwise variable selection
model1.back <- stepAIC(model1) # AIC 454.68, 8/15 predictors
model2.back <- stepAIC(model2) # AIC 317.25, 5/15 predictors
model3.back <- stepAIC(model3) # AIC 90.68, 8/15 predictors
model4.back <- stepAIC(model4) # AIC 687.63, 7/15 predictors
# old result for model 4 -- recheck? AIC 88.55, 11/15 predictors

# coefficients for best model following results of stepwise selection
summary(model1.back)
summary(model2.back)
summary(model3.back)
summary(model4.back)


### testing linear regression assumptions on best stepwise models
par(mfrow = c(2, 2))
plot(model1.back)

par(mfrow = c(2, 2))
plot(model2.back)

par(mfrow = c(2, 2))
plot(model3.back)

par(mfrow = c(2, 2))
plot(model4.back)

# scaling response variable of best stepwise models to attempt to improve residual plots
model1.scale <- lm(sqrt(peri_RBC) ~ gender + pre_hb + pre_inr + pre_ptt + repeat_status + 
                     exvivo_lung_perfusion + preoperative_ecls + intra_ecls_type, 
                   data = q1_dataclean)
par(mfrow = c(2, 2))
plot(model1.scale)

model2.scale <- lm(sqrt(peri_plasma) ~ gender + pre_ptt + repeat_status + 
                     exvivo_lung_perfusion + preoperative_ecls, 
                   data = q1_dataclean)
par(mfrow = c(2, 2))
plot(model2.scale)

model3.scale <- lm(sqrt(peri_platelets) ~ gender + age + pre_platelets + pre_inr + repeat_status + 
                          exvivo_lung_perfusion + preoperative_ecls + intra_ecls_type, 
                   data = q1_dataclean)
par(mfrow = c(2, 2))
plot(model3.scale)

model4.scale <- lm(sqrt(peri_cryoprecipitate) ~ gender + age + pre_platelets + repeat_status + 
                          exvivo_lung_perfusion + preoperative_ecls + intra_ecls_type, 
                   data = q1_dataclean)
par(mfrow = c(2, 2))
plot(model4.scale)




