library(tidyverse)
library(janitor)
library(lubridate)
library(mice)
library(corrplot)
library(psych)

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
  mutate(peri_RBC = intra_packed_cells + rbc_0_24hrs) %>%
  mutate(peri_plasma = intra_fresh_frozen_plasma + ffp_0_24hrs) %>%
  mutate(peri_platelets = intra_platelets + plt_0_24hrs) %>%
  mutate(peri_cryoprecipitate = intra_cryoprecipitate + cryo_0_24hrs) %>%
  mutate(las_status = case_when(las_score < 30 ~ "20-29",
                                las_score >= 30 & las_score < 40 ~ "30-39",
                                las_score >= 40 ~ "40+")) %>%
  mutate(gender = case_when(gender_male == T ~ "male",
                            gender_male == F ~ "female")) %>%
  mutate(transplant_type = case_when(type == "Single Left Lung" ~ "Single",
                                     type == "Single Right Lung" ~ "Single",
                                     type == "Bilateral" ~ "Double")) %>%
  mutate(repeat_status = case_when(first_lung_transplant == T ~ "First",
                                   redo_lung_transplant == T ~ "Redo")) %>%
  mutate(intra_ecls_type = case_when(ecls_ecmo == T & ecls_cpb == F ~ "ECMO",
                                     ecls_ecmo == F & ecls_cpb == T ~ "CPB",
                                     ecls_ecmo == F & ecls_cpb == F ~ "None")) %>% 
  mutate(death_count = death_date - or_date) %>% 
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
  select(peri_RBC, peri_plasma, peri_platelets, peri_cryoprecipitate, gender, 
         height, weight, age, bmi, transplant_reason, comorbidity_score, pre_hb,
         pre_hct, pre_platelets, pre_pt, pre_inr, pre_ptt, las_status, transplant_type, 
         repeat_status, exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type)


q2_data <- new_vars %>%
  select(peri_RBC, peri_plasma, peri_platelets, peri_cryoprecipitate, gender, 
         height, weight, age, bmi, transplant_reason, comorbidity_score, pre_hb,
         pre_hct, pre_platelets, pre_pt, pre_inr, pre_ptt, las_status, transplant_type, 
         repeat_status, exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type, 
         death_days, icu_los, hospital_los, alive_30days_yn, alive_90days_yn, alive_12mths_yn)


##### Descriptive Statistics #####
summary(q2_data)
  #1 missing in pre_ptt -> impute
  #12 missing in las_score -> impute
  #187 missing in pre_fibrinogen -> impute?
  #160 missing in death_date -> haven't died yet

# Numerical variables
q2_data$death_days <- as.numeric(q2_data$death_days)
num_vars <- names(dplyr::select_if(q2_data, is.numeric))
describe(q2_data[num_vars])

# Categorical variables
cat_vars <- names(dplyr::select_if(q2_data, is.character))
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
num_only <- selected_data %>% 
  dplyr::select(where(is.numeric))

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




