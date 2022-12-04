#question 1 - part 1 

library(tidyverse)
library(janitor)
library(lubridate)
library(mice)
library(corrplot)
library(psych)

#Import data
setwd("~/Documents/1877Final")

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
         pre_hct, pre_platelets, pre_pt, pre_inr, pre_ptt, las_score, transplant_type, 
         repeat_status, exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type)
summary(q1_data)

#take complete observations only - NAs in the las score
q1_datacomplete <- q1_data[complete.cases(q1_data),]

#convert to categorical after imputation
q1data <- q1_datacomplete %>%
  mutate(las_status = case_when(las_score < 30 ~ "20-29",
                              las_score >= 30 & las_score < 40 ~ "30-39",
                              las_score >= 40 ~ "40+"))%>%
  mutate(las_status = as.factor(las_status))%>%
  mutate(repeat_status = as.factor(repeat_status))%>%
  mutate(transplant_type = as.factor(transplant_type))%>%
  mutate(transplant_reason = as.factor(transplant_reason))%>%
  mutate(peri_RBC = as.numeric(peri_RBC))%>%
  mutate(peri_plasma = as.numeric(peri_plasma))%>%
  mutate(peri_platelets = as.numeric(peri_platelets))%>%
  mutate(peri_cryoprecipitate = as.numeric(peri_cryoprecipitate))%>%
  mutate(gender = as.factor(gender))%>%
  mutate(preoperative_ecls = as.factor(preoperative_ecls))%>%
  mutate(intra_ecls_type = as.factor(intra_ecls_type))
  

#view row and column count of new data frame
dim(no_outliers) 

#######################################################################
######MANOVA FOR DIFFERENCES BETWEEN TRANSFUSION OR NO TRANSFUSION#####
########################################################################

resman <- manova(cbind(height,weight, age, bmi, transplant_reason, 
                       comorbidity_score, pre_hb, pre_hct, pre_platelets, 
                       pre_pt, pre_inr, pre_ptt, las_status, repeat_status, 
                       exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type) ~ any_transfusion, data = q1data)
summary.aov(resman)

#significant variables: height, weight, age, comorbidity score, any transfusion, pre_hb, 
#pre_hct, las_status, repeat_status, intra_ecls_type

#Conducting t-test for continuous variables 
t.test(height ~ any_transfusion, data=q1data)
t.test(weight ~ any_transfusion, data=q1data)
t.test(age ~ any_transfusion, data=q1data)
t.test(pre_hb ~ any_transfusion, data=q1data)
t.test(pre_hct ~any_transfusion, data=q1data)
t.test(repeat_status ~ any_transfusion, data = q1data)
t.test(intra_ecls_type ~ any_transfusion, data = q1data)

#creating individidual boxplots  
#boxplot for height 
height <- ggboxplot(q1datattest, x = "any_transfusion", y = "height",
               xlab = "Transfusion Received", ylab = "Height in cm", color = "blue", palette = "jco",
               add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
height + stat_compare_means()
# Change method
height + stat_compare_means(method = "t.test")

#boxplot for weight
weight <- ggboxplot(q1datattest, x = "any_transfusion", y = "weight",
                    xlab = "Transfusion Received", ylab = "Weight", color = "black", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
weight + stat_compare_means()
# Change method
weight + stat_compare_means(method = "t.test")

#boxplot for age
age <- ggboxplot(q1datattest, x = "any_transfusion", y = "age",
                    xlab = "Transfusion Received", ylab = "Age", color = "purple", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
age + stat_compare_means()
# Change method
age + stat_compare_means(method = "t.test")

#boxplot for hb

pre_hb <- ggboxplot(q1datattest, x = "any_transfusion", y = "pre_hb",
                               xlab = "Transfusion Received", ylab = "Hemoglobin Levles Pre-Surgery", color = "dark orange", palette = "jco",
                               add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
pre_hb + stat_compare_means()
# Change method
pre_hb + stat_compare_means(method = "t.test")

#boxplot for hct
pre_hct <- ggboxplot(q1datattest, x = "any_transfusion", y = "pre_hct",
                    xlab = "Transfusion Received", ylab = "Hematocrit Levles Pre-Surgery", color = "dark red", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
pre_hct + stat_compare_means()
# Change method
pre_hct + stat_compare_means(method = "t.test")

#Fisher Exact test for categorical variables 
#######transplant reason################# - underlying diagnosis (not found to be sig for manova but still tried)
transplant_reason <- q1data$transplant_reason
any_transfusion <- q1data$any_transfusion

transplantframe <- data.frame(transplant_reason, any_transfusion)
transplantmatrix <- as.data.frame.matrix(table(transplantframe))

fisher.test(transplantmatrix)
pairwise_fisher_test(transplantmatrix, p.adjust.method = "fdr")

######comorbidity score########
comorbidity_score <- q1data$comorbidity_score

comorbframe <- data.frame(comorbidity_score, any_transfusion)
comorbmatrix <- as.data.frame.matrix(table(comorbframe))

fisher.test(comorbmatrix)
pairwise_fisher_test(comorbmatrix, p.adjust.method = "fdr")

###########las_status fisher test######
las_status <- q1data$las_status

lasframe <- data.frame(las_status, any_transfusion)
lasmatrix <- as.data.frame.matrix(table(lasframe))

fisher.test(lasmatrix)
pairwise_fisher_test(lasmatrix, p.adjust.method = "fdr")

#######repeat_status
repeat_status <- q1data$repeat_status

repeatframe <- data.frame(repeat_status, any_transfusion)
repeatmatrix <- as.data.frame.matrix(table(repeatframe))

fisher.test(repeatmatrix)
pairwise_fisher_test(repeatmatrix, p.adjust.method = "fdr")






















