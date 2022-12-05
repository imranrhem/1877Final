#### QUESTION 2: EFFECT OF TRANSFUSION ON PATIENT OUTCOMES #####
# Objectives: Survival Analysis and Patient Outcomes

library(tidyverse)
library(finalfit)
library(mice)
library(survival)
library(boot)
library(ggplot2)
library(survminer)


q2_data <- readRDS("q2_data.rds")

q2_data$death_days <- as.numeric(q2_data$death_days)

median(q2_data$peri_RBC) # 1
median(q2_data$peri_plasma) # 0
median(q2_data$peri_platelets) # 0 
median(q2_data$peri_cryoprecipitate) # 0 

q2_data <- q2_data %>%
  mutate(death_status = ifelse(!(is.na(death_days)), "Dead", "Alive")) %>%
  mutate(time = ifelse(!(is.na(death_days)), death_days, 365)) %>%
  mutate(high_RBC = ifelse(peri_RBC > median(peri_RBC), "High", "Low")) %>%
  mutate(high_plasma = ifelse(peri_plasma > median(peri_plasma), "High", "Low")) %>%
  mutate(high_platelets = ifelse(peri_platelets > median(peri_platelets), "High", "low")) %>%
  mutate(high_cryo = ifelse(peri_cryoprecipitate > median(peri_cryoprecipitate), "High", "Low"))

cat_vars <- names(dplyr::select_if(q2_data, is.character))
q2_data[cat_vars] <- lapply(q2_data[cat_vars], factor)

### Imputation ###

# Check if MCAR or MAR
las_imputation <- q2_data %>% 
  select(c(gender:comorbidity_score, las_status, transplant_type, repeat_status, intra_ecls_type))

las_imputation %>% missing_pairs(position = "fill")
  # related to transplant_reason
  # releated to gender (more missing in male)
  # related to transplant_type and repeat_status
  # Therefore MAR

### Mortality and Cox Models ###

### Survival analysis

survival_data <- q2_data %>%
  filter(time <= 365)

par(mfrow=c(2,2))

## Massive

# Create K-M Surival Curve
survival_data$massive_transfusion <- factor(survival_data$massive_transfusion)
massive_sf  <- survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = survival_data)
print(massive_sf)
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = survival_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Massive Transfusion", col = c('red', "blue"))
plot(survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = survival_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ massive_transfusion, data = survival_data) # p = 0.07

## RBC
rbc_sf  <- survfit(Surv(time, death_status == "Dead") ~ high_RBC, data = survival_data)
print(rbc_sf)
ggsurvplot(fit = rbc_sf, 
           data = survival_data,
           ####### Format Title #######
           title = "Survival Probability: RBC Transfusion",
           font.title = c(15, "black"),
           ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_classic will give a white background with no lines on the plot
             theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")), 
           ####### Format Axes #######
           xlab="Days From Surgery", # changes xlabel,
           ylab = "Survival Probability",
           font.x=c(15,"bold"), # changes x axis labels
           font.y=c(15,"bold"), # changes y axis labels
           font.xtickslab=c(13,"plain"), # changes the tick label on x axis
           font.ytickslab=c(13,"plain"),
           xlim=c(0,365),
           ####### Format Curve Lines #######
           palette = c("red","blue"),
           ####### Censor Details ########
           censor = TRUE, # logical value. If TRUE, censors will be drawn,
           censor.shape="|",
           censor.size = 5,
           ####### Confidence Intervals ########
           conf.int = F, # To Remove conf intervals use "FALSE"
           conf.int.fill = "purple", # fill color to be used for confidence interval
           surv.median.line = "hv", # allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal
           ######## Format Legend #######
           legend = "none", # If you'd prefer more space for your plot, consider removing the legend
           legend.title = "All Patients",
           legend.labs = c("RBC > 1 unit","RBC <= 1 unit"), # Change the Strata Legend
           ######## Risk Table #######
           risk.table = TRUE, # Adds Risk Table
           risk.table.height = 0.25 # Adjusts the height of the risk table (default is 0.25)
)

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ high_RBC, data = survival_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "RBC Transfusion", col = c('red', "blue"))
legend("topright", legend = c(">1 Unit RBC", "<= 1 Unit RBC"), lty = 1, col = c("red", "blue"))
plot(survfit(Surv(time, death_status == "Dead") ~ high_RBC, data = survival_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ high_RBC, data = survival_data) # p = 0.3

## Plasma
plasma_sf  <- survfit(Surv(time, death_status == "Dead") ~ high_plasma, data = survival_data)
print(plasma_sf)
ggsurvplot(fit = plasma_sf, 
           data = survival_data,
           ####### Format Title #######
           title = "Survival Probability: FFP Transfusion",
           font.title = c(15, "black"),
           ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_classic will give a white background with no lines on the plot
             theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")), 
           ####### Format Axes #######
           xlab="Days From Surgery", # changes xlabel,
           ylab = "Survival Probability",
           font.x=c(15,"bold"), # changes x axis labels
           font.y=c(15,"bold"), # changes y axis labels
           font.xtickslab=c(13,"plain"), # changes the tick label on x axis
           font.ytickslab=c(13,"plain"),
           xlim=c(0,365),
           ####### Format Curve Lines #######
           palette = c("red","blue"),
           ####### Censor Details ########
           censor = TRUE, # logical value. If TRUE, censors will be drawn,
           censor.shape="|",
           censor.size = 5,
           ####### Confidence Intervals ########
           conf.int = F, # To Remove conf intervals use "FALSE"
           conf.int.fill = "purple", # fill color to be used for confidence interval
           surv.median.line = "hv", # allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal
           ######## Format Legend #######
           legend = "none", # If you'd prefer more space for your plot, consider removing the legend
           legend.title = "All Patients",
           legend.labs = c("FFP > 0 units","FFP <= 0 units"), # Change the Strata Legend
           ######## Risk Table #######
           risk.table = TRUE, # Adds Risk Table
           risk.table.height = 0.25 # Adjusts the height of the risk table (default is 0.25)
)

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ high_RBC, data = survival_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "RBC Transfusion", col = c('red', "blue"))
legend("topright", legend = c(">1 Unit FFP", "<=0 Unit FFP"), lty = 1, col = c("red", "blue"))
plot(survfit(Surv(time, death_status == "Dead") ~ high_plasma, data = survival_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ high_plasma, data = survival_data)  # p = 0.5

## Platelets
platelets_sf  <- survfit(Surv(time, death_status == "Dead") ~ high_platelets, data = survival_data)
print(platelets_sf)
ggsurvplot(fit = platelets_sf, 
           data = survival_data,
           ####### Format Title #######
           title = "Survival Probability: Platelet Transfusion",
           font.title = c(15, "black"),
           ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_classic will give a white background with no lines on the plot
             theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")), 
           ####### Format Axes #######
           xlab="Days From Surgery", # changes xlabel,
           ylab = "Survival Probability",
           font.x=c(15,"bold"), # changes x axis labels
           font.y=c(15,"bold"), # changes y axis labels
           font.xtickslab=c(13,"plain"), # changes the tick label on x axis
           font.ytickslab=c(13,"plain"),
           xlim=c(0,365),
           ####### Format Curve Lines #######
           palette = c("red","blue"),
           ####### Censor Details ########
           censor = TRUE, # logical value. If TRUE, censors will be drawn,
           censor.shape="|",
           censor.size = 5,
           ####### Confidence Intervals ########
           conf.int = F, # To Remove conf intervals use "FALSE"
           conf.int.fill = "purple", # fill color to be used for confidence interval
           surv.median.line = "hv", # allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal
           ######## Format Legend #######
           legend = "none", # If you'd prefer more space for your plot, consider removing the legend
           legend.title = "All Patients",
           legend.labs = c("Platelets > 0 units","Platelets <= 0 units"), # Change the Strata Legend
           ######## Risk Table #######
           risk.table = TRUE, # Adds Risk Table
           risk.table.height = 0.25 # Adjusts the height of the risk table (default is 0.25)
)
title("platelets Transfusion")

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ high_platelets, data = survival_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Platelet Transfusion", col = c('red', "blue"))
legend("bottomright", legend = c(">1 Unit Platelets", "<=0 Unit Platelets"), lty = 1, col = c("red", "blue"))

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ high_platelets, data = survival_data)  # p = 0.9

## Cryo
cryo_sf  <- survfit(Surv(time, death_status == "Dead") ~ high_cryo, data = survival_data)
print(cryo_sf)
ggsurvplot(fit = platelets_sf, 
           data = survival_data,
           ####### Format Title #######
           title = "Survival Probability: Cryoprecipitate Transfusion",
           font.title = c(15, "black"),
           ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_classic will give a white background with no lines on the plot
             theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")), 
           ####### Format Axes #######
           xlab="Days From Surgery", # changes xlabel,
           ylab = "Survival Probability",
           font.x=c(15,"bold"), # changes x axis labels
           font.y=c(15,"bold"), # changes y axis labels
           font.xtickslab=c(13,"plain"), # changes the tick label on x axis
           font.ytickslab=c(13,"plain"),
           xlim=c(0,365),
           axes.offset = F,
           ####### Format Curve Lines #######
           palette = c("red","blue"),
           ####### Censor Details ########
           censor = TRUE, # logical value. If TRUE, censors will be drawn,
           censor.shape="|",
           censor.size = 5,
           ####### Confidence Intervals ########
           conf.int = F, # To Remove conf intervals use "FALSE"
           conf.int.fill = "purple", # fill color to be used for confidence interval
           surv.median.line = "hv", # allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal
           ######## Format Legend #######
           legend = "none", # If you'd prefer more space for your plot, consider removing the legend
           legend.title = "All Patients",
           legend.labs = c("Cryoprecipitate > 0 units","Cryoprecipitate <= 0 units"), # Change the Strata Legend
           ######## Risk Table #######
           risk.table = TRUE, # Adds Risk Table
           risk.table.height = 0.25 # Adjusts the height of the risk table (default is 0.25)
)

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ high_cryo, data = survival_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Cryoprecipitate Transfusion", col = c('red', "blue"))
legend("bottomright", legend = c(">1 Unit Cryoprecipitate", "<=0 Unit Cryoprecipitate"), lty = 1, col = c("red", "blue"))
plot(survfit(Surv(time, death_status == "Dead") ~ high_cryo, data = survival_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ high_cryo, data = survival_data)  # p = 0.2

### Proportional Cox Models

## Checking distributions of variables

# Numerical Variables
num_vars <- names(dplyr::select_if(survival_data, is.numeric))
num_vars
hist(survival_data$peri_RBC) # right skew
hist(log(survival_data$peri_RBC))
hist(survival_data$peri_plasma) # right skew
hist(survival_data$peri_platelets) # right skew
hist(survival_data$peri_cryoprecipitate) # right skew
hist(survival_data$height) # normal
hist(survival_data$weight) # normal
hist(survival_data$bmi) # bimodal
hist(survival_data$age) # left skew
hist(survival_data$comorbidity_score) # right skew

## High Cox
high_cox <- coxph(Surv(time, death_status == "Dead") ~ high_RBC + high_cryo + high_platelets + high_plasma + bmi + age + comorbidity_score + gender + transplant_reason + las_status + transplant_type + repeat_status, data = survival_data)
summary(high_cox)
cox.zph(high_cox)

ggforest(high_cox)

## Peri Cox
peri_cox  <- coxph(Surv(time, death_status == "Dead") ~ peri_RBC + peri_cryoprecipitate + peri_platelets + peri_plasma + bmi + age + comorbidity_score + gender + transplant_reason + las_status + transplant_type + repeat_status, data = survival_data)
summary(peri_cox)
cox.zph(peri_cox)

### Non-Mortality Patient Outcomes ###

# Hospital LOS
boxplot(q2_data$hospital_los ~ q2_data$high_RBC)
boxplot(q2_data$hospital_los ~ q2_data$high_plasma)
boxplot(q2_data$hospital_los ~ q2_data$high_platelets)
boxplot(q2_data$hospital_los ~ q2_data$high_cryo)

# ICU LOS 
boxplot(q2_data$icu_los ~ q2_data$high_RBC)
boxplot(q2_data$icu_los ~ q2_data$high_plasma)
boxplot(q2_data$icu_los ~ q2_data$high_platelets)
boxplot(q2_data$icu_los ~ q2_data$high_cryo)


