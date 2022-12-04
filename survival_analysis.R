#### QUESTION 2 #####
# Objectives: Survival Analysis and Patient Outcomes

library(tidyverse)
library(survival)
library(boot)

q2_data <- readRDS("q2_data.rds")

q2_data$death_days <- as.numeric(q2_data$death_days)

median(q2_data$peri_RBC) # 1

median(q2_data$peri_plasma) # 0
median(q2_data$peri_platelets) # 0 
median(q2_data$peri_cryoprecipitate) # 0 

q2_data <- q2_data %>%
  mutate(death_status = ifelse(!(is.na(death_days)), "Dead", "Alive")) %>%
  mutate(time = ifelse(!(is.na(death_days)), death_days, max(death_days, na.rm = T))) %>%
  mutate(high_RBC = ifelse(peri_RBC > mean(peri_RBC), "High", "Low")) %>%
  mutate(high_plasma = ifelse(peri_plasma > mean(peri_plasma), "High", "Low")) %>%
  mutate(high_platelets = ifelse(peri_platelets > mean(peri_platelets), "High", "low")) %>%
  mutate(high_cryo = ifelse(peri_cryoprecipitate > mean(peri_cryoprecipitate), "High", "Low"))

cat_vars <- names(dplyr::select_if(q2_data, is.character))
q2_data[cat_vars] <- lapply(q2_data[cat_vars], factor)

### Mortality and Cox Models ###

### Survival analysis

survival_data <- q2_data %>%
  filter(time <= 365)

par(mfrow=c(2,2))

## Massive

# Create K-M Surival Curve
q2_data$massive_transfusion <- factor(q2_data$massive_transfusion)
massive_sf  <- survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = q2_data)
print(massive_sf)
plot(survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = q2_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Massive Transfusion", col = c('red', "blue"))
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ massive_transfusion, data = q2_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ massive_transfusion, data = q2_data) # p = 0.07

## RBC

# Create K-M Surival Curve
rbc_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_RBC, data = q2_data)
print(rbc_sf)
plot(survfit(Surv(time, death_status == "Dead") ~ any_RBC, data = q2_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "RBC Transfusion", col = c('red', "blue"))
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ any_RBC, data = q2_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ any_RBC, data = q2_data) # p = 0.07

## Plasma
plasma_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_plasma, data = q2_data)
print(plasma_sf)
plot(survfit(Surv(time, death_status == "Dead") ~ any_plasma, data = q2_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "FFP Transfusion", col = c('red', "blue"))
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ any_plasma, data = q2_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ any_plasma, data = q2_data)  # p = 0.04

platelets_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data)
print(platelets_sf)
plot(platelets_sf, xlab = "Days from Surgery", ylab = "Survival", col = c("red", "blue"))
title("Platelet Transfusion"))

## Platelets
platelets_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data)
print(platelets_sf)
plot(survfit(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Platelet Transfusion", col = c('red', "blue"))
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))
title("platelets Transfusion")

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data)  # p = 0.7

platelets_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_platelets, data = q2_data)
print(platelets_sf)
plot(platelets_sf, xlab = "Days from Surgery", ylab = "Survival", col = c("red", "blue"))
title("Platelet Transfusion"))

## Cryo
cryo_sf  <- survfit(Surv(time, death_status == "Dead") ~ any_cryo, data = q2_data)
print(cryo_sf)
plot(survfit(Surv(time, death_status == "Dead") ~ any_cryo, data = q2_data), fun = "S", xlab = "Days From Transplant", ylab = "Survival", main = "Cryoprecipitate Transfusion", col = c('red', "blue"))
legend("topright", legend = c("Transfusion", "No Transfusion"), lty = 1, col = c("red", "blue"))
title("cryo Transfusion")

# Check Assumptions for log-rank test
plot(survfit(Surv(time, death_status == "Dead") ~ any_cryo, data = q2_data), fun = "cloglog")

# Log-rank test
survdiff(Surv(time, death_status == "Dead") ~ any_cryo, data = q2_data)  # p = 0.3

### Proportional Cox Models

## Checking distributions of variables

# Numerical Variables
num_vars <- names(dplyr::select_if(q2_data, is.numeric))
num_vars
hist(q2_data$peri_RBC) # right skew
hist(log(q2_data$peri_RBC))
hist(q2_data$peri_plasma) # right skew
hist(q2_data$peri_platelets) # right skew
hist(q2_data$peri_cryoprecipitate) # right skew
hist(q2_data$height) # normal
hist(q2_data$weight) # normal
hist(q2_data$bmi) # bimodal
hist(q2_data$age) # left skew
hist(q2_data$comorbidity_score) # right skew

## High Cox
high_cox <- coxph(Surv(time, death_status == "Dead") ~ high_RBC + high_cryo + high_platelets + high_plasma + bmi + age + comorbidity_score + gender + transplant_reason + las_status + transplant_type + repeat_status + intra_ecls_type, data = q2_data)
summary(high_cox)
cox.zph(high_cox)

## Peri Cox
peri_cox  <- coxph(Surv(time, death_status == "Dead") ~ peri_RBC + peri_cryoprecipitate + peri_platelets + peri_plasma + bmi + age + comorbidity_score + gender + transplant_reason + las_status + transplant_type + repeat_status + intra_ecls_type, data = q2_data)
summary(peri_cox)
cox.zph(peri_cox)

## FFP

## Platelets

## Cryoprecipitate


### Non-Mortality Patient Outcomes ###