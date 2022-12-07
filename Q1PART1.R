#converting to matrix 
q1datamatrix <- data.matrix(q1data)
mshapiro.test(cbind(x,y))
mshapiro_test(q1datamatrix)
#does not meet assumption for manova, however will still move forward

#conducting the manova model 
library(car)
resman <- manova(cbind(gender, height,weight, age, bmi, transplant_reason, 
                       comorbidity_score, pre_hb, pre_hct, pre_platelets, 
                       pre_pt, pre_ptt, las_status, repeat_status, 
                       exvivo_lung_perfusion, preoperative_ecls, intra_ecls_type) ~ any_transfusion, data = q1data)
summary.aov(resman)


#significant variables: height, weight, age, comorbidity score, any transfusion, pre_hb, 
#pre_hct, las_status, repeat_status, intra_ecls_type

#Conducting t-tests and wilcoxon tests for continuous variables 
shapiro.test(q1data$height)#normally distributed
shapiro.test(q1data$weight)#normally distributed
shapiro.test(q1data$age) #not normally distributed
shapiro.test(q1data$pre_hb) #not normally distributed
shapiro.test(q1data$pre_hct) #not normally distributed 
t.test(height ~ any_transfusion, data=q1data)
t.test(weight ~ any_transfusion, data=q1data)
q1data%>%
  cohens_d(height ~any_transfusion)
q1data%>%
  cohens_d(weight ~ any_transfusion)
wilcox.test(age ~ any_transfusion, data=q1data)
q1data %>%
  group_by(any_transfusion) %>%
  get_summary_stats(age, type = "median_iqr")

wilcox_effsize(age ~ any_transfusion, data=q1data)

wilcox.test(pre_hb ~ any_transfusion, data=q1data)

q1data %>%
  group_by(any_transfusion) %>%
  get_summary_stats(pre_hb, type = "median_iqr")

wilcox_effsize(pre_hb ~ any_transfusion, data=q1data)

#wilcox test for hct 
wilcox.test(pre_hct ~any_transfusion, data=q1data, var.equal = T)
q1data %>%
  group_by(any_transfusion) %>%
  get_summary_stats(pre_hct, type = "median_iqr")

wilcox_effsize(pre_hct~ any_transfusion, data=q1data)

#creating individidual boxplots  
#boxplot for height 
library(ggpubr)
height <- ggboxplot(q1data, x = "any_transfusion", y = "height",
                    xlab = "Transfusion Received", ylab = "Height in cm", color = "blue", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
height + stat_compare_means()
# Change method
height + stat_compare_means(method = "t.test")

#boxplot for weight
weight <- ggboxplot(q1data, x = "any_transfusion", y = "weight",
                    xlab = "Transfusion Received", ylab = "Weight", color = "black", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
weight + stat_compare_means()
# Change method
weight + stat_compare_means(method = "t.test")

#boxplot for age
age <- ggboxplot(q1data, x = "any_transfusion", y = "age",
                 xlab = "Transfusion Received", ylab = "Age", color = "red", palette = "jco",
                 add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
age + stat_compare_means()
# Change method
age + stat_compare_means(method = "wilcox.test")

#boxplot for hb

pre_hb <- ggboxplot(q1data, x = "any_transfusion", y = "pre_hb",
                    xlab = "Transfusion Received", ylab = "Hemoglobin Levles Pre-Surgery", color = "dark orange", palette = "jco",
                    add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
pre_hb + stat_compare_means()
# Change method
pre_hb + stat_compare_means(method = "wilcox.test")

#boxplot for hct
pre_hct <- ggboxplot(q1data, x = "any_transfusion", y = "pre_hct",
                     xlab = "Transfusion Received", ylab = "Hematocrit Levles Pre-Surgery", color = "dark orange", palette = "jco",
                     add = "jitter", font.label = list(size = 50, face = "plain"))
#  Add p-value
pre_hct + stat_compare_means()
# Change method
pre_hct + stat_compare_means(method = "wilcox.test")

#Fisher Exact test for categorical variables with pairwise comparisons 

######gender##########
gender <- q1data$gender
any_transfusion <- q1data$any_transfusion

genderframe <- data.frame(gender, any_transfusion)
gendermatrix <- as.data.frame.matrix(table(genderframe))

fisher.test(gendermatrix[,c(2,1)])

#######transplant reason################# - underlying diagnosis (not found to be sig for manova but still tried)
library(rstatix)

transplant_reason <- q1data$transplant_reason

transplantframe <- data.frame(transplant_reason, any_transfusion)
transplantmatrix <- as.data.frame.matrix(table(transplantframe))

fisher.test(transplantmatrix)
c <- pairwise_fisher_test(transplantmatrix, p.adjust.method = "fdr")
View(c)
#odds ratio for CF and COPD
transplantCFCOPDmatrix <- transplantmatrix[c(2,3), c(2,1)]
fisher.test(transplantCFCOPDmatrix)

#odds ratio for COPD and ILD 
transplantCOPDILDmatrix <- transplantmatrix[c(3,4), c(2,1)]
fisher.test(transplantCOPDILDmatrix)


#odds ratio for COPD and other 
transplantCOPDothermatrix <- transplantmatrix[c(3,6), c(2,1)]
fisher.test(transplantCOPDothermatrix)


######comorbidity score########
comorbidity_score <- q1data$comorbidity_score

comorbframe <- data.frame(comorbidity_score, any_transfusion)
comorbmatrix <- as.data.frame.matrix(table(comorbframe))

fisher.test(comorbmatrix)
pairwise_fisher_test(comorbmatrix, p.adjust.method = "fdr")

comorb02matrix <- comorbmatrix[c(1,3),c(2,1)]
fisher.test(comorb02matrix)

###########las_status fisher test######
las_status <- q1data$las_status

lasframe <- data.frame(las_status, any_transfusion)
lasmatrix <- as.data.frame.matrix(table(lasframe))

fisher.test(lasmatrix)
pairwise_fisher_test(lasmatrix, p.adjust.method = "fdr")

las30v40matrix <- lasmatrix[c(2,3),c(2,1)]
fisher.test(las30v40matrix)

#######repeat_status
repeat_status <- q1data$repeat_status

repeatframe <- data.frame(repeat_status, any_transfusion)
repeatmatrix <- as.data.frame.matrix(table(repeatframe))

fisher.test(repeatmatrix[,c(2,1)])
pairwise_fisher_test(repeatmatrix, p.adjust.method = "fdr")
# prop.test()
# 
# no group: 67
# yes group: 107
# 
# no group with first: 67/178 #0.3764045
# yes group with first: 102/107 #0.953271
# 
# propp <- prop.test(x=c(67,102),n=c(67,107))

#######intra_ecls_type
intra_ecls_type <- q1data$intra_ecls_type

eclsframe <- data.frame(intra_ecls_type, any_transfusion)
eclsmatrix <- as.data.frame.matrix(table(eclsframe))

fisher.test(eclsmatrix)
pairwise_fisher_test(eclsmatrix, p.adjust.method = "fdr")

eclsECMOnonematrix <- lasmatrix[c(1,3),c(2,1)]
fisher.test(eclsECMOnonematrix)