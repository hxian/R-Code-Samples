library(readxl)
library(ggplot2)
library(ggpubr)
library(dbplyr)
library(Hmisc)
library(psych)      # correlation p
library(writexl)      # export excel

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

library(XLConnect)
library(cutpointr)
library(report)
library(pROC)

library(tidyverse)
library(haven)
library(gtsummary)
library(gt)

# merge dataset: simoa (measured by yassine lab) and msd (measured by james)
# load datset

simoa_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/Summary Report ptau181 Plasma and AB42 CSF 2022_05_05.xlsx", sheet = "simoa_clean")
msd_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/Summary Report ptau181 Plasma and AB42 CSF 2022_05_05.xlsx", sheet = "msd_clean")

# merged_data = merge(simoa_data, msd_data, by="ID")
# dir = c("/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project")
# write_xlsx(merged_data, file.path(dir, "Merged_needs_edit.xlsx"))

# Manually delete duplicate entries for subjects 203708 and 1101827
View(merged_data)

# Make summary table
merged_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/Merged.xlsx")
View(merged_data)

tbl = merged_data %>% 
  select (Age, Sex, CDR, Diabetes, Plasma_pTau181_simoa, CSF_AB42_40_ratio_simoa, CSF_pTau181_MSD, CSF_AB42_40_ratio_MSD) %>% 
  tbl_summary (missing="ifany", statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} / {N} ({p}%)"), 
               digits = all_continuous() ~ 2, 
               label = list(Age ~ "Age at Blood Draw (years)", Sex ~ "Sex", CDR ~ "Clinical Dementia Rating", Diabetes ~ "Diabetes Status", 
                            Plasma_pTau181_simoa ~ "Plasma pTau181 (Simoa)", CSF_AB42_40_ratio_simoa ~ "CSF AB42:40 Ratio (Simoa)",
                            CSF_pTau181_MSD ~ "CSF pTau181 (MSD)", CSF_AB42_40_ratio_MSD ~ "CSF AB42:40 Ratio (MSD)")) %>% 
  add_n() %>% modify_header(label = "**Variable**") %>%  modify_caption("**Table 1. Basic Demographic Information and Summary Statistics**") %>% as_gt()

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl%>% gt::gtsave("summ_merge.html", path = dir)

# make boxplots for variables by CDR
# need new variable CDR_cat
merged_data$CDR_cat = NA
merged_data$CDR_cat[merged_data$CDR <= 0.5] = 1
merged_data$CDR_cat[merged_data$CDR > 0.5 & merged_data$CDR <= 1] = 2
merged_data$CDR_cat[merged_data$CDR > 1] = 3

# Boxplot: Plasma_pTau181_simoa by CDR_cat
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/boxplot_pla_ptau181_simoa_by_CDR.pdf", width = 5, height = 5)
ggplot(merged_data,aes(y = Plasma_pTau181_simoa, x = as.factor(CDR_cat), fill = as.factor(CDR_cat))) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_jitterdodge(), aes(fill = as.factor(CDR_cat))) + 
  theme_classic() +
  labs(x = "CDR Category", y = "Plasma P-Tau 181 (Simoa)") +
  scale_fill_discrete(name = "CDR", labels = c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1")) + 
  scale_x_discrete(labels=c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1"))
dev.off()

# Boxplot: CSF_AB42_40_ratio_simoa by CDR_cat
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/boxplot_csf_ab_ratio_simoa_by_CDR.pdf", width = 5, height = 5)
ggplot(merged_data,aes(y = CSF_AB42_40_ratio_simoa, x = as.factor(CDR_cat), fill = as.factor(CDR_cat))) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_jitterdodge(), aes(fill = as.factor(CDR_cat))) + 
  theme_classic() +
  labs(x = "CDR Category", y = "CSF AB42:AB40 Ratio (Simoa)") +
  scale_fill_discrete(name = "CDR", labels = c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1")) + 
  scale_x_discrete(labels=c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1"))
dev.off()

# Boxplot: CSF_pTau181_MSD by CDR_cat
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/boxplot_csf_ptau181_MSD_by_CDR.pdf", width = 5, height = 5)
ggplot(merged_data,aes(y = CSF_pTau181_MSD, x = as.factor(CDR_cat), fill = as.factor(CDR_cat))) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_jitterdodge(), aes(fill = as.factor(CDR_cat))) + 
  theme_classic() +
  labs(x = "CDR Category", y = "CSF P-Tau 181 (MSD)") +
  scale_fill_discrete(name = "CDR", labels = c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1")) + 
  scale_x_discrete(labels=c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1"))
dev.off()

# Boxplot: CSF_AB42_40_ratio_MSD by CDR_cat
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/boxplot_csf_ab_ratio_MSD_by_CDR.pdf", width = 5, height = 5)
ggplot(merged_data,aes(y = CSF_AB42_40_ratio_MSD, x = as.factor(CDR_cat), fill = as.factor(CDR_cat))) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_jitterdodge(), aes(fill = as.factor(CDR_cat))) + 
  theme_classic() +
  labs(x = "CDR Category", y = "CSF AB42:AB40 Ratio (MSD)") +
  scale_fill_discrete(name = "CDR", labels = c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1")) + 
  scale_x_discrete(labels=c("CDR <= 0.5", "0.5 < CDR <=1", "CDR > 1"))
dev.off()

# Boxplot Plasma pTau181 by Sex
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/boxplot_pla_ptau181_simoa_by_sex.pdf", width = 5, height = 5)
ggplot(merged_data,aes(y = Plasma_pTau181_simoa, x = as.factor(Sex), fill = as.factor(Sex))) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_jitterdodge(), aes(fill = as.factor(Sex))) + 
  theme_classic() +
  labs(x = "Sex", y = "Plasma P-Tau 181 (Simoa)") +
  stat_compare_means()
dev.off()

# Scatter plots
# Simoa Data: Correlate Plasma PTau181 with CSFAB42:40 Ratio
shapiro.test(merged_data$Plasma_pTau181_simoa)
shapiro.test(merged_data$CSF_AB42_40_ratio_simoa)

# use speatman because one of the variables is not normaly distributed
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_simoa.pdf", width = 5, height = 5)
ggplot(merged_data, aes(x=Plasma_pTau181_simoa, y=CSF_AB42_40_ratio_simoa)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  ylab("CSF AB42:AB40 Ratio (Simoa)") +
  xlab("Plasma P-Tau 181 (Simoa)") +
  stat_cor(method = "spearman")
dev.off()

# MSD Data: Correlate CSF PTau181 with CSFAB42:40 Ratio
shapiro.test(merged_data$CSF_pTau181_MSD)
shapiro.test(merged_data$CSF_AB42_40_ratio_MSD)

# use spearman because both variables not normally distributed
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_MSD.pdf", width = 5, height = 5)
ggplot(merged_data, aes(x=CSF_pTau181_MSD, y=CSF_AB42_40_ratio_MSD)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  ylab("CSF AB42:AB40 Ratio (MSD)") +
  xlab("CSF P-Tau 181 (MSD)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate MSD to Simoa, PTau181
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_simoa_MSD_ptau.pdf", width = 5, height = 5)
ggplot(merged_data, aes(x=Plasma_pTau181_simoa, y=CSF_pTau181_MSD)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  ylab("CSF P-Tau 181 (MSD)") +
  xlab("Plasma P-Tau 181 (Simoa)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate MSD to Simoa, CSF AB42:40
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_simoa_MSD_abratio.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=CSF_AB42_40_ratio_simoa, x=CSF_AB42_40_ratio_MSD)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("CSF AB42:AB40 Ratio (MSD)") +
  ylab("CSF AB42:AB40 Ratio (Simoa)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate CSF AB42 simoa with Plasma pTau181(Simoa)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_5.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=as.numeric(CSF_AB42_simoa), x=Plasma_pTau181_simoa)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("Plasma p-Tau181 (Simoa)") +
  ylab("CSF AB42 (Simoa)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate CSF AB42 simoa with CSF pTau181(MSD)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_6.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=as.numeric(CSF_AB42_simoa), x=CSF_pTau181_MSD)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("CSF p-Tau181 (MSD)") +
  ylab("CSF AB42 (Simoa)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate CSF AB42 MSD with Plasma pTau181(Simoa)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_7.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=as.numeric(CSF_AB42_MSD), x=Plasma_pTau181_simoa)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("Plasma p-Tau181 (Simoa)") +
  ylab("CSF AB42 (MSD)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate CSF AB42 MSD with CSF PTau181 (MSD)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_8.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=as.numeric(CSF_AB42_MSD), x=CSF_pTau181_MSD)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("CSF p-Tau181 (MSD)") +
  ylab("CSF AB42 (MSD)") +
  stat_cor(method = "spearman")
dev.off()

# Correlate Plasma pTau181 Simoa with age
shapiro.test(merged_data$Plasma_pTau181_simoa)
shapiro.test(merged_data$Age)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/results/scatter_9.pdf", width = 5, height = 5)
ggplot(merged_data, aes(y=Plasma_pTau181_simoa, x=Age)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE) +
  guides(fill = FALSE) +
  xlab("Age") +
  ylab("Plasma p-Tau181 (Simoa)") +
  stat_cor(method = "spearman")
dev.off()


# New Section: Cutpoint for predicting amyloid pathology
# Reference paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4774260/
# cutpoint using cutpoint r package, but first load data set with apoe genotype
df_all = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Simoa_msd project/Merged_apoe.xlsx")

# add new variable abeta_pos for abeta 42 pos/neg based on csf ab42 < 190 as positive (msd)
df_all$abeta_pos = ifelse(df_all$CSF_AB42_MSD < 190, 1, 0)
df_all$abeta_pos = as.factor(df_all$abeta_pos)
levels(df_all$abeta_pos) = c("Neg", "Pos")

# add new variable abeta_pos_2 for abeta pos/neg based on CSF Ab42/40 <0.09 as positive (msd)
df_all$abeta_pos_2 = ifelse(df_all$CSF_AB42_40_ratio_MSD < 0.09, 1, 0)
df_all$abeta_pos_2 = as.factor(df_all$abeta_pos_2)
levels(df_all$abeta_pos_2) = c("Neg", "Pos")

# add new variable abeta_pos_3 for abeta pos/neg based on CSF Ab42/40 <0.16 as positive (Quanterix simoa)
df_all$abeta_pos_3 = ifelse(df_all$CSF_AB42_40_ratio_simoa < 0.16, 1, 0)
df_all$abeta_pos_3 = as.factor(df_all$abeta_pos_3)
levels(df_all$abeta_pos_3) = c("Neg", "Pos")

# add new variable apoe_e4 based on apoe genotype
df_all$apoe_e4 = ifelse(df_all$apoe == "3/4" | df_all$apoe == "2/4" | df_all$apoe == "4/4", 1, 0)
df_all$apoe_e4 = as.factor(df_all$apoe_e4)
levels(df_all$apoe_e4) = c("Non E4", "E4")

# Using abeta_pos as outcome variable, find cutpoint for Plasma_pTau181_simoa
# make sure to filter out NA abeta values, df_all_1 is a copy of df_all without NA abtea_pos rows
df_all_1 = df_all[!is.na(df_all$abeta_pos),]
opt_cut <- cutpointr(df_all_1, Plasma_pTau181_simoa, abeta_pos)
summary(opt_cut)
plot(opt_cut)

# Now add age and apoe_e4 to model to improve fit, i.e. want increased AUC.
# plan:
# dicotaimize ptau181 based on cutoff from cutpointr
# run auc function with dicotimized ptau (using 1.233 cutoff) and abeta_pos to verify same AUC as cutpointr output
# glm binomial age, apoe, ptau181
# predict function to generate probabilities for each subject
# proc package auc function : abeta_pos (outcome) and predicted probabilities (x variable) (transformed variable, with age and apoe)
# Also use 2.1 for above steps as cutoff for ptau to compare performance
# then lastly, add two more abeta positive groups: MSD, CSF Ab42/40 <0.09 as positive, simoa CSF Ab42/40 <0.16 as positive
# Using the two additional abeta pos groups, repeat above steps

# Start here: outcome: abeta_pos, ptau cutoff: 1.233
# first we make dichotomized Plasma_pTau181_simoa using cutoff of 1.233 from above output
# then calculate auc
df_all_1$di_Plasma_pTau181_simoa = ifelse(df_all_1$Plasma_pTau181_simoa > 1.233, 1, 0)
roc(abeta_pos ~ di_Plasma_pTau181_simoa, data=df_all_1)
# results says:Area under the curve: 0.6336, which is the same as the out put from cutpointr

# next we run glm with di_Plasma_pTau181_simoa, age and apoe4
glm_model = glm(abeta_pos ~ di_Plasma_pTau181_simoa + Age + apoe_e4, binomial, df_all_1)

# using predict function with "response" option, we get predicted probabilities of all patients in df_all_1
predicted_probs = predict(glm_model, type = "response")
roc(abeta_pos ~ predicted_probs, data=df_all_1)
# above roc function returns area under curve as 0.7119, this is an improvement over 0.6336

# Now, use ptau cutoff at 2.1, repeat above calculations
df_all_1$di_Plasma_pTau181_simoa_2 = ifelse(df_all_1$Plasma_pTau181_simoa > 2.1, 1, 0)
glm_model = glm(abeta_pos ~ di_Plasma_pTau181_simoa_2 + Age + apoe_e4, binomial, df_all_1)
predicted_probs = predict(glm_model, type = "response")
roc(abeta_pos ~ predicted_probs, data=df_all_1)
# above roc function returns area under curve as 0.7227, this is an improvement over 0.6336

# Run model with just di_Plasma_pTau181 (cutoff at 2.1)
glm_model_2 = glm(abeta_pos ~ di_Plasma_pTau181_simoa_2, binomial, df_all_1)
predicted_probs_2 = predict(glm_model_2, type = "response")
roc(abeta_pos ~ predicted_probs_2, data=df_all_1)
