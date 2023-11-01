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

# Make datset to get summary stats for n=106
# Step 1, main with marged_all dataset, get new
main = read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation/106_PTinfo.xlsx")
View(main)

merged_all = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Merge:Ashley Project/Merged_all.xlsx")
View(merged_all)

merged_all = merged_all[, c("Sample_ID","age", "sex", "edu_yrs", "ethnicity")]

new = merge(main, merged_all, by="Sample_ID", all.x = TRUE)
View(new)

# step 2, Adding information from copy_of_YH
copy_of_YH = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Xulei Project/Lipid Analysis_volcano plot/new/Haotian/Copy of YH_CSF3_Batch23_Biomk.xlsx")
View(copy_of_YH)

copy_of_YH = copy_of_YH[,c("ID","edu_yrs", "sex", "age", "ethnicity")]
colnames(copy_of_YH)[1] = "Sample_ID"

new$age[is.na(new$age)] = copy_of_YH$age[match(new$Sample_ID,copy_of_YH$Sample_ID)][which(is.na(new$age))]
new$sex[is.na(new$sex)] = copy_of_YH$sex[match(new$Sample_ID,copy_of_YH$Sample_ID)][which(is.na(new$sex))]
new$edu_yrs[is.na(new$edu_yrs)] = copy_of_YH$edu_yrs[match(new$Sample_ID,copy_of_YH$Sample_ID)][which(is.na(new$edu_yrs))]
new$ethnicity[is.na(new$ethnicity)] = copy_of_YH$ethnicity[match(new$Sample_ID,copy_of_YH$Sample_ID)][which(is.na(new$ethnicity))]

# Step 3, Adding information from HMRI dataset
HMRI = read_excel("/Users/haotianxian/Documents/Dr. Yassine/Merge:Ashley Project/pass HMRI Full Data (HMRI) 20181105.xlsx")
View(HMRI)
HMRI = HMRI[,c("Sample_ID","Sex", "Age at Sample", "Ethnicity", "Edu Yrs")]
colnames(HMRI) = c("Sample_ID", "sex", "age", "ethnicity", "edu_yrs")

new$age[is.na(new$age)] = HMRI$age[match(new$Sample_ID,HMRI$Sample_ID)][which(is.na(new$age))]
new$sex[is.na(new$sex)] = HMRI$sex[match(new$Sample_ID,HMRI$Sample_ID)][which(is.na(new$sex))]
new$edu_yrs[is.na(new$edu_yrs)] = HMRI$edu_yrs[match(new$Sample_ID,HMRI$Sample_ID)][which(is.na(new$edu_yrs))]
new$ethnicity[is.na(new$ethnicity)] = HMRI$ethnicity[match(new$Sample_ID,HMRI$Sample_ID)][which(is.na(new$ethnicity))]

# Out put to manually add data from VCS and DHA Trial
write_xlsx(new, c("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation/summ_data.xlsx"))

# Load Data Back in
summ_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation/summ_data.xlsx")
View(summ_data)

# make summary table
summ_data[summ_data == "male"] <- "Male"
summ_data[summ_data == "female"] <- "Female"
summ_data[summ_data == "Hispanic"] <- "Hispanic or Latino"
summ_data[summ_data == "not Hispanic or Latinos"] <- "Not Hispanic or Latino"
summ_data[summ_data == "Hispanic or Latinos"] <- "Hispanic or Latino"

summ_data$sex = factor(summ_data$sex)
summ_data$ethnicity = as.factor(summ_data$ethnicity)
summ_data$CDR = as.factor(summ_data$CDR)
summ_data$Genotype = factor(summ_data$Genotype,
                      levels = c(23, 33, 24, 34, 44))
summ_data$Clinical_Status = as.factor(summ_data$Clinical_Status)

# make summary table
library(gtsummary)
tbl = summ_data %>% 
  select (age, sex, ethnicity, edu_yrs, CDR, Genotype, Clinical_Status, CSF_Ab42_pgperml, CSF_TotalTau_pgperml, P_Tau) %>% 
  tbl_summary (missing="ifany", statistic = list(all_continuous() ~ "{mean} ({min}, {max})", all_categorical() ~ "{n} / {N} ({p}%)"),
               digits = all_continuous() ~ 2, label = list(ethnicity ~ "Ethnicity", edu_yrs ~ "Education (Years)", 
              sex ~ "Gender", age ~ "Age", CDR ~ "CDR Score", Genotype ~ "APOE Genotype", Clinical_Status ~ "Clinical Status",
              CSF_Ab42_pgperml ~ "CSF AB42 (pg/mL)", CSF_TotalTau_pgperml ~ "CSF Total Tau (pg/mL)", P_Tau ~ "P-Tau (pg/mL)")) %>% modify_caption("**Demographic and Clinical Characteristics**") %>% add_n() %>% as_gt()
  
theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl %>% gt::gtsave("summ_106.docx", path = c("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation"))

# make summary table for age and male/female
tbl = summ_data %>% 
  select (age, sex) %>% 
  tbl_summary (missing="ifany", statistic = list(all_continuous() ~ "{mean} ({min}, {max})", all_categorical() ~ "{n} / {N} ({p}%)"),
               digits = all_continuous() ~ 2, label = list(sex ~ "Gender", age ~ "Age")) %>% modify_caption("**Demographic and Clinical Characteristics**") %>% add_n() %>% as_gt()

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl %>% gt::gtsave("summ_106_age_gender.html", path = c("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation"))

#### PART 2 #### Adding race data and making new summary table
# port data
summ_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/Haotian Analysis/data/summ_data.xlsx")
View(summ_data)

# clean data
summ_data[summ_data == "male"] <- "Male"
summ_data[summ_data == "female"] <- "Female"
summ_data[summ_data == "Hispanic"] <- "Hispanic or Latino"
summ_data[summ_data == "not Hispanic or Latinos"] <- "Not Hispanic or Latino"
summ_data[summ_data == "Hispanic or Latinos"] <- "Hispanic or Latino"

summ_data$sex = factor(summ_data$sex)
summ_data$ethnicity = as.factor(summ_data$ethnicity)
summ_data$CDR = as.factor(summ_data$CDR)
summ_data$Genotype = factor(summ_data$Genotype,
                            levels = c(23, 33, 24, 34, 44))
summ_data$Clinical_Status = as.factor(summ_data$Clinical_Status)

# make summary table

library(gtsummary)
tbl = summ_data %>% 
  select (age, sex, ethnicity, race, edu_yrs, CDR, Genotype, Clinical_Status, CSF_Ab42_pgperml, CSF_TotalTau_pgperml, P_Tau) %>% 
  tbl_summary (missing="ifany", statistic = list(all_continuous() ~ "{mean} ({min}, {max})", all_categorical() ~ "{n} / {N} ({p}%)"),
               digits = all_continuous() ~ 2, label = list(ethnicity ~ "Ethnicity", edu_yrs ~ "Education (Years)", 
                                                           sex ~ "Gender", age ~ "Age", CDR ~ "CDR Score", Genotype ~ "APOE Genotype", Clinical_Status ~ "Clinical Status",
                                                           CSF_Ab42_pgperml ~ "CSF AB42 (pg/mL)", CSF_TotalTau_pgperml ~ "CSF Total Tau (pg/mL)", P_Tau ~ "P-Tau (pg/mL)")) %>% modify_caption("**Demographic and Clinical Characteristics**") %>% add_n() %>% as_gt()

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl %>% gt::gtsave("summ_106_with_race.html", path = c("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/Haotian Analysis/results"))

###### Make Final Summary Table with Race, but not Ethnicity ####
# First port data and fix clinical status
summ_data = read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/Haotian Analysis/data/summ_data_with race.xlsx")
summ_data$Clinical_Status = as.character(summ_data$Clinical_Status)
summ_data$Clinical_Status[is.na(summ_data$Clinical_Status)] <- "NCI"
summ_data$Clinical_Status = factor(summ_data$Clinical_Status, levels = c("NCI", "MCI", "AD"))

# Now Clean Data
summ_data[summ_data == "male"] <- "Male"
summ_data[summ_data == "female"] <- "Female"

summ_data$sex = factor(summ_data$sex)
summ_data$ethnicity = as.factor(summ_data$ethnicity)
summ_data$CDR = as.factor(summ_data$CDR)
summ_data$Genotype = factor(summ_data$Genotype,
                            levels = c(23, 33, 24, 34, 44))

# Now Make Table
library(gtsummary)
tbl = summ_data %>% 
  select (age, sex, race, edu_yrs, CDR, Genotype, Clinical_Status, CSF_Ab42_pgperml, CSF_TotalTau_pgperml, P_Tau) %>% 
  tbl_summary (missing="ifany", statistic = list(all_continuous() ~ "{mean} ({min}, {max})", all_categorical() ~ "{n} / {N} ({p}%)"),
               digits = all_continuous() ~ 2, label = list(race ~ "Race", edu_yrs ~ "Education (Years)", 
                                                           sex ~ "Gender", age ~ "Age", CDR ~ "CDR Score", Genotype ~ "APOE Genotype", Clinical_Status ~ "Clinical Status",
                                                           CSF_Ab42_pgperml ~ "CSF AB42 (pg/mL)", CSF_TotalTau_pgperml ~ "CSF Total Tau (pg/mL)", P_Tau ~ "P-Tau (pg/mL)")) %>% modify_caption("**Demographic and Clinical Characteristics**") %>% add_n() %>% as_gt()

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl %>% gt::gtsave("summ_106_with_race_new.html", path = c("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/Haotian Analysis/results"))

