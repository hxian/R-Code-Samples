library(dplyr)
library(readxl)
##Total Glycosylation between CSF and plasma
df_csf <- read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/ashley analysis/data/csf_complete.xlsx", sheet = "Sheet 1")
df_pla <- read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/ashley analysis/data/pla_complete.xlsx", sheet = "Sheet 1")
df_complete <- read_excel("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/ashley analysis/data/csf_pla_complete.xlsx", sheet = "Sheet 1")
colnames(df_csf)
colnames(df_pla)
df_csf[['E2.site2glyc.perc']] <- df_csf$site2.E2.glyc.level * 100
df_csf[['E3.site2glyc.perc']] <- df_csf$site2.E3.glyc.level * 100
df_csf[['E4.site2glyc.perc']] <- df_csf$site2.E4.glyc.level * 100
df_csf[['totalsite2glyc']] <- (df_csf$E2.NHSS.NH.raw + df_csf$E3.NHSS.NH.raw + df_csf$E4.NHSS.NH.raw +
                                 df_csf$E2.NHSS.NHS.raw + df_csf$E3.NHSS.NHS.raw + df_csf$E4.NHSS.NHS.raw + 
                                 df_csf$E2.NHSS.NHSS.raw + df_csf$E3.NHSS.NHSS.raw + df_csf$E4.NHSS.NHSS.raw)/ 
  (df_csf$E2.raw + df_csf$E3.raw + df_csf$E4.raw + df_csf$E2.NHS.raw + df_csf$E3.NHS.raw +
     df_csf$E4.NHS.raw + df_csf$E2.NHSS.raw + df_csf$E3.NHSS.raw + df_csf$E4.NHSS.raw +
     df_csf$E2.NHSS.NH.raw + df_csf$E3.NHSS.NH.raw + df_csf$E4.NHSS.NH.raw +
     df_csf$E2.NHSS.NHS.raw + df_csf$E3.NHSS.NHS.raw + df_csf$E4.NHSS.NHS.raw + 
     df_csf$E2.NHSS.NHSS.raw + df_csf$E3.NHSS.NHSS.raw + df_csf$E4.NHSS.NHSS.raw) * 100
                                                                                                                  

df_pla[['E2.glyc.perc']] <- df_pla$E2.glyc.level * 100
df_pla[['E3.glyc.perc']] <- df_pla$E3.glyc.level * 100
df_pla[['E4.glyc.perc']] <- df_pla$E4.glyc.level * 100


colnames(df_complete)
shapiro.test(df_complete$E2.glyc.level.x)
shapiro.test(df_complete$E3.glyc.level.x)
shapiro.test(df_complete$E4.glyc.level.x)
shapiro.test(df_csf$totalsite2glyc)
shapiro.test(df_complete$total.glyc.perc.y) ##not normally distributed 
shapiro.test(df_complete$total.glyc.perc.x)
shapiro.test(df_complete$CSF_Ab42_pgperml) #not normally distributed
shapiro.test(df_complete$CSF_TotalTau_pgperml) #not normally distributed
shapiro.test(df_complete$P_Tau) #not normally distributed

library(ggpubr)
a<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_Ab42_pgperml", 
                add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
a
b<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 20, label.y =3000) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
b
c<- ggscatter(df_complete, x="total.glyc.perc.y", y="P_Tau", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =250) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
c
aa<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_Ab42_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
aa
bb<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
bb
cc<- ggscatter(df_complete, x="total.glyc.perc.x", y="P_Tau", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 67, label.y =200) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
cc
colnames(df_csf)
A <- ggscatter(df_csf, x="totalsite2glyc", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "CSF Ab42 (pg/mL)") 
A
B<- ggscatter(df_csf, x="totalsite2glyc", y="CSF_TotalTau_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 22, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF Total Tau (pg/mL)") 
B
C <- ggscatter(df_csf, x="totalsite2glyc", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =200) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF P-Tau (pg/mL)") 
C

##wide to long
csf_wide <- select(df_csf, Sample_ID, E2.site2glyc.perc, 
                   E3.site2glyc.perc, E4.site2glyc.perc)
colnames(csf_wide)[2] <- 'site2glyc.x'
colnames(csf_wide)[3] <- 'site2glyc.y'
colnames(csf_wide)[4] <- 'site2glyc.z'
str(csf_wide)
csf_wide <- as.data.frame(csf_wide)
csf_long<- reshape(data = csf_wide, 
                  idvar = "Sample_ID", 
                  varying = list(isoform=c(2:4)), 
                  direction="long", 
                  v.names = c("site2glyc"),
                  sep="_")
csf_long$time[csf_long$time == 1] <- "E2"
csf_long$time[csf_long$time == 2] <- "E3"
csf_long$time[csf_long$time == 3] <- "E4"

csf_long$time <- as.factor(csf_long$time)


library(ggplot2)
plot <-ggplot(data = csf_long, mapping = aes(x=time, y=site2glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
##Statistics
csf_long$site2glyc == ""
csf_long[csf_long$site2glyc == "", "site2glyc"]

sub_csf_long<-subset(csf_long, (!is.na(csf_long[,3]))) 

anov <- aov(data = sub_csf_long, site2glyc~time) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
x <- TukeyHSD(anov) ##Tukey to determine what group comparisons are significant
str(x)
print(x, digits=15)

##wide to long
pla_wide <- select(df_pla, Sample_ID, E2.glyc.perc, 
                   E3.glyc.perc, E4.glyc.perc)
colnames(pla_wide)[2] <- 'glyc.x'
colnames(pla_wide)[3] <- 'glyc.y'
colnames(pla_wide)[4] <- 'glyc.z'
str(pla_wide)
pla_wide <- as.data.frame(pla_wide)
pla_long<- reshape(data = pla_wide, 
                   idvar = "Sample_ID", 
                   varying = list(isoform=c(2:4)), 
                   direction="long", 
                   v.names = c("glyc"),
                   sep="_")
pla_long$time[pla_long$time == 1] <- "E2"
pla_long$time[pla_long$time == 2] <- "E3"
pla_long$time[pla_long$time == 3] <- "E4"

pla_long$time <- as.factor(pla_long$time)


plot <-ggplot(data = pla_long, mapping = aes(x=time, y=glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot

##Statistics
str(pla_long)
a <- lm(formula = glyc ~ time, data = pla_long) 
summary(a)

##wide to long
colnames(df_complete)
df_wide <- select(df_complete, Sample_ID, total.glyc.perc.x, 
                  total.glyc.perc.y)
str(df_wide)
df_wide <- as.data.frame(df_wide)
df_long<- reshape(data = df_wide, 
                   idvar = "Sample_ID", 
                   varying = list(matrix=c(2:3)), 
                   direction="long", 
                   v.names = c("total.glyc"),
                   sep="_")
df_long$time[df_long$time == 1] <- "CSF"
df_long$time[df_long$time == 2] <- "Plasma"

df_long$time <- as.factor(df_long$time)

plot <-ggplot(data = df_long, mapping = aes(x=time, y=total.glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  scale_x_discrete(limits=c("Plasma", "CSF")) + ##rearrange variables on the x-axis
  xlab("") + 
  ylab("Total Glycosylation %") +
  guides(fill = guide_legend(reverse = TRUE, title = NULL)) ##Remove legend title
plot
##Statistics
stat <- lm(formula = total.glyc ~ time, data = df_long)
summary(stat)
##ALL 106 CDR 
##CSF
df_complete$CDR[df_complete$CDR == 0] <- "0"
df_complete$CDR[df_complete$CDR == 0.5] <- "0.5"
df_complete$CDR[df_complete$CDR == 1] <- ">0.5"
df_complete$CDR[df_complete$CDR == 2] <- ">0.5"
df_complete$CDR[df_complete$CDR == 3] <- ">0.5"
df_complete$CDR <- as.factor(df_complete$CDR)
df_complete$CDR <- factor(df_complete$CDR, levels = c("0", "0.5", ">0.5"))
df_complete$CDR 
df_complete$Clinical_Status <- as.factor(df_complete$Clinical_Status)
str(df_complete)

full_csf_CDR <- ggplot(data=subset(df_complete, !is.na(CDR)), 
                       aes(x=CDR, y=total.glyc.perc.x, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  scale_x_discrete(limits=c("0", "0.5", ">0.5")) + ##rearrange variables on the x-axis
  labs(x= "CDR", y = "CSF Total Glycosylation %")  ##Remove legend title
full_csf_CDR

##Statistics
anov <- aov(data = df_complete, total.glyc.perc.x~CDR) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
TukeyHSD(anov) ##Tukey to determine what group comparisons are significant

colnames(df_csf)

# By Haotian: Make binary CDR variable for = 0 and > 0
df_csf$CDR_binary = ifelse(df_csf$CDR > 0, ">0", ifelse(!is.na(df_csf$CDR), "=0", df_csf$CDR))

# By Ashley: Make CDR categories of = 0, = 0.5 and > 0.5
df_csf$CDR[df_csf$CDR == 0] <- "0"
df_csf$CDR[df_csf$CDR == 0.5] <- "0.5"
df_csf$CDR[df_csf$CDR == 1] <- ">0.5"
df_csf$CDR[df_csf$CDR == 2] <- ">0.5"
df_csf$CDR[df_csf$CDR == 3] <- ">0.5"
df_csf$CDR <- as.factor(df_csf$CDR)
df_csf$CDR <- factor(df_csf$CDR, levels = c("0", "0.5", ">0.5"))
df_csf$CDR 

full_csfsite2_CDR <- ggplot(data=subset(df_csf, !is.na(CDR)), 
                       aes(x=CDR, y=totalsite2glyc, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "CDR", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csfsite2_CDR

# By Haotian: Now make the above figure broken down by ApoE4
df_csf$apoe_e4 = ifelse(df_csf$Genotype == 34, 1, 0)
df_csf$apoe_e4 = as.factor(df_csf$apoe_e4)

full_csfsite2_CDR <- ggplot(data=subset(df_csf, !is.na(CDR) & !is.na(apoe_e4)), 
                            aes(x=CDR, y=totalsite2glyc, fill=apoe_e4)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "CDR", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = "ApoE4")) 
full_csfsite2_CDR
View(df_csf)

# By Haotian: Now make the above figure broken down by genotype instead of apoe4, grouping CDR = 0 and CDR > 0 on the x-axis
df_csf$Genotype = as.factor(df_csf$Genotype)
full_csfsite2_CDR_binary <- ggplot(data=subset(df_csf, !is.na(CDR_binary) & !is.na(Genotype)), 
                            aes(x=CDR_binary, y=totalsite2glyc, fill=Genotype)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Genotype), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "CDR Binary", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = "Genotype")) 
full_csfsite2_CDR_binary


##Statistics
anov <- aov(data = df_csf, totalsite2glyc~CDR) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
TukeyHSD(anov) ##Tukey to determine what group comparisons are significant

##AD vs. MCI in CSF 
df_complete$Clinical_Status <- factor(df_complete$Clinical_Status, levels = c("MCI", "AD"))
df_csf$Clinical_Status <- factor(df_csf$Clinical_Status, levels = c("MCI", "AD"))

df_csf$Clinical_Status
aggregate(df_complete$total.glyc.perc.x, list(df_complete$Clinical_Status), mean)

full_csf_clin <- ggplot(data=subset(df_complete, !is.na(Clinical_Status)), 
                       aes(x=Clinical_Status, y=total.glyc.perc.x, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csf_clin
##Statistics n=31
t.test(df_complete$total.glyc.perc.x ~ df_complete$Clinical_Status)

full_csfsite2_clin <- ggplot(data=subset(df_csf, !is.na(Clinical_Status)), 
                            aes(x=Clinical_Status, y=totalsite2glyc, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csfsite2_clin

# By Haotian: draw above figure but add NCI (No Cogn. Imp.) to where clinical status is na, then broken down by Apoe4
df_csf$Clinical_Status = as.character(df_csf$Clinical_Status)
df_csf$Clinical_Status[is.na(df_csf$Clinical_Status)] <- "NCI"
df_csf$Clinical_Status = factor(df_csf$Clinical_Status, levels = c("NCI", "MCI", "AD"))


full_csfsite2_clin <- ggplot(data=subset(df_csf, !is.na(apoe_e4)), 
                             aes(x=Clinical_Status, y=totalsite2glyc, fill=apoe_e4)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "Clinical Status", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = "ApoE4")) 
full_csfsite2_clin

##Statistics n=31
aggregate(df_csf$totalsite2glyc, list(df_csf$Clinical_Status), mean)
t.test(df_csf$totalsite2glyc ~ df_csf$Clinical_Status)
##Plasma       
full_pla_CDR <- ggplot(data=subset(df_complete, !is.na(CDR)), 
                       aes(x=CDR, y=total.glyc.perc.y, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "CDR", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_pla_CDR

##Statistics
stat_pla <- lm(formula = total.glyc.perc.y ~ CDR, data = df_complete)
summary(stat_pla)
##AD vs. MCI in plasma
full_pla_clin <- ggplot(data=subset(df_complete, !is.na(Clinical_Status)), 
                        aes(x=Clinical_Status, y=total.glyc.perc.y, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_pla_clin
##Statistics
stat_pla <- lm(formula = total.glyc.perc.y ~ Clinical_Status, data = df_complete)      
summary(stat_pla) 

##CSF ab42 n=81
colnames(df_complete)
df_complete <- df_complete %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_complete$Ab42_status <- ifelse(df_complete$CSF_Ab42_pgperml <=190, 0, 1)
df_complete$Ab42_status[df_complete$Ab42_status == 0] <- "AB42+"
df_complete$Ab42_status[df_complete$Ab42_status == 1] <- "AB42-"
df_complete$Ab42_status <- as.factor(df_complete$Ab42_status)
df_complete$Ab42_status <- factor(df_complete$Ab42_status, levels = c("AB42-", "AB42+"))
plot <- ggplot(data=df_complete, 
                       aes(x=Ab42_status, y=total.glyc.perc.x, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = total.glyc.perc.x ~ Ab42_status, data = df_complete)      
summary(stat)
##plasma 
plot <- ggplot(data=df_complete, 
               aes(x=Ab42_status, y=total.glyc.perc.y, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = total.glyc.perc.y ~ Ab42_status, data = df_complete)      
summary(stat)

colnames(df_csf)
df_csf <- df_csf %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_csf$Ab42_status <- ifelse(df_csf$CSF_Ab42_pgperml <=190, 0, 1)
df_csf$Ab42_status[df_csf$Ab42_status == 0] <- "AB42+"
df_csf$Ab42_status[df_csf$Ab42_status == 1] <- "AB42-"
df_csf$Ab42_status <- as.factor(df_csf$Ab42_status)
df_csf$Ab42_status <- factor(df_csf$Ab42_status, levels = c("AB42-", "AB42+"))
plot <- ggplot(data=df_csf, 
               aes(x=Ab42_status, y=totalsite2glyc, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = totalsite2glyc ~ Ab42_status, data = df_csf)      
summary(stat)
       
## By Haotian: Make above figure but broken down by ApoE4
plot <- ggplot(data=df_csf, 
               aes(x=Ab42_status, y=totalsite2glyc, fill=apoe_e4)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "AB42 Status", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = "ApoE4"))
plot

##removing the 22 IDs that have already been published       
df_csf_new <- df_csf[-c(99,97,5,6,10,103,41,100,106,95,4,13,98,8,19,28,37,48,89,94,24,36),]
df_pla_new <- df_pla[-c(99,97,5,6,10,103,41,100,106,95,4,13,98,8,19,28,37,48,89,94,24,36),]
df_comp_new <- df_complete[-c(99,97,5,6,10,103,41,100,106,95,4,13,98,8,19,28,37,48,89,94,24,36),]
colnames(df_comp_new)
shapiro.test(df_comp_new$E2.glyc.level.x)
shapiro.test(df_comp_new$E3.glyc.level.x)
shapiro.test(df_comp_new$E4.glyc.level.x)
shapiro.test(df_csf_new$totalsite2glyc)
shapiro.test(df_comp_new$total.glyc.perc.y) ##not normally distributed 
shapiro.test(df_comp_new$total.glyc.perc.x)
shapiro.test(df_comp_new$CSF_Ab42_pgperml)
shapiro.test(df_comp_new$CSF_TotalTau_pgperml)
shapiro.test(df_comp_new$P_Tau)
a<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_Ab42_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
a
b<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 20, label.y =3000) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
b
c<- ggscatter(df_complete, x="total.glyc.perc.y", y="P_Tau", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =250) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
c
aa<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
aa
bb<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_TotalTau_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
bb
cc<- ggscatter(df_complete, x="total.glyc.perc.x", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 67, label.y =200) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
cc
colnames(df_csf)
A <- ggscatter(df_csf, x="totalsite2glyc", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "CSF Ab42 (pg/mL)") 
A
B<- ggscatter(df_csf, x="totalsite2glyc", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 22, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF Total Tau (pg/mL)") 
B
C <- ggscatter(df_csf, x="totalsite2glyc", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =200) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF P-Tau (pg/mL)") 
C

##wide to long
colnames(df_csf_new)
csf_wide <- select(df_csf_new, Sample_ID, E2.site2glyc.perc, 
                   E3.site2glyc.perc, E4.site2glyc.perc)
colnames(csf_wide)[2] <- 'site2glyc.x'
colnames(csf_wide)[3] <- 'site2glyc.y'
colnames(csf_wide)[4] <- 'site2glyc.z'
str(csf_wide)
csf_wide <- as.data.frame(csf_wide)
csf_long<- reshape(data = csf_wide, 
                   idvar = "Sample_ID", 
                   varying = list(isoform=c(2:4)), 
                   direction="long", 
                   v.names = c("site2glyc"),
                   sep="_")
csf_long$time[csf_long$time == 1] <- "E2"
csf_long$time[csf_long$time == 2] <- "E3"
csf_long$time[csf_long$time == 3] <- "E4"

csf_long$time <- as.factor(csf_long$time)
library(ggplot2)
plot <-ggplot(data = csf_long, mapping = aes(x=time, y=site2glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
##Statistics
csf_long$site2glyc == ""
csf_long[csf_long$site2glyc == "", "site2glyc"]

sub_csf_long<-subset(csf_long, (!is.na(csf_long[,3]))) 

anov <- aov(data = sub_csf_long, site2glyc~time) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
x <- TukeyHSD(anov) ##Tukey to determine what group comparisons are significant
str(x)
print(x, digits=15)
##wide to long
colnames(df_pla_new)
pla_wide <- select(df_pla_new, Sample_ID, E2.glyc.perc, 
                   E3.glyc.perc, E4.glyc.perc)
colnames(pla_wide)[2] <- 'glyc.x'
colnames(pla_wide)[3] <- 'glyc.y'
colnames(pla_wide)[4] <- 'glyc.z'
str(pla_wide)
pla_wide <- as.data.frame(pla_wide)
pla_long<- reshape(data = pla_wide, 
                   idvar = "Sample_ID", 
                   varying = list(isoform=c(2:4)), 
                   direction="long", 
                   v.names = c("glyc"),
                   sep="_")
pla_long$time[pla_long$time == 1] <- "E2"
pla_long$time[pla_long$time == 2] <- "E3"
pla_long$time[pla_long$time == 3] <- "E4"

pla_long$time <- as.factor(pla_long$time)
library(ggplot2)
plot <-ggplot(data = pla_long, mapping = aes(x=time, y=glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
##Statistics
str(pla_long)
a <- lm(formula = glyc ~ time, data = pla_long) 
summary(a)
##wide to long
colnames(df_comp_new)
df_wide <- select(df_comp_new, Sample_ID, total.glyc.perc.x, 
                  total.glyc.perc.y)
str(df_wide)
df_wide <- as.data.frame(df_wide)
df_long<- reshape(data = df_wide, 
                  idvar = "Sample_ID", 
                  varying = list(matrix=c(2:3)), 
                  direction="long", 
                  v.names = c("total.glyc"),
                  sep="_")
df_long$time[df_long$time == 1] <- "CSF"
df_long$time[df_long$time == 2] <- "Plasma"

df_long$time <- as.factor(df_long$time)
library(ggplot2)
plot <-ggplot(data = df_long, mapping = aes(x=time, y=total.glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  scale_x_discrete(limits=c("Plasma", "CSF")) + ##rearrange variables on the x-axis
  xlab("") + 
  ylab("Total Glycosylation %") +
  guides(fill = guide_legend(reverse = TRUE, title = NULL)) ##Remove legend title
plot
##Statistics
stat <- lm(formula = total.glyc ~ time, data = df_long)
summary(stat)
##84 CDR 
##CSF
df_comp_new$CDR[df_comp_new$CDR == 0] <- "0"
df_comp_new$CDR[df_comp_new$CDR == 0.5] <- "0.5"
df_comp_new$CDR[df_comp_new$CDR == 1] <- ">0.5"
df_comp_new$CDR[df_comp_new$CDR == 2] <- ">0.5"
df_comp_new$CDR[df_comp_new$CDR == 3] <- ">0.5"
df_comp_new$CDR <- as.factor(df_comp_new$CDR)
df_comp_new$CDR <- factor(df_comp_new$CDR, levels = c("0", "0.5", ">0.5"))

str(df_comp_new)
full_csf_CDR <- ggplot(data=subset(df_comp_new, !is.na(CDR)), 
                       aes(x=CDR, y=total.glyc.perc.x, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "CDR", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csf_CDR
##Statistics
anov <- aov(data = df_comp_new, total.glyc.perc.x~CDR) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
TukeyHSD(anov) ##Tukey to determine what group comparisons are significant


colnames(df_csf_new)
df_csf_new$CDR[df_csf_new$CDR == 0] <- "0"
df_csf_new$CDR[df_csf_new$CDR == 0.5] <- "0.5"
df_csf_new$CDR[df_csf_new$CDR == 1] <- ">0.5"
df_csf_new$CDR[df_csf_new$CDR == 2] <- ">0.5"
df_csf_new$CDR[df_csf_new$CDR == 3] <- ">0.5"
df_csf_new$CDR <- as.factor(df_csf_new$CDR)
df_csf_new$CDR <- factor(df_csf_new$CDR, levels = c("0", "0.5", ">0.5"))

full_csfsite2_CDR <- ggplot(data=subset(df_csf_new, !is.na(CDR)), 
                            aes(x=CDR, y=totalsite2glyc, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "CDR", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csfsite2_CDR
##Statistics
anov <- aov(data = df_csf_new, totalsite2glyc~CDR) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
TukeyHSD(anov) ##Tukey to determine what group comparisons are significant

##AD vs. MCI in CSF 
df_comp_new$Clinical_Status <- factor(df_comp_new$Clinical_Status, levels = c("MCI", "AD"))
df_comp_new$Clinical_Status <- factor(df_comp_new$Clinical_Status, levels = c("MCI", "AD"))
aggregate(df_comp_new$total.glyc.perc.x, list(df_comp_new$Clinical_Status), mean)
full_csf_clin <- ggplot(data=subset(df_comp_new, !is.na(Clinical_Status)), 
                        aes(x=Clinical_Status, y=total.glyc.perc.x, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csf_clin
##Statistics
t.test(df_comp_new$total.glyc.perc.x ~ df_comp_new$Clinical_Status)

full_csfsite2_clin <- ggplot(data=subset(df_csf_new, !is.na(Clinical_Status)), 
                             aes(x=Clinical_Status, y=totalsite2glyc, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_csfsite2_clin
##Statistics
t.test(df_csf_new$totalsite2glyc ~ df_csf_new$Clinical_Status)
##Plasma       
full_pla_CDR <- ggplot(data=subset(df_comp_new, !is.na(CDR)), 
                       aes(x=CDR, y=total.glyc.perc.y, fill=CDR)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=CDR), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "CDR", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_pla_CDR
##Statistics
stat_pla <- lm(formula = total.glyc.perc.y ~ CDR, data = df_comp_new)
summary(stat_pla)
##AD vs. MCI in CSF 
full_pla_clin <- ggplot(data=subset(df_comp_new, !is.na(Clinical_Status)), 
                        aes(x=Clinical_Status, y=total.glyc.perc.y, fill=Clinical_Status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Clinical_Status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
full_pla_clin
stat_pla <- lm(formula = total.glyc.perc.y ~ Clinical_Status, data = df_comp_new)
summary(stat_pla)

##CSF ab42 n=
colnames(df_comp_new)
df_comp_new <- df_comp_new %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_comp_new$Ab42_status <- ifelse(df_comp_new$CSF_Ab42_pgperml <=190, 0, 1)
df_comp_new$Ab42_status[df_comp_new$Ab42_status == 0] <- "AB42+"
df_comp_new$Ab42_status[df_comp_new$Ab42_status == 1] <- "AB42-"
df_comp_new$Ab42_status <- as.factor(df_comp_new$Ab42_status)
df_comp_new$Ab42_status <- factor(df_comp_new$Ab42_status, levels = c("AB42-", "AB42+"))
plot <- ggplot(data=df_comp_new, 
               aes(x=Ab42_status, y=total.glyc.perc.x, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = total.glyc.perc.x ~ Ab42_status, data = df_comp_new)      
summary(stat)
##plasma 
plot <- ggplot(data=df_comp_new, 
               aes(x=Ab42_status, y=total.glyc.perc.y, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = total.glyc.perc.y ~ Ab42_status, data = df_comp_new)      
summary(stat)

colnames(df_csf_new)
df_csf_new <- df_csf_new %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_csf_new$Ab42_status <- ifelse(df_csf_new$CSF_Ab42_pgperml <=190, 0, 1)
df_csf_new$Ab42_status[df_csf_new$Ab42_status == 0] <- "AB42+"
df_csf_new$Ab42_status[df_csf_new$Ab42_status == 1] <- "AB42-"
df_csf_new$Ab42_status <- as.factor(df_csf_new$Ab42_status)
df_csf_new$Ab42_status <- factor(df_csf_new$Ab42_status, levels = c("AB42-", "AB42+"))
plot <- ggplot(data=df_csf_new, 
               aes(x=Ab42_status, y=totalsite2glyc, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
plot
stat <- lm(formula = totalsite2glyc ~ Ab42_status, data = df_csf_new)      
summary(stat)

a<- ggscatter(df_comp_new, x="total.glyc.perc.y", y="CSF_Ab42_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
a
b<- ggscatter(df_comp_new, x="total.glyc.perc.y", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 20, label.y =3000) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
b
c<- ggscatter(df_comp_new, x="total.glyc.perc.y", y="P_Tau", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =250) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
c
aa<- ggscatter(df_comp_new, x="total.glyc.perc.x", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
aa
bb<- ggscatter(df_comp_new, x="total.glyc.perc.x", y="CSF_TotalTau_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
bb
cc<- ggscatter(df_comp_new, x="total.glyc.perc.x", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 67, label.y =200) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
cc
colnames(df_csf)
A <- ggscatter(df_csf_new, x="totalsite2glyc", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "CSF Ab42 (pg/mL)") 
A
B<- ggscatter(df_csf_new, x="totalsite2glyc", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 22, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF Total Tau (pg/mL)") 
B
C <- ggscatter(df_csf_new, x="totalsite2glyc", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =200) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF P-Tau (pg/mL)") 
C

# Haotian's New Analysis
# First Merge cognition data with working data set
# Load cognition dataset
CSF_no = read_xlsx('/Users/haotianxian/Documents/Dr. Yassine/HDL_Ashley/data/CSF Full Data_180.xlsx', sheet = "CSF_Ashley_no outliers")
CSF_no = subset(CSF_no, select = c("Sample_ID", "dom_imp012plus_adj", "Domains_01", "mem_avg_z", "attnexec_avg_z", "lang_avg_z", "glob_avg_z", "apoe_e4"))
View(CSF_no)

# Now Merge with CSF glyc data
cog_CSF = merge(df_csf, CSF_no, by = "Sample_ID")
View(cog_CSF)

# Now Merge with Plasma data
cog_pla = merge(df_pla, CSF_no, by ="Sample_ID")
View(cog_pla)

# Now Merge with Comeplte glyc data
cog_full = merge(df_complete, CSF_no, by = "Sample_ID")
View(cog_full)

# What are the ID's missing cognition data?
a = df_csf$Sample_ID
b = cog_CSF$Sample_ID
diff = setdiff(a, b)
diff = as.data.frame(diff)
write_xlsx(diff, path="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/ashley analysis/data/missing_cognition.xlsx")

# Make scatter for csf total glyc vs memory
A <- ggscatter(cog_full, x="total.glyc.perc.x", y="mem_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "Memory Average Z-Score") 
A

# Make scatter for csf secondary glyc vs memory
B <- ggscatter(cog_CSF, x="totalsite2glyc", y="mem_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Memory Average Z-Score") 
B
# Make scatter for pla total glyc vs memory
C <- ggscatter(cog_full, x="total.glyc.perc.y", y="mem_avg_z", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "Memory Average Z-Score") 
C

# Breaking down by apoe4, memory z score vs CSF secondary glyc
D <- ggscatter(subset(cog_CSF, !is.na(apoe_e4)), x="totalsite2glyc", y="mem_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Memory Average Z-Score") 
D
labels <- c("1" = "APOE4 +", "0" = "APOE4 -")
D + facet_grid(. ~ apoe_e4, labeller=labeller(apoe_e4 = labels))

# breaking down by apoe4, atten.exc z score vs CSF secondary glyc
E <- ggscatter(subset(cog_CSF, !is.na(apoe_e4)), x="totalsite2glyc", y="attnexec_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Attention Executive Function Average Z-Score") 
E
labels <- c("1" = "APOE4 +", "0" = "APOE4 -")
E + facet_grid(. ~ apoe_e4, labeller=labeller(apoe_e4 = labels))

# breaking down by apoe4, language z score vs CSF secondary glyc
F <- ggscatter(subset(cog_CSF, !is.na(apoe_e4)), x="totalsite2glyc", y="lang_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Language Average Z-Score") 
F
labels <- c("1" = "APOE4 +", "0" = "APOE4 -")
F + facet_grid(. ~ apoe_e4, labeller=labeller(apoe_e4 = labels))

# breaking down by apoe4, global z score vs CSF secondary glyc
G <- ggscatter(subset(cog_CSF, !is.na(apoe_e4)), x="totalsite2glyc", y="glob_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Global Average Z-Score") 
G
labels <- c("1" = "APOE4 +", "0" = "APOE4 -")
G + facet_grid(. ~ apoe_e4, labeller=labeller(apoe_e4 = labels))

# breaking down by apoe4, memory z score vs CSF secondary glyc but add AD patients and
# assign them memory z score of -2
# Now Merge with CSF glyc data
cog_CSF_new = merge(df_csf, CSF_no, by = "Sample_ID", all.x = TRUE)
View(cog_CSF_new)
# Now assign memeroy z scor of -2 if AD
cog_CSF_new$mem_avg_z = ifelse(is.na(cog_CSF_new$mem_avg_z) & cog_CSF_new$Clinical_Status == "AD", -2, cog_CSF_new$mem_avg_z)
# Now assign apoe4 status based on genotype
cog_CSF_new$apoe_e4 = ifelse(is.na(cog_CSF_new$apoe_e4) & cog_CSF_new$Genotype == 34, 1, cog_CSF_new$apoe_e4)
cog_CSF_new$apoe_e4 = ifelse(is.na(cog_CSF_new$apoe_e4) & cog_CSF_new$Genotype != 34, 0, cog_CSF_new$apoe_e4)

# Now draw new fig of memeory z score vs CSF secondary glyc
H <- ggscatter(subset(cog_CSF_new, !is.na(apoe_e4) & !(is.na(mem_avg_z))), x="totalsite2glyc", y="mem_avg_z", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "Memory Average Z-Score") 
H
labels <- c("1" = "APOE4 +", "0" = "APOE4 -")
H + facet_grid(. ~ apoe_e4, labeller=labeller(apoe_e4 = labels))

# Get UDS ID of the n=106 dataset
# First import HMRI Dataset
df_hmri <- read_excel("/Users/haotianxian/Documents/Dr. Yassine/Merge Ashley Project/pass HMRI Full Data (HMRI) 20181105.xlsx", sheet = "Sheet1")
df_translate = merge(df_hmri, df_pla, by = "Sample_ID")
View(df_translate)
write_xlsx(df_translate[,1:2], path="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/ashley analysis/data/UDS ID.xlsx")