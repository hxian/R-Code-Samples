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
library(dplyr)
library(XLConnect)
library(plotrix)

# Setting up Data
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

df_csf[['E2.glyc.perc']] <- df_csf$E2.glyc.level * 100
df_csf[['E3.glyc.perc']] <- df_csf$E3.glyc.level * 100
df_csf[['E4.glyc.perc']] <- df_csf$E4.glyc.level * 100

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

shapiro.test(df_csf$site2.E2.glyc.level)
shapiro.test(df_csf$site2.E3.glyc.level)
shapiro.test(df_csf$site2.E4.glyc.level)

shapiro.test(df_pla$E2.glyc.level)
shapiro.test(df_pla$E3.glyc.level)
shapiro.test(df_pla$E4.glyc.level)

# Now The Final Figures. By Haotian
# Figure 1A.
# CSF Secondary Glycosylation % by Genotype
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

View(csf_long)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig1.a.pdf", width = 5, height = 5)
plot <-ggplot(data = csf_long, mapping = aes(x=time, y=site2glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
my_comparisons <- list( c("E2", "E3"), c("E3", "E4"), c("E2", "E4"))
plot + theme_pubr() + 
  stat_compare_means(comparisons = my_comparisons) 
dev.off()

##Statistics
csf_long$site2glyc == ""
csf_long[csf_long$site2glyc == "", "site2glyc"]

sub_csf_long<-subset(csf_long, (!is.na(csf_long[,3]))) 

anov <- aov(data = sub_csf_long, site2glyc~time) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
x <- TukeyHSD(anov) ##Tukey to determine what group comparisons are significant
str(x)
print(x, digits=15)

View(sub_csf_long)

# grab means and se by isoform
mean(sub_csf_long$site2glyc[sub_csf_long$time == "E2"])
mean(sub_csf_long$site2glyc[sub_csf_long$time == "E3"])
mean(sub_csf_long$site2glyc[sub_csf_long$time == "E4"])

std_mean <- function(x) sd(x)/sqrt(length(x))
x=sub_csf_long$site2glyc[sub_csf_long$time == "E2"]
std_mean(x)

sub_csf_long %>%
  group_by(time) %>%
  summarise(mean=mean(site2glyc), se=std.error(site2glyc))

# grabbing means and sd
sub_csf_long %>%
  group_by(time) %>%
  summarise(mean=mean(site2glyc), sd=sd(site2glyc))

# Figure 1B. Plasma Total Glycosylation % by Genotype
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

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig1.b.pdf", width = 5, height = 5)
plot <-ggplot(data = pla_long, mapping = aes(x=time, y=glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
my_comparisons <- list(c("E2", "E4"))
plot + theme_pubr() + 
  stat_compare_means(comparisons = my_comparisons) 
dev.off()

##Statistics
str(pla_long)
a <- lm(formula = glyc ~ time, data = pla_long) 
summary(a)

# grabbing means and se
sub_pla_long<-subset(pla_long, (!is.na(pla_long[,3]))) 
sub_pla_long %>%
  group_by(time) %>%
  summarise(mean=mean(glyc), se=std.error(glyc))

# grabbing means and sd
sub_pla_long<-subset(pla_long, (!is.na(pla_long[,3]))) 
sub_pla_long %>%
  group_by(time) %>%
  summarise(mean=mean(glyc), sd=sd(glyc))

# Figure 1C.
#CSF Total Glycosylation % by Genotype
##wide to long
csf_wide_2 <- select(df_csf, Sample_ID, E2.glyc.perc, 
                   E3.glyc.perc, E4.glyc.perc)
colnames(csf_wide_2)[2] <- 'glyc.x'
colnames(csf_wide_2)[3] <- 'glyc.y'
colnames(csf_wide_2)[4] <- 'glyc.z'
str(csf_wide_2)
csf_wide_2 <- as.data.frame(csf_wide_2)

csf_long_2 <- reshape(data = csf_wide_2, 
                   idvar = "Sample_ID", 
                   varying = list(isoform=c(2:4)), 
                   direction="long", 
                   v.names = c("glyc"),
                   sep="_")

csf_long_2$time[csf_long_2$time == 1] <- "E2"
csf_long_2$time[csf_long_2$time == 2] <- "E3"
csf_long_2$time[csf_long_2$time == 3] <- "E4"

csf_long_2$time <- as.factor(csf_long_2$time)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig1.c.pdf", width = 5, height = 5)
plot <-ggplot(data = csf_long_2, mapping = aes(x=time, y=glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  xlab("") + 
  ylab("CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
my_comparisons <- list( c("E3", "E4"), c("E2", "E4"))
plot + stat_compare_means(comparisons = my_comparisons, method="t.test") + theme_pubr()
dev.off()

#Statistics
csf_long_2$glyc == ""
csf_long_2[csf_long_2$glyc == "", "glyc"]

sub_csf_long_2<-subset(csf_long_2, (!is.na(csf_long_2[,3]))) 

anov <- aov(data = sub_csf_long_2, glyc~time) ##one-way anova since I have more than 2 groups
summary(anov) ##p-value significant 
x <- TukeyHSD(anov) ##Tukey to determine what group comparisons are significant
str(x)
print(x, digits=15)

# grabbing means and se
sub_csf_long_2%>%
  group_by(time) %>%
  summarise(mean=mean(glyc), se=std.error(glyc))

# grabbing means and sd
sub_csf_long_2%>%
  group_by(time) %>%
  summarise(mean=mean(glyc), sd=sd(glyc))

#Figure 2.
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

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig2.pdf", width = 5, height = 5)
plot <-ggplot(data = df_long, mapping = aes(x=time, y=total.glyc, fill=time)) +
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  scale_x_discrete(limits=c("Plasma", "CSF")) + ##rearrange variables on the x-axis
  xlab("") + 
  ylab("Total Glycosylation %") +
  guides(fill = guide_legend(reverse = TRUE, title = NULL)) ##Remove legend title
my_comparisons <- list( c("Plasma","CSF"))
plot + stat_compare_means(comparisons = my_comparisons) + theme_pubr()
dev.off()
##Statistics
stat <- lm(formula = total.glyc ~ time, data = df_long)
summary(stat)

# grabbing mean and se
df_long%>%
  group_by(time) %>%
  summarise(mean=mean(total.glyc), se=std.error(total.glyc))

# grabbing mean and sd
df_long%>%
  group_by(time) %>%
  summarise(mean=mean(total.glyc), sd=sd(total.glyc))

# Figure 3.A.
colnames(df_complete)
df_complete <- df_complete %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_complete$Ab42_status <- ifelse(df_complete$CSF_Ab42_pgperml <=190, 0, 1)
df_complete$Ab42_status[df_complete$Ab42_status == 0] <- "AB42+"
df_complete$Ab42_status[df_complete$Ab42_status == 1] <- "AB42-"
df_complete$Ab42_status <- as.factor(df_complete$Ab42_status)
df_complete$Ab42_status <- factor(df_complete$Ab42_status, levels = c("AB42-", "AB42+"))

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig3a.pdf", width = 5, height = 5)
plot <- ggplot(data=df_complete, 
               aes(x=Ab42_status, y=total.glyc.perc.x, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
my_comparisons <- list( c("AB42-","AB42+"))
plot + stat_compare_means(comparisons = my_comparisons) + theme_pubr()
dev.off()

stat <- lm(formula = total.glyc.perc.x ~ Ab42_status, data = df_complete)      
summary(stat)

# grabbing mean and se
df_complete%>%
  group_by(Ab42_status) %>%
  summarise(mean=mean(total.glyc.perc.x), se=std.error(total.glyc.perc.x))

# grabbing mean and sd
df_complete%>%
  group_by(Ab42_status) %>%
  summarise(mean=mean(total.glyc.perc.x), sd=sd(total.glyc.perc.x))

# Figure 3.B.
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig3b.pdf", width = 5, height = 5)
plot <- ggplot(data=df_complete, 
               aes(x=Ab42_status, y=total.glyc.perc.y, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "Plasma Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title

plot + theme_pubr()
dev.off()

stat <- lm(formula = total.glyc.perc.y ~ Ab42_status, data = df_complete)      
summary(stat)

# grabbing mean and se
df_complete%>%
  group_by(Ab42_status) %>%
  summarise(mean=mean(total.glyc.perc.y), se=std.error(total.glyc.perc.y))

# Figure 3.C.
colnames(df_csf)
df_csf <- df_csf %>% 
  filter(!is.na(CSF_Ab42_pgperml))
df_csf$Ab42_status <- ifelse(df_csf$CSF_Ab42_pgperml <=190, 0, 1)
df_csf$Ab42_status[df_csf$Ab42_status == 0] <- "AB42+"
df_csf$Ab42_status[df_csf$Ab42_status == 1] <- "AB42-"
df_csf$Ab42_status <- as.factor(df_csf$Ab42_status)
df_csf$Ab42_status <- factor(df_csf$Ab42_status, levels = c("AB42-", "AB42+"))

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig3c.pdf", width = 5, height = 5)
plot <- ggplot(data=df_csf, 
               aes(x=Ab42_status, y=totalsite2glyc, fill=Ab42_status)) + 
  geom_boxplot (width = 0.75, size = 0.4, position = position_dodge(0.8), outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=Ab42_status), pch=21, colour="black", position = position_jitter(w=0.1, h=0)) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL)) ##Remove legend title
my_comparisons <- list( c("AB42-","AB42+"))
plot + stat_compare_means(comparisons = my_comparisons) + theme_pubr()
dev.off()

stat <- lm(formula = totalsite2glyc ~ Ab42_status, data = df_csf)      
summary(stat)

# grabbing mean and se
df_csf%>%
  group_by(Ab42_status) %>%
  summarise(mean=mean(totalsite2glyc), se=std.error(totalsite2glyc))

# grabbing mean and sd
df_csf%>%
  group_by(Ab42_status) %>%
  summarise(mean=mean(totalsite2glyc), sd=sd(totalsite2glyc))


# Figure 4 (Supplementary)
df_csf <- df_csf %>% 
  filter(!is.na(Genotype))

df_csf$apoe_e4 = ifelse(df_csf$Genotype == "34" | df_csf$Genotype == "24" | df_csf$Genotype == "44", 1, 0)
df_csf$apoe_e4 = as.factor(df_csf$apoe_e4)
levels(df_csf$apoe_e4) = c("Non E4", "E4")

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig4a.pdf", width = 5, height = 5)
plot <- ggplot(data=df_csf, 
               aes(x=Ab42_status, y=totalsite2glyc, fill=apoe_e4)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "", y = "CSF Secondary Glycosylation %") +
  guides(fill = guide_legend(title = NULL))
plot + stat_compare_means(aes(group = apoe_e4), label = "p.format", method="t.test") + theme_pubr() 
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig4b.pdf", width = 5, height = 5)
plot <- ggplot(data=df_csf, 
               aes(x=Ab42_status, y=total.glyc.perc, fill=apoe_e4)) + 
  geom_boxplot (outlier.shape = NA) +
  theme_classic() +
  geom_point(aes(color=apoe_e4), pch=21, colour="black", position = position_jitterdodge()) +
  labs(x= "", y = "CSF Total Glycosylation %") +
  guides(fill = guide_legend(title = NULL))
plot + theme_pubr() 
dev.off()

# export data set to check something
View(df_csf)
write_xlsx(subset(df_csf, select=c("Sample_ID", "total.glyc.perc", "totalsite2glyc","Genotype","apoe_e4", "CSF_Ab42_pgperml","Ab42_status")), path="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/df_csf.xlsx")

# grabbing mean and se for abeta42-
subset(df_csf, Ab42_status == "AB42-")%>%
  group_by(apoe_e4) %>%
  summarise(mean=mean(totalsite2glyc), se=std.error(totalsite2glyc))

# grabbing mean and sd for abeta42-
subset(df_csf, Ab42_status == "AB42-")%>%
  group_by(apoe_e4) %>%
  summarise(mean=mean(totalsite2glyc), sd=sd(totalsite2glyc))

# grabbing mean and se for abeta42+
subset(df_csf, Ab42_status == "AB42+")%>%
  group_by(apoe_e4) %>%
  summarise(mean=mean(totalsite2glyc), se=std.error(totalsite2glyc))

# grabbing mean and sd for abeta42+
subset(df_csf, Ab42_status == "AB42+")%>%
  group_by(apoe_e4) %>%
  summarise(mean=mean(totalsite2glyc), sd=sd(totalsite2glyc))

# Figure 5. Scatter plots
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.1.pdf", width = 5, height = 5)
a<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_Ab42_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
a
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.2.pdf", width = 5, height = 5)
b<- ggscatter(df_complete, x="total.glyc.perc.y", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 20, label.y =3000) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
b
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.3.pdf", width = 5, height = 5)
c<- ggscatter(df_complete, x="total.glyc.perc.y", y="P_Tau", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =250) +
  geom_point() + 
  labs(x = "Plasma Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
c
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.4.pdf", width = 5, height = 5)
aa<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Ab42 (pg/mL)") 
aa
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.5.pdf", width = 5, height = 5)
bb<- ggscatter(df_complete, x="total.glyc.perc.x", y="CSF_TotalTau_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 65, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF Total Tau (pg/mL)") 
bb
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.6.pdf", width = 5, height = 5)
cc<- ggscatter(df_complete, x="total.glyc.perc.x", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 67, label.y =200) +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF P-Tau (pg/mL)") 
cc
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.7.pdf", width = 5, height = 5)
A <- ggscatter(df_csf, x="totalsite2glyc", y="CSF_Ab42_pgperml", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =1300) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation % ", y = "CSF Ab42 (pg/mL)") 
A
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.8.pdf", width = 5, height = 5)
B<- ggscatter(df_csf, x="totalsite2glyc", y="CSF_TotalTau_pgperml", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 22, label.y =2500) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF Total Tau (pg/mL)") 
B
dev.off()

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig5.9.pdf", width = 5, height = 5)
C <- ggscatter(df_csf, x="totalsite2glyc", y="P_Tau", 
               add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman", label.x = 25, label.y =200) +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF P-Tau (pg/mL)") 
C
dev.off()

# Fig 6. HDL data
# Import HDL data
HDL_data = read_csv("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/DOBRIN HDL.CSV")
View(HDL_data)
# Filter by available HDL3 data
HDL_data = HDL_data[!is.na(HDL_data$`HDL 3 CSF`),]

# merge with df_csf, make sure n=22. if n !=22, then re-run lines 18-41 to reload csf_data (unfiltered by Ab42)
HDL_merged = merge(HDL_data, df_csf, by="Sample_ID")
View(HDL_merged)

# 6.A. scatter: CSF s-HDL-P vs. CSF Total Glycosylation %
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig6.a.pdf", width = 5, height = 5)
A<- ggscatter(HDL_merged, x="total.glyc.perc", y="HDL 3 CSF", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "CSF s-HDL-P") 
A
dev.off()

# 6.B. scatter: CSF s-HDL-P vs. CSF Secondary Glycosylation %
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig6.b.pdf", width = 5, height = 5)
B<- ggscatter(HDL_merged, x="totalsite2glyc", y="HDL 3 CSF", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "CSF s-HDL-P") 
B
dev.off()

# 6.C. scatter: Total CSF HDL-P vs. CSF Total Glycosylation %
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig6.c.pdf", width = 5, height = 5)
C<- ggscatter(HDL_merged, x="total.glyc.perc", y="Total CSF HDL", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Total Glycosylation %", y = "Total CSF HDL-P") 
C
dev.off()

# 6.d. scatter: Total CSF HDL-P vs. CSF Secondary Glycosylation %
pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig6.d.pdf", width = 5, height = 5)
D<- ggscatter(HDL_merged, x="totalsite2glyc", y="Total CSF HDL", 
              add = "reg.line", add.params = list(color = "blue")) +
  stat_cor(method = "spearman") +
  geom_point() + 
  labs(x = "CSF Secondary Glycosylation %", y = "Total CSF HDL-P") 
D + theme_pubr()
dev.off()

# Fig 7. heparin column experiments
#load data
heparin_data = read_xlsx("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/ApoE Gycosylation Paper Data_ Heparin Experiments - Sheet1.xlsx", sheet ="cleaned")
View(heparin_data)
heparin_data$Sialyation = as.factor(heparin_data$Sialyation)
levels(heparin_data$Sialyation) = c("Desialylated", "Sialyated")

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig7.d.pdf", width = 10, height = 5)
ggplot(heparin_data, aes(x = FN, y = APOE.conc, color =Sialyation, group = Sialyation)) +  
  geom_line() + geom_point() + theme_pubr() + guides(color = guide_legend(title = NULL)) +
  ylab("APOE Concentration (µg ml^1)") + 
  scale_x_continuous("Fraction Number", sec.axis = sec_axis(~ . * (1/20), name = "NaCl (M)", breaks = heparin_data$NaCl), breaks = heparin_data$FN)
dev.off()


# Fig 8. extra data rE3 triplicated runs
re3 = read_xlsx("/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/rE3 triplicate runs.xlsx", sheet ="cleaned")
re3$run = as.factor(re3$run)
levels(re3$run) = c("Run 1", "Run 2", "Run 3")
re3$value = as.numeric(re3$value)
re3$value = round(re3$value, 4)

pdf(file="/Users/haotianxian/Documents/Dr. Yassine/CSF Glycosylation (Dobrin)/making nice figures/results/fig8/fig8.pdf", width = 10, height = 5)

ggplot(re3, aes(x = Fraction, y = value, color =run, group = run)) +  
  geom_line() + geom_point() + theme_pubr() + guides(color = guide_legend(title = NULL)) +
  ylab("rE3") + 
  scale_x_continuous("Fraction Number", sec.axis = sec_axis(~ . * (1/20), name = "NaCl (M)", breaks = re3$NaCl), breaks = re3$Fraction)

dev.off()

# Making Table 2, getting the mean and SD Glycosylation % of each group: AD, NCI, MCI
# First Make those clinical status groups
df_complete$Clinical_Status = as.character(df_complete$Clinical_Status)
df_complete$Clinical_Status[is.na(df_complete$Clinical_Status)] <- "NCI"
df_complete$Clinical_Status = factor(df_complete$Clinical_Status, levels = c("NCI", "MCI", "AD"))

df_csf$Clinical_Status = as.character(df_csf$Clinical_Status)
df_csf$Clinical_Status[is.na(df_csf$Clinical_Status)] <- "NCI"
df_csf$Clinical_Status = factor(df_csf$Clinical_Status, levels = c("NCI", "MCI", "AD"))

# Plasma Total Gly % by Clinical Status mean(SD)
mean(df_complete$total.glyc.perc.y)

df_complete %>%
  group_by(Clinical_Status) %>%
  summarise(mean=mean(total.glyc.perc.y), sd=sd(total.glyc.perc.y))

# CSF Total Gly % by Clinal Statis
mean(df_complete$total.glyc.perc.x)

df_complete %>%
  group_by(Clinical_Status) %>%
  summarise(mean=mean(total.glyc.perc.x), sd=sd(total.glyc.perc.x))

# CSF Secondary Gly % by Clinical Status
mean(df_csf$totalsite2glyc)

df_csf %>%
  group_by(Clinical_Status) %>%
  summarise(mean=mean(totalsite2glyc), sd=sd(totalsite2glyc))
