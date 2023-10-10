# Load libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# Set working directories
setwd("your_directory)

# Single Trial analysis
# I loaded the data by just finding them in their folder
data <- read.table("IDED_theta250_single_trial_shepherd_repeat2_ID_ED.csv", sep = ',', header = TRUE)

#turn to factors
data$between <- factor(data$is_young, levels = c("1", "0"))
data$within  <- factor(data$within, levels = c("repeat", "  ID  ", "  ED  "))

# identify outliers
data %>%
  group_by(within, between) %>%
  identify_outliers(dv)
# there are no extreme outliers

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# data are normally distributed

# check for inhomogeneity of variances
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# there is no inhomogeneity of variances

# check homogeneit< of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)

############################
# Two-way mixed ANOVA test 
############################
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = c(within), effect.size = "pes"
)
res.aov
# there is a main effect of age group with p = 0.004 and pes=0.204

# pairwise comparison for main effects
pwc <- data %>%
  pairwise_t_test(dv ~between, paired = FALSE)
pwc

# get summary stats for significant effects
data %>%
  group_by(between) %>%
  get_summary_stats(dv, type = "mean_sd")

################################
# illustration of the results
################################
# only boxplots
  ggline(data, x = "within", y = "dv", color = "between",fill = "between",
         add = c("dotplot", "mean_ci"), add.params = list(dotsize = 0.5, fill = "between"), 
         position = position_dodge(width = 0.4))+
    labs(x = "Condition", y = "Slope", color = "Age group")

################################################
# one sample t-test per group per condition
################################################

t_young_repeat= t.test(data[data$between== "1" & data$within== "repeat", ]$dv, mu = 0)
t_young_repeat

t_young_ID= t.test(data[data$between== "1" & data$within== "  ID  ", ]$dv, mu = 0)
t_young_ID

t_young_ED= t.test(data[data$between== "1" & data$within== "  ED  ", ]$dv, mu = 0)
t_young_ED

t_old_repeat= t.test(data[data$between== "0" & data$within== "repeat", ]$dv, mu = 0)
t_old_repeat

t_old_ID= t.test(data[data$between== "0" & data$within== "  ID  ", ]$dv, mu = 0)
t_old_ID

t_old_ED= t.test(data[data$between== "0" & data$within== "  ED  ", ]$dv, mu = 0)
t_old_ED

# perform fdr correction on the values
# generate data frame with p values
library(fuzzySim)
pvalues <- data.frame(var = c("young_repeat", "young_ID", "young_ED", 
                              "old_repeat", "old_ID", "old_ED"), 
                      pval = c(t_young_repeat$p.value, t_young_ID$p.value, 
                               t_young_ED$p.value, t_old_repeat$p.value, 
                               t_old_ID$p.value, t_old_ED$p.value))
FDR(pvalues = pvalues, correction = "fdr")

# calculating cohens d for significant effects
data[data$between== "1" & data$within== "  ID  ", ] %>%
  cohens_d(dv ~ 1, mu = 0)
data[data$between== "1" & data$within== "  ED  ", ] %>%
  cohens_d(dv ~ 1, mu = 0)
data[data$between== "0" & data$within== "  ED  ", ] %>%
  cohens_d(dv ~ 1, mu = 0)

#############################
# show results in boxplot
############################
dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
dat_stat

# adding source functions for plotting
source("D:/Plot_functions.R")

pd <- position_dodge(0.4) # move them .25 to the left and right
singletrial_plot =ggplot(dat_stat, aes(x=within, y=dv, colour=between)) + 
  geom_dotplot(data = data, aes(x=within, y=dv, fill=between), binaxis='y', 
               stackdir='center',position = pd,
               dotsize = .5, colour = "grey", show.legend = FALSE)+
  geom_errorbar(aes(ymin=dv-ci, ymax=dv+ci), width=.1, position=pd, size = 1) +
  geom_line(aes(group=between),position=pd, size = 1) +
  geom_point(position=pd)+ 
  scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
  labs(x = "Condition", y = "z transformed Pi", color = "Age group")+
  ylim(c(-0.3, 0.750))+
  theme(legend.position="top", text = element_text(size = 10), 
        axis.line = element_line(size = 0.5))
singletrial_plot

# save plot as pdf
setwd("//linstore01/home/mdarna/PhD/A05-SFB1436/IDED_v1_Analysis/Output/5_StatFigures")
pdf("IDED_theta_single_trial_shepherd.pdf",         # File name
    width = 4, height = 3, # Width and height in inches
    bg = "white",    
    paper = "A4")          # Paper size
singletrial_plot
# Close the graphics device
dev.off() 
