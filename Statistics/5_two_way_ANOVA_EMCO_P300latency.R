# load libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# set directories
setwd("your_directory")

# P300 analysis - latency
# I loaded the data by just finding them in their folder
data <- read.table("stat_ERP_P300p_latency.csv", sep = ',', header = TRUE)

#turn to factors
data$between <- factor(data$between, levels = c("young", " old "))
data$within  <- factor(data$within, levels = c("repeat", "  ID  ", "  ED  "))

# get summary stats
data %>%
  group_by(between) %>%
  get_summary_stats(dv, type = "mean_sd")

# show boxplot
bxp <- ggboxplot(
  data, x = "between", y = "dv",
  color = "within", palette = "jco"
)
bxp

# identify outliers
data %>%
  group_by(within, between) %>%
  identify_outliers(dv)
# There are no extreme outliers

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# there is normal distribution

# check for homogeneity of variance, we want the p value to be higher than 0.05
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# we have homogeneity of variance

# check homogeneit< of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)

# Two-way mixed ANOVA test
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = c(within), effect.size = "pes"
)
get_anova_table(res.aov)

#########################
# pairwise tests
#########################
pwc <- data %>%
  pairwise_t_test(dv ~between, p.adjust.method = "bonferroni", paired = FALSE)
pwc

# get cohens d
d <- data %>%
  cohens_d(dv ~ between, paired = FALSE)
d

pwc <- data %>%
  group_by(within)%>%
  pairwise_t_test(dv ~between, p.adjust.method = "bonferroni", paired = FALSE)
pwc

#####################################
# Plotting
####################################
# adding source functions for plotting
source("D:/Plot_functions.R")

dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
dat_stat

pd <- position_dodge(0.4) # move them .25 to the left and right
latency_plot =ggplot(dat_stat, aes(x=within, y=dv*1000, colour=between)) + 
  geom_dotplot(data = data, aes(x=within, y=dv*1000, fill=between), binaxis='y', 
               stackdir='center',position = pd,
               dotsize = .5, colour = "grey", show.legend = FALSE)+
  geom_errorbar(aes(ymin=(dv-ci)*1000, ymax=(dv+ci)*1000), width=.1, position=pd, size = 1) +
  geom_line(aes(group=between),position=pd, size = 1) +
  geom_point(position=pd)+ 
  scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
  labs(x = "Condition", y = "P300 Latency (ms", color = "Age group")+
  ylim(c(300, 700))+
  theme(legend.position="top", text = element_text(size = 10), 
        axis.line = element_line(size = 0.5))
latency_plot

# save plot as pdf
setwd("your_directory")
pdf("IDED_P300_latency.pdf",         # File name
    width = 4, height = 3, # Width and height in inches
    bg = "white",    
    paper = "A4")          # Paper size
latency_plot
# Close the graphics device
dev.off() 

