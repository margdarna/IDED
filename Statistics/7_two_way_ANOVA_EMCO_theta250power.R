# Load libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# set working directory
setwd("your_directory")

# Theta analysis
# I loaded the data by just finding them in their folder
data <- read.table("stat_TFR_theta250.csv", sep = ',', header = TRUE)

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
# there are no extreme outliers

# check for inhomogeneity of variances
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# there is no inhomogeneity

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# data are fairly normally distributed

# check for homogeneity of variance, we want the p value to be higher than 0.05
data %>%
  group_by(between) %>%
  levene_test(dv ~ within)

# check homogeneit< of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)

############################
# Two-way mixed ANOVA test
############################
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = c(within), effect.size = "pes"
)
get_anova_table(res.aov)

# get summary stats
data %>%
  group_by(between, within) %>%
  get_summary_stats(dv, type = "mean_sd")

#pairwise comparisons for interaction
pwc <- data %>%
  group_by(between) %>%
  pairwise_t_test(dv ~ within, p.adjust.method = "holm", paired = TRUE)
pwc

pwc <- data %>%
  group_by(within) %>%
  pairwise_t_test(dv ~between, paired = FALSE)
pwc

# main effects
pwc <- data %>%
  pairwise_t_test(dv ~between, paired = FALSE)
pwc

pwc <- data %>%
  pairwise_t_test(dv ~ within, paired = TRUE)
pwc

# get cohens d
d <- data %>%
  group_by(between) %>%
  cohens_d(dv ~ within, paired = TRUE)
d

d <- data %>%
  group_by(within) %>%
  cohens_d(dv ~ between, paired = FALSE)
d

#################################
# illustration of the results
#################################

# only boxplots
  ggline(data, x = "within", y = "dv", color = "between", fill = "between",
         add = c("dotplot", "mean_ci"), add.params = list(dotsize = 0.5, fill = "between"), 
         position = position_dodge(width = 0.4))+
    labs(x = "Condition", y = "Frontocentral theta power (dB)", color = "Age group")

  # adding source functions for plotting
  source("D:/Plot_functions.R")
  
  dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
  dat_stat
  
  pd <- position_dodge(0.4) # move them .25 to the left and right
  theta_plot =ggplot(dat_stat, aes(x=within, y=dv, colour=between)) + 
    geom_dotplot(data = data, aes(x=within, y=dv, fill=between), binaxis='y', 
                 stackdir='center',position = pd,
                 dotsize = .5, colour = "grey", show.legend = FALSE)+
    geom_errorbar(aes(ymin=dv-ci, ymax=dv+ci), width=.1, position=pd, size = 1) +
    geom_line(aes(group=between),position=pd, size = 1) +
    geom_point(position=pd)+ 
    scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
    labs(x = "Condition", y = "Theta Power (dB)", color = "Age group")+
    ylim(c(-3, 6))+
    theme(legend.position="top", text = element_text(size = 10), 
          axis.line = element_line(size = 0.5))
  theta_plot
  
  # save plot as pdf
  setwd("your_directory")
  pdf("IDED_theta250_amplitude.pdf",         # File name
      width = 4, height = 3, # Width and height in inches
      bg = "white",    
      paper = "A4")          # Paper size
  theta_plot
  # Close the graphics device
  dev.off() 
  
