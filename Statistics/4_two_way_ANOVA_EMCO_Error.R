# load required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# set directories
setwd("your_directory")

######################
# Error analysis
######################
# I loaded the data by just finding them in their folder
data <- read.table("stat_error_repeat2_ID_ED.csv", sep = ',', header = TRUE)
#turn to factors
data$between <- factor(data$between, levels = c("young", " old "))
data$within  <- factor(data$within, levels = c("repeat", "  ID  ", "  ED  "))

# get summary stats
data %>%
  group_by(within) %>%
  get_summary_stats(dv, type = "mean_sd")

data %>%
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
# There is one extreme outlier, subj_num 26, for now they are not getting removed

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)
# the data is not normally distributed

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# the data are not normally distributed

# check for homogeneity of variance, we want the p value to be higher than 0.05
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# there is no homogeneity of variance within the repeat condition

# check homogeneity of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)
# there is no homogeneity of covariance

# Two-way mixed ANOVA test
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = c(within), effect.size = "pes"
)
get_anova_table(res.aov)

# main effects
pwc <- data %>%
  pairwise_t_test(dv ~ within, p.adjust.method = "holm", paired = TRUE)
pwc

# get cohens d
d <- data %>%
  cohens_d(dv ~ within, paired = TRUE)
d

################################################
# Two-way mixed robust ANOVA test
###############################################

# getting library
library(WRS2)
bwtrim(dv ~ between*within, id = subj_num, data = data)
# Perform post-hoc comparisons on the main effects:
# test the within effect using a bootstrap based approach for one-step M estimators
set.seed(123)
sppbb(dv ~ between*within, id = subj_num, data = data)
# testing with a different approach
# using post-hoc from repeated measure design 
with(data, rmmcp(y = dv, groups = within, block = subj_num))

# getting effect size
repeatData = data$dv[data$within == "repeat"]
IDData = data$dv[data$within == "  ID  "]
EDData = data$dv[data$within == "  ED  "]
set.seed(123)
dep.effect(repeatData, IDData)
dep.effect(repeatData, EDData)
dep.effect(IDData, EDData)

# illustration of the results
# only boxplots
  ggline(data, x = "within", y = "dv", color = "between", fill = "between",
         add = c("dotplot", "mean_ci"), add.params = list(dotsize = 0.5, fill = "between"), 
         position = position_dodge(width = 0.4))+
    labs(x = "Condition", y = "Error rate", color = "Age group")
# adding source functions for plotting
  source("D:/Plot_functions.R")

  
dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
dat_stat
  
pd <- position_dodge(0.4) # move them .25 to the left and right
error_plot =ggplot(dat_stat, aes(x=within, y=dv*100, colour=between)) + 
    geom_dotplot(data = data, aes(x=within, y=dv*100, fill=between), binaxis='y', 
                 stackdir='center',position = pd,
                 dotsize = .5, colour = "grey", show.legend = FALSE)+
    geom_errorbar(aes(ymin=(dv-ci)*100, ymax=(dv+ci)*100), width=.1, position=pd, size = 1) +
    geom_line(aes(group=between),position=pd, size = 1) +
    geom_point(position=pd)+ 
    scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
    labs(x = "Condition", y = "Error rate (%)", color = "Age group")+
    ylim(c(0, 20))+
    theme(legend.position="top", text = element_text(size = 10), 
          axis.line = element_line(size = 0.5))
error_plot
  
# save plot as pdf
setwd(""your_directory"")
pdf("IDED_error_repeat2_ID_ED.pdf",         # File name
    width = 4, height = 3, # Width and height in inches
      bg = "white",    
      paper = "A4")          # Paper size
error_plot
# Close the graphics device
dev.off() 
  
