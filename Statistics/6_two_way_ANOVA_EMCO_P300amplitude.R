library(tidyverse)
library(ggpubr)
library(rstatix)
#setwd("C:/Users/mdarna/Documents/PhD/A05-SFB1436/Output/4_Stats")
setwd("C:/Users/mdarna/Documents/PhD/A05-SFB1436/IDED_v1_Analysis/Output/repeat2_ascontrol/4_Stats")

# P300 amplitude analysis
# I loaded the data by just finding them in their folder
data <- read.table("stat_ERP_P300p.csv", sep = ',', header = TRUE)

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

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# data arw fairly normally distributed

# check for homogeneity of variance, we want the p value to be higher than 0.05
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# there is inhomogeneity of variance within the ED condition

# check homogeneit< of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)

# Two-way mixed ANOVA test
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = c(within), effect.size = "pes"
)
get_anova_table(res.aov)
# There is a significant main effect for within

# main effects
pwc <- data %>%
  pairwise_t_test(dv ~ within, p.adjust.method = "holm", paired = TRUE)
pwc


# get cohens d
d <- data %>%
  cohens_d(dv ~ within, paired = TRUE)
d


# illustration of the results
# only boxplots
ggline(data, x = "within", y = "dv", color = "between", fill = "between",
         add = c("dotplot", "mean_ci"), add.params = list(dotsize = 0.5, fill = "between"), 
         position = position_dodge(width = 0.4))+
    labs(x = "Condition", y = "P300 amplitude (uV)", color = "Age group")

# adding source functions for plotting
source("D:/Plot_functions.R")

dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
dat_stat

pd <- position_dodge(0.4) # move them .25 to the left and right
amplitude_plot =ggplot(dat_stat, aes(x=within, y=dv, colour=between)) + 
  geom_dotplot(data = data, aes(x=within, y=dv, fill=between), binaxis='y', 
               stackdir='center',position = pd,
               dotsize = .5, colour = "grey", show.legend = FALSE)+
  geom_errorbar(aes(ymin=dv-ci, ymax=dv+ci), width=.1, position=pd, size = 1) +
  geom_line(aes(group=between),position=pd, size = 1) +
  geom_point(position=pd)+ 
  scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
  labs(x = "Condition", y = "P300 amplitude (uV)", color = "Age group")+
  ylim(c(-10, 20))+
  theme(legend.position="top", text = element_text(size = 10), 
        axis.line = element_line(size = 0.5))
amplitude_plot

# save plot as pdf
setwd("//linstore01/home/mdarna/PhD/A05-SFB1436/IDED_v1_Analysis/Output/5_StatFigures")
pdf("IDED_P300_amplitude.pdf",         # File name
    width = 4, height = 3, # Width and height in inches
    bg = "white",    
    paper = "A4")          # Paper size
amplitude_plot
# Close the graphics device
dev.off() 
