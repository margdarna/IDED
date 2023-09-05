# Load required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# set working directory
setwd("//linstore01/home/mdarna/PhD/A05-SFB1436/IDED_v1_Analysis/Output/4_Stats")

# Reaction time ana
# I loaded the data by just finding them in their folder
data <- read.table("stat_RT_repeat_ID_ED.csv", sep = ',', header = TRUE)
#turn to factors
data$between <- factor(data$between, levels = c("young", " old "))
data$within  <- factor(data$within, levels = c("repeat", "  ID  ", "  ED  "))

# get summary stats
data %>%
  group_by(between) %>%
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
# there are no extreme outliers

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)
# data is normally distributed

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(within ~ between)
# data is normally distributed

# check for homogeneity of variance, we want the p value to be higher than 0.05
data %>%
  group_by(within) %>%
  levene_test(dv ~ between)
# p value is lower than 0.05 in the ID and ED condition, condition do not homogeneity of variance between age groups

# check homogeneity of covariances, we want the p value to be higher than 0.05
box_m(data[, "dv", drop = FALSE], data$within)
# there is homogeneity of covariance

# Two-way mixed ANOVA test
res.aov <- anova_test(
  data = data, dv = dv, wid = subj_num,
  between = between, within = within, effect.size = "pes"
)
res.aov
get_anova_table(res.aov)
# Main effect for:
# between, p< .001, pes = 0.612, F = 58.300
# within , p <.001, pes = 0.689, F = 81.942
# between:within, p= 0.064, pes = 0.086, F = 3.48

# pairwise comparisons for main effects
pwc_within <- data %>%
  pairwise_t_test(dv ~ within, p.adjust.method = "holm", paired = TRUE)
pwc_within

pwc_between <- data %>%
  pairwise_t_test(dv ~between, paired = FALSE)
pwc_between

# get cohens d
d_within <- data %>%
  cohens_d(dv ~ within, paired = TRUE)
d_within

d_between <- data %>%
  cohens_d(dv ~ between, paired = FALSE)
d_between





# investigating trend for interaction
pwc_within <- data %>%
  group_by(between)%>%
  pairwise_t_test(dv ~ within, p.adjust.method = "holm", paired = TRUE)
pwc_within

# calculate difference for every participant
data_wide = data %>% pivot_wider(names_from = within, values_from = dv)
data_wide$ID_repeat = data_wide$`  ID  ` - data_wide$`repeat`
data_wide$ED_repeat = data_wide$`  ED  ` - data_wide$`repeat`
data_wide$ED_ID = data_wide$`  ED  ` - data_wide$`  ID  `

data_long = data_wide %>% pivot_longer(cols=c('ID_repeat', 'ED_repeat', 'ED_ID'),
                    names_to='condition_dif',
                    values_to='dv')

# investigating trend for interaction
# perform ANOVA
res.aov <- anova_test(
  data = data_long, dv = dv, wid = subj_num,
  between = between, within = condition_dif, effect.size = "pes"
)
get_anova_table(res.aov)

pwc_between <- data_long %>%
  group_by(condition_dif)%>%
  pairwise_t_test(dv ~ between, p.adjust.method = "holm", paired = FALSE, , pool.sd=FALSE)
pwc_between

# get summary stats
summary_stats = data_long %>%
  group_by(condition_dif, between) %>%
  get_summary_stats(dv, type = "mean_sd")

# get cohens d
d_between <- data_long %>%
  group_by(condition_dif) %>%
  cohens_d(dv ~ between, paired = FALSE)
d_between




## Perform robust ANOVA
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

# checking the between effect
sppba(dv ~ between*within, id = subj_num , data = data)
# get effect size using robust cohen.s d
akp.effect(dv ~ between, data = data)

# illustration of the results
# only boxplots
ggline(data, x = "within", y = "dv", color = "between",
         add = c("dotplot", "mean_ci"), add.params = list(dotsize = 0.5), 
         position = position_dodge(width = 0.4))+
    labs(x = "Condition", y = "Mean Reaction Time (ms)", color = "Age group")

# adding source functions for plotting
source("D:/Plot_functions.R")

dat_stat <- summarySEwithin(data, measurevar="dv", betweenvars="between", withinvars="within")
dat_stat

# SVG graphics device
pd <- position_dodge(0.4) # move them .25 to the left and right
RT_plot =ggplot(dat_stat, aes(x=within, y=dv*1000, colour=between)) + 
  geom_dotplot(data = data, aes(x=within, y=dv*1000, fill=between), binaxis='y', 
               stackdir='center',position = position_dodge(0.4),
               dotsize = .5, colour = "grey", show.legend = FALSE)+
  geom_errorbar(aes(ymin=(dv-ci)*1000, ymax=(dv+ci)*1000), width=.1, position=pd, size = 1) +
  geom_line(aes(group=between),position=pd, size = 1) +
  geom_point(position=pd)+ 
  scale_fill_manual(values = c("grey", "grey"))+ theme_classic()+
  labs(x = "Condition", y = "Median Reaction Time (ms)", color = "Age group")+
  ylim(c(400, 1500))+
  theme(legend.position="top", text = element_text(size = 10), 
        axis.line = element_line(size = 0.5))
RT_plot

# save plot as pdf
setwd("//linstore01/home/mdarna/PhD/A05-SFB1436/IDED_v1_Analysis/Output/5_StatFigures")
pdf("IDED_RT_repeat2_ID_ED.pdf",         # File name
    width = 4, height = 3, # Width and height in inches
    bg = "white",    
    paper = "A4")          # Paper size
RT_plot
# Close the graphics device
dev.off() 
