library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(survminer)
library(survival)
library(ggplot2)

setwd("//linstore01/home/mdarna/PhD/A05-SFB1436/IDED_v1_Analysis/Output/4_Stats")

# Load in data
data <- read.table("ASST_survival.csv", sep = ',', header = TRUE)

data_young = data[data$Age == 1,]
data_old = data[data$Age == 2,]

# Table and Curves using log
fit <- survfit(Surv(ASSTD_surv, FAILED_ASSTD) ~ 1, data = data, conf.type = "log-log")
summary(fit)
ggsurvplot(fit, data = data)

# Comparison between two groups
fit_ASSTD <- survfit(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data, conf.type = "log-log")
summary(fit_ASSTD)
p = ggsurvplot(fit_ASSTD, data = data, conf.int = TRUE, xlab = "Stage number", title = "ASSTD",
           break.time.by = 1, xlim = c(0.5, 8.5))
p

survdiff(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data)

fit_ASSTG <- survfit(Surv(ASSTG_surv, FAILED_ASSTG) ~ Age, data = data, conf.type = "log-log")
ggsurvplot(fit_ASSTG, data = data, conf.int = TRUE, xlab = "Stage number", title = "ASSTG",
           break.time.by = 1, xlim = c(0.5, 8))

survdiff(Surv(ASSTG_surv, FAILED_ASSTG) ~ Age, data = data)

## assess the difference in trials needed to finish the ED
data_all <- read.table("stat_ASST_error_trials_to_criterion.csv", sep = ',', header = TRUE)
data_EDS <- data_all[data_all$within == " EDS ",]
# exclude NaNs
data_nonan <-data_EDS[data_EDS$dv !='NaN',]
# exclude participants that did not complete the stage
data <-data_nonan[data_nonan$dv !=50,]

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
  group_by(between) %>%
  identify_outliers(dv)
# there are two extreme outliers in the young group

# check if normally distributed, we want the p value to be higher than 0.05
data %>%
  group_by(within, between) %>%
  shapiro_test(dv)
# the data is not normally distributed

# draw ggplot
ggqqplot(data, "dv", ggtheme = theme_bw()) +
  facet_grid(1~ between)
# data is normally distributed

data_young = data[data$between == "young",]
data_old = data[data$between == " old ",]

# performing t test
t.test(data_young$dv, data_old$dv, var.equal = FALSE)
