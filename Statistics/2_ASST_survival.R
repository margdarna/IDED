# load required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(survminer)
library(survival)
library(ggplot2)

# set directories
setwd("your_directory")

# Load data
data <- read.table("ASST_survival.csv", sep = ',', header = TRUE)

# separate by age group
data_young = data[data$Age == 1,]
data_old = data[data$Age == 2,]

# Comparison between two groups
fit_ASSTD <- survfit(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data, conf.type = "log-log")
summary(fit_ASSTD)

survdiff(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data)


# plot results
p = ggsurvplot(fit_ASSTD, data = data, conf.int = TRUE, xlab = "Stage number", title = "ASSTD",
           break.time.by = 1, xlim = c(0.5, 8.5))
p
