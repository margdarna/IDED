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

# Comparison between two groups
fit_ASSTD <- survfit(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data, conf.type = "log-log")
summary(fit_ASSTD)
p = ggsurvplot(fit_ASSTD, data = data, conf.int = TRUE, xlab = "Stage number", title = "ASSTD",
           break.time.by = 1, xlim = c(0.5, 8.5))
p

survdiff(Surv(ASSTD_surv, FAILED_ASSTD) ~ Age, data = data)
