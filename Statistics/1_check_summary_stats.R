# check the differences between summary statistics characteristics between the two groups

# Load required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# set working directory, here change accordingly
setwd("your_directory")

# I loaded the data by just finding them in their folder
data_all <- read.table("Protocol.csv", sep = ',', header = TRUE)

# ignore excluded participants
data <- data_all[data_all$Excluded == 0,]

#################################################
# check if we have a statistical difference in age
#################################################
# test t-test assumptions:
# Shapiro-Wilk normality test for young age
with(data, shapiro.test(age[age_cohort == "young"]))# p = 0.33
# Shapiro-Wilk normality test for Women's weights
with(data, shapiro.test(age[age_cohort == "old"]))# p = 0.12
# data is normally distributed

# check for homogeneity of variance
res.ftest <- var.test(age ~ age_cohort, data = data)
res.ftest
# there is no equality of variances

# performing t-test with unequal variances
res <- t.test(age ~ age_cohort, data = data, var.equal = FALSE)
res

###########################################################
# check if we have statistical difference in education years
###########################################################

# transform Education in numeric values
data$Education <- as.numeric(data$Education) 
# test t-test assumptions:
# Shapiro-Wilk normality test for young age
with(data, shapiro.test(Education[age_cohort == "young"]))# p = 0.02
# Shapiro-Wilk normality test for Women's weights
with(data, shapiro.test(Education[age_cohort == "old"]))# p = 0.07
# data is not normally distributed

#using non-parametric wilcoxon-rank test

# performing t-test with unequal variances
res <- wilcox.test(Education ~ age_cohort, data, alternative = "two.sided",
                   exact = FALSE)
res
# p = 0.7474
# W = 172


###############################################################
# check if we have a statistical difference in MWTB performance
###############################################################
# test t-test assumptions:
# Shapiro-Wilk normality test for young age
with(data, shapiro.test(MWTB[age_cohort == "young"]))# p = 0.66
# Shapiro-Wilk normality test for Women's weights
with(data, shapiro.test(MWTB[age_cohort == "old"]))# p = 0.021
# data of old participants is not normally distributed

#using non-parametric wilcoxon-rank test

# performing t-test with unequal variances
res <- wilcox.test(MWTB ~ age_cohort, data, alternative = "two.sided",
                   exact = FALSE)
res
# p = 0.0008
# W = 294

######################################################
# check if we have a difference in gender distribution
######################################################
# Create a data frame from the main data set.
str_data = data.frame(data$Gender,data$age_cohort)

# Create a contingency table with the needed variables.           
gender_data = table(data$Gender,data$age_cohort) 

chisq.test(gender_data)
# p = 0.92
# X-sqared = 0.01008

# As the p-value is greater than the .05, we conclude that the gender is independent 
# of the age cohort.
