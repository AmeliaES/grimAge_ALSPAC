# A few error/warnings came up when calculating DNAm ages
# let's check them here

setwd("/Volumes/cmvm/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/")
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
# ----------------------------------------------------
output <- read.csv("Output/scores_15.csv")

# ----------------------------------------------------
head(output)

sum(output$Comment == "") # only 85 without a warning

unique(output$Comment)

# " Warning: meanMethBySample is <0.25 MAJOR WARNING:Predicted gender does not agree with column Female."

# ----------------------------------------------------
# Plot the meanXchromosome values against Female,
plot(as.factor(output$Female),output$meanXchromosome, ylab = "meanXchromosome" , xlab = "Sex")
# Outliers of one sex is not near the mean of the other sex, so all looks fine here.

# ----------------------------------------------------
# How many have the mean < 0.25 warning?
sum(str_detect(output$Comment, "<0.25")) #103
nrow(output) #253



