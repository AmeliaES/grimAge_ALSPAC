# A few error/warnings came up when calculating DNAm ages
# let's check them here

setwd("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/")
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(ggplot2)
library(gridExtra)
# ----------------------------------------------------
ages <- c("0", "7", "15", "17")
outputs <- lapply(ages, function(age){
  x <- read.csv(paste0("Output/scores_", age, ".csv"))
  x$Female <- as.factor(x$Female)
  return(x)
  })

# ----------------------------------------------------
# check what warnings we got
lapply(outputs, function(x) unique(x$Comment))

# " Warning: meanMethBySample is <0.25 MAJOR WARNING:Predicted gender does not agree with column Female."
# ----------------------------------------------------
# Plot the meanXchromosome values against Female for all ages
# Check whether outliers of one sex is near the mean of the other sex
labels <- paste0("Age ", ages, " years")
plots <- lapply(1:4, function(i){
          ggplot(outputs[[i]], aes(x = Female, y = meanXchromosome)) +
            geom_boxplot() +
            xlab("Sex") +
            ggtitle(labels[i])
          })

do.call("grid.arrange", c(plots, ncol=2))

# ----------------------------------------------------




