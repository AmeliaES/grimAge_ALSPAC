# Calculate DNAm age acceleration

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
# Let's calculate age acceleration for GrimAge at each time point, adjusting for cell estimates
# Adjust for batch effects?

lapply(1:4, function(i){
  outputs[[i]]$DNAmGrimAgeAdjWBCAge <- residuals(lm(DNAmGrimAge ~ Age * Female + CD8T + CD4T + NK + Bcell + Mono + Gran, data = outputs[[i]]))
  write.csv(outputs[[i]], paste0("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Output/scores_acceleration_", ages[i] ,".csv"))
  head(outputs[[i]])
})

