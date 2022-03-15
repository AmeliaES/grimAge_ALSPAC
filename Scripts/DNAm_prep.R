# Prep DNAm for upload to GrimAge calculator (https://dnamage.genetics.ucla.edu/home)

# -------------------------------------------------------
# Ensure you are connected to GenScotDepression on DataStore 
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)

# Check working directory
getwd()
# -------------------------------------------------------
# Load in DNAm data (R object is called "betas")
load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj") # takes a wee while as file is large
# Convert matrix to dataframe object
betas <- as.data.frame(betas)
# Make a back up copy just in case we make a silly mistake when wrangling 
bkup <- betas

# -------------------------------------------------------
# Load in subset of CpGs used to make GrimAge (to make the file smaller)
# system("mkdir Input")
system("
       cd Input
       wget https://horvath.genetics.ucla.edu/html/dnamage/datMiniAnnotation.csv
       ls
       cd ../
       pwd
       ")

datMiniAnnotation <- read.csv("Input/datMiniAnnotation.csv")

head(datMiniAnnotation)

# -------------------------------------------------------
# Subset betas for the time points we want using samplesheet references

processBetas <- function(time){
# Get corresponding time points for betas, ie. was DNAm sampled at birth, age 7, 15 or 17 years?
load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/samplesheet/data.Robj") # R object is called "samplesheet"

# Merge these two cols together to get unique individual ID
samplesheet <- samplesheet %>%
  unite("Subject", c(cidB3421, QLET))

# Get DNAm sample names for the time point of interest
samplesheet <- samplesheet %>%
  filter(!str_detect(Subject, "M")) %>% # this removes DNAm for mothers
  filter(str_detect(time_code, time)) %>%
  filter(is.na(duplicate.rm))  # remove duplicates

# Subset to get betas of individuals at a specific time point
betas <- betas %>%
  dplyr::select(samplesheet$Sample_Name)

# additionally we want to rename these cols by the "Subject" col in samplesheet (not the DNAm sample name)
# reorder samplesheet subject to same order as betas cols
samplesheet <- samplesheet[match(colnames(betas), samplesheet$Sample_Name),]
# test that worked:
# sapply(1:ncol(betas), function(i) colnames(betas)[i] == samplesheet$Sample_Name[i]) %>% sum()
# rename cols
betas %>%
  setnames(old = colnames(betas), new = samplesheet$Subject)

betas$ProbeID <- rownames(betas)

betas <- betas %>%
  dplyr::select(ProbeID, everything())

return (betas)
}

# -------------------------------------------------------
# Run function over DNAm samples from birth, age 7, 15 and 17 years
timePoints <- c("cord", "F7", "F17",  "TF3")
betasList <- lapply(timePoints , processBetas)
head(betasList[[1]][,1:10])
# -------------------------------------------------------
# Subset CpGs used to make GrimAge (to make the file smaller)

betasListSub <- lapply(betasList, function(x){
                  x[row.names(x) %in% datMiniAnnotation$Name,]
                })
names(betasListSub) <- timePoints

# -------------------------------------------------------
# See how many CpGs are in each dataframe

lapply(betasList, nrow) # 482,855 CpGs at each time point

lapply(betasListSub, nrow) # 26,973 CpGs in our final dataset

nrow(datMiniAnnotation) # 28,587 CpGs from the GrimAge calculator

# -------------------------------------------------------
# Save a  csv for each time point
write.csv(betasListSub$cord, "/Volumes/ALSPAC/users/amelia/Data/betas_0.csv", row.names = FALSE)
write.csv(betasListSub$F7, "/Volumes/ALSPAC/users/amelia/Data/betas_7.csv", row.names = FALSE)
write.csv(betasListSub$TF3, "/Volumes/ALSPAC/users/amelia/Data/betas_15.csv", row.names = FALSE)
write.csv(betasListSub$F17, "/Volumes/ALSPAC/users/amelia/Data/betas_17.csv", row.names = FALSE)

# -------------------------------------------------------
# Compress csv files
system("zip /Volumes/ALSPAC/users/amelia/Data/betas_0.zip /Volumes/ALSPAC/users/amelia/Data/betas_0.csv
       rm /Volumes/ALSPAC/users/amelia/Data/betas_0.csv")

system("zip /Volumes/ALSPAC/users/amelia/Data/betas_7.zip /Volumes/ALSPAC/users/amelia/Data/betas_7.csv
       rm /Volumes/ALSPAC/users/amelia/Data/betas_7.csv")

system("zip /Volumes/ALSPAC/users/amelia/Data/betas_15.zip /Volumes/ALSPAC/users/amelia/Data/betas_15.csv
       rm /Volumes/ALSPAC/users/amelia/Data/betas_15.csv")

system("zip /Volumes/ALSPAC/users/amelia/Data/betas_17.zip /Volumes/ALSPAC/users/amelia/Data/betas_17.csv
       rm /Volumes/ALSPAC/users/amelia/Data/betas_17.csv")





