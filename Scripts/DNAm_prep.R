# Prep DNAm for upload to DNAm age calculator (https://dnamage.genetics.ucla.edu/home)
# -------------------------------------------------------
# Create compressed csv files with ALL cPGs listed in "datMiniAnnotation3.csv" ie. 30084 rows (fill rows with NA for missing CpGs)
# Create accompanying annotation file with details on Age, Female for each participant 
# Participant = column of beta file, but row of annotation file

# -------------------------------------------------------
# Ensure you are connected to "ALSPAC" on "DataStore"
# Check R is pointing to the conda env:
R.home()
.libPaths()

# Load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(haven)

# Check working directory
getwd()

########################################################
# 1. Create cleaned beta files for each DNAm measurement time point
########################################################
# Load in DNAm data (R object is called "betas")
load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj") # takes a wee while as file is large
# Convert matrix to dataframe object
betas <- as.data.frame(betas)
# Make a back up copy just in case we make a silly mistake when wrangling 
# (though this uses more memory so you might not want to do this)
bkup <- betas

# -------------------------------------------------------
# Load in subset of CpGs used to make DNAm age (to make the file smaller)
# This is a list of 30,084 CpGs listed in "datMiniAnnotation3.csv" found from the Hovarth Calulation website under the "New Methylation Age Calculator" (https://dnamage.genetics.ucla.edu/new) - Advanced Analysis.
system("mkdir Input")
system("
       wget -P Input/ https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/datMiniAnnotation3.csv
       ls
       ")

datMiniAnnotation <- read.csv("Input/datMiniAnnotation3.csv")

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
# Subset CpGs to those used to make DNAm age (to make the file smaller)
# Insert NA for any CpGs (rows) we don't have data for

betasListSub <- lapply(betasList, function(x){
                  x[row.names(x) %in% datMiniAnnotation$Name,]
                })
names(betasListSub) <- timePoints

# -------------------------------------------------------
# See how many CpGs (rows) are in each dataframe
nrow(datMiniAnnotation) # XX CpGs from the DNAm age calculator

lapply(betasList, nrow) 

lapply(betasListSub, nrow) # Should match nrows of datMiniAnnotation

# -------------------------------------------------------
# Check how many participants (columns) are in each dataframe



# -------------------------------------------------------
# Save a  csv for betas at each time point
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

########################################################
# 2. Create an annotation file for each respective beta file
########################################################
# Make sure that the rows (samples) in the sample annotation file have the same order as the columns (samples) in the methylation file.
# col names = Subject, Age, Female
# age in years;  female = 1, male = 0

# -------------------------------------------------------
# Read in data to find age and sex of each participant
otherData <- read_dta("/Volumes/ALSPAC/data/B3421_Whalley_04Nov2021.dta")

library(haven)
## Some useful haven functions:
# Get the labels for the col names
vars <- var_label(data)
# Get labels of values for each variable
val_labels(data)

vars[str_detect(names(vars), "FKMS1030")]
vars[str_detect(vars, "age")]
# -------------------------------------------------------
# Get age and sex for each participant

# -------------------------------------------------------
# Check participants (ie. rows) in same order as column names of beta files (betasListSub)


# -------------------------------------------------------
# Save annotation files for each respective age (eg. ALSPAC/users/amelia/Data/annotation_0.csv)




