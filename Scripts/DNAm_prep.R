# Prep DNAm for upload to DNAm age calculator (https://dnamage.genetics.ucla.edu/home)
# -------------------------------------------------------
# Create compressed csv files with ALL cPGs listed in "datMiniAnnotation3.csv" ie. 30084 rows (fill rows with NA for missing CpGs)
# Create accompanying annotation file with details on Age, Female for each participant 
# Participant = column of beta file, but row of annotation file

# -------------------------------------------------------
# Ensure you are connected to "ALSPAC" on "DataStore" (smb://cmvm.datastore.ed.ac.uk/cmvm/scs/groups/ALSPAC)
# Check R is pointing to the conda env:
R.home()
.libPaths()

# Load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)

# Check working directory
getwd()

########################################################
# 1. Create cleaned beta files for each DNAm measurement time point
########################################################
# Load in DNAm data (R object is called "betas")
load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj") # takes a wee while as file is large
# Or on Eddie:
#load("/exports/eddie/scratch/s1211670/ALSPAC/data.Robj")
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
       ls Input/
       ")

datMiniAnnotation <- read.csv("Input/datMiniAnnotation3.csv")

head(datMiniAnnotation)

# -------------------------------------------------------
# Subset betas for the time points we want using samplesheet references

processBetas <- function(time){
# Get corresponding time points for betas, ie. was DNAm sampled at birth, age 7, 15 or 17 years?
load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/samplesheet/data.Robj") # R object is called "samplesheet"
# On eddie:
# load("/exports/eddie/scratch/s1211670/ALSPAC/sample_sheet_data.Robj")

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

# Check how many of these CpGs are in ALSPAC:
lapply(betasList, function(x){
  nrow(x[x$ProbeID %in% datMiniAnnotation$Name,])
})

# = 28463

# This is less than the number of CpGs supplied in datMiniAnnotation, 
# therefore we add NA to rows for these missing CpGs:

betasListSub <- lapply(betasList, function(x){
                  newDF <- x[x$ProbeID %in% datMiniAnnotation$Name,]
                  missing <- datMiniAnnotation$Name[!datMiniAnnotation$Name %in% x$ProbeID]
                  missingDF <- data.frame("ProbeID" = missing)
                  idDF <- data.frame(matrix(NA, nrow = length(missing), ncol = ncol(x)-1))
                  colnames(idDF) <- colnames(x)[-1]
                  missingDF <- cbind(missingDF, idDF)
                  newDF <- rbind(newDF, missingDF)
                  return(newDF)
                })
names(betasListSub) <- timePoints


# -------------------------------------------------------
# See how many CpGs (rows) are in each dataframe
nrow(datMiniAnnotation) # 30084 CpGs from the DNAm age calculator

lapply(betasList, nrow) # 482855

lapply(betasListSub, nrow) # 30084 Should match nrows of datMiniAnnotation
# :-) yes it does!
# -------------------------------------------------------
# Check how many participants (columns) are in each dataframe

lapply(betasListSub, ncol)
# $cord
# [1] 906

# $F7
# [1] 971

# $F17
# [1] 718

# $TF3
# [1] 254

# -------------------------------------------------------
# Save a  csv for betas at each time point
write.csv(betasListSub$cord, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_0.csv", row.names = FALSE)
write.csv(betasListSub$F7, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_7.csv", row.names = FALSE)
write.csv(betasListSub$TF3, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_15.csv", row.names = FALSE)
write.csv(betasListSub$F17, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_17.csv", row.names = FALSE)

# On eddie:
# write.csv(betasListSub$cord, "/exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_0.csv", row.names = FALSE)
# write.csv(betasListSub$F7, "/exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_7.csv", row.names = FALSE)
# write.csv(betasListSub$TF3, "/exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_15.csv", row.names = FALSE)
# write.csv(betasListSub$F17, "/exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_17.csv", row.names = FALSE)


# -------------------------------------------------------
# Compress csv files
system("zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_0.zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_0.csv
       rm /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_0.csv")

system("zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_7.zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_7.csv
       rm /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_7.csv")

system("zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_15.zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_15.csv
       rm /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_15.csv")

system("zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_17.zip /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_17.csv
       rm /Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/betas_17.csv")

# On eddie:
# system("zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_0.zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_0.csv
#        rm /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_0.csv")

# system("zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_7.zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_7.csv
#        rm /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_7.csv")

# system("zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_15.zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_15.csv
#        rm /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_15.csv")

# system("zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_17.zip /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_17.csv
#        rm /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/betas_17.csv")


########################################################
# 2. Create an annotation file for each respective beta file
########################################################
# Make sure that the rows (samples) in the sample annotation file have the same order as the columns (samples) in the methylation file.
# col names = Subject, Age, Female
# age in years;  female = 1, male = 0

# -------------------------------------------------------
# Get age (at each time point) and sex for each participant from samplesheet
annotationFunc <- function(time){
  # Get corresponding time points for betas, ie. was DNAm sampled at birth, age 7, 15 or 17 years?
  load("/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/samplesheet/data.Robj") # R object is called "samplesheet"
  # On eddie:
  # load("/exports/eddie/scratch/s1211670/ALSPAC/sample_sheet_data.Robj")

  # Merge these two cols together to get unique individual ID
  samplesheet <- samplesheet %>%
    unite("Subject", c(cidB3421, QLET))

  # Get DNAm sample names for the time point of interest
  samplesheet <- samplesheet %>%
    filter(!str_detect(Subject, "M")) %>% # this removes DNAm for mothers
    filter(str_detect(time_code, time)) %>%
    filter(is.na(duplicate.rm)) %>% # remove duplicates
    select(c(Subject, age, Sex))

  # Test the order of participants is the same as in the betas file
  test_match_order <- function(x,y) {
  if (all(x==y)) print('Perfect match in same order')
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
  }

  test_match_order(colnames(betasListSub[[time]][-1]),samplesheet$Subject)

  # Change column names 
  colnames(samplesheet)[2:3] <- c("Age", "Female")

  # Recode Female column
  samplesheet$Female <- recode(samplesheet$Female, "F" = 1, "M" = 0)

  return(samplesheet)
}

# -------------------------------------------------------
timePoints <- c("cord", "F7", "F17",  "TF3")
annotationsList <- lapply(timePoints , annotationFunc)
names(annotationsList) <- timePoints

# Change age in cord from NA to 0
annotationsList$cord$Age <- 0

# Check
lapply(annotationsList, head)

# -------------------------------------------------------
# Save annotation files for each respective age (eg. ALSPAC/users/amelia/Data/annotation_0.csv)

write.csv(annotationsList$cord, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/annotation_0.csv", row.names = FALSE)
write.csv(annotationsList$F7, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/annotation_7.csv", row.names = FALSE)
write.csv(annotationsList$TF3, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/annotation_5.csv", row.names = FALSE)
write.csv(annotationsList$F17, "/Volumes/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/annotation_7.csv", row.names = FALSE)

# On eddie:
# write.csv(annotationsList$cord, "/exports/eddie/scratch/s1211670/ALSPAC/annotation/annotation_0.csv", row.names = FALSE)
# write.csv(annotationsList$F7, "/exports/eddie/scratch/s1211670/ALSPAC/annotation/annotation_7.csv", row.names = FALSE)
# write.csv(annotationsList$TF3, "/exports/eddie/scratch/s1211670/ALSPAC/annotation/annotation_15.csv", row.names = FALSE)
# write.csv(annotationsList$F17, "/exports/eddie/scratch/s1211670/ALSPAC/annotation/annotation_17.csv", row.names = FALSE)


