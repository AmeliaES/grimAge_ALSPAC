# Calculate GrimAge for ALSPAC participants using the [online GrimAge calculator](https://dnamage.genetics.ucla.edu/home)   


Follow instructions from the GrimAge calculator tutorial here: https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/TUTORIALonlineCalculator.pdf


### 1. Make conda environment
```
cd /Users/aes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Studies/grimAge_ALSPAC

# Check the dependencies
cat environment.yml

# Create the conda environment (use --force if overwriting existing evironment eg. to add more dependencies)
conda env create --prefix ./env --file environment.yml --force
```

### 2. Activate conda environment and open R studio
```
# Activate the environment and start R studio
conda activate ./env

# Include an "&" after RStudio so even if the terminal closes instance to RStudio remains open
/Applications/RStudio.app/Contents/MacOS/RStudio & 

# Once in R studio check R.home() and .libPaths() point to the conda environment

```

### 3. Prep files for GrimAge calculator
Source `Scripts/DNAm_prep.R` to create compressed csv files to upload to the GrimAge calculator.
Briefly, this script:
1. Reads in DNAm betas from "ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj"
2. Subsets CpGs to those used by the GrimAge calculator (https://horvath.genetics.ucla.edu/html/dnamage/datMiniAnnotation.csv). There is 28,587 CpGs in datMiniAnnotation.csv. ALSPAC contains 26,973 of these CpGs, this is how many CpGs are in the final saved csv files from which GrimAge is calculated from.
3. Saves zip compressed csv files for each DNAm collection time point (birth, age 7, 15 and 17 years).   

These csv files are saved on DataStore at:
 `/ALSPAC/users/amelia/Data/`   

### 4. Upload csv files to GrimAge online calculator
Create an account to submit methylation data (https://dnamage.genetics.ucla.edu/submit). Upload zip compressed csv files. **Check the box which says "Normalize Data"** - see tutorial doc (link at top) for more info on this. We will calculate age acceleration using additional covariates (cell count estimates) so we will not include an annotation file as an optional additional upload.


