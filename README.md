# Calculate DNA methylation age for ALSPAC participants

## Background information
*Insert info about what DNA methylation age is*

Follow instructions from the DNA methylation age calculator (Horvath, 2013) tutorial here: https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/TUTORIALonlineCalculator.pdf. I reccomend reading through this and Horvath (2013) before running scripts.

DNA methylation ages calculated include:
- pan tissue clock (Horvath 2013)
- blood based epigenetic clocks (Hannum et al., 2013)
- skin and blood clock (Horvath et al., 2018)
- estimators of phenotypic age (Levine et al., 2017)
- mortality risk DNAmGrimAge (Lu et al., 2019a)
- estimator of telomere length: DNAmTL (Lu et al., 2019b)
- epigenetic estimators of plasma proteins e.g. DNAmPAI1 (Lu et al., 2019a)
-  DNAm estimate of smoking pack years(DNAmPACKYRS) (Lu et al., 2019a)

ALSPAC participants (n = ~900) had DNA methylation measured in umbilical cord blood and peripheral blood from age 7 year, and age 15 or 17 years clinic visits using Illumina Infinium HumanMethylation450 BeadChip (450 K) array (Relton et al., 2015). Final sample sizes for DNA methylation ages are:
- age 0 years: n = 906
- age 7 years: n = 971
- age 15 years: n = 718
- age 17 years: n = 254

## Calculate DNA methylation age and age acceleration in ALSPAC
### 1. Make conda environment
```
# Change directory to the location you want to calculate DNAm ages in and have cloned this repo.
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

### 3. Prep files for DNA methylation age calculator
Source `Scripts/DNAm_prep.R` to create compressed csv files to upload to the DNA methylation age calculator. You may need to change relative paths in this script to wherever the ALSPAC server is mapped to on your computer. *If it takes too long on your local computer please see `Scripts/eddie.md` for running this script on Eddie.*
Briefly, this script:
1. Reads in DNAm betas from "ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj"
2. Subsets CpGs to those used by the DNA methylation age calculator listed in "datMiniAnnotation3.csv" (https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/datMiniAnnotation3.csv). There is 30,084 CpGs in datMiniAnnotation3.csv. ALSPAC contains 28,463 of these CpGs. For the remaining CpGs with no data in ALSPAC have rows filled in with NA, as specified by the tutorial.
3. Saves zip compressed csv files for each DNAm collection time point (birth, age 7, 15 and 17 years). In these files each row is a CpG and each column is a participant.
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/betas_0.zip
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/betas_7.zip
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/betas_15.zip
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/betas_17.zip
    
4. Creates annotation files (.csv) for each respective beta file. These contain column names: Subject, Age, Female for each participant. Each row is a participant.
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/annotation_0.csv
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/annotation_7.csv
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/annotation_15.csv
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Input/annotation_17.csv

### 4. Upload csv files to DNA methylation age online calculator
Create an account to submit methylation data (https://dnamage.genetics.ucla.edu/new). Upload zip compressed csv files. **Check the boxes which say "Normalize Data" and "Advanced Analysis"** - see tutorial doc (link at top) for more info on this. Upload corresponding age annotation file. An email is sent when analysis has run with a file for DNA methylation ages and age acceleration values, see the tutorial document for details on output.
I saved the output to:
    - ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Output/

### 5. Calculate age acceleration using additional covariates (cell count estimates)
TODO - make script for this
`AgeAccelerationResidual=residuals(lm(DNAmAge~Age))`
Output saved in:
`ALSPAC/data/genomics/B3421/methylation/B3421/grimage_ALSPAC/Output/`

## References
- Horvath S (2013) DNA methylation age of human tissues and cell types. Genome Biol
14(10):R115 PMID: 24138928

- Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, Klotzle B, Bibikova M, Fan
JB, Gao Y, Deconde R, Chen M, Rajapakse I, Friend S, Ideker T, Zhang K (2013) Genomewide methylation profiles reveal quantitative views of human aging rates. Mol Cell. 2013
Jan 24;49(2):359-67.

- Horvath S, Oshima J, Martin GM, Lu AT, Quach A, Cohen H, Felton S, Matsuyama M,
Lowe D, Kabacik S, Wilson JG, Reiner AP, Maierhofer A, Flunkert J, Aviv A, Hou L,
Baccarelli AA, Li Y, Stewart JD, Whitsel EA, Ferrucci L, Matsuyama S, Raj K. (2018)
Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome
and ex vivo studies. Aging (Albany NY). 2018 Jul 26;10(7):1758-1775. doi:
10.18632/aging.101508. PMID: 30048243 PMCID: PMC6075434

- Levine ME, Lu AT, Quach A, Chen BH, Assimes TL, Bandinelli S, Hou L, Baccarelli AA,
Stewart JD, Li Y, Whitsel EA, Wilson JG, Reiner AP, Aviv A, Lohman K, Liu Y, Ferrucci
L, Horvath S. An epigenetic biomarker of aging for lifespan and healthspan. Aging (Albany
NY). 2018 Apr 18;10(4):573-591. doi: 10.18632/aging.101414. PMID: 29676998 PMCID:
PMC5940111

- Lu AT, Quach A, Wilson JG, Reiner AP, Aviv A, Raj K, Hou L, Baccarelli AA, Li Y,
Stewart JD, Whitsel EA, Assimes TL, Ferrucci L, Horvath S. (2019). DNA methylation
GrimAge strongly predicts lifespan and healthspan. Aging (Albany NY). 2019 Jan
21;11(2):303-327. doi: 10.18632/aging.101684. PMID: 30669119 PMCID: PMC6366976

- Lu AT, Seeboth A, Tsai PC, Sun D, Quach A, Reiner AP, Kooperberg C, Ferrucci L, Hou
L, Baccarelli AA, Li Y, Harris SE, Corley J, Taylor A, Deary IJ, Stewart JD, Whitsel EA,
Assimes TL, Chen W, Li S, Mangino M, Bell JT, Wilson JG, Aviv A, Marioni RE, Raj K
Horvath S (2019) DNA methylation-based estimator of telomere length. Aging (Albany
NY). 2019 Aug 18;11(16):5895-5923. doi: 10.18632/aging.102173. PMID: 31422385
PMCID: PMC6738410

- Caroline L Relton, Tom Gaunt, Wendy McArdle, Karen Ho, Aparna Duggirala, Hashem Shihab, Geoff Woodward, Oliver Lyttleton, David M Evans, Wolf Reik, Yu-Lee Paul, Gabriella Ficz, Susan E Ozanne, Anil Wipat, Keith Flanagan, Allyson Lister, Bastiaan T Heijmans, Susan M Ring, George Davey Smith, Data Resource Profile: Accessible Resource for Integrated Epigenomic Studies (ARIES), International Journal of Epidemiology, Volume 44, Issue 4, August 2015, Pages 1181–1190, https://doi.org/10.1093/ije/dyv072

- Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH,
Wiencke JK, Kelsey KT (2012) DNA methylation arrays as surrogate measures of cell
mixture distribution. BMC Bioinformatics 2012, 13:86 doi:10.1186/1471-2105-13-86



