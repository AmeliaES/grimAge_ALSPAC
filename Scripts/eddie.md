## Run scripts on EDDIE
It was taking too long to load in the betas file over the VPN connection so decided to run the DNAm_prep.R script on Eddie.

1. Stage in methylation data to scratch
```
qlogin -q staging
mkdir /exports/eddie/scratch/s1211670/ALSPAC/
mkdir /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas
mkdir /exports/eddie/scratch/s1211670/ALSPAC/annotation

# Copy DNAm betas from datastore into scratch space
cp /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/betas/data.Robj /exports/eddie/scratch/s1211670/ALSPAC/

# Do the same with the sample sheet
cp /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/samplesheet/data.Robj /exports/eddie/scratch/s1211670/ALSPAC/sample_sheet_data.Robj

# And other raw data
cp /exports/cmvm/datastore/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta /exports/eddie/scratch/s1211670/ALSPAC/

# Check
ls /exports/eddie/scratch/s1211670/ALSPAC/
```

2. Run `DNAm_prep.R` in R ensuring paths to `/exports/eddie/scratch/s1211670/ALSPAC/` are specified.
```
qlogin -l h_vmem=64G -P sbms_CrossImmune
. /etc/profile.d/modules.sh
module load anaconda
cd /exports/igmm/eddie/GenScotDepression/amelia/grimAge_ALSPAC
source activate ./env
R
```

3. Stage outputs from scratch to datastore
```
mkdir /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input

# Copy betas 
cp /exports/eddie/scratch/s1211670/ALSPAC/QC_Betas/* /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/

# Copy annotation files
cp /exports/eddie/scratch/s1211670/ALSPAC/annotation/* /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/

# check
ls /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/methylation/B3421/grimAge_ALSPAC/Input/
```

4. Now return to the original `README.md` to upload files to the calculator.



