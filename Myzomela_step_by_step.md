# Step by step bioinformatics guide for Sonya's Myzomela project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## STEP 1: Log in to ARC move into Myzomela project directory
Let's log into Arcus-htc (replace OSS with your single sign on).
```
ssh OSS@oscgate.arc.ox.ac.uk
ssh arcus-htc

```

Move into data storage for project zool-zost
```
cd /data/zool-zost/Myzomela
```

## STEP 2: Download reference genome and build reference dictionary



In our first assessment of hybridization across the genome we will calculate a single D-statistic for the genome. This will give us an idea of the overall extent of introgression between the two species. For this we will use the program [DSuite](https://github.com/millanek/Dsuite) which allows for the fast calculation of the D-statistic from SNP data in VCF format.

First, set path to DSuite and specify the VCF file
```
DSUITE=~/Dropbox/DPhil/BIN/Dsuite/Build/Dsuite
VCF=../VCF/ZFified_NorfolkHybridization_Zlatv1_Biallelic_NoIndels_MinQC20_MinDP4_MaxMiss0.5.vcf.gz
```
To run the analysis we will need to create a tab delimited file ("Pops.txt") that assigns samples to populations. That file looks like this, where "Outgroup" is used to assigns individuals to the outgroup (in this case we use the Reunion grey white-eye a.k.a *Z.borbonicus*) and remaining samples are assigned to their corresponding populations:
