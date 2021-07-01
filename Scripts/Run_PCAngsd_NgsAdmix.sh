#!/bin/bash
##########################################################################################################
# CONDUCT PRINCIPAL COMPONENT ANALYSIS AND ESTIMATE ADMIXTURE PROPORTIONS FROM GLs
# Requires the following:
# 1. Genotype likelihood file in .beagle format (produced by ANGSD)
# 2. bam.list previously used as input for ANGSD
# 3. PCAngsd - http://www.popgen.dk/software/index.php/PCAngsd
# 4. NgsAdmix - http://www.popgen.dk/software/index.php/NgsAdmix

# A. Sendell-Price, July 2021
##########################################################################################################

#-- SETTING UP ENVIRONMENT ----------------------------
#Define path to .beagle file
BEAGLE=/home/zoo/Myzomela/Lcass.v1_WholeGenome.beagle.gz

#Define prefix for output
OUT=Myzomela_Indv81_SNPs5M

#Define path to PCAngsd
PCAngsd=/home/zoo/BIN/pcangsd/pcangsd.py

#Note on NgsAdmix
#NgsAdmix is included as part of the angsd installation on Nesoi
#so we dont need to define the path for this software


#-- ESTIMATE COVARIANCE MATRIX USING PCANGSD -----------

#make directory for PCA and move into it
mkdir PCA
cd PCA

#Calculate covariance matrix
#Note: By default the maximum number of iterations performed is 100, this has been
#increased to 10000 using "-iter" flag to ensure there are sufficient itetations
#to converge.
python $PCAngsd \
-beagle $BEAGLE \
-out $OUT \
-iter 10000

#Once complete PCANGSD will write two file:
# .log - a log file providing details of run
# .cov - a covariance matrix (.cov) which can be used
#        to perfom PCA in R

#move out of PCA directory
cd ../

#-- ESTIMATE ADMIXTURE PROPORTIONS USING NGSADMIX -------

#make directory for Admixture analysis and move into it
mkdir ADMIX
cd ADMIX

#Using a for loop run NGSadmix using different k values
for K in 2 3 4 5 6 7 8 9 10
do

NGSadmix \
-likes $BEAGLE \
-K $K \
-o ${OUT}_K${K} -P 10

done

#move out of ADMIX directory
cd ../


#-- EXTRACT SAMPLE IDs FROM BAM.LIST -------------------

#Create a list of sample names using bam.list using awk and sed
cat bam.list | awk -F/ '{print $NF}' | sed 's/.sorted.bam//' > samplesIDs.txt
#1             #2                      #3                    #4
#1 - print bam.list to screen (sample_bams/ABY15.sorted.bam)
#2 - using awk field seperator flag (-F) seperate by "/" and print final chunk (ABY15.sorted.bam)
#3 - using sed substitute flag substitute ".sorted.bam" with "" (ABY15)
#4 - write to file "sampleIDs.txt"
