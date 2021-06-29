#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=0-47:1
#SBATCH --time=120:00:00
#SBATCH --job-name=ANGSD_BEAGLE
#SBATCH --output=ANGSD_BEAGLE_%A_%a.log
#SBATCH --error=ANGSD_BEAGLE_%A_%a.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

##########################################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION WITH ANGSD v.0.921 AND IMPUTE AND PHASE GENOTYPES USING BEAGLEv.4
# Will output the following files per scaffold:
# 1. VCF with genotype probability (GP) and genotype likelihood (GL) fields
# 2. MAFs
# 3. Genotype likelihoods (GLs) in beagle format
# 4. 2 x VCFs, one with imputed genotypes and one with imputed/phased genotypes

# A. Sendell-Price, June 2021
##########################################################################################################

#Load required modules
module load angsd/0.921
module load java/1.8.0

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed using samtools index and reference file indexed using samtools faidx 
REF=/data/zool-zir/Myzomela/Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
BAMs=bam.list

#Set minumum number individuals needed for ANGSD to output GLs and GPs for a site
#(half total number is sensible)
MIN_INDV=40

# Create direcories for output (ignored (with warning) if already exists)
mkdir MAFs VCFs BEAGLEs BEAGLE_LOGs
mkdir VCFs/RAW
mkdir VCFs/IMPUTED
mkdir VCFs/PHASED

# Use slurm array task ID to get list of scaffold names
SCAFFOLDS=scaffold_lists/scaffold.list.$(printf '%02d' $SLURM_ARRAY_TASK_ID)
END=$(cat $SCAFFOLDS | wc -l)

# For each line within scaffold list do the following
for i in $(eval echo "{1..$END}")
do

# Extact scaffold ID
REGION_2_LOOKUP=$(cat $SCAFFOLDS | head -n $i | tail -n 1)

# Estimate GLs and GPs using ANGSD
# Note: No attempt has been made to call an individuals genotype, instead these will be imputed later
angsd -b $BAMs -ref $REF \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -minMapQ 20 -minQ 20 \
-GL 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -doVcf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6 -minInd $MIN_INDV \
-r $REGION_2_LOOKUP -out Lcass.v1_${REGION_2_LOOKUP}

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes reads with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -baq = perform BAQ computation (1: same as in SAMtools)
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -GL = calculate genotype likelihoods (1: using SAMtools model)
# -doMajorMinor = infer major and minor alleles (1: from GLs)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doVcf = output a VCF file (1: yes)
# -doGlf = ouput genotype likelihoods file (4: beagle likelihood format)
# -minMaf = minumum minor allele frequency tolerated
# -SNP_pval = significance threshold for determining true polymorphism
# -minInd = only output sites with information for atleast [int] samples (half is sensible)
#           Although we will impute missing data, this means we will only impute sites where we
#			have a good idea what is going on.
# ==============================

# Use beagle4.1 to impute and phase genotypes
java -Xmx15000m -jar beagle.27Jan18.7e1.jar \
gl=Lcass.v1_${REGION_2_LOOKUP}.vcf.gz out=Lcass.v1_${REGION_2_LOOKUP}.imputed

java -Xmx15000m -jar beagle.27Jan18.7e1.jar \
gt=Lcass.v1_${REGION_2_LOOKUP}.imputed.vcf.gz out=Lcass.v1_${REGION_2_LOOKUP}.imputed.phased

# Perform some tidying up
mv Lcass.v1_${REGION_2_LOOKUP}.mafs.gz MAFs/
mv Lcass.v1_${REGION_2_LOOKUP}.vcf.gz VCFs/RAW/
mv Lcass.v1_${REGION_2_LOOKUP}.beagle.gz BEAGLEs/
mv Lcass.v1_${REGION_2_LOOKUP}.imputed.vcf.gz VCFs/PHASED/
mv Lcass.v1_${REGION_2_LOOKUP}.imputed.phased.vcf.gz VCFs/IMPUTED/
mv Lcass.v1_${REGION_2_LOOKUP}.imputed*.log BEAGLE_LOGs/
rm Lcass.v1_${REGION_2_LOOKUP}.arg

done
