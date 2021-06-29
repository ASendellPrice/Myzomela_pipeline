#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=0-47:1
#SBATCH --time=120:00:00
#SBATCH --job-name=Merge_Beagles
#SBATCH --output=Merge_Beagles.log
#SBATCH --error=Merge_Beagles.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

######################################################################################
# MERGE SINGLE SCAFFOLD BEAGLES / VCFS INTO SINGLE WHOLE GENOME BEAGLES / VCFS
# A. Sendell-Price, June 2021
######################################################################################

#STEP 1: Merge BEAGLE files
#Create list of beagle files to merge
ls BEAGLEs/*.beagle.gz > beagle.list

#Extract BEAGLE header and start file
BEAGLE1=$(head -n 1 beagle.list)
zcat $BEAGLE1 | head -n 1 > Lcass.v1_WholeGenome.beagle

#Extract content from each scaffold bagle and append to Lcass.v1_WholeGenome.beagle
for SCAFFOLD in $(cat beagle.list)
do
	#tail -n +2 means print from second line (missing out the header)
	zcat $SCAFFOLD | tail -n +2 >> Lcass.v1_WholeGenome.beagle
done

#Compress the file and remove beagle list
zcat Lcass.v1_WholeGenome.beagle
rm beagle.list


#STEP 2: Merge VCF files
#Create list of VCF files to merge
ls VCFs/PHASED/*.vcf.gz > vcf.list

#Extract BEAGLE header and start file
VCF1=$(head -n 1 vcf.list)
zcat $VCF1 | grep "#" > Lcass.v1_WholeGenome.imputed.phased.vcf

#Extract content from each scaffold bagle and append to Lcass.v1_WholeGenome.beagle
for SCAFFOLD in $(cat vcf.list)
do
	#"grep "#" -v" is a reverse pattern search where we omit header lines (which start with a hash)
	zcat $SCAFFOLD | grep "#" -v >> Lcass.v1_WholeGenome.imputed.phased.vcf
done

#Compress the file and remove vcf list
zcat Lcass.v1_WholeGenome.imputed.phased.vcf
rm vcf.list
