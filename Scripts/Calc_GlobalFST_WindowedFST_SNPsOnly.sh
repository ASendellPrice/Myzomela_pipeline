#!/bin/bash

###############################################################################################################
# CALCULATE FST FROM SFS USING ANGSD AND REALSFS
# INPUT: Population bam lists plus list of psuedochromosomes
# OUTPUT: snp by snp, global, and windowed fst estimates per chromosome
# A.Sendell-Price, Aug 2021
###############################################################################################################

#Load angsd environment
conda activate ANGSD

#Define window and step size to use
WINDOW_SIZE=50000
STEP_SIZE=50000

#Set population bam lists and reference assembly
#You will need one bam list per population
POP1_BAMs=pop1.bam.list
POP2_BAMs=pop2.bam.list
REF=/data/Users/Sonya_Myzomela/Lcass_2_Tgutt_ZW/Lcass_2_Tgutt_ZW.fasta

for CHROM in $(cat /data/Users/Sonya_Myzomela/chr.short.names | grep -v "W")
do
	
  #Create directory for chromosome and move into it
	mkdir ${CHROM}
	cd ${CHROM}

	#Combined bam lists
	cat ../${POP1_BAMs} > Pop1.Pop2.bam.list
	cat ../${POP2_BAMs} >> Pop1.Pop2.bam.list

	#Run angsd to find variable sites
	FILTERS="-remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 30 -minQ 20 -minInd 10 -baq 1"
	OPT="-doMaf 1 -doMajorMinor 1 -doPost 2 -GL 2 -minMaf 0.05 -SNP_pval 1e-6 "
	angsd -bam Pop1.Pop2.bam.list -ref $REF $FILTERS $OPT -r $CHROM -out ${CHROM}.pop1.pop2.snps

	#Create list of variable sites by extracting columns 1 and 2 from the ".mafs.gz" file
	zcat ${CHROM}.pop1.pop2.snps.mafs.gz | tail -n +2 | cut -f 1,2 > ${CHROM}.pop1.pop2.snps

	#Index list of variable sites
	angsd sites index ${CHROM}.pop1.pop2.snps

	#Calculate the allele frequency likelihoods for pop1 and pop2 using angsd
	#Note: As we dont have an ancestral state reference we will instead supply the reference assembly (this means we will 
	#need to "fold" the 2DSFS when we generate it - i.e. use the minor allele in lieu of the derived allele).
	FILTERS="-remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 30 -minQ 20 -baq 1"
	OPT="-dosaf 1 -gl 2"
	angsd -bam pop1.bam.list -ref $REF -anc $REF $FILTERS $OPT -r $CHROM -sites ${CHROM}.pop1.pop2.snps -out ${CHROM}.pop1
	angsd -bam pop2.bam.list -ref $REF -anc $REF $FILTERS $OPT -r $CHROM -sites ${CHROM}.pop1.pop2.snps -out ${CHROM}.pop2

	#Calculate the 2D site frequency spectrum
	#-fold 1: tells realSFS to output the "folded" site frequency spectrum
	realSFS -P 20 ${CHROM}.pop1.saf.idx ${CHROM}.pop2.saf.idx -fold 1 > ${CHROM}.pop1.pop2.folded.sfs

	#Calculate FST on a per site basis
	realSFS fst index ${CHROM}.pop1.saf.idx ${CHROM}.pop2.saf.idx -sfs ${CHROM}.pop1.pop2.folded.sfs -fstout ${CHROM}.pop1.pop2

	#Get the global FST estimate for chromosome and save output to file
	realSFS fst stats ${CHROM}.pop1.pop2.fst.idx > ${CHROM}.pop1.pop2.global.fst

	#Estimate FST in a sliding window
	#(-type 2 sets the left most position of the first window to 1 i.e. the first bp of the chromosome)
	realSFS fst stats2 ${CHROM}.pop1.pop2.fst.idx -win $WINDOW_SIZE -step $STEP_SIZE -type 2 > ${CHROM}.pop1.pop2.fst.size${WINDOW_SIZE}_step${STEP_SIZE}

done

