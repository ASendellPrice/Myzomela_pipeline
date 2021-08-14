#!/bin/bash

###############################################################################################################
# CALCULATE FST FROM SFS USING ANGSD AND REALSFS
# INPUT: Population bam lists plus list of psuedochromosomes
# OUTPUT: snp by snp, global, and windowed fst estimates per chromosome
# A.Sendell-Price, Aug 2021
###############################################################################################################

#Define paths to angsd and realSFS
ANGSD=/data/Users/BIN/angsd/angsd
REAL_SFS=/data/Users/BIN/angsd/misc/realSFS

#Define window and step size to use
WINDOW_SIZE=50000
STEP_SIZE=50000

#Define path to list of chromosomes (needs to match chrom names in reference assembly)
CHROM_LIST=/data/Users/Sonya_Myzomela/Ref_Genome/PseudoChroms.txt

#For each line in CHROM_LIST do the following . . .
for LINE in {1..33}
do
	#Get long chromosome name (e.g. "PseudoCM012113.1_Taeniopygia_guttata_isolate_Black17_chromosome_Z,_whole_genome_shotgun_sequence")
	CHROM=$(head -n $LINE $CHROM_LIST | tail -n 1)

	#Get short chromosome name 
	#(e.g. "Z" rather than "PseudoCM012113.1_Taeniopygia_guttata_isolate_Black17_chromosome_Z,_whole_genome_shotgun_sequence")
	CHROM_SHORT=$(echo $CHROM | cut -d "," -f 1 | cut -d "_" -f 7)

	#Make a directory for the chromosome and move into it
	mkdir chr${CHROM_SHORT}
	cd chr${CHROM_SHORT}

	#Calculate the allele frequency likelihoods for pop1 and pop2 using angsd
	#Note: As we dont have an ancestral state reference we will instead supply the reference assembly (this means we will 
	#need to "fold" the 2DSFS when we generate it - i.e. use the minor allele in lieu of the derived allele).
	REF=/data/Users/Sonya_Myzomela/Ref_Genome/Lcass_2_Tgut_pseudochromosomes.fasta.gz
	$ANGSD -bam ../pop1.bam.list -ref $REF -anc $REF -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop1.chr${CHROM_SHORT}
	$ANGSD -bam ../pop2.bam.list -ref $REF -anc $REF -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop2.chr${CHROM_SHORT}

	#Calculate the 2D site frequency spectrum
	#-fold 1: 			tells realSFS to output the "folded" site frequency spectrum
	#-m 0: tells		tells realSFS to use the standard EM algorithm instead of the accelerated EM algorithm
	#-maxIter 50000: 	sets the maximum number of EM iterations to 50,000 instead of the default 100 (to provide sufficient iterations for covergence)
	#-tole 1e-6: 		when the difference in successive likelihood values in the EM algorithm gets below this value the optimisation will stop.
	$REAL_SFS pop1.chr${CHROM_SHORT}.saf.idx pop2.chr${CHROM_SHORT}.saf.idx -fold 1 -m 0 -maxIter 50000 -tole 1e-4 -P 20 > pop1.pop2.chr${CHROM_SHORT}.folded.sfs
	
	#Calculate FST on a per site basis
	$REAL_SFS fst index pop1.chr${CHROM_SHORT}.saf.idx pop2.chr${CHROM_SHORT}.saf.idx -sfs pop1.pop2.chr${CHROM_SHORT}.folded.sfs -fstout pop1.pop2.chr${CHROM_SHORT}

	#Get the global FST estimate for chromosome and save output to file
	$REAL_SFS fst stats pop1.pop2.chr${CHROM_SHORT}.fst.idx > chr${CHROM_SHORT}.global.fst

	#Estimate FST in a sliding window
	#(-type 2 sets the left most position of the first window to 1 i.e. the first bp of the chromosome)
	$REAL_SFS fst stats2 pop1.pop2.chr${CHROM_SHORT}.fst.idx -win $WINDOW_SIZE -step $STEP_SIZE -type 2 > pop1.pop2.chr${CHROM_SHORT}.fst.size${WINDOW_SIZE}_step${STEP_SIZE}

	cd ../
done
