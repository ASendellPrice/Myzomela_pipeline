#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=62:00:00
#SBATCH --array=1-150:1
#SBATCH --job-name=Filter2Bam_Pipeline
#SBATCH --output=Filter2Bam.%A_%a.out
#SBATCH --error=Filter2Bam.%A_%a.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

#######################################################################################################
# SET SAMPLE NAME AND DIRECTORY
# This pipeline requires the number of jobs to match the number of lines
# in file sample.list.txt (update #SBATCH --array=1-65:1)
#######################################################################################################

# STEP 1:
# Specify file containing samples to be processed
Sample_List=samples.txt

# STEP 2:
# Use slurm array task ID to alocate sample name and directory
SAMPLE_NAME=$(head -n $SLURM_ARRAY_TASK_ID $Sample_List | tail -n 1)
SAMPLE_DIRECTORY=raw_data/${SAMPLE_NAME}
#Create directories for output
mkdir filtered_reads
mkdir filtered_reads/${SAMPLE_NAME}
mkdir fastp_QC_reports
mkdir sample_bams

# STEP 3:
#move into sample directory
cd $SAMPLE_DIRECTORY


#######################################################################################################
# CONDUCT FILTERING OF WGS RAW READS
# We will use FastP - a tool designed to provide fast all-in-one preprocessing for FastQ files
# https://github.com/OpenGene/fastp
#######################################################################################################

# STEP 1:
# Define path to Fastp:
FASTP=/data/Users/Sonya_Myzomela/BIN/fastp

# STEP 2:
# Set up for loop to conduct filtering for each read pair
for ReadPair in $(ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4 -d'_')
do

  #Use Fastp to conduct automated filtering of fastq files
  #Note: based on initial test we will trim the first 10bp from start of each read
  $FASTP \
  -i ${ReadPair}_1.fq.gz \
  -o ../../filtered_reads/${SAMPLE_NAME}/Filtered_${ReadPair}_1.fq.gz \
  -I ${ReadPair}_2.fq.gz \
  -O ../../filtered_reads/${SAMPLE_NAME}/Filtered_${ReadPair}_2.fq.gz \
  --trim_front1 10 \
  --trim_front2 10 \
  --detect_adapter_for_pe

  #Remove .json file as not needed
  rm fastp.json

  #Rename QC report and move to fastp qc report folder
  mv fastp.html ../../fastp_QC_reports/${ReadPair}.html

done


#######################################################################################################
# MAP FILTERED READS TO REFERENCE ASSEMBLY
#######################################################################################################

# STEP 1:
# Define path to BWA and reference assembly
# Note: If not already done, will need to index reference assembly ($BWA index $REF)

BWA=/data/Users/Sonya_Myzomela/BIN/bwa/bwa
REF=/data/Users/Sonya_Myzomela/Lcass_2_Tgutt_ZW/Lcass_2_Tgutt_ZW.fasta

# STEP 2:
# For each pair of reads conduct mapping using BWA MEM
for ReadPair in `ls Filtered_${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5 -d'_'`
do
  	$BWA mem $REF \
  	${ReadPair}_1.fq.gz \
  	${ReadPair}_2.fq.gz	\
  	-R "@RG\tID:${ReadPair}\tSM:${SAMPLE_NAME}" \
  	| samtools view -bS - \
  	| samtools sort - > ${ReadPair}.bam
done

#######################################################################################################
# MERGE SAMPLE BAMS INTO A SINGLE FILE
# Currently we have a bam file for each read pair per sample. We will now merge these to create a
# single bam file per sample
#######################################################################################################

# Count number of bam files per sample
BAM_Count=$(ls Filtered_${SAMPLE_NAME}*.bam | wc -l)

# If number of bams is greater than 1 then merge bams into a single file using samtools merge
# else, rename single bam file
if [ $BAM_Count -gt 1 ]
then
  samtools merge ${SAMPLE_NAME}.bam Filtered_${SAMPLE_NAME}*.bam -f
  rm Filtered_${SAMPLE_NAME}*.bam
else
  mv Filtered_${SAMPLE_NAME}*.bam ${SAMPLE_NAME}.bam
fi

#Sort bam file and remove unsorted version
samtools sort ${SAMPLE_NAME}.bam > ${SAMPLE_NAME}.sorted.bam
rm ${SAMPLE_NAME}.bam

#Move and index sorted bam file to bam directory
mv ${SAMPLE_NAME}.sorted.bam /data/zool-zir/Myzomela/sample_bams/${SAMPLE_NAME}.sorted.bam
samtools index /data/zool-zir/Myzomela/sample_bams/${SAMPLE_NAME}.sorted.bam
