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
FASTP=/data/zool-zir/Myzomela/BIN/fastp

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
# Load required modules in ARC
module load java/1.8.0
module load samtools
module load bowtie2

# STEP 2:
# Specify "plate" name for sequencing run
PLATE=Novogene_WGS_2020

# STEP 3:
# Specify path to reference genome (cornetti assembly) including file prefix "ZOLAv0"
GENOME_DB=/data/zool-zir/Myzomela/Ref_Genome/L_cass

# STEP 4:
# Move into directory containing sample's filtered reads
cd /data/zool-zir/Myzomela/filtered_reads/${SAMPLE_NAME}

# STEP 5:
# For each pair of reads conduct mapping using bowtie2
for ReadPair in $(ls Filtered_${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5 -d'_')
do
  	bowtie2 -x $GENOME_DB \
    -1 ${ReadPair}_1.fq.gz \
    -2 ${ReadPair}_2.fq.gz \
    --rg-id $ReadPair --rg SM:$SAMPLE_NAME --rg LB:$PLATE \
    --rg PU:$PLATE --rg PL:illumina 2> ${ReadPair}.log  | \
    samtools view -bhS - | \
    samtools sort - > ${ReadPair}.bam
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
