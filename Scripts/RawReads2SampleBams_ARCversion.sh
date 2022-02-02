#!/bin/bash
#SBATCH --clusters=all
#SBATCH --nodes=1
#SBATCH --array=1-81:1
#SBATCH --time=2-00:00:00 
#SBATCH --job-name=Bam2Fasta
#SBATCH --partition=long
#SBATCH --output=Bam2Fasta_%a.log
#SBATCH --error=Bam2Fasta_%a.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

# Submit like so:
# source RawReads2SampleBams.sh SAMPLE_LIST REF_PATH

SAMPLE_NAME=$(head -n $SLURM_ARRAY_TASK_ID $1 | tail -n 1)

#######################################################################################################
# SET RAW DATA PATH AND CREATE DIRECTORIES
#######################################################################################################

#Set path to directory containing raw reads for sample
SAMPLE_DIRECTORY=Non-zosterops/${SAMPLE_NAME}
mkdir filtered_reads_myzomela
mkdir filtered_reads_myzomela/${SAMPLE_NAME}
mkdir fastp_QC_reports_myzomela
mkdir sample_bams_myzomela

#move into sample directory
cd $SAMPLE_DIRECTORY

#######################################################################################################
# CONDUCT FILTERING OF WGS RAW READS
# We will use FastP - a tool designed to provide fast all-in-one preprocessing for FastQ files
# https://github.com/OpenGene/fastp
#######################################################################################################
 
# Define path to Fastp:
FASTP=/data/zool-zost/BIN/fastp

# Set up for loop to conduct filtering for each read pair
for ReadPair in $(ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4 -d'_')
do
  #Use Fastp to conduct automated filtering of fastq files
  #Note: based on initial test we will trim the first 10bp from start of each read
  $FASTP \
  -i ${ReadPair}_1.fq.gz \
  -o ../../filtered_reads_myzomela/${SAMPLE_NAME}/Filtered_${ReadPair}_1.fq.gz \
  -I ${ReadPair}_2.fq.gz \
  -O ../../filtered_reads_myzomela/${SAMPLE_NAME}/Filtered_${ReadPair}_2.fq.gz \
  --trim_front1 10 \
  --trim_front2 10 \
  --detect_adapter_for_pe

  #Remove .json file as not needed
  rm fastp.json

  #Rename QC report and move to fastp qc report folder
  mv fastp.html ../../fastp_QC_reports_myzomela/${ReadPair}.html

done

#######################################################################################################
# MAP FILTERED READS TO REFERENCE ASSEMBLY
#######################################################################################################

#Move to directory containing filtered reads
cd ../../filtered_reads_myzomela/${SAMPLE_NAME}/
  
# Define path to BWA and reference assembly
# Note: If not already done, will need to index reference assembly ($BWA index $REF)
BWA=/data/zool-zost/BIN/bwa/bwa
REF=$2
  
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

######################################################################################################
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
mv ${SAMPLE_NAME}.sorted.bam ../../sample_bams_myzomela/${SAMPLE_NAME}.sorted.bam
samtools index ../../sample_bams_myzomela/${SAMPLE_NAME}.sorted.bam

