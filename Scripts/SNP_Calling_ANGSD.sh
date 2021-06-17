#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=0-47:1
#SBATCH --time=120:00:00
#SBATCH --job-name=ANGSD_SNP_Calling
#SBATCH --output=ANGSD_SNP_Calling.log
#SBATCH --error=ANGSD_SNP_Calling.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sonya.clegg@zoo.ox.ac.uk

# Load ANGSD
module load angsd/0.921

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed (using samtools index) 
REF=/data/zool-zir/Myzomela/Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
BAMs=bam.list

#Use slurm array task ID to get list of scaffold names
SCAFFOLDS=scaffold_lists/scaffold.list.$(printf '%02d' $SLURM_ARRAY_TASK_ID)
END=$(cat $SCAFFOLDS | wc -l)

for i in $(eval echo "{1..$END}")
do 

REGION_2_LOOKUP=$(cat $SCAFFOLDS | head -n $i | tail -n 1)

#Call SNPs using ANGSD
angsd -b $BAMs -ref $REF \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -geno_minDepth 5 -checkBamHeaders 0 \
-doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 4 -doVcf 1 -minMaf 0.01 -SNP_pval 1e-6 \
-r $REGION_2_LOOKUP -out Lcass.v1_${REGION_2_LOOKUP}

mv Lcass.v1_${REGION_2_LOOKUP}.geno.gz GENOs/
mv Lcass.v1_${REGION_2_LOOKUP}.mafs.gz MAFs/
mv Lcass.v1_${REGION_2_LOOKUP}.vcf.gz VCFs/
rm Lcass.v1_${REGION_2_LOOKUP}.arg

done
