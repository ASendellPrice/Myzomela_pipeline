# Step by step bioinformatics guide for Sonya's *Myzomela* project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## Creating pseudo-chromosome assembly based on synteny with T. guttata assembly.

In our previous work we had used the program [**Satsuma**](http://bioinformatics.oxfordjournals.org/content/26/9/1145.long) to infer the mean position of Z. lateralis scaffolds (and their orientation) relative to chromosomes of the T. guttata assembly. This was rather long-winded and involved complicated R-scripts to translate Z. lateralis positions within a VCF/Beagle file into T. guttata chromosomal positions. For large VCF/Beagle files this creates a complete headache as is difficult to allocate the memory required to do this . . . 

However, the developers of Satsuma have recently introduced a new tool: the Chromosembler! Which makes this process a whole lot simpler. Chromosembler takes two fasta files: 
1. query sequence - the helmeted honeyeater reference (scaffold-level assembly)
2. target sequence - the zebra finch reference (chromosome level-assembly)

and outputs two fasta files:
1. A super-scaffold assembly where query scaffolds have been joined into larger and fewer pieces
2. A pseudo-chromosome assembly where super scaffolds have been combined into "chromosomes" and gaps filled with Ns.

The outputted pseudo-chromosome assembly can then be used as the reference when mapping filtered sequencing reads (as in Step 4)

### Let's try running Chromosembler on Nesoi

STEP 1: Log into Nesoi and move to the data directory '/data/Users/'

STEP 2: Create a directory for your work (something like 'Sonya_Myzomela') and move into it

STEP 3: Create a directory for Satsuma work (something like 'Satsuma_Chromosembler_Run') and move into it

STEP 4: Create a new tmux session called 'Satsuma'
```
tmux new -s Satsuma
```

STEP 5: Download the query and target genomes and then unzip them
```
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz

gunzip GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz
gunzip GCA_008360975.1_HeHo_1.0_genomic.fna.gz
```

STEP 6: Set path to Chromosemble which I have installed on Nesoi
```
CHROMOSEMBLE=/data/Users/BIN/satsuma2/bin/Chromosemble
export SATSUMA2_PATH=/data/Users/BIN/satsuma2/bin/
```

STEP 7: Run Chromosemble
```
$CHROMOSEMBLE \
-t GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
-q GCA_008360975.1_HeHo_1.0_genomic.fna \
-o Lcass_2_Tgutt_ZW
```
* -t = target sequence (zebra finch assembly)
* -q = query sequence (Honeyeater assembly)
* -o = name for the output directory

STEP 8: Detatch from session

This will take a few hours to run, so detatch from tmux session using cntrl+b followed by d

## Editing the chromosome names in the pseudochromosome assembly

Chromosemble will output a psuedoassembly in fasta format. This file consists of a header line specifying the chromosome (or scaffold) name followed by the sequence. Below is a screenshot showing the first few lines of the pseudoassembly.

![](https://github.com/ASendellPrice/Myzomela_pipeline/blob/main/head_assembly.png)

As you can see the pseudochromosome names are really long. Such long chromosome names can conflict with programmes used in downstream analyses, so we will shorten these so that they are formatted like so: "chr1".

STEP 1: Reconnect to Nesoi, and reconnect to the tmux session

STEP 2: Move into the directory containing output from Chromosemble and rename pseudochromosomes.fasta to something sensible
```
cd Lcass_2_Tgutt_ZW
mv pseudochromosomes.fasta Lcass_2_Tgutt_ZW.fasta
```

STEP 3: Using the cat, grep and sed commands extract the chromosome names from the fasta file and save these to a text file we will call "long.chrom.names"
```
cat Lcass_2_Tgutt_ZW.fasta | grep "chromosome" | sed 's/>//g' > long.chrom.names
```
![](https://github.com/ASendellPrice/Myzomela_pipeline/blob/main/chromname_extraction.png)

STEP 4: For each long chromosome name in our list "long.chrom.names" replace the long chromosome name in the fasta file with a shortened version. This will be implemented in a for loop.

```
for CHROM in $(cat long.chrom.names)
do
    #Extract chromosome name / number from long name
    CHROM_SHORT=$(echo $CHROM | cut -d "," -f 1 | cut -d "_" -f 8)
    #Use sed to replace long chrom name with short chrom name 
    #This may take several minutes
    sed -i "s/$CHROM/chr$CHROM_SHORT/g" Lcass_2_Tgutt_ZW.fasta
done
```

![](https://github.com/ASendellPrice/Myzomela_pipeline/blob/main/create_short_chrom_name.png)

Let's check this has worked by extracting the header lines from the fasta file again. Short chromosome names should be printed to your screen.
```
cat Lcass_2_Tgutt_ZW.fasta | grep "chr"
```

STEP 6: Close tmux session (control +b followed by x)


# Filter and align myzomela reads to pseudoassembly

Our raw sequencing reads need to pass through a number of bioinformatic steps before we can call SNPs. The main steps are as follows:

* Filter raw sequencing reads to remove any sequencing adapter content/duplicates and remove bases of low quality - conducted using [fastp](https://github.com/OpenGene/fastp).
* Map filtered sequencing reads to our reference genome - conducted using [Burrows-Wheeler Aligner (BWA)](https://github.com/lh3/bwa).
* As most of the samples have been sequenced across multiple lanes we have multiple reads pairs and therefore multiple bams per sample. Where samples have multiple bams we will merge them using the aptly named "merge" command of [samtools](http://satsuma.sourceforge.net). Merged bams will then be indexed using the samtools "index" command.

All three of these steps are conducted using the pipeline (RawReads2SampleBams.sh) that I have developed for our Novogene sequencing data. I have added this script to the "scripts" directory in the Myzomela folder.

**Note:** This script can also be viewed on [gitHub](https://github.com/ASendellPrice/Myzomela/blob/main/Scripts/RawReads2SampleBams.sh).

This single script takes care of all the filtering, alignment and merging steps, but does have a few assumptions.

### ASSUMPTION 1:
Raw sequencing reads for an individual are saved in a directory called "sample id" and this directory sits within a parent directory called "raw_data" like so:
```
#Let's preview the structure of directory "raw_data" . . .
tree raw_data
```
```
raw_data/
├── ABY76
│   ├── ABY76_FDSW202513070-1r_HHGCJDSXY_L2_1.fq.gz
│   ├── ABY76_FDSW202513070-1r_HHGCJDSXY_L2_2.fq.gz
│   └── MD5.txt
├── AMB1
│   ├── AMB1_FDSW202513071-1r_HHGCJDSXY_L2_1.fq.gz
│   ├── AMB1_FDSW202513071-1r_HHGCJDSXY_L2_2.fq.gz
│   └── MD5.txt
└── AMB2
    ├── AMB2_FDSW202513073-1r_HCKW5DSXY_L3_1.fq.gz
    ├── AMB2_FDSW202513073-1r_HCKW5DSXY_L3_2.fq.gz
    ├── AMB2_FDSW202513073-1r_HHG3GDSXY_L3_1.fq.gz
    ├── AMB2_FDSW202513073-1r_HHG3GDSXY_L3_2.fq.gz
    └── MD5.txt
```

**NOTE:** As we have used paired-end sequencing, we have two files per lane of sequencing (Forward: "_1.fq.gz" and Reverse: "_2.fq.gz"). Together these are called a **read pair**. Here is some explanation of what each read file name means:

```
1    2                3         4  5 6
AMB2_FDSW202513073-1r_HHG3GDSXY_L3_1.fq.gz
```
1. AMB2 = Name of sample
2. FDSW202513073-1r = reference no. unique to Illumina HighSeq machine
3. HHG3GDSXY = reference no. of flow cell
4. L3 = Lane number on flow cell
5. Read direction (1 = forward, 2 = reverse)
6. fq.gz = compressed fastq format

**ALSO NOTE:** Not all samples have the same number of read pairs present in their raw_data directories. Here ABY76 and AMB1 both have a single pair of sequencing reads, whereas AMB2 has two. This is because AMB2 has been sequenced across multiple (2) sequencing lanes (likely this sample was of low concentration and required twice the sequencing effort to produce the data quantity we requested).

### ASSUMPTION 2:
The working directory contains the text file "samples.txt" which lists the samples we are going to work with. This can be quickly generated like so:
```
ls raw_data > samples.txt
```

### ASSUMPTION 3:
Our reference genome has been indexed using "bwa index". This can be done by running the following command from the directory containing the psuedo-chromosome assembly:
```
/data/Users/Sonya_Myzomela/BIN/bwa/bwa index Lcass_2_Tgutt_ZW.fasta
```
Note: this can take a while to index, so best to do from a tmux session

### RUNNING THE SCRIPT:

If the three assumptions are met we can run the script like so:
```
source Scripts/RawReads2SampleBams.sh samples.txt /data/Users/Sonya_Myzomela/Lcass_2_Tgutt_ZW/Lcass_2_Tgutt_ZW.fasta
```


## ESTIMATING FST FROM GENOTYPE LIKELIHOODS
The following outlines how you can estimate fst (windowed and global) between a pair of populations. Note: You will need to decide which samples will form the populations you wish to contrast, these could be geographical comparisons e.g. Tanna vs. Mainland, or could be based on sex e.g. Vanuatu males vs Vanuatu females. This step by step tutorial is based on the FST pages of ANGSD, see [link](http://www.popgen.dk/angsd/index.php/Fst).

STEP 1: Log into nesoi and create a new detachable session (we have done this before so I will leave it to you to figure out)

STEP 2: Move into the following directory '/data/Users/Sonya_Myzomela/' and create a new directory for the fst analysis and move into it (call it something informative). 
```
mkdir Fst_pop1_vs_pop2
cd Fst_pop1_vs_pop2
```

STEP 3: Create a list of bam files for each sample within a given population. I have transfered bam files from arc to nesoi, they are stored here: '/data/Users/Sonya_Myzomela/sample_bams'. 

1. You can get a list of all the sample bams available like so:
```
ls /data/Users/Sonya_Myzomela/sample_bams_ARC/*.bam
```

2. Create your list file for population 1 using nano
```
nano pop1.bam.list
```

3. Now do the same for population 2

STEP 4: Run script "Calc_GlobalFST_WindowedFST.sh" which will do the following for each chromosome:
1. Calculate the allele frequency likelihoods for pop1 and pop2 using angsd
2. Calculate the 2D "folded" site frequency spectrum
3. Calculate FST on a per site basis
4. Calculate global (whole chromosome) estimate of FST
5. Calculate FST in windows (currently set to non-overlapping 50kb windows)
View script online [here](https://github.com/ASendellPrice/Myzomela_pipeline/blob/main/Scripts/Calc_GlobalFST_WindowedFST.sh).

```
source ../Scripts/Calc_GlobalFST_WindowedFST.sh
```

##STOP HERE!









## STEP 7: Estimating genotype likelihoods (GLs) and imputing genotypes (GTs)
A challenge of working with low coverage sequencing data is that we cannot be 100% certain of the genotypes (GTs) of our samples. A solution to this challenge is to use a probabilistic measurement of the genotypes in the form of genotype likelihoods (GLs) and/or genotype probabilities (PLs). There are an increasing number of tools that take genotype likelihoods/probabilities as input for downstream analyses, however where traditional genotype inputs are reruired GTs can be imputed from GLs/PLs based on local linkage patterns. 

In this step we will estimate genotype likelihoods using the software [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) and impute sample genotypes using [BEAGLE4](https://faculty.washington.edu/browning/beagle/b4_0.html). All of this is performed using the script "GL_Estimation_GT_Imputation.sh".

Before we can submit the script we will need to create a list of sample bam files. From '/data/zool-zir/Myzomela' type the following commands:

```
#Create list of sample bams
ls sample_bams/*.bam > bam.list
```

Let's look at the script 'GL_Estimation_GT_Imputation.sh' using the command less (used up and down arrows to scroll through script and hit 'q' key to close). As you will see the script is fully annotated and should (hopefully) be understandable. To speed things up processing of the 1000+ scaffolds is conducted across a multiple job array.
```
less GL_Estimation_GT_Imputation.sh
```

We can submit this script with the simple command:
```
sbatch GL_Estimation_GT_Imputation.sh
```









## STEP 11: Conducting PCA and admixture analyses
To get an idea of the genetic structuring of the *Myzomela* samples we will conduct two different analyses:
1. Principal Component Analysis using [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd)
2. Admixture using [NgsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix)

The computational parts of these analyses are conducted using the script **Run_PCAngsd_NgsAdmix.sh**. Let's have a look at this file using the less command
```
less Myzomela_pipeline/Scripts/Run_PCAngsd_NgsAdmix.sh
```
We can run the script like so:
```
source Myzomela_pipeline/Scripts/Run_PCAngsd_NgsAdmix.sh
```
After a number of hours the script will produce the following files:
**Within PCA directory:**
1. Myzomela_Indv81_SNPs5M.cov - covariance matrix which can be used to conduct PCA in R
2. Myzomela_Indv81_SNPs5M.log - log file from PCAngsd that summarises run

**Within ADMIX directory:**
One of each per K (K2 as an example)
1. Myzomela_Indv81_SNPs5M_K2.filter 
2. Myzomela_Indv81_SNPs5M_K2.fopt.gz - an estimate of the allele frequency in each of the 3 assumed ancestral populations. There is a line for each locus.
3. Myzomela_Indv81_SNPs5M_K2.log - log file from PCAngsd that summarises run
4. Myzomela_Indv81_SNPs5M_K2.qopt - contains an estimate of the individual's ancestry proportion (admixture) from each of the three assumed ancestral populations for all individuals. There is a line for each individual.

**Bonus file:**
1. sampleIDs.txt - file listing samples within beagle file

These files can be transfered from Nesoi to your laptop via Filezilla

### Using R to conduct PCA from PCAngsd covariance matrix
```r
#Read in covariance matrix for chunk
C <- as.matrix(read.table("Myzomela_Indv81_SNPs5M.cov"))
  
#Compute eigenvalues and eigenvectors
e <- eigen(C)
  
#Extract PCs 1 and 2 
PC1_PC2 <- as.data.frame(e$vectors[,1:2])
colnames(PC1_PC2) <- c("PC1","PC2")

write.table(PC1_PC2, "PCs.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#Load sample names
sample_ids <- read.delim("samplesIDs.txt", header = FALSE)
colnames(sample_ids) <- "sample.id"

#Append sample ID column to PC dataframe
PC1_PC2 <- cbind(sample_ids, PC1_PC2)

#You can also add other info to each sample e.g. population, species/subspecies etc.
#To do so create simple text files containing info for each sample.
```



