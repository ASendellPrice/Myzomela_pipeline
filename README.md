# Step by step bioinformatics guide for Sonya's *Myzomela* project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## STEP 1: Creating pseudo-chromosome assembly based on synteny with T. guttata assembly.

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

# <---- this is where we stopped

STEP 9: Shorten chromsome names in pseudochromosome assembly fasta file

The pseudochromosome assembly outputted by Chromosemble will have really long chromosome names e.g. "PseudoCM012081.1_Taeniopygia_guttata_isolate_Black17_chromosome_1,_whole_genome_shotgun_sequence", we will shorten these so that they are formatted like so: "chr1"

#Extract chromosome names from pseudochromosome assembly fasta file and save these to a text file called "long.chrom.names"
```
cat *.fna | grep ">" | awk 'g/>//'s > long.chrom.names
#   1          2           3        4 5 
```
1. read file with cat command and pipe output directly to the grep command
2. Using grep command search for lines containing ">" - as each chromosome header in the fasta file is initiated with ">" by searching for this symbol we can pull out all the lines containing chromosome headers, e.g (">PseudoCM012081.1_Taeniopygia_guttata_isolate_Black17_chromosome_1,_whole_genome_shotgun_sequence") 

4. Pipe these to the awk command
5. 

3. search for lines containing ">" using the grep command - each chromosome in the fasta file starts with a header line e.g. (>PseudoCM012081.1_Taeniopygia_guttata_isolate_Black17_chromosome_1,_whole_genome_shotgun_sequence) so we can use ">" to seach for these.
4. using awk replace ">" with "" (nothing)  






## STEP 2: Transfer psuedochromosome assembly to ARC
Due to memory demands we will need to conduct the read processing (filtering and mapping) using the ARC HTC cluster
**Note: You will need to replace SSO with your username, you will also be prompted for your ARC password**
```
rsync -S -av Lcass_2_Tgutt_ZW SSO@htc-login.arc.ox.ac.uk:/data/zool-zir/Myzomela/
```





## STEP 1: Log in to ARC move into *Myzomela* project directory
Let's log into Arcus-htc (replace SSO with your single sign on):
```
ssh SSO@htc-login.arc.ox.ac.uk
```

Move into Myzomela directory in "/data/zool-zir/":
```
cd /data/zool-zir/Myzomela
```

## STEP 2: Download reference genome and build reference dictionary
Although there are currently no publically available genomes for *Myzomela* species there are three available for species in the wider *Meliphagidae family*. These include:
* *Lichenostomus cassidix* (helmeted honeyeater) [GCA_008360975.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_008360975.1/)
* *Grantiella picta* (painted honeyeater) [GCA_013397955.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_013397955.1)
* *Anthochaera phrygia* (Regent honeyeater) [GCA_009430485.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_009430485.1)

All of these are scaffold level assemblies i.e. not assembled to chromosome level. This isn't too bad as we can find the position of scaffolds relative to a chromosome level assembly later using [Satsuma](http://satsuma.sourceforge.net). More about that later.

We will be using the *Lichenostomus cassidix* as a reference genome as this is the least fragmented of the genome available, having the lowest number of scaffolds (1,912). In contrast *Grantiella picta* and *Anthochaera phrygia* have 15k and 265k scaffolds respectively.   

So now we need to download the assembly from NCBI using command "wget" and save the file in a new directory called "Ref_Genome".

```
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0//GCA_008360975.1_HeHo_1.0_genomic.fna.gz -P Ref_Genome/
```

The downloaded assembly has been compressed using gzip, so let's unzip this
```
gunzip Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
```
When it comes to aligning sequencing reads to the reference genome we will need to have created an index for our reference genome. In fact we need to create two types of reference, one that can be used by bowtie2 (used when aligning sequencing reads to reference genome) and one which can be used by samtools (used when sorting bam files a.k.a aligned reads).

**Indexing** a genome is analogous to indexing a book. If you want to know on which page a certain word appears or a chapter begins, it is much more efficient/faster to look it up in a pre-built index than going through every page of the book until you found it. Same goes for alignments. Indices allow the aligner (Bowtie2) to narrow down the potential origin of a query sequence within the genome, saving both time and memory.

Create bowtie2 index:
This will produce a number of files with prefix "L_cass" which stands for *Lichenostomus cassidix*.
```
module load bowtie2
bowtie2-build Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna Ref_Genome/L_cass
```

Create genome index using samtools [faidx](http://www.htslib.org/doc/faidx.html):
This will produce two files with suffixes ".fai" and ".gzi".
Irritatingly, requires the assembly to be re-compressed (this time using bgzip).
```
module load samtools
bgzip Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna
samtools faidx Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
```

## STEP 3: Create sample list
Before we dive into the bioinformatics pipeline, we will need to create a list of samples we want to use in the planned analyses. Some analyses (e.g. introgression, SFS, phylogeny building) require use of an outgroup, so if you plan to conduct these you will need to identify publically available WGS data from a suitable outgroup (or we can use other *Meliphagidae* samples included in the Novogene dataset).

The sample list you create needs to be a simple text file (called samples.txt) with one sample ID per line. Like below:

```
ABY15
ABY76
AMB1
```
This can be created with nano (a command line text editor). To close nano, press control + x, hit y for "yes" and hit enter (to confirm you want to save).

```
nano samples.txt
```

## STEP 4: Filter raw sequencing reads, map filtered reads to reference genome, and create a single bam file per sample

Our raw sequencing reads need to pass through a number of bioinformatic steps before we can call SNPs. The main steps are as follows:

* Filter raw sequencing reads to remove any sequencing adapter content/duplicates and remove bases of low quality - conducted using [fastp](https://github.com/OpenGene/fastp).
* Map filtered sequencing reads to a reference genome of choice (in our case *Lichenostomus cassidix*) - conducted using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
* As most of the samples have been sequenced across multiple lanes we have multiple reads pairs and therefore multiple bams per sample. Where samples have multiple bams we will merge them using the aptly named "merge" command of [samtools](http://satsuma.sourceforge.net). Merged bams will then be indexed using the samtools "index" command.

All three of these steps are conducted using the pipeline (RawReads2SampleBams.sh) that I have developed for our Novogene sequencing data. This script is saved in a shared folder in project "zool-zir" and can be copied to our working directory like so:

```
#Create a directory where we will store scripts
mkdir Scripts
#Copy "RawReads2SampleBams.sh" from shared folder to "Scripts" directory
cp /data/zool-zir/Bioinformatics_pipelines/RawReads2SampleBams.sh Scripts/
```

**Note:** This script can also be viewed on [gitHub](https://github.com/ASendellPrice/Myzomela/blob/main/Scripts/RawReads2SampleBams.sh).

This single script will launch a job array (1 job per sample) that takes care of all the filtering, alignment and merging steps. This script has a couple of assumptions.

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
The working directory contains the text file "samples.txt" which lists the samples we are going to work with. We created this during step 3.

### ASSUMPTION 3:
The reference genome we want to use is stored in directory "Ref_Genome" and has been indexed with both samtools and bowtie2 (this was done in step 2).

### ASSUMPTION 4:
Line four of RawReads2SampleBams.sh is updated so that job array range ("--array=1-150:1") matches the number of samples in "samples.txt". This needs to be updated to "--array=1-65:1" (as there are 65 samples in "samples.txt" that we want to process). This can be updated using nano. Can you remember how to do this?

```
#This is what the first nine lines of the script look like:
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=62:00:00
#SBATCH --array=1-150:1
#SBATCH --job-name=Filter2Bam_Pipeline
#SBATCH --output=Filter2Bam.%A_%a.out
#SBATCH --error=Filter2Bam.%A_%a.error
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@zoo.ox.ac.uk

```

**Note:** Also need to update line nine so that notification emails are sent to the correct email address.

## STEP 5: Submitting the script

Running the pipeline is simply a case of submitting the following command on arcus-htc:

```
sbatch Scripts/RawReads2SampleBams.sh
```

## STEP 6: Downloading fastp QC reports
Open a new terminal and from a directory on your laptop where you want to save the fastp reports type the following command (changing SSO to your username):

```
scp -r SSO@arcus-b.arc.ox.ac.uk:/data/zool-zir/Myzomela/fastp_QC_reports/* ./
```

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

## STEP 8: Merging single scaffold files

To be added ...

## STEP 9: Logging into Nesoi and migrating datasets

Now that the big bioinformatics work has been done we can now move from ARC to Nesoi. Like ARC we will access the server using a secure shell. Type the following command (replace SSO with your single sign-on), after which you will be prompted for your password (different from ARC)
```
ssh SSO@zoo-nesoi.zoo.ox.ac.uk
```

Welcome to Nesoi goddesses of islands!! You will now be in your home directory "/home/zoo/" followed by your SSO. You can check this using the print working directory command (pwd)
```
pwd
```

To keep things tidy we will now create a directory for the Myzomela project and move into it
```
mkdir Myzomela
cd Myzomela
```

Copying whole genome beagle and VCF files from ARC (you will be prompted for your ARC password)
```
scp -r SSO@arcus-b.arc.ox.ac.uk:/data/zool-zir/Myzomela/EstimateGLsImputeGTs/Lcass.v1_WholeGenome.* ./
```

Copy bam.list previously used by ANGSD (again, you will be prompted for your ARC password)
```
scp -r SSO@arcus-b.arc.ox.ac.uk:/data/zool-zir/Myzomela/EstimateGLsImputeGTs/bam.list ./
```

Copy scripts from GitHub
```
git clone https://github.com/ASendellPrice/Myzomela_pipeline
```

## STEP 10: Launch a detatchable shell using tmux
Unlike ARC where we submit jobs to a job scheduler, on Nesoi we run jobs interactively in real time. This is great as means we dont have to wait for our analyses to leave a queue before they starts. However, if your laptop disconnects from the server the analysis will stop as ssh requires a constant connection. This can be overcome by launching a detatchable session once we are logged into the server using the command tmux (terminal multiplexer).

Launch a new tmux session called "mysession" (you can call it anything you like)
```
tmux new -s mysession
```

To detatch from the session press control b, immediately followed by d. You can then reattach to the session using:
```
tmux a -t mysession
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
ls /data/Users/Sonya_Myzomela/sample_bams/*.bam
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
