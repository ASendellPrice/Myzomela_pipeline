# Step by step bioinformatics guide for Sonya's *Myzomela* project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## STEP 1: Log in to ARC move into *Myzomela* project directory
Let's log into Arcus-htc (replace OSS with your single sign on):
```
ssh OSS@oscgate.arc.ox.ac.uk
ssh arcus-htc
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

All three of these steps are conducted using the pipeline (RawReads2SampleBams.sh) that I have developed for our Novogene sequencing data. This single script will launch a job array (1 job per sample) that takes care of all the filtering, alignment and merging steps. This script has a couple of assumptions.

* raw sequencing reads for an individual are saved in a directory called "sample id" and this directory sits within a directory called "raw_data" like so:
```
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

* the text file "samples.txt" is present (we created this during the step above).
* the reference genome is stored in directory "Ref_Genome" and has been indexed with both samtools and bowtie2 (this was done in step 2).
* line 4 of RawReads2SampleBams.sh is updated so that the job array range matches the number of samples in
"samples.txt"

```
head -n 10 Scripts/RawReads2SampleBams.sh
```

```
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
Note: also need to updae line nine so that notification emails are sent to the correct email address. 

CHANGE PERMISSIONS --   chmod -R a+rwx Myzomela
