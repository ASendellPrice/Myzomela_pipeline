# Step by step bioinformatics guide for Sonya's *Myzomela* project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## STEP 1: Log in to ARC move into *Myzomela* project directory
Let's log into Arcus-htc (replace OSS with your single sign on):
```
ssh OSS@oscgate.arc.ox.ac.uk
ssh arcus-htc

```

Move into Myzomela directory in "/data/zool-zost/":
```
cd /data/zool-zost/Myzomela
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

The downloaded assembly has been compressed using gzip. Irritatingly, some of the downstream steps require the assembly to be compressed using bgzip instead. So now we need to decompress the file (using gunzip) and compress again using bgzip.

```
gunzip Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
module load samtools
bgzip Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna
```

When it comes to aligning sequencing reads to the reference genome we will need to have created an index for our reference genome.

**Indexing** a genome is analogous to indexing a book. If you want to know on which page a certain word appears or a chapter begins, it is much more efficient/faster to look it up in a pre-built index than going through every page of the book until you found it. Same goes for alignments. Indices allow the aligner (Bowtie2) to narrow down the potential origin of a query sequence within the genome, saving both time and memory.

Create genome index using samtools [faidx](http://www.htslib.org/doc/faidx.html):
This will produce two files with suffixes ".fai" and ".gzi"
```
samtools faidx Ref_Genome/GCA_008360975.1_HeHo_1.0_genomic.fna.gz
```

## STEP 3: Create sample list
Before we dive into the bioinformatics pipeline, we will need to create a list of samples we want to use in the planned analyses. Some analyses (e.g. introgression, SFS, phylogeny building) require use of an outgroup, so if you plan to conduct these you will need to identify publically available WGS data from a suitable outgroup (or we can use other *Meliphagidae* samples included in the Novogene dataset).

The sample list you create needs to be a simple text file (called samples.txt) with one sample ID per line. Like below:

Note: make sure there are no spaces or tabs in the file.
```
cat samples.txt
```

```
ABY15
ABY76
AMB1
```

## STEP 4: Filter raw sequencing data for samples.
