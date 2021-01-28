# Step by step bioinformatics guide for Sonya's *Myzomela* project: From rawdata to pop gen analyses
A. Sendell-Price, Jan 2021.

## STEP 1: Log in to ARC move into Myzomela project directory
Let's log into Arcus-htc (replace OSS with your single sign on).
```
ssh OSS@oscgate.arc.ox.ac.uk
ssh arcus-htc

```

Move into data storage for project zool-zost
```
cd /data/zool-zost/Myzomela
```

## STEP 2: Download reference genome and build reference dictionary
Although there are currently no publically available genomes for *Myzomela* species there are three available for species in the wider *Meliphagidae family*. These include:
* *Lichenostomus cassidix* (helmeted honeyeater) [GCA_008360975.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_008360975.1/)
* *Grantiella picta* (painted honeyeater) [GCA_013397955.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_013397955.1)
* *Anthochaera phrygia* (Regent honeyeater) [GCA_009430485.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_009430485.1)

All of these are scaffold level assemblies i.e. not assembled to chromosome level. This isn't too bad as we can find the position of scaffolds relative to a chromosome level assembly later using [Satsuma](http://satsuma.sourceforge.net). More about that later. We will be using the Lichenostomus cassidix as a reference genome as this has the lowest number of scaffolds (1,912).

So now we need to download the assembly from NCBI using command "wget" and save the file in a new directory called "Ref_Genome".

```
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0//GCA_008360975.1_HeHo_1.0_genomic.fna.gz -P Ref_Genome/

```
