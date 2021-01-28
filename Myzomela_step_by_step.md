# Step by step bioinformatics guide for Sonya's Myzomela project: From rawdata to pop gen analyses
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

All of these are scaffold level assembly i.e. not assembled to chromosome level. This isn't too bad as we can find the position of scaffolds relative to a chromosome level assembly later using [Satsuma](http://satsuma.sourceforge.net). More about that later. We will be using the Lichenostomus cassidix as a reference genome as this has the lowest number of scaffolds (1,912).


In our first assessment of hybridization across the genome we will calculate a single D-statistic for the genome. This will give us an idea of the overall extent of introgression between the two species. For this we will use the program [DSuite](https://github.com/millanek/Dsuite) which allows for the fast calculation of the D-statistic from SNP data in VCF format.


```
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/360/975/GCA_008360975.1_HeHo_1.0//GCA_008360975.1_HeHo_1.0_genomic.fna.gz -P Ref_Genome/


```
To run the analysis we will need to create a tab delimited file ("Pops.txt") that assigns samples to populations. That file looks like this, where "Outgroup" is used to assigns individuals to the outgroup (in this case we use the Reunion grey white-eye a.k.a *Z.borbonicus*) and remaining samples are assigned to their corresponding populations:
