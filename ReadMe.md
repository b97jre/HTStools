# HTStools 


This repository contains years of java files that is used to parse and handle various functions when handling high throughput data. Unfortunately it is not annotated well. 



To run different functions you can call the HTStools.jar file. You also need to have the jar folder in same location as the HTStools.jar file since this contains files that HTStools need to run many of the functions. 


I will add examples that I have there based on usage.

## Commands used for extracting information regarding phased vcf files.

### Turn a unphased vcf file into a phased on using another phased vcf file

Phasing the files using the counts from the unphased RNA vcf files and matching that with the phased info from another vcf file 

```bash
java -Xmx7g -jar HTStools.jar -p databases -phaseVCFfile 
-R <reference.fa> \
-VCF  <reference.phased.vcf> \
-unphased <sample.vcf> \
-readPhased 
```
Output will be a phased sample vcf file `sample.phased.vcf` 



### Extracting one vcf file for one sample and only including the heterozygous sites 

```bash
java -Xmx7g -jar HTStools.jar -p databases -getPhasedVCFinfo \
-R REFERENCE.fa \
-BED REFERENCE_GENES.bed \
-VCF Unified.output.raw.snps.indels.DNAreadPhased.vcf \
-sample SAMPLE_NAME 
```
Output will be a phased sample vcf file `SAMPLE_NAME.heterozygous.Unified.output.raw.snps.indels.DNAreadPhased.vcf` with all heterozygous sites only.



### Phasing the RNA files using the counts from the unphased RNA vcf files and matching that with the phased info from the DNA read-phased information 

```bash
java -Xmx7g -jar HTStools.jar -p databases -phaseVCFfile \
-R <reference.fa> \
-VCF  <reference.phased.vcf>\
-unphased <sample.vcf>\
 -readPhased 
```
Output will be a phased sample vcf file `sample.phased.vcf` 




### Turn phased file into skelly format.

Using phased vcf files and bed files containing the location of the genes on the chromosomes to create correct output for ASE analysis using skelly using a inhouse java script. 
Bed files to identify the locations of the genes and the phased vcf files to identify the longest stretch of phased heterozygous SNPs for each gene.

Default is that the file is fully phased by heritage. If the file is phased by reads this has to be flagged using the read phase flag `-readPhased`

```bash
java -Xmx7g -jar HTStools.jar -p databases -SkellyFormat \
-readPhased -R <reference.fa> \
-BED <reference.gene.bed> \
-VCF  <sample.phased.vcf> 
```


### Generate haplogenomes from phased vcf file and reference

Generate two haplo genomes based on genome reference sequences and phased DNA sample 
in vcf file


```bash
java -Xmx20g -jar HTStools.jar -p databases -printPhasedGenome\
-R REFERENCE.fa \
-phasedVCF ReadPhased.vcf \
-sample SAMPLE_NAME
```
output is  `REFERENCE.sample.phased.mother.fa` and `REFERENCE.sample.phased.father.fa`


### Get phased vcf file based on mpileup file and phased vcf file 

 Use DNA phased heterozygous vcf file with RNA mpileupe files to get the RNA phased vcf files

java -Xmx20g -jar HTStools.jar -p databases -parseMpileUpToVCF 
-R REFERENCE.fa \
-phasedVCF SAMPLE_NAME.heterozygous.Unified.output.raw.snps.indels.DNAreadPhased.vcf \
-mpileupFile Aligned.out.sorted.F_256.mpileup \
-sample SAMPLE_NAME \
-mother sample_SAMPLE_NAME_phased_Mother \
-father sample_SAMPLE_NAME_phased_Father 

This will generate the phased RNA vcf file with the name 
`Aligned.out.sorted.F_256.mpileup.SAMPLE_NAME.heterozygous.Unified.output.raw.snps.indels.DNAreadPhased.vcf"


