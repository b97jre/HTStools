# HTStools 


This repository contains years of java files that is used to parse and handle various functions when handling high throughput data. Unfortunately it is not annotated well. 



To run different functions you can call the HTStools.jar file. You also need to have the jar folder in same location as the HTStools.jar file since this contains files that HTStools need to run many of the functions. 


I will add examples that I have there based on usage.


## Handling phased vcf files and making them possible to run using the Skelly approach. 



### Turn a unphased vcf file into a phased on using a phased vcf file


Phasing the RNA files using the counts from the unphased RNA vcf files and matching that with the phased info from the DNA read-phased information 

```bash
java -Xmx7g -jar HTStools.jar -p databases -phaseVCFfile -R <reference.fa> -VCF  <reference.phased.vcf> -unphased <sample.vcf> -readPhased 
```
Output will be a phased sample vcf file `sample.phased.vcf` 

### Turn phased file into skelly format.

Using phased vcf files and bed files containing the location of the genes on the chromosomes to create correct output for ASE analysis using skelly using a inhouse java script. 
Bed files to identify the locations of the genes and the phased vcf files to identify the longest stretch of phased heterozygous SNPs for each gene.

Default is that the file is fully phased by heritage. If the file is phased by reads this has to be flagged using the read phase flag `-readPhased`

```bash
java -Xmx7g -jar HTStools.jar -p databases -SkellyFormat -readPhased -R <reference.fa> -BED <reference.gene.bed> -VCF  <sample.phased.vcf> 
```


