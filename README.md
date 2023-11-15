# MELFUtilities
Several script useful to analyse or explore genomic data

# Download
&nbsp;git clone https://github.com/Yuma248/MELFUtilities
        
# Scripts
## ExtSeq_dDocent
This scritp extracts sequences from a dDocent reference fasta file, according to infromation from a VCF file with your filtered SNPs.

Usage:  
ExtSeq.pl  
&emsp;-inR \<path to reference file fron dDocent, with contig sequences\>  
&emsp;-inV \<path to VCF file with SNP\>  
&emsp;-out \<output fasta file name\>  

Example:  
ExtSeq.pl -inR /dDocent/reference.fasta -inV /dDocent/Filtered_SNPs.vcf -out myseq.fasta  

## ExtractBQFRLMRM  
This script will extract the best quality SNP, first SNP, less missing data SNP and a random SNP per contig. It is usefull when SNPs were call de novo to reduce linkage, by selecting one SNP per contig.  

Usage:  
ExtractBQFRLMRM.pl  
&emsp;-inf \<input vcf_file\>  
&emsp;-bq \<output best quality SNP, default Best.vcf\>  
&emsp;-frt \<First SNP ouput, default First.vcf\>  
&emsp;-rndm \<random SNP output, default Random.vcf\>  
&emsp;-lm \<output SNP with lest missing data, defaul Lmissin.vcf\>  

For example:  
ExtractBQFRLMRM.pl -inf /my/input.vcf -bq /my/bestSNP.vcf -frt /my/firstSNP.vcf -rndm /my/randomSNP.vcf -lm /my/lessmissingSNP.vcf  

## UpFigshare  
This script will upload a file to a figshare repository, it needs the whole path of the file to be uploaded, a URL from figshare, and a token to access the repository or project. If the file is not required to be uploaded in a specific project or item, you can use the general URL https://api.figshare.com/v2/account/articles. But if you need to put in a specific project or item (folder) you will require the project ID (using the -p option) and/or item ID (usign the -f option). If you already upload one file in the project and item, you can integrate the item ID in the URL and add just /files at the end, see the example.  

Usage:  
UpFigshare.sh  
&emsp;-i \<path to input file\>  
&emsp;-u \<URL of figshare account/project/item, required if you have specific project and item ready, default https://api.figshare.com/v2/account/articles\>  
&emsp;-t \<Access token for the figshare account\>  
&emsp;-p \<project ID, this is required just if you require to upload the file in specific project\>  
&emsp;-f \<item ID, this is required just if you need to upload the file in specific item (folder)\>  
&emsp;-n \<item name, this is required if you are creating a new item (folder) and you want a specific name, default NEW UPLOAD\>  

For example:  
UpFigshare.sh -i /sanpper/genome/genomev2.fasta -u https://api.figshare.com/v2/account/articles/ -t 75050303931z87ab7c72038ab9eaf02d853766bd8f7cc695f390d5b9cdeda1fd230c462a7cb7c7a67ef7c507f27f3fde647c4145667664533374d54609ef477874c4aa11 -p 24892

## basic_filtering_dart
This script performs basic filtering on a dart CSV file and output a table with the remaining SNPs after each filtering step, the basic stats (#samples, #loc. #polymorphic loc, observerd heterozygozity, espected heterozygozity, and Fis inbreeding coefficient) after filtering, a FST plot, a PCA plot and a Admixture/Structure plot.   

To use this fuction in R you will need the to install and load the next R libraries:
 
library(dartR)  
library(tidyverse)  
library(vegan)  
library(related)  
library(igraph)  
library(reshape2)  
library(ggplot2)  
library(HardyWeinberg)  
library(ggpubr)  
library(gridExtra)  
library(LEA)  
  
And either copy and save the fuction in R or use the command  
  
source("path/to/MELFUtilties/basic_filtering_dart.R")  


You will have to upload your dart data and samples metadata (sample pop map or pop info file) in R, using commands like: 

dartData <- gl.read.dart(filename="Report_Sps0000_SNP_mapping.csv",ind.metafile="SamplesAllMetadata.csv")

meta <- read.csv("SamplesAllMetadata.csv", header=TRUE)

 
Usage:  
basic_filter(dartData, "output prefix", maxmisi = 50, mincalL = 0.80, minrep = 0.99, minmaf = 0.03, secd = TRUE, HWEF = TRUE, depthr = c(5,75), npopsHE = 5, maxsim = 0.85, PDFplots = TRUE)  
&emsp;dartD    - DArT genotype file    
&emsp;name  - prefix for output files (string)  
&emsp;maxmisi  - maximum proportion of missing data per individual   
&emsp;mincalL  - minimum call rate per locus, inverse of maximum missing data por locus  
&emsp;minrep  - minimum reproducibility   
&emsp;minmaf  - minimum minor allele frequency  
&emsp;maxsim  - max reletedness to remove duplicates  
&emsp;secd  - remove secundaries (TRUE, FALSE)  
&emsp;depthr  - min and maximum depth of coverage  
&emsp;HW  - filter for HWE (TRUE, FALSE)  
&emsp;npopsHE  - Maximum number of population where a locus can be out of HWE  
&emsp;PDFplots  - If you want to create a PDF with basic results (filter steps table, diversity table, Fst heat map, PCA)   


Example:   
filterdata <- basic_filter(YumadartData, "Yuma", maxmisi = 50, mincalL = 0.80, minrep = 0.99, minmaf = 0.03, secd = TRUE, HWEF = TRUE, depthr = c(5,75), npopsHE = 5, maxsim = 0.85, PDFplots = TRUE)
 
 
 
 







