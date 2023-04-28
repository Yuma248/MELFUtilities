# MELFUtilities
Several script useful to analyse or explore genomic data

# Download
&nbsp;git clone https://github.com/Yuma248/MELFUtilities
        
# Scripts
## ExtSeq_dDocent
This scritp extracts sequences from a dDocent reference fasta file, according to infromation from a VCF file with your filtered SNPs.

Usage:  
ExtSeq.pl  
&emsp;-inR <path to reference file fron dDocent, with contig sequences>  
&emsp;-inV <path to VCF file with SNP>  
&emsp;-out <output fasta file name>  

Example:  
ExtSeq.pl -inR /dDocent/reference.fasta -inV /dDocent/Filtered_SNPs.vcf -out myseq.fasta  

## ExtractBQFRLMRM  
This script will extract the best quality SNP, first SNP, less missing data SNP and a random SNP per contig. It is usefull when SNPs were call de novo to reduce linkage, by selecting one SNP per contig.  

Usage:  
ExtractBQFRLMRM.pl  
&emsp;-inf <input vcf_file>  
&emsp;-bq <output best quality SNP, default Best.vcf>  
&emsp;-frt <First SNP ouput, default First.vcf>  
&emsp;-rndm <random SNP output, default Random.vcf>  
&emsp;-lm <output SNP with lest missing data, defaul Lmissin.vcf>  

For example:  
ExtractBQFRLMRM.pl -inf /my/input.vcf -bq /my/bestSNP.vcf -frt /my/firstSNP.vcf -rndm /my/randomSNP.vcf -lm /my/lessmissingSNP.vcf  
