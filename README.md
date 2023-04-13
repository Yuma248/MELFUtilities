# MELFUtilities
Several script useful to analyses or explore genomic data

# Download
        git clone https://github.com/Yuma248/MELFUtilities
        
# ExtSeq_dDocent
This scritp extracts sequences from a dDocent reference fasta file, according to infromation from a VCF file with your filtered SNPs.

Usage: 
ExtSeq.pl
        -inR <path to reference file fron dDocent, with contig sequences>
        -inV <path to VCF file with SNP>
        -out <output fasta file name>

Example:
ExtSeq.pl -inR /dDocent/reference.fasta -inV /dDocent/Filtered_SNPs.vcf -out myseq.fasta
# ExtracctBQFRLMRM 
This script will extract the best quality SNP, first SNP, less missing data SNP and a random SNP per contig. It is usefull when SNPs were call de novo to reduce linkage, by selecting one SNP per contig.

Usage:
ExtractBQFRLMRM.pl
        -inf <input vcf_file>
        -bq <output best quality SNP, default Best.vcf>
        -frt <First SNP ouput, default First.vcf>
        -rndm <random SNP output, default Random.vcf>
        -lm <output SNP with lest missing data, defaul Lmissin.vcf>

For example:
ExtractBQFRLMRM.pl -inf /my/input.vcf -bq /my/bestSNP.vcf -frt /my/firstSNP.vcf -rndm /my/randomSNP.vcf -lm /my/lessmissingSNP.vcf
