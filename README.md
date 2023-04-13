# MELFUtilities
Several script useful to analyses or explore genomic data

# Download
        git clone https://github.com/Yuma248/MELFUtilities
        
#ExtSeq_dDocent
This scritp extracts sequences from a dDocent reference fasta file, according to infromation from a VCF file with your filtered SNPs.

Usage: ExtSeq.pl
        -inR <path to reference file fron dDocent, with contig sequences>
        -inV <path to VCF file with SNP>
        -out <output fasta file name>

Example:
ExtSeq.pl -inR /dDocent/reference.fasta -inV /dDocent/Filtered_SNPs.vcf -out myseq.fasta
