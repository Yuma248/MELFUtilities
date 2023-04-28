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

## UpFigshare.sh  
This script will upload a file to a figshare repository, it needs the whole path of the file to be uploaded, a URL from figshare, and a token to access the repository or project. If the file is not required to be uploaded in a specific project or item, you can use the general URL https://api.figshare.com/v2/account/articles. But if you need to put in a specific project or item you will require the ID for the project and/or item. If you already upload one file in the project an item, you can integrate the ID in the URL and add just /files at the end, see the example.  

Usage:  
UploadFSY.sh  
&emsp;-i <path to input file>  
&emsp;-u <URL of figshare account/project/item, required if you have specific project and item ready, default https://api.figshare.com/v2/account/articles>  
&emsp;-t <Access token for the figsahre account>  
&emsp;-p <project ID, this is required just if you require to upload in specific project>  
&emsp;-f <item ID, this is required just if you need to upload in specific item>  
&emsp;-n <item name, this is required if you are creating a new item, default NEW UPLOAD>  

For example:  
UploadFSY.sh -i /sanpper/genome/genomev2.fasta -u https://api.figshare.com/v2/account/articles/22707369/files -t 75050303931z87ab7c72038ab9eaf02d853766bd8f7cc695f390d5b9cdeda1fd230c462a7cb7c7a67ef7c507f27f3fde647c4145667664533374d54609ef477874c4aa11 

