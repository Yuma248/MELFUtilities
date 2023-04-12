#!/usr/bin/perl
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-inR$/){$inputR=shift @ARGV;}
elsif ($_=~ /^-inV$/){$inputVCF=shift @ARGV;}
elsif ($_=~ /^-out$/){$output=shift @ARGV;}
}

if (not defined ($inputR && $inputVCF && $output)) {print "\nThis scritp extracts sequences from a dDocent reference fasta file, according to infromation from a VCF file with your filtered SNPs.\n\nUsage: ExtSeq.pl\n\t-inR <path to reference file fron dDocent, with contig sequences>\n\t-inV <path to VCF file with SNP>\n\t-out <output fasta file name>\n\nExample:\nExtSeq.pl -inR /dDocent/reference.fasta -inV /dDocent/Filtered_SNPs.vcf -out myseq.fasta"; exit;}
open (IN1, "<$inputVCF") or die "Can not open inputfile $!\n";
print "Checking VCF file: $inputVCF \n";
my @lines1=<IN1>;
close IN1;
chomp @lines1;
print "Extracting sequences from reference file: $inputR\n";
our $count=0;
foreach $line1 (@lines1){
next unless ($line1 =~ /^dDocent/);
my @loci = split /\s/, $line1, 2;
$count=$count+1;
print "This is the name of the locus number $count: $loci[0]\n";
my $cmd = "cat $inputR | grep -A 1 \"$loci[0]\$\" >> $output";
system ($cmd);
}
print " A total of $count loci were imported\n";
