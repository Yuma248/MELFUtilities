#!/usr/bin/perl
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-inf$/){$input=shift @ARGV;}
elsif ($_=~ /^-bq$/){$outputBQ=shift @ARGV;}
elsif ($_=~ /^-frt$/){$outputFirst=shift @ARGV;}
elsif ($_=~ /^-rndm$/){$outputRndm=shift @ARGV;}
elsif ($_=~ /^-lm$/){$outputLM=shift @ARGV;}
}

if (not defined ($input)) {print "\nThis script will extract the best quality SNP, first SNP, less missing data SNP and a random SNP per contig. It is usefull when SNPs are call de novo to reduce linkage, by selecting one SNP per contig.\n\nUsage: extBQFRLMRM.pl -inf <input vcf_file> -bq <output best quality SNP,default Best.vcf> -frt <First SNP ouput, default First.vcf> -rndm <random SNP output, default Random.vcf> -lm <output SNP with lest missing data, default Lmisin.vcf>\n\nFor example:\nExtractBQFRLMRM.pl -inf /my/input.vcf -bq /my/bestSNP.vcf -frt /my/firstSNP.vcf -rndm /my/randomSNP.vcf -lm /my/lessmissingSNP.vcf\n\n"; exit;}
if (not defined ($outputBQ)) {$outputBQ="Best.vcf"}; 
if (not defined ($outputFirst)) {$outputFirst="First.vcf"}; 
if (not defined ($outputRndm)) {$outputRndm="Random,vcf"}; 
if (not defined ($outputLM)) {$outputLM="Lmissin.vcf"}; 
open (OUT1, ">>$outputBQ") or die "Can nor open Bestquality output $!\n";
open (OUT2, ">>$outputFirst")or die "Can nor open First output $!\n";
open (OUT3, ">>$outputRndm")or die "Can nor open Random output $!\n";
open (OUT4, ">>$outputLM")or die "Can nor open LessMissing output $!\n";
open (IN, "<$input") or die "Can nor open input $!\n";
our $loci="E000";
our $qualy=0;
our $missing=0;
our $BQline="";
our $Fline="";
our $LMline="";
our @Rline=();
our $uniSNP=0;
{First:while (<IN>){
	chomp($_);
	our $line = $_;
	next if ($line=~ /^##/);
	if ($line=~ /^#CHR/){
		print OUT1 "$line\n";
		print OUT2 "$line\n";
		print OUT3 "$line\n";
		print OUT4 "$line\n";
	}
	elsif ($line=~ /^d/){
		my $new=$line;
		my @data=split(/\t/, $line);
		my $locus=$data[0];
		#		my @locus=split(/_/, $data[0]);
		my @ext=split(/\;/,$data[7]);
		my @miss=split (/\=/,$ext[4]);
		my @geno1=split(//,$data[3]);
		my @geno2=split(//,$data[4]);
		my $pb1=scalar(@geno1);
		my $pb2=scalar(@geno2);
		if ($data[4]=~/,/ or $pb1!=$pb2 ){
			next First;
	}
	our $data3="";
	our $data4="";
		if ($pb1>1){
			for (my $count=0;$count<=$pb1;$count++){
				if ($geno1[$count] ne $geno2[$count] and $geno1[$count] ne "N" and $geno2[$count] ne "N"){
				$data3=$geno1[$count];
				$data4=$geno2[$count];
				}
				}
				$new=~ s/$data[3]/$data3/;
				$new=~ s/$data[4]/$data4/;
				}
		if ($loci ne $locus){
			if ($loci ne "E000"){
				$uniSNP=$uniSNP+1;
				print OUT1 "$BQline\n";
				print OUT2 "$Fline\n";
				print OUT4 "$LMline\n";
				my $rn=int(rand(scalar(@Rline)-1));
				print OUT3 "$Rline[$rn]\n";
	}
			$loci=$locus;
			$qualy=$data[5];
			$missing=$miss[1];
			$Fline=$new;
			$BQline=$new;
			$LMline=$new;
			@Rline= ();
			push (@Rline, $new);
	}
		elsif ($loci eq $locus){
			push (@Rline, $new);
			if ($data[5]>$qualy){
				$BQline=$new;
				$qualy=$data[5];
	}
	if ($miss[1]>$missing){
				$LMline=$new;
				$missing=$miss[1];
	}
	}
	}
	}
	print OUT1 "$BQline\n";
	print OUT2 "$Fline\n";
	my $rn=int(rand(scalar(@Rline)-1));
	print OUT3 "$Rline[$rn]\n";
	$uniSNP=$uniSNP+1;
	print "$uniSNP\n";
close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
}
