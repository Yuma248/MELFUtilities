#!/usr/bin/perl
########Original escript from Brauer to create a gINLAnd inputfile ###########
########Script to create 2 population SFS from a Genepop file###########
########Genepop file should have extention .gen #######################
########just run perl genepop_to_SFS.pl folder/to/Genpopfile.gen#######


use strict;
use warnings;


use Cwd;
my $workingdirectory = getcwd();

#Global Constants
my $Line_1 = 1;
#the genotypes are coded 001, 002, 003, and 004 (A, T, C, and G) and each individual has two values
my @Valid_Genotypes = ("001", "002", "003", "004");#changed  ("010", "011", "012", "013")to suit dDocent pipeline converting vcf to stru to genepop 010, 011, 012, 013
my $NO_DATA = "000";
@ARGV ||  die "Error: No arguments provided.\n\tUsage: perl genepop_to_SFS.pl folder/to/Genepopfile.gen\n";
my $input_dirname = $ARGV[0];;
my $output_dirname = "${workingdirectory}/results";
print "$input_dirname\n";
sub Is_Valid_Genotype {
   my $n = scalar(@_);
   my $result = "FALSE";
   if ($n != 1) {
     print "too many paramters passed in\n";
   } else {
     foreach my $item (@_){
       foreach my $value (@Valid_Genotypes) {
         if ($item eq $value) {
           $result = "TRUE"; } } } }
  return $result;
}
unless(-e $output_dirname or mkdir $output_dirname) { die "Unable to create $output_dirname directory"; }
opendir ( DIR, $input_dirname ) || die "Error in opening dir $input_dirname\n";
my @files = grep {/.*\.gen/} readdir(DIR);

my $num_args = $#files + 1;
if ($num_args <= 0 ) { die "no arguments found, you must provide at least on input file name";}

foreach my $infile (@files) {
  # Create and open the required files
  my $filename=$infile;
  $filename=~s/.gen//; #Stripped file name to create a matched out file name

  my $IN_DATA_FILE_NAME = "$input_dirname/$infile";
  open(IN_DATA_FILE, $IN_DATA_FILE_NAME) || die "can't open in file $IN_DATA_FILE_NAME, aborting attempt";

  # Variables in scope per input file.
  my $Individual_Count = 0;
  my $line_Count = 0;
  my $MARKER_DONE = "FALSE";
  my $MARKER_COUNT = 0;
  my $Population_Count = 0;
  my @Markers = ();
  my @Population_Marker_Count = ();
  my @Population_Allele_Count = ();
  my @Total_Marker_Count = ();
  my @Total_Allele_Count = ();
  my $NO_DATA = "000";
  my @Population_Key;
  my @Individual_Marker;
  my $OK = "FALSE";


  # Process the data from the current input file
  #
  while (<IN_DATA_FILE>) {
   chomp;
   $line_Count++;
   # Find the list of Markers for this file, ignoring the files first line.
   if ($line_Count != $Line_1) {
     if ("$MARKER_DONE" eq "FALSE") {
      if ("$_" ne "Pop") {
          $MARKER_COUNT++;
          my $MARKER_Index = $MARKER_COUNT - 1;
          $Markers[$MARKER_Index] = $_;
                  #print "$MARKER_Index\n";
      } else {
          # we've got the marker list
          $MARKER_DONE = "TRUE";
          # Create the first populations key_set
          @Population_Key = ($NO_DATA) x $MARKER_COUNT;
                    @Individual_Marker = ($NO_DATA) x $MARKER_COUNT;
          push @Population_Marker_Count, [ (0) x $MARKER_COUNT ] for (0);
          push @Population_Allele_Count, [ (0) x $MARKER_COUNT ] for (0);
                  push @Total_Marker_Count, [ (0) x $MARKER_COUNT ] for (0);
          push @Total_Allele_Count, [ (0) x $MARKER_COUNT ] for (0); }
   } else {
   # At this point we've pulled off the marker info and need to start to construct the population data.
     if ("$_" eq "Pop") {
        $Population_Count++;
        $Individual_Count = 0;
        push @Population_Marker_Count, [ (0) x $MARKER_COUNT ] for ($Population_Count);
        push @Population_Allele_Count, [ (0) x $MARKER_COUNT ] for ($Population_Count);
      } else {
         $OK = "FALSE";
         my @Individual_Alleles;
         # take the current line and break it up into tokens.
         my @tokens = split('\s', $_);
         #ignore the first 2 tokens (the individual identifier and the comma.)
          for (my $i = 2; $i <@tokens; $i++) {
            # If that's the case we need to dump the data as we don't know how it fits.
            $Individual_Alleles[$i - 2] = $tokens[$i];
            if (($i - 1) == $MARKER_COUNT) {$OK = "TRUE"; } else {$OK = "FALSE";}
          }
          # check that this has the correct number of datapoints - ignore it if not.
          if ($OK eq "TRUE") {
              # If this is the first individual of the first population in the file, it contains the keys, so save them.
              # SABR TBD: need to handle the case where the first individual in the population had bad data, but the rest don't
              # so also check if the stored key if invalid and try and find one if it's not.
              if (($Individual_Count == 0) && ($Population_Count == 0)) {
                 # find the first genometype and then start the count.
                 # loop through the individuals alleles
                 # split them into single genometypes and count for the first one
                 for (my $m = 0; $m <@Markers; $m++) {
                      my @Genome = (substr( $Individual_Alleles[$m], 0, 3 ), substr( $Individual_Alleles[$m], 3, 3 ));
                      $Population_Key[$m] = $Genome[0];
                 }
              }
              #process each individual (including the one the key came from) to count the matches
              for (my $m = 0; $m <@Markers; $m++) {
                 #check if we have a valid population key, if not see if this individual can provide one before progressing with the processing
                 if (Is_Valid_Genotype($Population_Key[$m]) ne "TRUE") {
                    my @Genome = (substr( $Individual_Alleles[$m], 0, 3 ), substr( $Individual_Alleles[$m], 3, 3 ));
                      $Population_Key[$m] = $Genome[0];
                 }
                 #Now continue on with the count.
                 if (Is_Valid_Genotype($Population_Key[$m]) eq "TRUE") {
                   $Individual_Marker[$m] = $Individual_Alleles[$m];

                   my @value = (substr( $Individual_Alleles[$m], 0, 3 ), substr( $Individual_Alleles[$m], 3, 3 ));
                   for (my $v = 0; $v <@value; $v++) {
                     if ($value[$v] eq $Population_Key[$m]) {
                        my $temp_marker_count = $Population_Marker_Count[$Population_Count][$m];
                        $Population_Marker_Count[$Population_Count][$m] = $temp_marker_count + 1;
                                                $Total_Marker_Count[0][$m]=$Total_Marker_Count[0][$m]+1;
                     }
                     if ($value[$v] ne $NO_DATA) {
                        my $temp_allele_count = $Population_Allele_Count[$Population_Count][$m];
                        $Population_Allele_Count[$Population_Count][$m] = $temp_allele_count + 1;
                                                $Total_Allele_Count[0][$m]=$Total_Allele_Count[0][$m]+1;
                     }
                   }
                 }
               }
              # at this point, the individuals counts have been added to the populations counts
             $Individual_Count++;
                #print "$Individual_Count\n";
            } else {
               #skip this one, but record the event.
               print "Dropped a line on input during parse of population $Population_Count, Individual $Individual_Count due to incorrect number of data points\n";
            }
     }
   } # end of the poplution search
 }#exclude the first line
} # end of file read

  # END THE DATA FILE PROCESSING

  close(IN_DATA_FILE);
  #write to the outputfile
  # MARKERS
  #foreach (my $i = 0; $i <@Markers; $i++) {
  #print "$Markers[$i]\t";
    #if ($i != $#Markers) {
      #print  "$Markers[$i] 1\t";
      #print  "$Markers[$i] 2\t";
    #} else {
      #print "$Markers[$i]\n";
   #   print OUTFILE2 "$Markers[$i]\n";
    #}
  #}
  my @MAXIMUS=();
  #print "$MARKER_COUNT\n";
  foreach (my $population = 0; $population < $Population_Count +1; $population++) {
  my $ma1=0;
  my $ma2=0;
  $MAXIMUS[$population]=0;
                foreach (my $population2 = 0; $population2 < $Population_Count +1; $population2++) {
                my @COUNTV=();
                next if ($population>=$population2);
                $MAXIMUS[$population2]=0;
            for (my $count = 0; $count <$MARKER_COUNT; $count++) {
                my $alelle1=$Total_Marker_Count[0][$count]+0;
                my $alelle2=($Total_Allele_Count[0][$count]-$alelle1);
                my $totalna=$Total_Allele_Count[0][$count];
                my $totpop1=$Population_Allele_Count[$population][$count];
                my $totpop2=$Population_Allele_Count[$population2][$count];
                #print "$population\t$population2\t$totalna\t$totpop1\t$totpop2\t$alelle1\t$alelle2\n";
                if ($alelle1<=$alelle2){
                $ma1=$Population_Marker_Count[$population][$count]+0;
                $ma2=$Population_Marker_Count[$population2][$count]+0;
                unless (defined $COUNTV[$ma1][$ma2]){$COUNTV[$ma1][$ma2]=0;}
                $COUNTV[$ma1][$ma2]=$COUNTV[$ma1][$ma2]+1;}
                else {
                $ma1=$Population_Allele_Count[$population][$count]-$Population_Marker_Count[$population][$count];
                $ma2=$Population_Allele_Count[$population2][$count]-$Population_Marker_Count[$population2][$count];
                unless (defined $COUNTV[$ma1][$ma2]){$COUNTV[$ma1][$ma2]=0;}
                $COUNTV[$ma1][$ma2]=$COUNTV[$ma1][$ma2]+1;}
                #print "$ma1\t$ma2\t$COUNTV[$ma1][$ma2]\t$population\t$population2\n";
                if ($ma1>$MAXIMUS[$population]){$MAXIMUS[$population]=$ma1;}
                if ($ma2>$MAXIMUS[$population2]){$MAXIMUS[$population2]=$ma2;}
                }
        my $OUTFILESFSF = "$output_dirname/$filename-SFSf_jointMAFpop$population2-$population";
        open (OUTFILE3, '>', $OUTFILESFSF) || die "Can't open the output file. exiting";
        print OUTFILE3 "1 observations\n";
        my @DSFS0=(0..$MAXIMUS[$population]);
        my @DSFS1=(0..$MAXIMUS[$population2]);
        my $DSFSV="";
        my $DSFSD="";
        foreach my $D0 (@DSFS0){$DSFSD=$DSFSD."\t"."d0_".$D0;}
        print OUTFILE3 "$DSFSD\n";
        foreach my $D1 (@DSFS1){
        #print "$numale\n";
        my $DSFSV="d0_$D1";
        foreach my $D0 (@DSFS0){
        if (defined $COUNTV[$D0][$D1]){
        $DSFSV=$DSFSV."\t".$COUNTV[$D0][$D1];
        }else {
        $COUNTV[$D0][$D1]=0;
        $DSFSV=$DSFSV."\t".$COUNTV[$D0][$D1];
        }
        }
        print OUTFILE3 "$DSFSV\n";
        }
        close OUTFILE3;
   }
   }
  }
closedir(DIR);
