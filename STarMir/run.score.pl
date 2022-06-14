#!/usr/bin/perl

# usage: >perl run.score.pl <working directory> <species> 
#
# calculate score.  Concatenate all of the input files and then calculate the score for each discovered
# miRNA:target pair.  Score is added to all input files.  Almighty kludge as is all the score project.
#
# Note: needed to do score calculation after the R script as it requires the combination of all six logit files to calculate. separate script would add impossible complexity.
# A complexity occurs in that the script has to handle the situation where a logit file is missing, as one is not produced if there are no binding sites for that case.
#
# 2022-4-21 (war)
# removed all reference to Carp module for distribution version
#
use warnings;
use strict;

use File::Basename;
#use Carp::Assert; # use Carp for debugging  no Carp to shut off debugging
my $debug = 0;

use Cwd;

my $cwd = cwd;

print STDERR "working directory is $cwd\n";

print  "\nEntering $0\n" if $debug;
print  "Arguments:\n\t" if $debug;
print  join("\n\t", @ARGV) . "\n" if $debug;

# holders for arguments to script
my $workDir; # directory containing logit files, generally the directory the top level script is called in
my $species;

# Globals
# Constants

# efficacy hash (using both final eff sets)
my %EFF1 =
( '8mer' => 0.0198010483,
  '7mer' => 0.0132759554,
  '6mer' => 0.0065819796,
  'seedless' => 0.0042824619);

my %EFF2 =
( '8mer' => 0.0296250981,
  '7mer' => 0.0129698398,
  '6mer' => 0.0065819796,
  'seedless' => 0.0055265616);

# seed type mapping
# note: no GU allowed, sites containing GU pairs labeled seedless
my %SEEDKIND =
( '8mer' => '8mer',
  '7mer-m8' => '7mer',
  '7mer-A1' => '7mer',
  '6mer' => '6mer',
  'offset-6mer' => '6mer',
  'other' => 'seedless');

# input fields
my $TARGETCOL = 1;
my $MIRNACOL = 2;
my $TARGETLEN = 3;
my $SEEDTYPE = 6;
my $LOGITPROBSEEDLESS = 36;
my $LOGITPROBSEED = 39;

# hashes for  =accumulating values
my %eff1 = (); # eff using first trial coeffs
my %eff2 = (); # eff using second trial coeffs

# collect arguments
if (scalar @ARGV >= 2) {
    $workDir= $ARGV[0];
    $species = $ARGV[1];
}
else {
    die("usage: >perl run.score.pl <working directory> <species> \n");
} # end if else

my @types = ("seed","seedless");
my @components = ("3pUTR", "CDS", "5pUTR");


foreach my $type (@types) {
	foreach my $component (@components) {

	    print  "In main loop, type=$type, component=$component\n" if $debug;

		my $logitBase = "Output-bsite-$type-$component-$species.txt";
		my $logitFile = "${workDir}/$logitBase";

		# determine if logit site file exits.  It will not be created
		#   there are no binding sites for a given component/type
		#   best way I could figure to signal this script there was no 
		#   data

		if ( -e $logitFile) {

		    print  "logit file: $logitFile\n" if $debug; 
		    open(INFILE, $logitFile) or next;

		    <INFILE>; # discard header
		    while(my $inline = <INFILE>) {
			chomp $inline;
			$inline =~ s/\"//g; # dammed R quotes mess everything up

			my @richard = split "\t", $inline;

			# find the pair find the site type calculate the eff
			my $key = $richard[$TARGETCOL] . "#" . $richard[$MIRNACOL];

			print  'key is : ' . $key  . "\n" if $debug;
			
			# we have to use the type component of the filename to determine whether
			# we are dealing with the seed or seedless case.  This is because
			# we require no wobble pairs, and RNA hybrid does not.
			# to directly test for wobble pairs we have to go back to the 
			# RNAhybrid output.
			my $seedType = $SEEDKIND{$richard[$SEEDTYPE]};
			if ( ($type eq 'seed') 
			    and ($seedType ne 'seedless') ) {
			    # seed case

			    $eff1{$key} = 0 if !defined $eff1{$key};
			    $eff1{$key} += $richard[$LOGITPROBSEED] * $EFF1{$seedType};
			    $eff2{$key} = 0 if !defined $eff2{$key};
			    $eff2{$key} += $richard[$LOGITPROBSEED] * $EFF2{$seedType};

    
			}
			else {
			    # seedless case

			    $eff1{$key} = 0 if !defined $eff1{$key};
			    $eff1{$key} += $richard[$LOGITPROBSEEDLESS] * $EFF1{'seedless'};
			    $eff2{$key} = 0 if !defined $eff2{$key};
			    $eff2{$key} += $richard[$LOGITPROBSEEDLESS] * $EFF2{'seedless'};

			}
		    } # end while lines in input
		    close INFILE;
}
		else {
		    print  "$logitFile not found\n" if $debug;
		} # end if else binding site file found
	} # end foreach component
} # end foreach type

# at this point we have populated the hashes
# I reopen and reread the files to avoid messing up the sorts, or anything else the earlier
# scripts might have done.  Again an almighty kludge. Preservation of promenance trumps style.

# again for each case
foreach my $type (@types) {
	foreach my $component (@components) {

		my $logitBase = "Output-bsite-$type-$component-$species.txt";
		my $logitFile = "${workDir}/$logitBase";
		my $outBase = "Final-$type-$component-$species.txt";
		my $outFile = "${workDir}/$outBase";

		# determine if logit site file exits.  It will not be created
		#   there are no binding sites for a given component/type
		#   best way I could figure to signal this script there was no 
		#   data

		if ( -e $logitFile) {

		    print  "logit file: $logitFile\n" if $debug; 
		    open(INFILE, $logitFile) or next; # next for robustness.  Have allowed for missing file above.
		    open (OUTFILE, ">$outFile") or die(); 
		    my $inline = <INFILE>;
		    chomp $inline; # get header line		    
		    print OUTFILE $inline . "\t" . 'Score' . "\n"; # add score column

		    while(my $inline = <INFILE>) {
			chomp $inline;
			$inline =~ s/\"//g; # dammed R quotes mess everything up		    
			my @richard = split "\t", $inline;
			print join "\n" ,@richard . "\n" if $debug;

			my $key = $richard[$TARGETCOL] . "#" . $richard[$MIRNACOL];

			# choice of eff1 hash is arbitrary.  The difference between the values in the two hashes is not signifigant.
			print 'key is : ' . $key  . ':' . $eff1{$key} . "\n" if $debug;
			
			print OUTFILE $inline . "\t" . $eff1{$key} . "\n";
		    } # end while lines in input
		    close INFILE;
		    close OUTFILE;
		}
		else {
		    print  "$logitFile not found\n" if $debug;
		} # end if else binding site file found
	} # end for each component
} # end foreach type

close OUTFILE;
print  "Exiting $0\n" if $debug;

exit 0;


