#! /bin/perl

# AddCols.pl -- William Rennie 2011-8-10
#
# usage: AddCols.pl <input file name> <output file name>
#
# This script does the final processing to prepare data
# for logistic regression processing.
# it works on a single file, adding two additional columns
# via table lookup in the Appropriate enrichment file.
### these columns no longer added 3/22/2012 ###
#
# It also selects values for Acc and AU based on the 
# appropriate window size for training species and
# state.
### Values are no longer eliminated.  All window sizes are passed on to the
# R script, which resolves wsat is used. ###
#
# Some headers are remapped to reflect the requirements
# of the logistic regression R script.
#
# Certain fields are eliminated as the features are no longer
# needed.
#
# Generally, this script is used to finally prepare
# the input files for the R script.
#
# If the input file contains no lines save the header, do not 
#  create an output file.
# 
# 3/22/2012 (war)
# major modifications. Removed all the enrichment elements of the
# code and altered the column headers to conform to new R script requirements.
# reduced the size of the script by so much, its practically absurd.
#
# 3/30/2012 (war)
# changes in the R script now render the entire mapping of column names
# moot.  Leaving the mapping array in place in anticipation of future changes.
# Its like trying to pick up jello with a pitchfork sometimes.
#
# 2012-4-21 (war)
# Modified to remove references to Carp for distribution.
#
use strict;
use warnings;

use Getopt::Std;
use File::Basename;
# use Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals

my $InFName; # Input Filename
my $OutFName;# Output Filename 

my @InputLines; # array of lines from the input file.

# data are (<name in input file>, <name in output file>, <input col> <output col?>)
# 0 do not use in output, 1 use in output.
# output column name is the name required by the R script.
#
my @colMapping = (
    ['SiteID', 'SiteID', 0, 1],
    ['Transcript','Transcript', 1, 1],
    ['miRNA', 'miRNA', 2, 1],

    ['Ghybrid','Ghybrid', 3, 1],

#    ['Site_range', 'site_position', 4, 1],
#    ['Seed_range', 'seed_position', 5, 1],

    ['Site_range', 'Site_range', 4, 1],
    ['Seed_range', 'Seed_range', 5, 1],

    ['Seed_kind', 'Seed_kind', 6, 1],

    ['GU_pair', 'GU_pair', 7, 0],
    ['tA1', 'tA1', 8, 0],
    ['region12_17','region12_17', 9, 1],
    ['region4-15', 'region5-15', 10, 0],
    ['2-8_range', '2-8_range',11, 0],
    ['9-11range', '9-11range', 12, 0],
    ['12+range', '12+range', 13, 0],
    ['SiteAcc', 'SiteAcc', 14, 1],

    ['Old_Seed_be', 'Old_Seed_be', 15, 0],
    ['Old_Seed_en', 'Old_Seed_en', 16, 0],
    ['Old_Seed_Acc', 'Old_Seed_Acc', 17, 0],

    ['New_Seed_be', 'New_Seed_be', 18, 0],
    ['New_Seed_en', 'New_Seed_en', 19, 0],
    ['SeedAcc', 'SeedAcc', 20, 1],

    ['UpsAcc_5nt', 'UpsAcc_5nt', 21, 1],
    ['UpAU_5nt', 'UpAU_5nt', 22, 1],
    ['DwsAcc_5nt', 'DwsAcc_5nt', 23, 1],
    ['DwAU_5nt', 'DwAU_5nt', 24, 1],

    ['UpsAcc_10nt', 'UpsAcc_10nt', 25, 1],
    ['UpAU_10nt', 'UpAU_10nt', 26, 1],
    ['DwsAcc_10nt', 'DwsAcc_10nt', 27, 1],
    ['DwAU_10nt', 'DwAU_10nt', 28, 1],

    ['UpsAcc_15nt', 'UpsAcc_15nt', 29, 1],
    ['UpAU_15nt', 'UpAU_15nt', 30, 1],
    ['DwsAcc_15nt', 'DwsAcc_15nt', 31, 1],
    ['DwAU_15nt', 'DwAU_15nt', 32, 1],

    ['UpsAcc_20nt', 'UpsAcc_20nt', 33, 1],
    ['UpAU_20nt', 'UpAU_20nt', 34, 1],
    ['DwsAcc_20nt', 'DwsAcc_20nt', 35, 1],
    ['DwAU_20nt', 'DwAU_20nt', 36, 1],

    ['UpsAcc_25nt', 'UpsAcc_25nt', 37, 1],
    ['UpAU_25nt', 'UpAU_25nt', 38, 1],
    ['DwsAcc_25nt', 'DwsAcc_25nt', 39, 1],
    ['DwAU_25nt', 'DwAU_25nt', 40, 1],

    ['UpsAcc_30nt', 'UpsAcc_30nt', 41, 1],
    ['UpAU_30nt', 'UpAU_30nt', 42, 1],
    ['DwsAcc_30nt', 'DwsAcc_30nt', 43, 1],
    ['DwAU_30nt', 'DwAU_30nt', 44, 1],


    ['Gtotal', 'Gtotal', 45, 1],
    ['Gnucl', 'Gnucl', 46, 1],
    ['UTR_length', 'UTR_length', 20, 1]
    );

# System State
my $SEED;   # Either Seed or seedless 
my $REGION; # Either 5UTR, CDS, or 3UTR
my $trainSpecies; # currently either human or mouse (8/19/2011)
                  # contains model that will be used for the prediction. 

if (scalar @ARGV == 2) {
    #$EnrichBase = $ARGV[0]; # not used
    $InFName= $ARGV[0];
    $trainSpecies = $ARGV[1];
} else {
    my $bsnme = basename($0);
    die("usage: $bsnme <Input File> <Training Species>\n");
}

open(INFILE, $InFName)
    or die("Could not open $InFName for input\n");

# Determine if the input file contains data or is just a header line.

@InputLines = <INFILE>;
close INFILE;

if (scalar (@InputLines) < 2 ) {
    # only header line present
    print STDERR "No data in input file $InFName output file not created\n";
    print STDERR "Exiting $0\n" if $debug;
    exit 0;
} # end if

my $bsename = basename($InFName);
$bsename =~ m/^([Ss].+)-([CU].+)-SiteFeatures-En-Seed-.+\.out$/;  # extract state from filename
($SEED, $REGION) = ($1, $2);

print STDERR "State is:\n\tseed = $SEED \n\tregion = $REGION \n\ttrainspecies = $trainSpecies\n" if $debug;

# Output file name format is fixed.  The R script
# That does logistic regression obtains information
# by parsing the filename of the file passed to it.
#

# creating output filename
my $tmpseed = lc($SEED);
my $tmpreg = $REGION;
$tmpreg =~ s/UTR3/3pUTR/;
$tmpreg =~ s/UTR5/5pUTR/;
my $tmptrain = lc($trainSpecies);

$OutFName = dirname($InFName) . "/" . "bsite-$tmpseed-$tmpreg-$tmptrain.txt";
open(OUTFILE,">$OutFName")
    or die("Could not open $OutFName for output\n");


# open table input file.

# throw away header line in input
my $input = shift @InputLines;

# print headerline
my @outList = ();
foreach my $afield (@colMapping) {
    if ($afield->[3] == 1) {
	push @outList, "\"$afield->[1]\""; # print col name 
    }
} # end foreach
print OUTFILE join("\t", @outList) . "\n";

# For each line after the header in the input file.

while (my $inline = shift @InputLines)
{
    chomp $inline;
    my @fields = split(':', $inline);
    my @outList = ();

    for (my $i = 0; $i < scalar @colMapping; $i++) {
	if ($colMapping[$i]->[3] == 1) {
	    push @outList, $fields[$i];
	}
    } # end for
    print OUTFILE join("\t", @outList) . "\n";;

} # end while <INLINE>

	    
print STDERR "Exiting $0\n" if $debug;

exit 0;


