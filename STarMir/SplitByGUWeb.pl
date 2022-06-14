#!/usr/bin/perl

# SplitByGU.pl by William Rennie for Wadsworth Center 2011-08-9
#
# divide files into categories based on existence of seed regions.
# criteria is exactly 0 G-U pairings
#
# usage $ ./SplitByGU <Input File>
#
# Written for website
#
# For naming to work, this script should be called first
# followed by SplitByRegionWeb.pl
# 
use warnings;
use strict;

use File::Basename;
# no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n", @ARGV) . "\n" if $debug;

# Globals
my $InFname;        # SiteFeatures-En-Seed file
my $GUFIELD = 7;    # Column containing the count of GU
                    #   pairings


if (scalar @ARGV == 1) {
        $InFname = $ARGV[0];
} else {
    my $bsname = basename($0);
    die("usage:>perl $bsname <Sitefeatures-En-Seed file>\n");
}

# Open input file
open (SF, "<$InFname") 
    or die "Could not open $InFname for read.";

# Open output files
my $outSeedFname = dirname($InFname) . "/Seed-" . basename($InFname);
open (SEED,    ">$outSeedFname")     
    or die "Could not open $outSeedFname for write.";
my $outGuFname = dirname($InFname) . "/seedless-" . basename($InFname);
open (SEEDLESS,     ">$outGuFname")      
    or die "Could not open $outGuFname for write.";

# read in header
my $headerLine = <SF>;

# print headerline to each output file

print SEED $headerLine;
print SEEDLESS  $headerLine;

my $line = 1;
while (<SF>) {
    # do not chomp.
    my ($gu) = (split(/:/))[$GUFIELD];
    
    # Case 1:  GU count is exactly zero
    if ($gu ==  0) {
	print SEED $_ ;

    }
    # Case 2, no seed, any other value of GU
    else {
	print SEEDLESS $_;
    } # end else if

} # end while

close SF;
close SEED;
close SEEDLESS;

print STDERR "Exiting $0\n" if $debug;

exit 0;

