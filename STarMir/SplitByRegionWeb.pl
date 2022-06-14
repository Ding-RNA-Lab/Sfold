#!/usr/bin/perl -w

# SplitByRegion.pl by Adam Wolenc for Wadsworth Center 2011-02-16
#
# divide all SiteFeatures-En-Seed- files into regions following the list files
#
# usage $ ./SplitByRegion.pl
#
# 2011-03-02 Adam.
#   Cols 3 and 4 are site range and seed range in absolute positions.
#   Moved relative positons of site and seed to two columns at the end.
#   One last column indicates region name.
#   Region name when site straddles boundary depends on where the site begins.
# 2011-8-4 (war) Modified for web
#          -- Added arguments for start and end of coding region
#          -- output to only three files
#             for UTR5, CDs and UTR3
#          -- Print headers for columns 
#             Taken from input file.
#


use strict;
use warnings;

use File::Basename;
# no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals
my $InFname;    # SiteFeatures-En-Seed file
                # Also includes a header line on Website
my $CDSSTART;   # Start of coding region
my $CDSEND;     # End of coding region
my $HEADERLINE; 

if (scalar @ARGV == 3) {
    $CDSSTART =  $ARGV[0];
    $CDSEND = $ARGV[1];
    $InFname = $ARGV[2];
} else {
    my $bsname = basename($0);
    die("usage:>perl $bsname <start of coding region> <end of coding region> <Sitefeatures-En-Seed file>\n");
}

# Open input file
open (SF, "<$InFname") 
    or die "Could not open $InFname for read.";

# Open output files
my $utr5Outfile = dirname($InFname). "/UTR5-" . basename($InFname);
open (UTR5,    ">$utr5Outfile")     
    or die "Could not open $utr5Outfile for write.";

my $cdsOutfile = dirname($InFname) . "/CDS-" . basename($InFname);
open (CDS,     ">$cdsOutfile")      
    or die "Could not open $cdsOutfile for write.";

my $utr3Outfile = dirname($InFname) . "/UTR3-" . basename($InFname);
open (UTR3,    ">$utr3Outfile")     
    or die "Could not open $utr3Outfile for write.";

# read in header
$HEADERLINE = <SF>;
chomp $HEADERLINE;
#Add the new fields
$HEADERLINE .= "\tAdj_Site_Pos\tAdj_Seed_Pos\tRegion";

# print headerline to each output file

print UTR5 "$HEADERLINE\n";
print CDS  "$HEADERLINE\n";
print UTR3 "$HEADERLINE\n";

my $line = 1;
while (<SF>) {
    # do not chomp.
    (my $seqname, my $site) = (split(/:/))[1,4];
    (my $st, my $en) = split(/-/, $site);
    
    # implicit assumption: cds region length is greater than or equal to site length.
    #assert($CDSEND - $CDSSTART >= 0) if DEBUG;
    #assert($en - $st >= 0) if DEBUG;
    
    # if site is longer than coding region, skip
#     if(($CDSEND - $CDSSTART) < ($en - $st))
#     {
# 	print
# 	    "Warning: " . $InFname . ":" . $line . " ~ CDS[" . $CDSSTART . ", " .
# 	    $CDSEND . "] < " . "site[" . $st . ", " . $en . "] \n";
# 	next;
#     } # end if site is larger than coding region
    
    # Case 1:  entire site is contained in UTR5 region
    if ($en <  $CDSSTART) {
	my $adjusted = subtract($_, 0);  # relative pos is same as absolute pos
	print UTR5 $adjusted . ":UTR5\n";

    } # end if UTR5

    # Case 2 site overlaps UTR5 and CD region
    elsif ($st <  $CDSSTART && $en >= $CDSSTART) {
	my $adjusted = subtract($_, 0);  # relative pos is same as absolute pos
	print UTR5 $adjusted . ":UTR5\n";
    } # end elseif overlap UTR5, CDS

    # Case 3 site contained in coding region
    elsif ($st >= $CDSSTART && $en <= $CDSEND) {
	my $adjusted = subtract($_, $CDSSTART - 1 );  # adjust to relative pos
	print CDS $adjusted . ":CDS\n";
    } # end elseif CDS
    
    # Case 4 site overlaps coding region and UTR3
    elsif ($st <= $CDSEND && $en >  $CDSEND) {
      my $adjusted = subtract($_, $CDSSTART -1);  # relative pos is same as absolute pos
	print CDS $adjusted . ":CDS\n";
    } # end elseif overlap CDS, UTR3

    # Case 5 site is contained in UTR3 region
    elsif ($st >  $CDSEND) {
      my $adjusted = subtract($_, $CDSEND);  # adjust to relative pos
      print UTR3 $adjusted . ":UTR3\n";
    } # end elseif UTR3

    # Technically its impossible to get here.
    else {
	print STDERR "Error: I couldn't figure out what region $site belongs to within the $seqname sequence.\n";
    } # end else
    $line++;
} # end while

close SF;

close UTR5;
close CDS;
close UTR3;

print STDERR "Exiting $0\n" if $debug;

exit 0;


sub subtract {
    my $record = shift;
    chomp $record;  # we will be adding records to the end.
    my $offset = shift;
    my @rec = split(/:/, $record);
    my $site_range = subtract_from_range($rec[4], $offset);
    my $seed_range = subtract_from_range($rec[5], $offset);
    return join(":", @rec, $site_range, $seed_range);
} # end sub subtract

sub subtract_from_range {
    my $range = shift;
    my $offset = shift;
    (my $be, my $en) = split(/-/, $range);
    $be -= $offset;
    $en -= $offset;
    return "$be-$en";
} # end sub subtract_from_range
