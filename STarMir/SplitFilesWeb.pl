#!/usr/bin/perl 

# SplitFiles.pl -- William Rennie 2011-8-12
#
# divide SiteFeatures-En-Seed- files into six files
#   required by logit.R
# Add additional data files needed by logit to
# The split files.
#
# runs SplitbyGU.pl
# runs SplitbyRegionWeb.pl
#
# usage $ ./SplitFilesWeb.pl <path to scripts> <path to starmir_params> <SiteFeatures-En-Seed file>  
#                             <coding region start> <coding region end> 
#

use strict;
use warnings;

# no Carp::Assert;
use File::Basename;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals
my $InFname;    # SiteFeatures-En-Seed file
                # Also includes a header line on Website
my $CDSSTART;   # Start of coding region
my $CDSEND;     # End of coding region
my $TRAINSPECIES; # Training Species
my $BASEDIR;    # Base directory for scripts
my $PARAMDIR;    # starmir_params directory

if (scalar @ARGV == 6) {
    $BASEDIR = $ARGV[0];
    $PARAMDIR = $ARGV[1];
    $InFname = $ARGV[2];
    $CDSSTART = $ARGV[3];
    $CDSEND = $ARGV[4];
    $TRAINSPECIES = $ARGV[5];
} else {
    my $bsname = basename($0);
    die("usage $bsname <Base directory for scripts> <starmir param directory> <SiteFeatures-En-Seed file>  <coding region start> <coding region end> <training species>");
}


# Break the input file by
# region 5UTR|CDS|3UTR

# assert($CDSSTART >=0 && $CDSEND >=0) if DEBUG;
# assert($CDSSTART <= $CDSEND) if DEBUG;

my @args = ('perl',"$BASEDIR/SplitByRegionWeb.pl", $CDSSTART, $CDSEND, $InFname);

system(@args);

die("@args: Failed to execute\n") if($? == -1);
die("@args: Exited with value $?\n") if ($? > 0);

# Now Split each of the created files by seed, seedless

@args = ('perl',"$BASEDIR/SplitByGUWeb.pl"
	 , dirname($InFname) . "/" . "UTR3-" . basename($InFname));

system(@args);

die("@args: Failed to execute\n") if($? == -1);
die("@args: Exited with value $?\n") if ($? > 0);

@args = ('perl',"$BASEDIR/SplitByGUWeb.pl"
	 , dirname($InFname) . "/" . "CDS-" . basename($InFname));

system(@args);

die("@args: Failed to execute\n") if($? == -1);
die("@args: Exited with value $?\n") if ($? > 0);

@args = ('perl',"$BASEDIR/SplitByGUWeb.pl"
	 , dirname($InFname) . "/" . "UTR5-" . basename($InFname));

system(@args);

die("@args: Failed to execute\n") if($? == -1);
die("@args: Exited with value $?\n") if ($? > 0);

#Add the additional data columns to the six created files

foreach my $type ('Seed', 'seedless') {
    foreach my $region ('UTR5', 'CDS', 'UTR3') {

 	@args = ('perl', "$BASEDIR/AddCols.pl"
		 , dirname($InFname) . "/" . "$type-$region-" . basename($InFname)
		 , $TRAINSPECIES);
	
	print STDERR "Executing:\n" if $debug;
	print STDERR join('|',@args). "\n" if $debug;

	system(@args);

 	die("@args: Failed to execute\n") if($? == -1);
 	die("@args: Exited with value $?\n") if ($? > 0);
	
    } # end foreach region
} # end foreach type

print STDERR "Exiting $0\n" if $debug;

exit 0;
