#!/usr/bin/perl
#usage: >perl run.logit.pl <species> <trainSpecies> <cdsBegin> <cdsEnd>
# The script is used for computation of conservation score and relative percentage 
# of siteBegin in each component, then based on the features computed, functional 
# probabilities (linear and nonlinear) of each binding site are computed using 
# a logistic prediction model.

use warnings;
#use strict;

use Carp::Assert;

print STDERR "\nEntering $0\n" if DEBUG;
print STDERR "Arguments:\n\t" if DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if DEBUG;

### DEBUG
exit 0;

$0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;
my $basePath = $1;
my $Rargs = "--slave --no-save --no-restore --no-environ --silent --args";
my $RexePath = "/usr/local/amd64/bin/";
$RexePath = "" if !-e "${RexePath}R";

die("usage:>perl run.logit.pl <species> <trainSpecies> <cdsBegin> <cdsEnd>\n") unless scalar @ARGV == 4;
my ($species, $trainSpecies, $cdsBegin, $cdsEnd) = @ARGV;

#the file of conservation score for species "mouse" and "human"
$conserFile = "${basePath}logitCoef/Refseq-conserScore-$species.txt";


my $wdir = `pwd`;

@types = ("seed","seedless");
@components = ("3pUTR","CDS","5pUTR");
foreach $type (@types) {
    foreach $component (@components) {
	
	$bindsiteBase = "bsite-$type-$component-$trainSpecies.txt";
	$bindsiteFile = "${basePath}$bindsiteBase";
	$logitCoefFile = "${basePath}logitCoef/logReg-coef-$type-$component-$trainSpecies.txt";
	
	
	my $exitStatus = system("${RexePath}R $Rargs bindsiteFile=$bindsiteFile outputFile=Prob.$bindsiteBase conserFile=$conserFile logitCoefFile=$logitCoefFile cdsBegin=$cdsBegin cdsEnd=$cdsEnd component=$component type=$type < ${basePath}cnsvScore_prediction.R");
	
    } # end foreach component
} # end foreach type

print STDERR "Exiting $0\n" if DEBUG;

exit 0;
