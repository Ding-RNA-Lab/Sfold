#!/usr/bin/perl -w

# 2010-12-30 Adam Wolenc for Wadsworth Center.
# Make binding_sites files from Res-Hyb files for input to disruptEn program
# This task was once performed by Match_NE_Hyb.pl, but it didn't produce all
# the sites for some reason, and it depended on CalcNE output. No one knows why.
# Now, running disruptEn does not require that CalcNE be run prior.
# usage ./Bsites.pl <ResHyb_File> <length_offset> <Output_file>

# 2011-02-02 Adam. Adapted for web server.

# 10/18/2011 (war)
# Added standard tracking
#

# 2022-4-14 (war)
# removed calls to the Carp::Assert package
#
use strict;
use warnings;

#no Carp::Assert;
my $debug = 0;

use File::Basename;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

my $ResHybFname;
my $CDlen;
my $BsiteOutFname;
if (scalar @ARGV == 3) {
  $ResHybFname = $ARGV[0];
  $CDlen = $ARGV[1];
  $BsiteOutFname = $ARGV[2];
} else {
  die("usage: ./Bsites.pl <ResHyb_File> <length_offset> <Output_file>\n");
}

open(IN,"<$ResHybFname") || die "Unable to open the $ResHybFname file to read";
my %site;
while (my $line = <IN>) {
  chomp $line;
  my @rec = split(/:/, $line);
  my $genename = $rec[0];
  $genename =~ s/-3pUTR//;
  my $pos = $rec[5];
  my ($st, $en) = split(/-/, $pos);
  $st += $CDlen;
  $en += $CDlen;
  my $key = sprintf "%08d-%08d", $st, $en;
  $site{$key} = "$st\t$en\n";  # prevent duplicates.
}
close IN;

#output in order with no duplicates.
open(OUT,">$BsiteOutFname") || die "Unable to open the $BsiteOutFname file to write";
for my $key (sort keys %site) {
  print OUT $site{$key};
}
close OUT;

print STDERR "Exiting $0\n" if $debug;

exit 0;
