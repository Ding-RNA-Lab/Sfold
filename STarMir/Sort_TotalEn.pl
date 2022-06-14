#!/usr/bin/perl -w

# Sort_TotalEn.pl
# a simple sort, to
# insure a relationship between site id numbers
# and the position of the site
#
# TODO:
# probably could be replaced with a carefully crafted
# call to unix sort
#
# 2022-4-14 (war)
# modified to remove calls to the carp module

use warnings;
use strict;

#no Carp::Assert;
use File::Basename;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals
my $InFname1; # output of Get_TotalEn.pl columns remapped
              # for web version
my $OutFname; # temp as we recopy the file on its old name in the 
              # next processing step. 

# Column Mapping for Energy file $InFname1
my $POS = 7;
my $MIRNAME = 4;


if (scalar @ARGV == 2) {
  $InFname1 = $ARGV[0];
  $OutFname = $ARGV[1];

} else {
    my $bsname = basename($0);
  die("usage: $bsname <Ener_file> <Output_file>\n");
}


open(INFILE,"<$InFname1") || die "Unable to open the $InFname1 file to read";

my @inRecords = ();;
# Read the input file into an array
while(<INFILE>) {
    chomp;
    push @inRecords, [split ':'];
} # end while

close INFILE;


# sort the records
my @sortRecords = sort {
     $a->[$MIRNAME] cmp $b->[$MIRNAME]
	 ||
    ((split('-',$a->[$POS]))[0]) 
	<=> 
    ((split('-', $b->[$POS]))[0])

} @inRecords;

open(OUTFILE,">$OutFname") || die "Unable to open the $OutFname file to write";

# print out sorted records
foreach my $arecord (@sortRecords) {
    print OUTFILE join(':', @{$arecord}) . "\n";
} # end foreach

close OUTFILE;


print STDERR "Exiting $0\n" if $debug;

exit 0;
