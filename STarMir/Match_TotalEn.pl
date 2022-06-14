#!/usr/bin/perl -w

# Match_TotalEn.pl
#
# a script to match TotalEner values in the first file
# to the sites in the second file
# usage > perl Match_TotalEn.pl <Ener_file> <Site_file>  
# -----------------------------------------------------

# 2010-12-29 Adam. strict, cleanup, warnings, etc.
# 2011-01-03 Adam. produce a warning for each difference in set
# 2011-8-2 (war) adaptation for web site
#    Mainly, allowing caller to specifiy output file name. to
#    give some flexibility and consistency in web scripting.
#
# 2011-8-4 (war) Web adaptation
#   remapped columns on energy input file
#   to reflect difference in web version of
#   CalcNE
#
# 2022-4-14 (war)
# removed calls to the Carp module

use warnings;
use strict;

#no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals
my $InFname1; # output of Get_TotalEn.pl columns remapped
              # for web version
my $InFname2; # SiteFeatures-  GetSiteSeqFeaturesMIT.pl
my $OutFname;
  # $OutFname =~ s/SiteFeatures/SiteFeatures-En/;

# Column Mapping for Energy file $InFname1
my $TOTALEN = 0;
my $DISEN = 1;
my $NUCLEN = 8;
my $POS = 7;
my $SEQNAM = 2;
my $SEQLEN = 3;
my $HYBEN = 6;

# Column Mapping for SiteFeaturesFile $InFname2
my $SEQNAM2 = 0;
my $SEQLEN2 = 1;
my $HYBEN2 = 3;
my $BEG = 4;
my $EN = 5;


if (scalar @ARGV == 3) {
  $InFname1 = $ARGV[0];
  $InFname2 = $ARGV[1];
  $OutFname = $ARGV[2];

} else {
  die("usage:>perl Match_TotalEn.pl <Ener_file> <Site_file> <Output_file>\n");
}


open(IN1,"<$InFname1") || die "Unable to open the $InFname1 file to read";
open(IN2,"<$InFname2") || die "Unable to open the $InFname2 file to read";
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";

my %G = ();
my $siteId = 1;  # site id corresponding to line in TotalEn . . . 
while (<IN1>) {
  #   0  NM_175116-3pUTR
  #   1  706
  #   2  let-7a
  #   3  -6.549
  #   4  -3.718
  #   5  -20.4
  #   6  8-37
  #   7  GACUA,UUGAU_8,12
  #   8  GC,UG_14,15
  #   9  ACU,UGG_17,19
  #  10  ACUGCU,UGAUGG_23,28
  #  11  UCA,AGU_35,37
  chomp;
  my ($seqname, $seqlen, $totalEn, $nuclEn, $hybEn, $pos) = (split(/:/))[$SEQNAM,$SEQLEN,$TOTALEN,$NUCLEN,$HYBEN,$POS];
  # $seqname =~ s/-3pUTR$//;
  my ($beg, $end) = (split(/-/, $pos))[0,1] ;

  #assert($beg >= 1);  # non-positive position
  #assert($end >= $beg);

  my $key = join(":", $seqname, $seqlen, $hybEn, $beg, $end);

  #assert(!defined $G{$key});  # enforce unique keys

  $G{$key}->{'total'} = $totalEn;
  $G{$key}->{'nucl'}  = $nuclEn;
  # Add generated site ID number. This number links subsequent files to the 
  # conformation in sitelisting.out
  $G{$key}->{'SiteId'} = $siteId++;
  $G{$key}->{'used'}  = 0;
}
close IN1;

my $numline = 0;
while (<IN2>) {
  #   0  NM_175116
  #   1  706
  #   2  let-7a
  #   3  -20.4
  #   4  8
  #   5  37
  #   6  0.495
  #   7  23
  #   8  26
  #   9  0.791
  #  10  35
  #  11  37
  #  12  0.496
  #  ...
  #  36  0.633

  # Output:
  #   0  1  
  #   0  NM_175116
  #   1  706
  #   2  let-7a
  #   3  -20.4
  #   4  8
  #   5  37
  #   6  0.495
  #   7  23
  #   8  26
  #   9  0.791
  #  10  35
  #  11  37
  #  12  0.496
  # ...
  #  36  0.633
  #  37  -6.549
  #  38  -3.718
  chomp;
  my ($seqname, $seqlen, $hybEn, $beg, $end) = (split(/:/))[$SEQNAM2,$SEQLEN2,$HYBEN2,$BEG,$EN];

  #assert($beg >= 1);  # non-positive position
  #assert($end >= $beg);

  my $key = join(":", $seqname, $seqlen, $hybEn, $beg, $end);
  if (defined $G{$key}) {
    print OUT $G{$key}->{'SiteId'}, ':', $_, $G{$key}->{'total'}, ":", $G{$key}->{'nucl'}, "\n";
    $G{$key}->{'used'} = 1;
    ++$numline;
  } else {
    print STDERR "$key appears in SiteFeat file but not in TotalEn file.\n";
  }
}

close IN2;
close OUT;

for my $key (keys %G) {
  if (!($G{$key}->{'used'})) {
    print STDERR "$key appears in TotalEn file but not in SiteFeat file.\n";
  }
}

print "Matched $numline lines\n" if $debug;

print STDERR "Exiting $0\n" if $debug;

exit 0;
