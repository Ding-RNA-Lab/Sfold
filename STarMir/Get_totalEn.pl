#!/usr/bin/perl -w
#This script is used by the STarMir module.
# matching disrupt ener. files to the hybrid ener. file 
# of the possible binding sites
#
# 2011-02-04 Adam. Total rewrite to correct faults in matching.
#  This script is now a functional combination of the scripts Match_disEn_Hyb and Match_TotalEn
#  part of the Prediction Steps protocol.

# usage  >perl Get_totalEn.pl <Hyb_Ener_File> <Pos_Adj> <bsites_File> <disen_File> <output_file>
#
# 2011-12-8 (war)
# Final web cleanup. 
#
use Cwd qw/ abs_path /;
use warnings;

#no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# get the directory in which this script is located
my $base_dir = `dirname $0`;
chomp $base_dir;
$base_dir = abs_path($base_dir);
# load parameters contained in option file
require "$base_dir/starmir_param.pl";
#-------------------------------------------

my $skip_filter = 0;
if (scalar @ARGV >= 5) {
  $FilNuclFname = $ARGV[0];
  $Pos_Adj = $ARGV[1];
  $BsiteFname = $ARGV[2];
  $DisEnFname = $ARGV[3];
  $OutFname = $ARGV[4];
  $skip_filter = $ARGV[5] if defined ($ARGV[5]);  # if true, skip both the init energy filter
                                                  #  and the overlap filter.
} else {
  die(
		"args: " . (scalar @ARGV) . ": [@ARGV]\n" .
		"usage: >perl Get_totalEn.pl <Hyb_Ener_File> <Pos_Adj> <bsites_File> <disen_File> <output_file>\n"
	); 
}

# load disEn data. Key is position range.
my %disEnDB;
open(BSITES, "<$BsiteFname") or die("Could not open $BsiteFname for read.");
open(DISENS, "<$DisEnFname") or die("Could not open $DisEnFname for read.");
while (my $bsite = <BSITES>) {
  my $disen = <DISENS>;
  die("Disen and Bsites files should have the same number of records.") if ($disen eq "");
  # modification because disen files SOMETIMES have a leading space in records which causes
  # the wrong value to be taken 7-10-2012 (war)
  $disen =~ s/^\s+//;
  my $disruptEn = (split(/\s+/, $disen))[1];
  my ($beg, $end) = split(/\s+/, $bsite);
  $beg -= $Pos_Adj;
  $end -= $Pos_Adj;
  my $key = sprintf "%08d-%08d", $beg, $end;
  $disEnDB{$key} = $disruptEn;
}
close BSITES;
close DISENS;

# associate disen data with each record.
# load all records and store in hash of arrays with key of total en so that the best non-overlapping
# records may be selected.
my %recordsByTotalEn;
my $range_ubound = 0;
open(FILNUC,"<$FilNuclFname") || die "Unable to open the $FilNuclFname file to read";
while (<FILNUC>) {
  # example data
  #  0	NM_133985_CD300
  #  1	476
  #  2	miR-1907
  #  3	22
  #  4	-20.6
  #  5	276-298
  #  6	A  GG    UGUG    AAG   C A
  #  7	 GC  CCAG    CCUC    GC C 
  #  8	 UG  GGUC    GGAG    CG G 
  #  9	   GA    UA      ACGA  A  
  # 10	4-276,298
  # 11	-1.396
  chomp;
  tr/\r//d;

  my $site = {};
  $site->{'full_record'} = $_;
  my @rec = split(/:/);

  my $range = $rec[5];
  my ($beg, $end) = split(/-/, $range);
  $range_ubound = $end if ($end > $range_ubound);
  $site->{'beg'} = $beg;
  $site->{'end'} = $end;

  $site->{'hybEn'} = $rec[4];
  $site->{'nuclEn'} = $rec[11];

  my $key = sprintf "%08d-%08d", $beg, $end;
  if (defined $disEnDB{$key}) {
    $site->{'disEn'} = $disEnDB{$key};
    $site->{'totalEn'} = $site->{'hybEn'} + $site->{'disEn'};  # output from disruptEn program is negated
    # (equivalent to the more logically clear expression (($hybEn) - (-$disEn)) )

    push @{$recordsByTotalEn{$site->{'totalEn'}}}, $site;

  } else {
    die("Could not find $key in the disEn results.");
  }
}
close FILNUC;

print "range_ubound=$range_ubound\n"; 
# output only the best overlapping records, where best is defined as largest totalen  
my @occupied;
for (my $i = 1; $i <= $range_ubound; ++$i) {
  $occupied[$i] = 0;
}
# walk through the records in order of totalen, and only output those whose ranges
# are totally unoccupied. Maintain occupied array.
# if skip_filter is true, output unconditionally even if there is overlap.
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";
for my $sitesKey (sort { $a <=> $b } keys %recordsByTotalEn) {
  for my $site (@{$recordsByTotalEn{$sitesKey}}) {

    if ($skip_filter || $site->{'nuclEn'} + $NE_threshold < 0) {
      my $accept = 1;
      for (my $i = $site->{'beg'}; $i <= $site->{'end'}; ++$i) {
        if ($occupied[$i] == 1) {
          $accept = 0;
          next;
        }
      }
  
      if ($skip_filter || $accept == 1) {
        # Example output (convoluted legacy format)
        #  0	-23.500
        #  1	-6.100
        #  2	NM_133985_CD300
        #  3	476
        #  4	miR-1907
        #  5	22
        #  6	-29.6
        #  7	205-243
        #  8	-7.007
        #  9	4.09
        # 10	A      CCUC   GAGU   GUUGGUAGAUAC       A
        # 11	 GCCUCC    AGA    UCC            GCUGCUC 
        # 12	 UGGAGG    UCU    AGG            CGACGAG 
        # 13	                     AGA                 
        # 14	4-205,243
        my @rec = split(/:/, $site->{'full_record'});
        printf OUT "%.3f:%.3f:", $site->{'totalEn'}, -($site->{'disEn'});
        print  OUT join(":", @rec[0..5,11]) . ":";
        printf OUT "%.3f:", ($site->{'nuclEn'} + $NE_threshold);
        print  OUT join(":", @rec[6..10]) . ":\n";

        for (my $i = $site->{'beg'}; $i <= $site->{'end'}; ++$i) {
          $occupied[$i] = 1;
        }
      }  # else it conflicts with an already accepted site
    }  # else its nucleation energy does not overcome initiation energy.
  }
}

close OUT;

print STDERR "Exiting $0\n" if $debug;

exit 0;
