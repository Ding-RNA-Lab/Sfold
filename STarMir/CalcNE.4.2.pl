#!/usr/bin/perl -w

# usage: % perl CalcNE.4.2.pl <AlignResult_filename.txt> <list_filename> <-Option>

# CalcNE.3.0.pl  by Dang Long 10/27/04
# modified 1/30/06
# to calculate nucl.energy for WU-Blast ResHyb based on UTR struct. output of sfold
# sequences_ UTRs added with CD nts;
# 
# 3/01/05
# New modification: 
# - calculate the N.E. for each 3-bp window, and chose the minimum energy of the windows (as MNE2)
# - no mismatch in the nucleation sites. 
# - the list file has the CD ranges -> only calculte NE of the ResHyb after that ranges
# - calculate the N.E. for the longest "open" window, and chose the minimum energy of the windows (as BNE)
# defined values for free energies from Turner website version 3.0 
# Stacking energies (base pairs)

# 2010-03-19 (Adam): CalcNE.4.1.pl  optimized for speed.

# 2010-12-22 Adam code audit. Cleanup, warnings etc.
#  Strip out the 'new modification' above, which uses 3nt windows. This is
#  not used by the rest of the process, and has lead to confusion.
# 2010-12-28 to 2010-12-31 Adam adapted script to accept Res-*-Hyb-1.out instead of *-hit_1.out
#  Set the script to ignore all helixes shorter than 4nt (was 2).
#  Directly produce Fil-0-...csv file, obviating FilterBNE.pl
#  Made $mdata structure one level deeper, grouped by conformation (RNAhybrid record).
#   (This revealed a design flaw whereby nucleation energy leaked across conformations.)
#  Produce Hyb-Fil-0 file instead of Fil-0-...-hit_1_NE.csv file. This repairs the leakage problem.
#  An earlier script produces bsites file independently. No need for Match_NE_Hyb now.
#  Removed the SNE and BNE features as they do not make sense in the new paradigm.
#  MNE is now simply NE. Removed 0-0 seed range placeholder.
#  Original code had implicit assumption that the first open position would
#   always be the start of the longest open range. In rare cases, this
#   was not true.
# 2011-01-03 (Adam) This script to implements the following model for Nucl En
# 1. Possible miRNA binding sites are miRNA-mRNA complementary matches of at
#     least four contiguous matches anywhere in the miRNA. Each such site is
#     refered to as a "Helix" in group discussions (but not in the publication)
# 2. For each such helix, identify the particular single-stranded 4nt block 
#     within the helix that would form the most stable 4-base-pair duplex with 
#     the miRNA.
# 3. Identify open blocks in the target site using sample of 1000 structures
#     from Sfold. Consider each structure separately. Calculate the average
#     nucleation potential over the whole set of structures in the sample.

use strict;
use warnings;

# no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

#globals

#SE . see stack.dat
#    5' --> 3'       5' --> 3'
#       01     equiv    23
#       32      <=>     10
#    3' <-- 5'       3' <-- 5'
our %SE = (
  'AAUU' , -0.9, 'ACGU' , -2.2,
  'AGCU' , -2.1, 'AGUU' , -0.6,
  'AUAU' , -1.1, 'AUGU' , -1.4,

  'CAUG' , -2.1, 'CCGG' , -3.3,
  'CGCG' , -2.4, 'CGUG' , -1.4,
  'CUGG' , -2.1,

  'GAUC' , -2.4, 'GAUU' , -1.3,
  'GCGC' , -3.4, 'GCGU' , -2.5,
  'GGUC' , -1.5,
  'GGUU' , -0.5,
  'GUGU' ,  1.3,

  'UAUA' , -1.3, 'UAUG' , -1.0,
  'UGUG' ,  0.3);
#compute mirror
for my $key (keys %SE) {
  my $alt = substr($key, 2, 2) . substr($key, 0, 2);
  $SE{$alt} = $SE{$key};
}

our %dsfe_memo;
# sizes of a nucleation site
use constant NUCL_LEN => 4;
# number of structures to be processed
use constant NSTRUCT => 1000;

# -------------------------------------------------------

#2011-02-02 Adam. Adapted to web server.
my $HybOutFile;
my $Pos_Adj;
my $NuclFname;
my $SFold_outdir;
if (scalar @ARGV == 4) {
    $HybOutFile = $ARGV[0];
    $Pos_Adj = $ARGV[1];
    $NuclFname = $ARGV[2];
    $SFold_outdir = $ARGV[3];
} else {
   die("usage: %perl CalcNE.4.2.pl <HybOutFile> <Pos_Adj> <NuclFname> <SFold_outdir>\n");
}

# read the "ResHyb" file into the memory 
print "Reading ResHyb file into memory.\n";
open(RESHYB, "<$HybOutFile") || die "Unable to open the $HybOutFile file to read";

#2010-03-05 (Adam): Use hash table rather than array to improve lookup time.
# O(n + m) rather than O(n * m)

#mdata refers to hashmap where:
#  key is UTR name, value is hashmap where:
#    key is full ResHyb record, value is array where:
#      each element is a hashmap describing an interesting helix.
my $mdata = {};

while (my $line = <RESHYB>) {
  #   0  NM_024451-3pUTR
  #   1  1288
  #   2  let-7a
  #   3  22
  #   4  -22.6
  #   5  702-725
  #   6  CUGUGC,GAUAUG_702,707
  #   7  CCU,GGA_711,713
  #   8  UUGCCUC,GAUGGAG_719,725

  #2010-03-19 (Adam): Hot spot.  Optimized
  chomp $line;
  my @record = split(/:/, $line);
  #fields 6 and beyond contain information about one helix each
  for (my $helix_idx = 6; $helix_idx < scalar @record - 4; ++$helix_idx) {
    next if ($record[$helix_idx] eq "");
    #parse data about this helix
    my $data = {};
    $data->{'initial'} = join(":", @record[0..5, ($#record-3)..$#record]);  # initial part of each record, reproduced verbatim. 

    my ($helix_seq, $helix_pos) = split(/_/, $record[$helix_idx]);
    ($data->{'utr_st'}, $data->{'utr_end'})  = split(/,/, $helix_pos);
    ($data->{'utr_hit'}, $data->{'mir_hit'}) = split(/,/, $helix_seq);
    $data->{'utr_hit'} = reverse($data->{'utr_hit'});
    $data->{'mir_hit'} = reverse($data->{'mir_hit'});

    push @{$mdata->{$line}}, $data;
  }
} # end while $line =<RESHYB>
close RESHYB;
print "Reading ResHyb file into memory. Done.\n";

# reading the list and struct. files and calc. N.E.
print "Reading list and struct. files and calc. N.E.\n";
open(OUTFILE, ">$NuclFname")|| die "Unable to open the $NuclFname file to write";
my $fname = "$SFold_outdir/bp.out";
if (!-e $fname||-z $fname) {
  print "$fname does not exist\n";
}

print "Loading structures from $fname.\n";
#is_paired array of NSTRUCT hashes
# each hash contains one element for each position that is paired
# key is position, value is 1.  Unpaired positions are undefined.
my @is_paired;
for (my $i = NSTRUCT - 1; $i >= 0;--$i) {
  $is_paired[$i] = {};
}
open(BPDOTOUT, "<$fname") || die "Unable to open file $fname";
my $j = -1;
while (my $line = <BPDOTOUT>) {
  #2010-03-19 (Adam): hot spot.  optimized expressions in this block.
  if ($line =~ /tructure/o) {
    ++$j;
    next;
  } else  {
    $line =~ /([0-9]+)\s+([0-9]+)/o;
    $is_paired[$j]->{$1} = 1;
    $is_paired[$j]->{$2} = 1;
  }
}
close BPDOTOUT;

print "Calculating NE\n";
#2010-03-05 (Adam): O(n) lookups changed to to O(1) via hashtable.
for my $sitekey ( keys %$mdata ) {
  # build a helix string (this is an auxiliary task that is convenient to perform here)
  my $helix_list_str;  # one field for each helix, composed as we iterate through them, output verbatim
  for my $od ( @{$mdata->{$sitekey}} ) {
    # append to helix_list_str
    $helix_list_str .= ":"
      . reverse($od->{'utr_hit'}) . "," . reverse($od->{'mir_hit'}) . "_"
      . $od->{'utr_st'}  . "," . $od->{'utr_end'};
  }

  # calculate the average of best NE for each structure
  my $NE_total = 0.0;
  # for each structure
  for my $structure_num (0..(NSTRUCT-1)) {
    my $NE = 0.0;  # Within structure, the best NE, to take part in average.
                   # If there is no open block of 4nt, NE is default.
                   # Default is 0.0 (per Ye in 2011-01-04 email).
     # for each helix in the rnahybrid diagram
    for my $od ( @{$mdata->{$sitekey}} ) {
      if (length($od->{'mir_hit'}) >= NUCL_LEN) {
        # find most favorable block in this helix for this structure
        my $utr_st  = $od->{'utr_st'}  + $Pos_Adj;
        my $utr_end = $od->{'utr_end'} + $Pos_Adj;
        for (my $offset = $utr_st;
                $offset <= $utr_end - NUCL_LEN + 1;
              ++$offset) {
         # entire block must be open
          my $block_is_fully_open = 1;
          for (my $pos_on_cd300 = $offset;
                  $pos_on_cd300 < $offset + NUCL_LEN;
                ++$pos_on_cd300) {
            if (defined $is_paired[$structure_num]->{$pos_on_cd300}) {
              # block is not open. stop looking.
              $block_is_fully_open = 0;
              last;
            }
          }
          if ($block_is_fully_open) {
            my $dG = dsfe(substr($od->{'utr_hit'}, $offset - $utr_st, NUCL_LEN),
                          substr($od->{'mir_hit'}, $offset - $utr_st, NUCL_LEN));
            # update best NE if this one is better.
            $NE = $dG if ($dG < $NE);
          }
        }
      }  # else too short
    } # end for each helix
    # maintain sum
    $NE_total += $NE;
  } # end for each structure
  # average best NE for all structures
  my $avg_best_NE = $NE_total / NSTRUCT;

  # output one line for each conformation (previously, this script output one line per helix)
  my $site_id = (@{$mdata->{$sitekey}})[0]->{'initial'};
  my $dumb_range = "4-" 
    . (@{$mdata->{$sitekey}})[0]->{'utr_st'} . "," 
    . (@{$mdata->{$sitekey}})[$#{$mdata->{$sitekey}}]->{'utr_end'};
  printf OUTFILE "$site_id:$dumb_range:%4.3f\n", $avg_best_NE;
} # end while <GENELIST>

close OUTFILE;
print "Reading list and struct. files and calc. N.E. Done.\n";

print STDERR "Exiting $0\n" if $debug;

exit 0;

#2010-03-05 (Adam): Memoize.
sub dsfe {
  (my $be3str, my $be5str) = @_;  # the UTR strand from 3' end and the mir strand from 5' end
  if (defined $dsfe_memo{"$be3str$be5str"}) {
    return $dsfe_memo{"$be3str$be5str"};
  }
  my $fe = 0;
  my $seq_len = length($be3str);
  my @strbe3 = split(//, $be3str);
  my @strbe5 = split(//, $be5str);
  for (my $i = 0; $i < $seq_len - 1; ++$i) {
    my $dinucl = join("", $strbe5[$i], $strbe5[$i+1], $strbe3[$i+1], $strbe3[$i]);
    print "$dinucl is not defined.\n" if (!defined $SE{$dinucl});
    $fe += $SE{$dinucl};
  }
  $dsfe_memo{"$be3str$be5str"} = $fe;
  return $fe;
}

