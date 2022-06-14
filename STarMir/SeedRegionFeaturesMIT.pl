#!/usr/bin/perl

# SeedRegionFeaturesMIT.pl by Adam Wolenc for Wadsworth Center 2010-07-29
# Re-write of SeedRegionFeatures.pl to match the seen conformation labeling system
# used by MIT.  Basically, only Poss 2-8 matter for determining seed. tA1
# is also relevant. 

# Modified 2010-12-10 Adam.  Added 4 output columns.
# 5    Seeded: Start and end Pos of the seed helix. Positions are on the mRNA target
#      Seedless: Start and end Pos of the 1st helix. Positions are on the mRNA target
# 6 (Primary name)
# 7 (GU pair)
# 8 (tA1)
# 9 (region 12-17)
# 10(region 4-15)
# 11   Pos on target aligning with 1st miRNA Pos from  2 up   to  8 that is paired -DASH-
#      Pos on target aligning with 1st miRNA Pos from  8 down to  2 that is paired.
# 12   Pos on target aligning with 1st miRNA Pos from  9 up   to 11 that is paired -DASH-
#      Pos on target aligning with 1st miRNA Pos from 11 down to  9 that is paired.
# 13   Pos on target aligning with 1st miRNA Pos from 12 up to miRNA end that is paired -DASH-
#      Pos on target aligning with 1st miRNA Pos from miRNA end down to 12 that is paired.
# In columns 11 to 13, if there are no pairs in the range, report 0.

# parsing a RNAhybrid result file (compact format) to get features of seed regions
# usage  >perl SeedRegionFeaturesMIT.pl <Input_File> <Output_File> 
#
# 7/29/11 (war)
# Modified for web

# 5/21/2012 (ccliu)
# fixed a bug for initial value of $curr_tar_pos. For few cases of binding sites, the first or last nt may be paired, original code only shifts one for $curr_tar_pos, but if both ends are paired, the shift should be two.
#
use warnings;
use strict;

#no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

use constant DEBUGGING => 0;  # set this to 1 to pause at every record with seed for manual QC

my $InFname;
my $OutFname;
if (scalar @ARGV == 2) {
    $InFname = $ARGV[0];
    $OutFname= $ARGV[1];
}
else {
   die("usage: >perl SeedRegionFeaturesMIT.pl <Input_File> <Output_File>\n");
}

open(IN,"<$InFname") || die "Unable to open the $InFname file to read";
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";

print OUT "Target name:Target length:miRNA name:Hybrid Energy:Site Pos:";
print OUT "Seed_or_1st_h_range:Primary name:GU pair:tA1:region 12-17:region 4-15";
print OUT "2-8_range:9-11_range:12-end_range\n";

my $record_num=0;
while (<IN>) {
  # 0  X66610-3pUTR
  # 1  317
  # 2  hsa-let-7a
  # 3  22
  # 4  -22.0
  # 5  0.488194
  # 6  36
  # 7  C    GUC        UAGA    C
  # 8   ACUA    ACCUACU    CUCA 
  # 9   UGAU    UGGAUGA    GAGU 
  #10  U    AUGU       UG       

  chomp;
  ++$record_num;
  my @rec = split(/:/);
  my $tar_name = $rec[0];
  my $tar_len = $rec[1];
  my $mir_name = $rec[2];
  my $mir_len = $rec[3];
  my $h_ener = $rec[4];
  my $rnah_be_tar = $rec[6];  # often this has an extra base that is not actually part of the site.

  my $tar_misses = lc(reverse($rec[7]));
  my $tar_matches = lc(reverse($rec[8]));
  my $mir_matches = lc(reverse($rec[9]));
  my $mir_misses = lc(reverse($rec[10]));

  my @tar_misses_c = split(//, $tar_misses);
  my @tar_matches_c = split(//, $tar_matches);
  my @mir_matches_c = split(//, $mir_matches);
  my @mir_misses_c = split(//, $mir_misses);

  my $tar_either = merge($tar_matches, $tar_misses);
  my $tar_either_no_space = $tar_either;
  $tar_either_no_space =~ s/ //g;
  my @tar_either_c = split(//, $tar_either);

  my $mir_either = merge($mir_matches, $mir_misses);
  my $mir_either_no_space = $mir_either;
  $mir_either_no_space =~ s/ //g;
  my @mir_either_c = split(//, $mir_either);

  my $conformation="other";

  # map all mir positions to positions on the diagram.
  my $i = 0;
  my $count = 0;
  if ($mir_len != length ($mir_either_no_space)) {
    print "Error. mir length ($mir_len) does not match diagram (". length ($mir_either_no_space) .") in line " . ($record_num+1) . "\n";
    next;
  }
  my @pos_map;
  my @pos_on_tar;
  for (my $i = 1; $i <= $mir_len; ++$i){
    $pos_map[$i] = -1;
    $pos_on_tar[$i] = 1;
  }
  # without a target length, the best we can do is a relative position
  my $curr_tar_pos = 1;
  # special case where first or last base in mRNA is paired, offset shifts by one
  $tar_misses =~ /[^ ]/;
  my $last_tar_miss = $-[0];
  $tar_matches =~ /[^ ]/;
  my $last_tar_match = $-[0];

  reverse($tar_misses) =~ /[^ ]/;
  my $first_tar_miss = $-[0];
  reverse($tar_matches)=~ /[^ ]/;
  my $first_tar_match = $-[0];

#  if ($last_tar_miss > $last_tar_match || $first_tar_miss > $first_tar_match) {
 #   $curr_tar_pos = 0;
 # }  
#if each end of mRNA is paired, offset shifts by one for each time. (by ccliu)
  if ($last_tar_miss > $last_tar_match) {
    $curr_tar_pos--;
  }
  if ($first_tar_miss > $first_tar_match) {
    $curr_tar_pos--;
  }
  my $lbound = 1;
  my $ubound = 1;
  while ($pos_map[$mir_len] == -1 && $i < length $mir_either) {
    if ($tar_either_c[$i] ne " ") {
      --$curr_tar_pos;
    }
    if ($mir_either_c[$i] ne " ") {
      ++$count;
      $pos_map[$count] = $i;
      if ($tar_matches_c[$i] ne " ") {
       $pos_on_tar[$count] = $curr_tar_pos;
        if ($ubound == 1) {
          $ubound = $curr_tar_pos;
        }
        $lbound = $curr_tar_pos;
      }
    }
    ++$i;
  }
  # calculate site length
  my $site_len = $ubound - $lbound + 1;
  # calculate site range. Translate relative positions to absolute positions on target.
  # also, calculate the range of the first helix.
  my $be_tar = -1;  # undefined
  my $end_tar = -1;  # undefined
  my $helix1_be_tar = -1;  # undefined.
  my $helix1_end_tar = -1;  # undefined.
  my $mir_gap_detected = 0;  # true after encountering first bulge on mirna
  for (my $i=1; $i<=$mir_len; ++$i) {
    if ($pos_on_tar[$i] == 1) {
      $pos_on_tar[$i] = 0;
      if ($helix1_be_tar != -1) {
        $mir_gap_detected = 1;
      }
    } else {
      # translate position relative to end of target to an absolute position on the target.
      $pos_on_tar[$i] += ($rnah_be_tar + $site_len +1);
      # maintain site begin and end values
      if ($end_tar == -1) {
        $end_tar = $pos_on_tar[$i];
      }
      $be_tar = $pos_on_tar[$i];

      # maintain 1st helix begin and end values
      if ($helix1_end_tar == -1) {
        # helix 1 is just beginning. We know the end pos of it
        $helix1_end_tar = $pos_on_tar[$i];
        $helix1_be_tar = $pos_on_tar[$i];
      } else {
        # update the beginning of first helix only if there have been no gaps on either side.
        if ($pos_on_tar[$i] == $helix1_be_tar - 1 && $mir_gap_detected == 0) {
          $helix1_be_tar = $pos_on_tar[$i];
        }
      }
    }
  }
  # assert($be_tar >= 1) if DEBUG;
  # assert($end_tar >= $be_tar) if DEBUG;
  # assert($helix1_be_tar >= 1) if DEBUG;
  # assert($helix1_end_tar >= $helix1_be_tar) if DEBUG;
  # assert($helix1_be_tar >= $be_tar) if DEBUG;
  # assert($helix1_end_tar <= $end_tar)if DEBUG;

  # calculate special positions on target for Chaochun
  my $beg_e = 2;
  my $beg_s = 8;
  my $mid_e = 9;
  my $mid_s = 11;
  my $end_e = 12;
  my $end_s = $mir_len;
  ++$beg_e while ($pos_on_tar[$beg_e] == 0 && $beg_e < 8);
  --$beg_s while ($pos_on_tar[$beg_s] == 0 && $beg_s > 2);
  ++$mid_e while ($pos_on_tar[$mid_e] == 0 && $mid_e < 11);
  --$mid_s while ($pos_on_tar[$mid_s] == 0 && $mid_s > 9);
  ++$end_e while ($pos_on_tar[$end_e] == 0 && $end_e < $mir_len);
  --$end_s while ($pos_on_tar[$end_s] == 0 && $end_s > 12);
  $beg_e = $pos_on_tar[$beg_e];
  $beg_s = $pos_on_tar[$beg_s];
  $mid_e = $pos_on_tar[$mid_e];
  $mid_s = $pos_on_tar[$mid_s];
  $end_e = $pos_on_tar[$end_e];
  $end_s = $pos_on_tar[$end_s];
  
  next if ($pos_map[1] == -1);  #hopeless.
  # Look for A in Pos 1 of target. It can be either miss or match.
  my $tA1 = (substr($tar_either, $pos_map[1], 1) eq "a") ? 1 : 0;

  # gather all possible seed annotations.
  my %possible_annotation;
  # walk through the diagram to determine seed len and pos
  for (my $gu_allowed = 0; $gu_allowed <= 7; ++$gu_allowed) {
    for (my $start_pos = 2; $start_pos <= 8; ++$start_pos) {

      my $mirlen = 0;  # to start
      my $gu_remaining = $gu_allowed;  # to start
      my $GU_count = 0;

      for (my $i = $pos_map[$start_pos]; $i <= $pos_map[8] && $i != -1; ++$i) {

        my $tar_match_ch = $tar_matches_c[$i];
        my $mir_match_ch = $mir_matches_c[$i];

        #there are 3 possibilities. 0 means space, 1 means base
        #tar_miss
        # mir_match
        #  mir_miss
        #010-WC = match, mir len increments
        #010-GU, some GU's left = match, mir len increments, gus left decrments
        #010-GU, no GU's left = essentially a bulge in both, exit loop
        #?0? = mismatch or bulge. exit loop.

        if ($mir_match_ch ne " ") {
          if (($tar_match_ch eq "g" && $mir_match_ch eq "u") ||
              ($tar_match_ch eq "u" && $mir_match_ch eq "g")) {
            #010-GU
            if ($gu_remaining) {
              #010-GU, some GU's left
              ++$mirlen;
              ++$GU_count;
              --$gu_remaining;
            } else {
              #010-GU, no GU's left
              last;
            }
          } else {
            #010-WC
            ++$mirlen;
          }
        } else {
          #?0?
          last;
        }

        my $end_pos = $start_pos + $mirlen - 1;
        my $label = "$mirlen-mer $start_pos-$end_pos $GU_count GU";
        #store this as a possible annotation. Value is the range on the target.
        my $seed_be_tar  = $pos_on_tar[$start_pos + $mirlen - 1];
        my $seed_end_tar = $pos_on_tar[$start_pos];
        $possible_annotation{$label} = $seed_be_tar."-".$seed_end_tar;

      }
    }
  }

  # walk through the diagram to determine if a 4-mer exists in 12-17
  # start counting at 1st match, stop at 1st bulge after 1st match.
  # gather all possible seed annotations.
  my %possible_twelve;
  # walk through the diagram to determine seed len and pos
  for (my $start_pos = 12; $start_pos <= 17; ++$start_pos) {

    my $mirlen = 0;  #to start

    for (my $i = $pos_map[$start_pos]; $i <= $pos_map[17] && $i != -1; ++$i) {

      my $tar_match_ch = $tar_matches_c[$i];
      my $mir_match_ch = $mir_matches_c[$i];

      #there are 3 possibilities. 0 means space, 1 means base
      #tar_miss
      # mir_match
      #  mir_miss
      #010-WC = match, mir len increments
      #010-GU, essentially a bulge in both, exit loop
      #?0? = mismatch or bulge. exit loop.

      if ($mir_match_ch ne " " && 
        !(($tar_match_ch eq "g" && $mir_match_ch eq "u") ||
          ($tar_match_ch eq "u" && $mir_match_ch eq "g"))) {
        #010-WC
        ++$mirlen;
      } else {
        #?0? or 010-GU
        last;
      }

      my $end_pos = $start_pos + $mirlen - 1;
      my $label = "$mirlen-mer $start_pos-$end_pos";
      $possible_twelve{$label} = $mirlen;
    }
  }

  my $comp_12_17 = 0;
  if (defined $possible_twelve{"4-mer 12-15"}
   || defined $possible_twelve{"4-mer 13-16"}
   || defined $possible_twelve{"4-mer 14-17"}) {
    $comp_12_17=1;
  }

  # close-up of mirna Poss 4-15 only:
  my $four_len = 0;
  my $comp_4_15 = 0;
  if ($pos_map[15] != -1) {
    my $len_4_15 = $pos_map[15] - $pos_map[4] + 1;
    my $tar_misses_4_15 = substr($tar_misses, $pos_map[4], $len_4_15);
    my $tar_matches_4_15 = substr($tar_matches, $pos_map[4], $len_4_15);
    my $mir_matches_4_15 = substr($mir_matches, $pos_map[4], $len_4_15);
    my $mir_misses_4_15 = substr($mir_misses, $pos_map[4], $len_4_15);

    my $four_start = 0;
    my $four_end = 0;
    my $pos_on_mir = 4;  # to start

    # walk through the diagram to determine if 4-15 are WC
    # start counting at 4 if it is a match. Stop at any gap.
    for (my $i = 0; $i < $len_4_15; ++$i) {

      my $tar_miss_ch = substr($tar_matches_4_15, $i, 1);
      my $tar_match_ch = substr($tar_matches_4_15, $i, 1);
      my $mir_match_ch = substr($mir_matches_4_15, $i, 1);
      my $mir_miss_ch = substr($mir_misses_4_15, $i, 1);

      #there are 2 possibilities. 0 means space, 1 means base
      #tar_miss
      # mir_match
      #  mir_miss
      #010-WC = match, mir pos increments
      #other = any exception precludes the 4-15 conformation
  
      if ($mir_match_ch ne " " && 
        !(($tar_match_ch eq "g" && $mir_match_ch eq "u") ||
          ($tar_match_ch eq "u" && $mir_match_ch eq "g"))) {
        #case 010-WC
        #if we haven't started to track the seed, begin now
        $four_start = $pos_on_mir if ($four_start==0);
        #if we are already tracking the seed, continue now
        $four_end = $pos_on_mir;
        ++$pos_on_mir;
      } else {
        last;
      }
    }
    $four_len = ($four_end - $four_start + 1) if ($four_start != 0);
    $comp_4_15=1 if ($four_len >= 12);
  }

  # look for primary conformations
  my $GU_count = -1;
  my $gu_allowed = 0;
  my $seed_range_on_tar = "$helix1_be_tar-$helix1_end_tar";  # default, seedless case.
  while ($conformation eq "other" && $gu_allowed <= 9) {
    if ($tA1 == 1 && $possible_annotation{"7-mer 2-8 $gu_allowed GU"}) {
      $seed_range_on_tar = $possible_annotation{"7-mer 2-8 $gu_allowed GU"};
      $conformation = "8mer";
      $GU_count=$gu_allowed;
    } elsif ($tA1 == 0 && $possible_annotation{"7-mer 2-8 $gu_allowed GU"}) {
      $seed_range_on_tar = $possible_annotation{"7-mer 2-8 $gu_allowed GU"};
      $conformation = "7mer-m8";
      $GU_count=$gu_allowed;
    } elsif ($tA1 == 1 && $possible_annotation{"6-mer 2-7 $gu_allowed GU"}) {
      $seed_range_on_tar = $possible_annotation{"6-mer 2-7 $gu_allowed GU"};
      $conformation = "7mer-A1";
      $GU_count=$gu_allowed;
    } elsif ($tA1 == 0 && $possible_annotation{"6-mer 2-7 $gu_allowed GU"}) {
      $seed_range_on_tar = $possible_annotation{"6-mer 2-7 $gu_allowed GU"};
      $conformation = "6mer";
      $GU_count=$gu_allowed;
    } elsif ($possible_annotation{"6-mer 3-8 $gu_allowed GU"}) {
      $seed_range_on_tar = $possible_annotation{"6-mer 3-8 $gu_allowed GU"};
      $conformation = "offset-6mer";
      $GU_count=$gu_allowed;
    }
    ++$gu_allowed;
  }

  if (DEBUGGING) {
    print $rnah_be_tar;
    print $rec[7] . "\n";
    print "" . (' ' x length $rnah_be_tar) . $rec[8] . "\n";
    print "" . (' ' x length $rnah_be_tar) . $rec[9] . "\n";
    print "" . (' ' x length $rnah_be_tar) . $rec[10] . "\n";
    print     "" . join(":", $tar_name, $tar_len, $mir_name, $h_ener, "$be_tar-$end_tar",
                        $seed_range_on_tar, $conformation, $GU_count, $tA1, $comp_12_17, $comp_4_15,
                        "$beg_s-$beg_e", "$mid_s-$mid_e", "$end_s-$end_e") . "\n";
  }

  print OUT "" . join(":", $tar_name, $tar_len, $mir_name, $h_ener, "$be_tar-$end_tar",
                      $seed_range_on_tar, $conformation, $GU_count, $tA1, $comp_12_17, $comp_4_15,
                      "$beg_s-$beg_e", "$mid_s-$mid_e", "$end_s-$end_e") . "\n";

  if (DEBUGGING) {
    for (my $p = 1; $p < scalar @pos_map; ++$p) {
      print "$p\t" . $pos_map[$p] . "\t" . $pos_on_tar[$p] . "\n";
    }
    print "--------------\n";
  
    print "" . <STDIN>;
  }
}
close IN;
close OUT;

print STDERR "Exiting $0\n" if $debug;

exit 0;

sub merge {
  # fill in blank spaces in $line1 with characters from
  #  $line2, resulting in a single, merged sequence.
  # line1 and line2 must be the same length for this to
  #  make any sense.
  my $line1 = shift;
  my $line2 = shift;
  my $ret = "";
  for (my $i = 0; $i < length $line1; ++$i) {
    if (substr($line1, $i, 1) ne " ") {
      $ret.=substr($line1, $i, 1);
    } else {
      $ret.=substr($line2, $i, 1);
    }
  }
  return $ret;
}
