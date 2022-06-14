#!/usr/bin/perl -w

# GetSiteSeqFeatures.pl
# Website version
#
#
# 2010-12-15 modification of GetSiteSeqFeatures.pl for dealing with MIT seeds
# Using files produced by ParseHybFile_2.pl and SeedRegionFeaturesMIT.pl
# Append the following 3 fields
# Site Accessibility (using positions from the output of SeedRegionFeaturesMIT.pl)
# Seed Accessibility (as originally calculated by GetSiteSeqFeatures.pl)
# Seed Accessibility (using positions from the output of SeedRegionFeaturesMIT.pl)
# For each of 6 window sizes (5, 10, 15, 20, 25, 30), append the following 4 fields
#   Upstream openness
#   Upstream AU %
#   Downstream openness
#   Downstream AU %
# Openness is defined as average of a range of relevant values in the sstrand.out

# 2011-02-11 Adam. Adapted script for full-length mRNA. no 30 nt of CDS

# 2011-8-11 war
# Modifications to fit into Website
#  - take the associated sstrand filename in
#    as an argument.
#
# 2022-4-14 (war)
# modified for distribution
# removed calls to Carp::Assert package

# Output Key
#     0  NM_019683  \
#     1  966         | Site ID taken from Sites- input file
#     2  let-7a      |
#     3  -18.8      / 
#     4  547        \
#     5  577         > Site begin, end, and openness based on Seed- input file.
#     6  0.281      /
#     7  573        \
#     8  576         > OLD METHOD Seed begin, end, and openness based on Sites- input file.
#     9  0.418      /   (best 4nt within positions 1 to 10)
#    10  571        \
#    11  577         > NEW METHOD Seed begin, end, and openness based on Seed- input file.
#    12  0.242      /   (matches MIT seed helix, or first helix if no MIT seed)
#    13  0.360      \
#    14  0.200       | Upstream Openness, Upstream AU, Downstream Openness, Downstream AU
#    15  0.900       | (5nt window)
#    16  0.000      /
#    17  0.283      \
#    18  0.400       |
#    19  0.556       | (10nt window)
#    20  0.200      /
#    21  0.437      \
#    22  0.533       |
#    23  0.594       | (15nt window)
#    24  0.400      /
#    25  0.386      \
#    26  0.500       |
#    27  0.667       | (20nt window)
#    28  0.450      /
#    29  0.431      \
#    30  0.520       |
#    31  0.604       | (25nt window)
#    32  0.440      /
#    33  0.432      \
#    34  0.500       |
#    35  0.618       | (30nt window)
#    36  0.433      /

use warnings;
use strict;

#no Carp::Assert;
my $debug = 0;

use File::Basename;

use constant DEBUGGING => 0;  # 0 disables debug mode

print STDERR "Entering $0\n" if $debug ;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

my $USAGEMESG = "usage: >$0 <ParseHyb_Input_File> <SeedRegionMIT_Input_File> <sstrand.out file> <output file>\n";

my $PHFFname = $ARGV[0] or die($USAGEMESG);
my $SRMFname = $ARGV[1] or die($USAGEMESG);
my $SSFname = $ARGV[2] or die($USAGEMESG);
my $OutFname = $ARGV[3] or die($USAGEMESG);


# # seq. accessibility files' directory prefix and file postfix
# my $FolderPre = "";
# my $StrFpost = "output/sstrand.out";
# my $PathToSfoldOutput = ("/local/biorna/bxm19/Human/Human-Full/");
# my @SubfolderList = ("Hsa_FULL_1", "Hsa_FULL_2","Hsa_FULL_3","Hsa_FULL_4");
# my $Folder_post = "/";

 my %Address;
# # load sequence lengths and locations of sstrand.out files.
# for my $subfolder (@SubfolderList)  {
#   my $LFname = join("",$PathToSfoldOutput,$subfolder,".txt");
#   open(LFILE,"<$LFname") || die "Unable to open the $LFname file to read";
#   <LFILE>;  # skip the first line
#   while (<LFILE>) {
#     chomp;
#     my $seqname = (split(/\t/))[0];
#     assert(!defined $Address{$seqname});  # duplicates not allowed.
#     $Address{$seqname} = join("", $PathToSfoldOutput, $subfolder, "/", $seqname, $Folder_post, $StrFpost);
#   }
#   close LFILE;
# }

# Produce one record of output for each record of input.
open(PHF,"<$PHFFname") or die "Unable to open the $PHFFname file to read";
open(SRM,"<$SRMFname") or die "Unable to open the $SRMFname file to read";
open(OUT,">$OutFname") or die "Unable to open the $OutFname file to write";

my $old_name = "";
my @bases;
my @acc;

<SRM>;  # SRM has a header row. PHF does not.

while (<PHF>) {
  #Example PHF record:
  #   0  NM_019683
  #   1  966
  #   2  miR-21
  #   3  22
  #   4  -15.4
  #   5  n,G_253-AG_22,21
  #   6  m,GAUAUU_254,259-UUGUAG_20,15
  #   7  n,U_260-UC_14,13
  #   8  m,UUUGGU_261,266-AGACUA_12,7
  #   9  n,_267-UU_6,5
  #  10  m,GCUG_267,270-CGAU_4,1
  #  11  n,A_271-_0
  # ...as many fields as there are diagram features.
  my @phf_rec = split(/:/);

  my $srm_line = <SRM>;
  die ("$PHFFname and $SRMFname do not have an equal number of lines") if (!$srm_line);
  my @srm_rec = split(/:/, $srm_line);
  #Example SRM record:
  #   0  NM_019683
  #   1  966
  #   2  let-7a
  #   3  -18.8
  #   4  547-577
  #   5  571-577
  #   6  7mer-m8
  #   7  1
  #   8  0
  #   9  0
  #  10  0
  #  11  571-577
  #  12  556-570
  #  13  547-555
  my ($MIT_site_be, $MIT_site_en) = split(/-/, $srm_rec[4]);

  #assert($MIT_site_be >= 1);  # non-positive position
  #assert($MIT_site_en >= $MIT_site_be);

  my ($MIT_seed_be, $MIT_seed_en) = split(/-/, $srm_rec[5]);

  #assert($MIT_seed_be >= 1);  # non-positive position
  #assert($MIT_seed_en >= $MIT_seed_be);
  #assert($MIT_seed_be >= $MIT_site_be);
  #assert($MIT_seed_en <= $MIT_seed_en);

  my $seqname = $phf_rec[0];

  # # load sstrand values for this sequence.
  # if (!defined $Address{$seqname}) {
  #   print "$seqname is not in gene list\n";
  #   next;
  # }

  # don't bother loading if this is the one currently loaded.
  if ($seqname ne $old_name) {
      $old_name = $seqname;

      # if (-e $Address{$seqname}) {
      # 	  open(SS,"<$Address{$seqname}");
      # } else { 
      # 	  print "Error: with $seqname no $Address{$seqname} file\n";
      # 	  next;
      # }
      
      open(SS, $SSFname) or die("Could not open $SSFname for input\n");

      # TODO(Adam): what to do about upstream openness in positions less than 30 ?
      my $array_idx = 0;
      while (my $line = <SS>) {
	  # skip nothing.
	  ($bases[$array_idx], $acc[$array_idx]) = (split(/\s+/,$line))[2,4];
	  ++$array_idx;
      }
      close SS;

      # add 30 open A's to this.
      for (my $a = 0; $a < 30; ++$a) {
	  $bases[$array_idx] = "A";
	  $acc[$array_idx] = 1.0;
	  ++$array_idx;
      }
  } # end if new sequence name

  # output what we have so far.
  print OUT $seqname,":",$phf_rec[1],":",$phf_rec[2],":",$phf_rec[4],":";

  # compute site accessibility (new way does not include dangling bases.)
  my $beg = $MIT_site_be;
  my $end = $MIT_site_en;
  
  #assert($beg > 0);
  #assert($end >= $beg);

  my $Site_Acc = 0;
  print "\t" if (DEBUGGING);
  for (my $i = $MIT_site_be-1; $i <= $MIT_site_en-1; ++$i) {
    $Site_Acc += $acc[$i];
    print $bases[$i] if (DEBUGGING);
  }
  $Site_Acc /= ($end - $beg + 1);

  print "\t" if (DEBUGGING);
  # compute seed accessibility (the new way)
  my $MIT_Seed_Acc = 0;
  for (my $i= $MIT_seed_be-1; $i<= $MIT_seed_en-1; ++$i) {
    $MIT_Seed_Acc += $acc[$i];
    print $bases[$i] if (DEBUGGING);
  }
  $MIT_Seed_Acc /= ($MIT_seed_en - $MIT_seed_be + 1);
  print "\n" if (DEBUGGING);

  # compute seed accessibility (the old way)
  # Adam deduced function from Dang code: For each helix of length 4 or more
  # which has at least 4 nt that align with miRNA position 10 or less, compute
  # the average single-stranded probability of each 4-nucleotide subsequence
  # for which all positions align with miRNA position 10 or less. Return the
  # maximum of all averages computed in this way.
  # Also track the s_be and s_end of the winning subsequence.
  my $Seed_Acc = 0;
  my $Acc4mer  = 0;
  my $s_be     = 0;
  my $s_end    = 0;
  # for each input field from 5 (the first feature to end
  for (my $i= 5; $i < scalar @phf_rec; ++$i) {
    # only look at m records (helixes), not n records (bulges).
    next if ($phf_rec[$i] =~ /n/);

    # an m record is like    m,GAUAUU_254,259-UUGUAG_20,15
    # get the length of the feature
    next if ($phf_rec[$i] !~ /-/);
    my $temp2 = (split(/-/, $phf_rec[$i]))[1];  # UUGUAG_20,15
    next if ($temp2 !~ /_/);
    my $temp3 = (split(/_/, $temp2))[1];  #        20,15
    next if ($temp3 !~ /,/);
    my ($m_pos1, $m_pos2) = (split(/,/, $temp3))[0,1];
    # skip this if it is less than 3 or more than 10
    next if ((($m_pos1 - $m_pos2) < 3) || (($m_pos2 + 3) > 10));

    $temp2 =  (split(/-/, $phf_rec[$i]))[0];  # m,GAUAUU_254,259
    next if ($temp2 !~ /_/);
    $temp3 = (split(/_/, $temp2))[1];  #          254,259
    next if ($temp3 !~ /,/);
    my ($pos1, $pos2) = (split(/,/, $temp3))[0,1];

    my $k = $m_pos1;
    # 2010-12-08 Adam Repaired off-by-one error. (@acc is 0 bound not 1 bound.)
    for (my $j = $pos1-1; $j <= ($pos2-3)-1; ++$j) {
      if ($k <= 10) {
        $Acc4mer = ($acc[$j] + $acc[$j+1] + $acc[$j+2] + $acc[$j+3]) / 4;
        if ($Acc4mer > $Seed_Acc) {
		      $Seed_Acc = $Acc4mer;
          $s_be = $j;
          $s_end = $j + 3;
	      }
	    }
      --$k;
	  }
  }

  # incremental print
  print OUT $beg,":",$end,":";
  printf OUT "%4.3f:", $Site_Acc;
  print OUT $s_be,":",$s_end,":";
  printf OUT "%4.3f:", $Seed_Acc;
  print OUT $MIT_seed_be,":",$MIT_seed_en,":";
  printf OUT "%4.3f:", $MIT_Seed_Acc;
  
  # compute and print openness and au content figures
  # 2010-12-09 Adam. New edge behavior. Use up to 30nt from the CDS region, and
  #  up to 30 imaginary A's, which are open, after the end of the 3UTR.
  #  WARNinG: Window size greater than 30 has undefined behavior.
  my $step = 5;
  my $num_windows = 6;
  #assert($step*$num_windows <= 30);  # script will require modification if this assertion is false.
  for (my $win = $step; $win <= $step*$num_windows; $win += $step) {
    my $upO  = 0;
    my $dwO  = 0;
    my $upAU = 0;
    my $dwAU = 0;

    #assert($win > 0);
    #assert($win <= 30);
    #assert($beg > 0);
    #assert($end > 0);
  
    # calculate upstream openness and upstream AU
    #     upstreamBINDINGSITEdownstream
    #     |      ||         ||        |      |
    #  0  up_s   up_e        dn_s     dn_e   en
    #             up_e+1    dn_s-1 
    my $up_s = ($beg-$win) -1;
    my $up_e = ($beg-1)    -1;
    my $dn_s = ($end+1)    -1;
    my $dn_e = ($end+$win) -1;

    my $adjusted_win = $win;  # account for the fact that sometimes there are not 30 upsteam bases...
    if ($up_s < 0) {
      # undefined, and possible.
      # we really don't know what to do in this case.
      $up_s = 0;
      $up_e = 0 if ($up_e < 0);
      $adjusted_win = $up_e - $up_s + 1;
    }

    #assert(($up_e - $up_s + 1) == $adjusted_win);
    #assert(($dn_e - $dn_s + 1) == $win);
    #assert($up_s >= 0);
    #assert($up_s <= $up_e);
    #assert($up_e < $dn_s);
    #assert($dn_s <= $dn_e);
    #assert($dn_e <= scalar @acc);
    
    for (my $j= $up_s; $j <= $up_e; ++$j) {
      #assert($j >= 0);
      #assert($j < scalar @acc);
      $upO += $acc[$j];
      if ($bases[$j] =~ /A|U/) {
        ++$upAU;
      }
    }
    $upO  /= $adjusted_win;
    $upAU /= $adjusted_win;

    # calculate downstream openness and downstream AU
    for (my $j= $dn_s ; $j <= $dn_e; ++$j) {
    
	#assert($j >= 0);
	#assert($j < scalar @acc);
      
	$dwO += $acc[$j];
	if ($bases[$j] =~ /A|U/) {
	    ++$dwAU;
	}
    }
    $dwO  /= $win;
    $dwAU /= $win;
    
    # incremental print
    printf OUT "%4.3f:%4.3f:%4.3f:%4.3f:",$upO,$upAU,$dwO,$dwAU;
  }
  print OUT "\n";
} # end while lines in PHF file

close OUT;
close PHF;
if (<SRM>) {
  print "ERROR: $PHFFname and $SRMFname do not have an equal number of lines\n";
}
close SRM;

print STDERR "Exiting $0\n" if $debug;

exit 0;
