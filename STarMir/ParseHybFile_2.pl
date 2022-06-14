#!/usr/bin/perl 

# ParseHybFile_2.pl
#
# parsing a RNAhybrid result file (compact format) 
# usage  >perl ParseHybFile_2.pl <Input_File>
#
# 7/29/11 (war)
# Created web version
# Following changes:
#    Output file name now taken as a parameter
#
# 2010-12-08 Adam made this code readable.
# Refactored repeated code, improved variable naming, and added many comments
# Each output field represents a helix or a bulge.
# n means this is field describes a bulge.  m means this is field describes a helix.
# There is no difference between n and m field other than this.
# Format is :
#    [n or m],[ information about the target side ]-[ information about the mirna side ]:
# Information about the target or mirna follows the same format:
#    [ sequence fragment on the target or mirna ]_[ range on target or mirna ]
# In the case of bulges, one sequence or the other may be blank.
# If any sequence is 2 nt or greater, the beginning and ending positions are shown,
# separated by a comma.
# If any sequence is 0 or 1 nt, only the beginning position is shown, without a comma.

use strict;
use warnings;

#no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

my $InFname = shift @ARGV or die("usage: >perl ParseHybFile_2.pl <Input File> <Output File>\n");
my $OutFname = shift @ARGV or die("usage: >perl ParseHybFile_2.pl <Input File> <Output File>\n");


open(IN,"<$InFname") || die "Unable to open the $InFname file to read";
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";

while (<IN>) {
  # 0 NM_019683-3pUTR
  # 1 966
  # 2 miR-21
  # 3 22
  # 4 -15.4
  # 5 1.000000
  # 6 253
  # 7 G      U             A
  # 8  GAUAUU  UUUGGU  GCUG
  # 9  UUGUAG  AGACUA  CGAU
  #10 AG      UC      UU    
  chomp;
  my @rec = split(/:/);

  my $sbj_name = $rec[0];
  my $sbj_len  = $rec[1];

  # if subject_name is like "XX|1234", remove the XX|
  if ($sbj_name =~ /\|/) {
    $sbj_name = (split(/\|/, $sbj_name))[1];
  }

  my $q_name   = $rec[2];
  my $q_len    = $rec[3];
  my $hyb_ener = $rec[4];
  my $be_tar   = $rec[6];

  # output what we have so far. this is amended later.
  print OUT $sbj_name,":",$sbj_len,":",$q_name,":",$q_len,":",$hyb_ener,":";

  my @tar_mm = split(//, $rec[7]);
  my @tar_m = split(//, $rec[8]);
  my @mir_m = split(//, $rec[9]);
  my @mir_mm = split(//, $rec[10]);

  my $tar_fragment= "";
  my $mir_fragment = "";

  # if first charater of mir matches is blank, start in mode m
  # otherwise, start in mode n.
  my $curr_mode = 'n';
  $curr_mode = "m" if ($mir_m[0] =~ /\w/);
  # by default, previous mode is the current mode.
  my $prev_mode = $curr_mode;

  my $tar_frag_start = $be_tar;
  my $tar_frag_end = $be_tar;
  my $mir_frag_start = $q_len;
  my $mir_frag_end = $q_len;

  # For each position from 0 to length of diagram (all 4 elements of diagram have the same length)
  # compose helixes and output encoded helix information when appropriate.
  for (my $i= 0; $i < scalar @mir_m; ++$i) {
    if ($mir_m[$i] =~ /\w/) {
      $curr_mode = "m";
    } else {
      $curr_mode = "n";
    }

    # if mode has changed, do output and reset tar_fragment and/or mir_fragment.
    # Otherwise, accumulate more tar_fragment and mir_fragment
    if ($curr_mode ne $prev_mode) {
      # do output and reset
      do_output($prev_mode, $tar_fragment, $mir_fragment, \$tar_frag_start, \$mir_frag_start);
      $tar_fragment = "";
      $mir_fragment = "";
      $prev_mode = $curr_mode;
    } else {
      # Update mode for next iteration
      if ($mir_m[$i] =~ /\w/) { 
        $curr_mode = "m";
      } else {  
        $curr_mode = "n";
      }
    }
    # either way, add to tar_fragment and/or mir_fragment
    if ($mir_m[$i] =~ /\w/) { 
      $tar_fragment .= $tar_m[$i];
      $mir_fragment .= $mir_m[$i];
    } else {  
      if ($tar_mm[$i] =~ /\w/) {
        $tar_fragment .= $tar_mm[$i];
      }
      if ($mir_mm[$i] =~ /\w/) {
        $mir_fragment .= $mir_mm[$i];
      }
    }

    # In special case of the last character in the diagram, always do one last output.
    if ($i == $#mir_m) {
      do_output($prev_mode, $tar_fragment, $mir_fragment, \$tar_frag_start, \$mir_frag_start);
    }
  }

  print OUT "\n";
}

close IN;
close OUT;

print STDERR "Exiting $0\n" if $debug;

exit 0;

sub do_output {
  my $prev_mode = shift;
  my $tar_fragment = shift;
  my $mir_fragment = shift;
  my $tar_frag_start = shift;
  my $mir_frag_start = shift;
  my $tar_frag_end = $$tar_frag_start + length($tar_fragment)-1;
  my $mir_frag_end = $$mir_frag_start - length($mir_fragment)+1;
  # attempts have been made to reproduce the original, subtle behavior of output without altering it.
  if ($tar_frag_end > $$tar_frag_start) {
    print OUT $prev_mode,",",$tar_fragment,"_",$$tar_frag_start,",",$tar_frag_end,"-";
  } else {  
    print OUT $prev_mode,",",$tar_fragment,"_",$$tar_frag_start,"-";
  }
  if ($mir_frag_end < $$mir_frag_start) {
    print OUT $mir_fragment,"_",$$mir_frag_start,",",$mir_frag_end,":";
  } else {  
    print OUT $mir_fragment,"_",$$mir_frag_start,":";
  }
  #assert($$tar_frag_start > 0);
  $$tar_frag_start = $tar_frag_end+1;
  $$mir_frag_start = $mir_frag_end-1;
}
