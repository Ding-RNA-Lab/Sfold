#!/usr/bin/perl -w

use warnings;
use strict;

use Cwd qw/ abs_path /;

# get the directory in which this script is located
$0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;


my $base_dir = abs_path($1);

#require "$base_dir/starmir_param.pl";

# This script is used by the STarMir module.
# This script takes TotalEn-Hyb-Fil-mirna.out from Dang's script

#die(" Usage: $0 <STarMir_output>\n\n") if $#ARGV != 0;
(my $file, my $siteNumber) = @ARGV;

my $target_name = "";
my $miRNA_name = "";

my $dG_total ;
my $dG_disrupt;
my $dG_hybrid;
my $dG_N;
my $miRNA_len;
my $tar_startpos;
my $tar_endpos;
my $tarconf1;
my $tarconf2;
my $tarconf3;
my $tarconf4;

my $mirnaid;

open(INFILE, $file) or die("file $file not found");
my $nsites = 0;
while (<INFILE>) {
  $nsites++;
  chomp;
  next if /^\s*$/;
  if ($nsites != $siteNumber) {
      next;
  }
  else{
  my(@e) = split(':');

  $target_name = $e[2] if !$target_name;
  # $miRNA_name = $e[4] if !$miRNA_name;
  # reading in this information even though its
  # not currently used to provide flexibility

  my $dG_total = $e[0];
  my $dG_disrupt = $e[1];
  my $dG_hybrid = $e[6];
  my $dG_N = $e[9];

  my $miRNA_len = $e[5];

  $e[7] =~ /(\d+)-(\d+)/;
  $tar_startpos =  $1;
  $tar_endpos = $2;
  $tarconf1 =  $e[10];
  $tarconf2 = $e[11];
  $tarconf3 = $e[12];
  $tarconf4 = $e[13];
  $mirnaid = $e[4];
  last;
  } # end else
 
} # end while



my $i = $nsites;

  # create png files for web display
  #my $cmd = "$base_dir/plot_Hybrid_Web.R $i $target_name $mirnaid[$i-1] \"$tar_startpos[$i-1]-$tar_endpos[$i-1]\" \"$tarconf1[$i-1]\" \"$tarconf2[$i-1]\" \"$tarconf3[$i-1]\" \"$tarconf4[$i-1]\" \"$dG_hybrid[$i-1]\" > /dev/null";

  #system($cmd);

  # create downloadable pdf files
  my $pdfcmd = "$base_dir/plot_Hybrid_PDF.R $i $target_name $mirnaid \"$tar_startpos-$tar_endpos\" \"$tarconf1\" \"$tarconf2\" \"$tarconf3\" \"$tarconf4\" \"$dG_hybrid\" > /dev/null";

  system($pdfcmd);

#} # end for each site

exit 0;
