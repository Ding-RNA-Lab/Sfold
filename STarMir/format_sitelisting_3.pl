#!/usr/bin/perl -w

# format_sitelisting_3.pl
# Cloned from format_sitelisting_2.pl
#
# 10/4/2011 (war)
# removed all output but the ascii conformation
# 
# 12//17/2012 (war)
# added code to make phg images for hybrid conformation.  Images named
# for site, and written to the same directory as the site listing.
#
use warnings;
use strict;

use Cwd qw/ abs_path /;

# get the directory in which this script is located
$0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;


my $base_dir = abs_path($1);

require "$base_dir/starmir_param.pl";

# This script is used by the STarMir module.
# This script takes TotalEn-Hyb-Fil-mirna.out from Dang's script
# as input and produce the text output file for the web that
# contains the structural conformation

#die(" Usage: $0 <STarMir_output>\n\n") if $#ARGV != 0;

my $target_name = "";
my $miRNA_name = "";

my @dG_total = ();
my @dG_disrupt = ();
my @dG_hybrid = ();
my @dG_N = ();
my @miRNA_len = ();
my @tar_startpos = ();
my @tar_endpos = ();
my @tarconf1 = ();
my @tarconf2 = ();
my @tarconf3 = ();
my @tarconf4 = ();

my @mirnaid = ();
my $nsites = 0;
while (<>) {
    chomp;
  next if /^\s*$/;

  my(@e) = split(':');

  $target_name = $e[2] if !$target_name;
  # $miRNA_name = $e[4] if !$miRNA_name;
  # reading in this information even though its
  # not currently used to provide flexibility

  push @dG_total, $e[0];
  push @dG_disrupt, $e[1];
  push @dG_hybrid, $e[6];
  push @dG_N, $e[9];

  push @miRNA_len, $e[5];

  $e[7] =~ /(\d+)-(\d+)/;
  push @tar_startpos, $1;
  push @tar_endpos, $2;
  push @tarconf1, $e[10];
  push @tarconf2, $e[11];
  push @tarconf3, $e[12];
  push @tarconf4, $e[13];
  push @mirnaid, $e[4];

  $nsites++;
} # end while

die("No target site was identified based on the current energy thresholds.\n") if !$nsites;

print
	"Target name   = $target_name\n",
	"\n";

for (my $i=1; $i<=$nsites; $i++) {
  my $bp = $tarconf2[$i-1];
  $bp =~ s/[^\s]/|/g;
  my $my_spos = sprintf("%5d", $tar_startpos[$i-1]);
  my $my_len = sprintf("%5d", $miRNA_len[$i-1]);

  print <<EndOfOutput;
Site $i -- $mirnaid[$i-1]

    5'->3'        $tarconf1[$i-1]
    Target $my_spos  $tarconf2[$i-1]  $tar_endpos[$i-1]
                  $bp
    miRNA  $my_len  $tarconf3[$i-1]  1
    3'->5'        $tarconf4[$i-1]
 

EndOfOutput

  my @args; # arguments for the call to the R script

  # create png files for web display
  #my $cmd = "$base_dir/plot_Hybrid_Web.R $i $target_name $mirnaid[$i-1] \"$tar_startpos[$i-1]-$tar_endpos[$i-1]\" \"$tarconf1[$i-1]\" \"$tarconf2[$i-1]\" \"$tarconf3[$i-1]\" \"$tarconf4[$i-1]\" \"$dG_hybrid[$i-1]\" > /dev/null";

  #system($cmd);

  # create downloadable pdf files
  my $pdfcmd = "$base_dir/plot_Hybrid_PDF.R $i $target_name $mirnaid[$i-1] \"$tar_startpos[$i-1]-$tar_endpos[$i-1]\" \"$tarconf1[$i-1]\" \"$tarconf2[$i-1]\" \"$tarconf3[$i-1]\" \"$tarconf4[$i-1]\" \"$dG_hybrid[$i-1]\" > /dev/null";

  system($pdfcmd);

} # end for each site

exit 0;
