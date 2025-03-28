#!/usr/bin/perl -w

#
# This script takes a target mRNA, identifies ribozyme cleavage
# sites, and calculates energies for each of the sites.
#

use Getopt::Std;
use Cwd qw/ cwd /;
use Cwd qw/abs_path/;

use strict;
use warnings;

# the binary for the comand line version should also be located in Sfold-main/bin
# this ensures that we can deduce the location of all the
# helper scripts and executables.
my($bindir) = `dirname $0`;
chomp($bindir);
require "$bindir/parse_lib.pl";
my($sfold_root) = abs_path ("$bindir/../");

# customize locations of external programs/directories
my($loc_dG_hybrid) = "$bindir/rz_dghybrid.pl";
my($loc_dG_switch) = "$bindir/rz_dgswitch.pl";
my($loc_dG_disrupt) = "$bindir/disruptEn";

# the fact that sfold root was hard coded worries
# me.  In the sfold code base this is supposed to
# only happen in config files.

my($sfold_param) = join("", $sfold_root, "/param/");
die " Error: unable to locate folder '$sfold_param'"
  if (!-e $sfold_param);
my($sfold_bin) = join("", $sfold_root, "/bin/");
die " Error: unable to locate folder '$sfold_bin'"
  if (!-e $sfold_bin);

# added an option for an output file as it was suggested that
# most users will not be able to handle output redirection.

my $opt_string = "a:b:f:m:s:t:o:";
my %opt = ();
my  $tgseq = "";
my($tgseqlen);
my($h3len) = 11;
my($h1len) = 9;
my($triplet) = "GUC";
my(@sites) = ();
my(@startpos) = ();
my($sitelen) = 21;
# number of matched triplets at the 5' end of target that do not
# have long-enough arms to be included as a target site
my($num_dropped_sites_5p) = 0;

my($Rz_core) = "CUGAUGAGUCCGUGAGGACGAA";
my(@Rz_seq, @Rz_5p, @Rz_3p);

my(@dG_hybrid) = ();
my(@dG_switch) = ();
my(@dG_disrupt) = ();

&getopts ($opt_string, \%opt) || &usage;
if (!defined($opt{m}) || !defined($opt{s}) || !defined($opt{f}) || !defined($opt{o})) {
  &usage;
}

my($bpfile) = $opt{s};
my($fefile) = $opt{f};
my($outfile)= $opt{o};
# read target sequence
undef $/;
open(SEQFILE, "$opt{m}") || die "Unable to open '$opt{m}'";
my($input) = <SEQFILE>;
$tgseq = &parseFASTA($input);
$tgseq =~ tr/T/U/;
$tgseqlen = length($tgseq);
close(SEQFILE);
$/ = "\n";

# Helix III length (including the first 2 base pairs of the triplet)
$h3len = $opt{a} if (defined $opt{a});
die " Error: invalid h3len" if ($h3len !~ /^\d+$/);

# Helix I length
$h1len = $opt{b} if (defined $opt{b});
die " Error: invalid h1len" if ($h1len !~ /^\d+$/);

# total length of target site
$sitelen = $h3len + $h1len + 1;

# Preferred cleavage triplet for hammerhead ribozyme
$triplet = $opt{t} if (defined $opt{t});
die " Error: invalid triplet" if ($triplet !~ /^[ACGU]{1}U[ACU]{1}$/);

# find the number of matched triplets we have to skip before 
# we have one with full arms. this helps identify the "site ID"
# reported in the final output
for (my $i=0; $i<$h3len-2; $i++) {
  $num_dropped_sites_5p++ if (substr($tgseq, $i, 3) eq "$triplet");
}

# identify hits in the target sequence and calculate energies
my($list_pos) = "";
my($list_len1) = "";
my($list_len2) = "";

for (my $i=0; $i<=($tgseqlen-$sitelen); $i++) {
  my($mysite) = substr($tgseq, $i, $sitelen);
  if ( substr($mysite, $h3len-2, 1) eq substr($triplet, 0, 1) &&
       substr($mysite, $h3len-1, 1) eq substr($triplet, 1, 1) &&
       substr($mysite, $h3len, 1) eq substr($triplet, 2, 1) ) {
    push @startpos, $i+1;
    push @sites, $mysite;

    # determine Rz sequence
    # reverse complement the two arms
    my($my_Rz_3p, $my_Rz_5p);
    $my_Rz_3p = reverse(substr($mysite, 0, $h3len));
    $my_Rz_3p =~ tr/AUCG/UAGC/;
    $my_Rz_5p = reverse(substr($mysite, $h3len+1, $h1len));
    $my_Rz_5p =~ tr/AUCG/UAGC/;

    # construct the Rz sequence
    my($my_Rz_seq) = join("", $my_Rz_5p, $Rz_core, $my_Rz_3p);

    push @Rz_seq, $my_Rz_seq;
    push @Rz_3p, $my_Rz_3p;
    push @Rz_5p, $my_Rz_5p;

    push @dG_switch, &compute_dG_switch($i+1, $my_Rz_seq, $my_Rz_5p, $my_Rz_3p);

    $list_pos .= ($i+1) . ",";
    $list_len1 .= $h3len . ",";
    $list_len2 .= $h1len . ",";
  }
}

$list_pos =~ s/,$//;
$list_len1 =~ s/,$//;
$list_len2 =~ s/,$//;

open(DG_HYBRID, "$loc_dG_hybrid -m $opt{m} -s $list_pos -a $list_len1 -b $list_len2 |")
  || die " Error: unable to run $loc_dG_hybrid";
while (<DG_HYBRID>) {
  my(@e) = split;
  push @dG_hybrid, $e[1];
}
close(DG_HYBRID);


# temporary file name
my($randfn) = &randname();

# write target site starting positions to temp file for disruptEn
$list_pos =~ s/,/\n/g;
open(TMPFILE, ">$randfn") || die " Error: unable to create temp file";
# 2011-02-04 Adam. Design relic.
for my $one_start_pos (split(/\n/, $list_pos)) {
  print TMPFILE "$one_start_pos\t" . ($one_start_pos + $sitelen - 1) . "\n";
}
close(TMPFILE);

open(DG_DISRUPT, 
  "SFOLDBIN=$sfold_bin/ $loc_dG_disrupt -m $opt{m} -s $bpfile -p $randfn -f $fefile |")
  || die " Error: unable to run $loc_dG_disrupt";
while (<DG_DISRUPT>) {
  my(@e) = split;
  push @dG_disrupt, $e[1];
}
close(DG_DISRUPT);

#unlink $randfn;

# now sum energies to get dG_total and print output
    open(OUTFILE, ">", $outfile) or die("could not open $outfile for writing\n");
    
print OUTFILE <<EndOfOutput;
~~~~~~~~~~~~~~~Output for hammerhead ribozyme cleavage site energies~~~~~~~~~~~~~~~~~

Column 1: site ID (only those sites with long enough arms will be shown)
Column 2: target starting position
Column 3: target ending position
Column 4: dG_total = dG_hybrid - dG_switch - dG_disrupt (kcal/mol)
Column 5: dG_hybrid (kcal/mol)
Column 6: dG_switch (kcal/mol)
Column 7: dG_disrupt (kcal/mol)
Column 8: target site sequence (5p --> 3p)
Column 9: ribozyme sequence (5p --> 3p)

--------------------------------------------------------------------------------

EndOfOutput

for (my $i=0; $i<=$#startpos; $i++) {
  printf(OUTFILE "%3d %5d %5d %7.2f %7.2f %7.2f %7.2f %s %s\n", $i+1+$num_dropped_sites_5p, 
    $startpos[$i], $startpos[$i]+$sitelen-1,
    $dG_hybrid[$i]-$dG_switch[$i]-$dG_disrupt[$i], $dG_hybrid[$i],
    $dG_switch[$i], $dG_disrupt[$i], $sites[$i], $Rz_seq[$i]);
}

close(OUTFILE);
exit;


# fold the Rz seq and determine the dG_switch
sub compute_dG_switch {
  my($pos, $my_Rz_seq, $my_Rz_5p, $my_Rz_3p) = @_;
  my($Rz_ec_fname) = join("", "Rz_ecentroid_", $pos, ".ct");

  # temporary file name
  my($tmpdirname) = &randname();

  my($cur_dir) = &cwd();
  $cur_dir .= "/" if ($cur_dir !~ /\/$/);

  # prepare directory structures for running Sfold;
  # you don't have to do all these if you are using the latest version
  # of Sfold, this is needed here only because the version of Sfold on
  # the web server is not of the latest...
  mkdir "$tmpdirname" || die " Error: unable to create temp directory";
  chdir "$tmpdirname" || die " Error: unable to change to directory '$tmpdirname'";
  mkdir "output" || die " Error: unable to create Sfold temp output directory";
  `ln -sf $sfold_param param`;
  chdir "output" || die " Error: unable to change to Sfold temp output dir";

  # create fasta file for Rz sequence
  `echo "> Rz sequence" > $tmpdirname.fa`;
  `echo "$my_Rz_seq" >> $tmpdirname.fa`;

  # execute Sfold on the Rz sequence
  `$sfold_bin/sfold $tmpdirname.fa &> sfold.out`;
  die " Error: Sfold execution failed in compute_dG_switch()" if ($?);

  # Change back to original directory
  chdir "$cur_dir" || die " Error: unable to change to directory '$cur_dir'";

  # run rz_dgswitch.pl to determine dG_switch
  my($dG) = `$loc_dG_switch $tmpdirname/output/fe.out $my_Rz_5p $my_Rz_3p | awk '{print \$3}'`;
  chomp $dG;

  # delete entire temp directory
  `rm -rf $tmpdirname`;

  return $dG;
}

# print help messages
sub usage
{
    print STDERR "\n";
    print STDERR "Usage: rz_dghybrid.pl -m target.fasta -a h3len -b h1len -t triplet\n";
    print STDERR "         -s bp.out -f fe.out -o <output file>\n";
    print STDERR "  where h3len is the length of Helix III at target sites\n";
    print STDERR "        h1len is the length of Helix I at target sites\n";
    print STDERR "        triplet is one of the NUH triplets\n";
    print STDERR "        bp.out is the Sfold file containing sampled structures\n";
    print STDERR "        fe.out is the Sfold file containing structure free energies\n";
    print STDERR "        <output file>, is the name of the file the energies should be written to.";
    print STDERR "\n";
    exit;
}

# generate a random temp file name
sub randname {
  my($lim) = 10000;
  my($fn) = "tmp";
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);

  return $fn;
}
