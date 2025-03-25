#!/usr/bin/perl -w

#
# Globals
#
use vars qw/ %opt /;

require File::Spec;
use Cwd qw/ cwd abs_path /;
use warnings;

# program paths
$curdir = &dir_full_path(".");
$bindir = &dir_full_path($0);
$sfoldbin = &dir_full_path($0);
$prefilter = join("", $bindir, "/filter.soligo.disrupten.pl");
$disruptEn = join("", $bindir, "/disruptEn");

# user options
$inituprob = 0.5;
$minuprob = 0.5;
$minsites = 10;
$maxsites = 20;

# default command-line options
$newdir = "";
$bpfile = "./bp.out";
$fefile = "./fe.out";
$seqfile = "";
$outfile = "STDOUT";
$apfile = "./sstrand.out";
$sitelen = 20;

&usage() if ($#ARGV == -1);
&init();

if ($opt{d}) {
  $newdir = &dir_full_path($opt{d});
  die "Error: unable to locate directory '$newdir'!" if (!-e $newdir);

  chdir "$newdir" || die "Error: unable to change to directory '$newdir'!";
  print STDERR "Changed to directory $newdir\n";
}

$bpfile = $opt{b} if ($opt{b});
$fefile = $opt{f} if ($opt{f});

if (!$opt{i}) {
  die "Error: you must specify an input sequence file!";
} else {
  $seqfile = $opt{i};
}

$outfile = $opt{o} if ($opt{o});
$apfile = $opt{s} if ($opt{s});

if ($opt{l}) {
  die "Error: invalid target site length specified!"
    if ($opt{l} !~ /^\d+$/ || $opt{l} <= 0);

  $sitelen = $opt{l};
}


die "Error: unable to locate file '$bpfile'!" if (!-e $bpfile);
die "Error: unable to locate file '$fefile'!" if (!-e $fefile);
die "Error: unable to locate file '$seqfile'!" if (!-e $seqfile);
if ($outfile ne "STDOUT") {
  open($OUTFH, ">$outfile") || die "Error: unable to create file '$outfile'!";
} else {
  $OUTFH = STDOUT;
}
die "Error: unable to locate file '$apfile'!" if (!-e $apfile);

#
# OK, we are ready to start the calculation
#

$randfn = randname();
$uprob = $inituprob;
@denergies = ();

my($num) = `$prefilter -c -g -l $sitelen -p -s $apfile -t $uprob | wc -l`;
chomp($num);

while ($num < $minsites) {
  # threshold too high perhaps, let's lower it down a bit at a time
  last if ($uprob-0.1 < $minuprob);
  $uprob -= 0.1;

  $num = `$prefilter -c -g -l $sitelen -p -s $apfile -t $uprob | wc -l`;
  chomp($num);
}

`$prefilter -c -g -l $sitelen -p -s $apfile -t $uprob | sort -n -k 3 -r | head -$maxsites | sort -n -k 1 > $randfn`;
# 2011-02-04 Adam. The version of disruptEn that accepts a start pos and site len is obsolete. I replaced
#  it with a version that takes start and end pos, in a file.
`gawk 'BEGIN {OFS="\t"} {print \$1, \$1 + $sitelen - 1}' $randfn > $randfn.bsites`;
open (RUNDE, 
  "SFOLDBIN=$sfoldbin/ $disruptEn -m $seqfile -s $bpfile -p $randfn.bsites -f $fefile |")
  || die "Error: unable to run $disruptEn";
while (<RUNDE>) {
  my(@e) = split;
  push @denergies, $e[1];
}
close(RUNDE);

print $OUTFH <<EndOfHeader;
~~~~~~~~~~~~~~~~Filtered output for design of antisense oligos~~~~~~~~~~~~~~~~~

Column 1: starting target position
Column 2: ending target position
Column 3: target sequence (5p --> 3p)
Column 4: antisense oligo (5p --> 3p)
Column 5: GC content
Column 6: average unpaired probability for target site nucleotides
Column 7: binding site disruption energy (kcal/mol)

FILTER CRITERIA: ("<=": less than or equal to)
                 (">=": greater than or equal to)

 A) 40% <= GC % <= 60%;
 B) No GGGG in the target sequence;
 C) Average unpaired probability for target site nucleotides >= $uprob;
 D) For each peak in the accessibility profile that is above the threshold
    probability of $uprob, all sites targeted to this same peak are
    ranked by their average unpaired probability (the higher the better) and
    at most n sites are selected for each peak, where n is determined by
    max([width of peak/site length], 2);
 E) Among sites satisfying criteria A-D, the top $maxsites unique ones with
    the highest average unpaired probability are listed.

NOTE:
 i) The average unpaired probability is used in filter criteria C, D and E to
    cut down the number of reported sites in order to make the disruption
    energy calculation manageable on our web servers.
--------------------------------------------------------------------------------

EndOfHeader

#NOTE:
# i) In order to make the disruption energy calculation manageable on our web
#    servers, the average unpaired probability is used in filter criteria C,
#    D and E to cut down the number of reported sites. However, the average
#    unpaired probability of a site has not been proven to be correlated to
#    its disruption energy.


$nsites = 0;
open(SFILE, "$randfn") || die "Error: unable to open file '$randfn'";
while (<SFILE>) {
  my($spos, $sense, $prob, $gcp) = split;
  my($antisense) = scalar reverse "$sense";
  $antisense =~ tr/ATCGUN/TAGCAN/;

  printf $OUTFH "%4d %4d %s %s %6s %6.3f %7.1f\n", $spos, $spos+$sitelen-1,
    $sense, $antisense, $gcp, $prob, $denergies[$nsites];

  $nsites++;
}
close(SFILE);

#unlink $randfn;
#unlink "$randfn.bsites";

#
# End calculation.
#

close($OUTFH);

# change back to original directory
if ($opt{d}) {
  chdir "$curdir" || die "Error: unable to change back to directory '$curdir'!"
}

exit;


#
# Return the absolute path of the input directory (or the absolute
# path of the directory containing the input file)
#
sub dir_full_path()
{
  my($in) = @_;

  if (!-d $in) {
    my($volume, $dir, $file) = File::Spec->splitpath($in);
    if ($dir eq "") {
      $in = ".";
    } else {
      $in = $dir;
    }
  }

  return abs_path($in);
}


#
# Command line options processing
#
sub init()
{
    use Getopt::Std;
    my $opt_string = 'b:d:f:hi:l:o:s:';
    getopts( "$opt_string", \%opt ) or usage();

    usage() if ($opt{h});
}


# return the larger of the two input values
#
sub max() {
  my($x, $y) = @_;

  if ($x > $y) {
    return $x;
  } else {
    return $y;
  }
}


sub randname {
  my($lim) = 10000;
  my($fn) = "tmp";
  $fn .= sprintf "%d", (rand($lim) + 1);
  $fn .= sprintf "%d", (rand($lim) + 1);
  $fn .= sprintf "%d", (rand($lim) + 1);

  return $fn;
}


sub round() {
  my($num) = @_;

  return int($num + 0.5 * ($num <=> 0));
}


#
# Message about this program and how to use it
#
sub usage()
{
  print STDERR << "EOF";

Usage: $0 [-b file] [-d dir] [-f file] [-h] -i file [-l length]
       [-o file] [-s file]

 -b file      : file containing sampled structures (bp.out from Sfold)
                [default = $bpfile]
 -d dir       : change to this directory before reading from input files and
                writing to output files. If specified, all relative paths of
                input and output files will be resolved with respect to this
                new directory
 -f file      : file containing free energies of all sampled structures (fe.out
                from Sfold) [default = $fefile]
 -h           : print this help message
 -i file      : sequence file in FASTA format
 -l length    : target site length [default = $sitelen]
 -o file      : name of file for output [default = standard output]
 -s file      : file containing accessibility profile (sstrand.out from Sfold)
                [default = $apfile]

Example: $0 -i ./seq.fa -l 20 -o ./disrupten.out

EOF

  exit;
}
