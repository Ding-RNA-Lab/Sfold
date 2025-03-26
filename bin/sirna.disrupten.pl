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
#$sfoldbin = &dir_full_path("/nfs/sfold/bin");
$sfoldbin = &dir_full_path($0);
$prefilter = join("", $bindir, "/filter.sirna.disrupten.pl");
$disruptEn = join("", $bindir, "/disruptEn");

# user options
$inituprob = 0.5;
$minuprob = 0.5;
$minsites = 0;
$maxsites = 20;

# default command-line options
$newdir = "";
$bpfile = "./bp.out";
$fefile = "./fe.out";
$seqfile = "";
$outfile = "STDOUT";
$apfile = "./sstrand.out";
$sifile = "./sirna.out";
$sitelen = 19;

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
$sifile = $opt{r} if ($opt{r});

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
die "Error: unable to locate file '$sifile'!" if (!-e $sifile);

#
# OK, we are ready to start the calculation
#

$randfn = randname();
$uprob = $inituprob;
@denergies = ();

my($num) = `$prefilter -c -d -g -i $sifile -l $sitelen -p -s $apfile -t $uprob | wc -l`;
chomp($num);

while ($num < $minsites) {
  # threshold too high perhaps, let's lower it down a bit at a time
  last if ($uprob-0.1 < $minuprob);
  $uprob -= 0.1;

  $num = `$prefilter -c -d -g -i $sifile -l $sitelen -p -s $apfile -t $uprob | wc -l`;
  chomp($num);
}

# 2011-02-04 Adam. The version of disruptEn that accepts a start pos and site len is obsolete. I replaced
#  it with a version that takes start and end pos, in a file.
`$prefilter -c -d -g -i $sifile -l $sitelen -p -s $apfile -t $uprob | sort -n -k 3 -r | head -$maxsites | sort -n -k 1 > $randfn`;
`gawk 'BEGIN {OFS="\t"} {print \$1, \$1 + $sitelen - 1}' $randfn > $randfn.bsites`;

print "SFOLDBIN=$sfoldbin/ $disruptEn -m $seqfile -s $bpfile -p $randfn.bsites -f $fefile\n";
#print "SFOLDBIN=/home/williamrennie/Desktop/SfoldRun1/Sfold-main/bin/ $disruptEn -m $seqfile -s $bpfile -p $randfn.bsites -f $fefile\n";
open (RUNDE, 
  "SFOLDBIN=$sfoldbin/ $disruptEn -m $seqfile -s $bpfile -p $randfn.bsites -f $fefile |")
  || die "Error: unable to run $disruptEn";
  
while (<RUNDE>) {	
  my(@e) = split;
  push @denergies, $e[1];
}
close(RUNDE);

print $OUTFH <<EndOfHeader;
~~~~~~~~~~~~~~~~Filtered output for siRNAs with disruption energy~~~~~~~~~~~~~~~~~

Column 1: starting target position
Column 2: ending target position
Column 3: sense siRNA (5p --> 3p)
Column 4: antisense siRNA (5p --> 3p)
Column 5: siRNA GC content
Column 6: differential stability of siRNA duplex ends (DSSE, in kcal/mol)
Column 7: average unpaired probability for target site nucleotides
Column 8: binding site disruption energy (kcal/mol)

FILTER CRITERIA: ("<=": less than or equal to)
                 (">=": greater than or equal to)
                 ( ">": greater than)

 A) 30% <= GC % <= 70%;
 B) Exclusion of target sequence with at least one of AAAA,
    CCCC, GGGG, or UUUU;
 C) DSSE > -1 kcal/mol (asymmetry rule);
 D) Average unpaired probability for target site nucleotides >= $uprob;
 E) For each peak in the accessibility profile that is above the threshold
    probability of $uprob, all sites targeted to this same peak are
    ranked by their average unpaired probability (the higher the better) and
    at most n sites are selected for each peak, where n is determined by
    max([width of peak/site length], 2);
 F) Among sites satisfying criteria A-E, the top $maxsites unique ones with
    the highest average unpaired probability are listed.

NOTE:
 i) The average unpaired probability is used in filter criteria D, E and F to
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
  my($spos, $sense, $prob, $gcp, $dsse) = split;
  my($antisense) = scalar reverse "$sense";
  $antisense =~ tr/ATCGUN/UAGCAN/;

  printf $OUTFH "%4d %4d %sTT %sTT %6s %5.1f %6.3f %7.1f\n", $spos, $spos+$sitelen-1,
    $sense, $antisense, $gcp, $dsse, $prob, $denergies[$nsites];

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
    my $opt_string = 'b:d:f:hi:l:o:r:s:';
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
       [-o file] [-r file] [-s file]

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
 -r file      : Sfold output file sirna.out [default = $sifile]
 -s file      : file containing accessibility profile (sstrand.out from Sfold)
                [default = $apfile]

Example: $0 -i ./seq.fa -l 19 -o ./disrupten.out

EOF

  exit;
}
