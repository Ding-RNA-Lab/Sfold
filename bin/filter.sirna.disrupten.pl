#!/usr/bin/perl -w

#
# Globals
#
#use vars qw/ %opt /;
use warnings;
#use strict;

my $DEBUG=0;
#
# Command line options processing
#
use Getopt::Std;
my $opt_string = 'cdghi:l:ps:t:o:';
getopts( 'cdghi:l:ps:t:o:', \my %opt );

#print (STDERR $opt{l});
#print (STDERR $opt{o});
#print (STDERR $opt{s});

my $apfile = "";
my $sifile = "";
my $sitelen = 19;
my $threshold = 0.000;
my $gcp = 0;
my $nx4 = 0;
my $dsseopt = 0;
my $combinesites = 0;
my $output_file = "";

#&usage() if ($#ARGV == -1);


if (!$opt{o}) {
  die "Error: you must specify an output file!";
} 
else {
 print $opt{o};
   $output_file = $opt{o};
}
open(OUTPUTFILE, '>', $output_file) or die("could not open $output_file for output.\n");

if (!$opt{i}) {
  die "Error: you must specify Sfold output file sirna.out!";
} else {
  $sifile = $opt{i};
}
die "Error: unable to locate Sfold output file sirna.out!" if (!-e $sifile);

if (!$opt{s}) {
  die "Error: you must specify an input file containing accessibility profile!";
} else {
  $apfile = $opt{s};
}
die "Error: unable to locate input file with accessibility profile!" if (!-e $apfile);

if ($opt{l}) {
  die "Error: invalid target site length specified!"
    if ($opt{l} !~ /^\d+$/ || $opt{l} <= 0);

  $sitelen = $opt{l};
}

if ($opt{t}) {
  die "Error: invalid threshold probability specified!"
    if ($opt{t} !~ /^[\d\.]+$/ || $opt{t} < 0 || $opt{t} > 1);

  $threshold = $opt{t};
}

if ($opt{p}) {
  $gcp = 1;
}

if ($opt{g}) {
  $nx4 = 1;
}

if ($opt{c}) {
  $combinesites = 1;
}

if ($opt{d}) {
  $dsseopt = 1;
}

# parse sirna.out to read in the DSSE values
open(SIRNA, "$sifile") || die "Error: unable to open file '$sifile'";
while (<SIRNA>) {
  # skip the headers
  last if (/---------------------/);
}

my @dsse = ();
# the first window in sirna.out starts at position 3.
# so here we set windows at positions 1 and 2 to some
# small values so that they won't get through the
# filters
$dsse[1] = $dsse[2] = -999;
my $idx = -1;
while (<SIRNA>) {
  chomp;
  next if (/^\s*$/);

  if (/\s*(\d+)\-/) {
    # line 1 of record
    $idx = $1;
  } else {
    my(@e) = split;
    if ($#e == 9) {
      # line 2 of record
      if ($idx == -1) {
        die " Error: sirna.out not in the right format";
      } else {
        $dsse[$idx] = $e[6];
        $idx = -1;
      }
    }
  }
}
close(SIRNA);


open(SSTRAND, "$apfile") || die "Error: unable to open file '$apfile'";
# file format for sstrand.out
# we need the 4th column for accessibility profile for w=1
#    1 U A  0.992 0.988 0.142
@uprob = ();
$seqlen = 0;
$seq = "";
while (<SSTRAND>) {
  chomp;

  die "Error: file '$apfile' not in the right format"
    if (!/^\s*(\d+)\s+(\w)\s+\w\s+([\d\.]+)\s+/);
#    if (!/^\s*(\d+)\s+(\w)\s+\w\s+([\d\.]+)\s+[\d\.]+\s+[\d\.]+$/);

  $seqlen++;
  die "Error: missing line in file '$apfile'" if ($seqlen != $1);

  $seq .= $2;
  $uprob[$1] = $3;
}
close(SSTRAND);

my $nwin = $seqlen - $sitelen + 1;
$seq =~ tr/atcgun/ATCGUN/;
$seq =~ tr/T/U/;

die " Error: lengths of sstrand.out and sirna.out do not match!"
  if ($nwin != $#dsse);

if ($DEBUG) {
  print STDERR "\n";
  print STDERR " Sequence length = $seqlen\n";
  print STDERR " Number of $sitelen nt-long windows = $nwin\n";
}

my @spos = ();
my %sites = ();
my %avgprob = ();
my %percent_gc = ();
my %site_dsse = ();
for (my $i=1; $i<=$nwin; $i++) {
  my($thiswin) = substr($seq, $i-1, $sitelen);

  if ($nx4) {
    # no GGGG,AAAA,UUUU,CCCC in the target sequence
    next if ($thiswin =~ /GGGG/ || $thiswin =~ /AAAA/ ||
             $thiswin =~ /UUUU/ || $thiswin =~ /CCCC/);
  }

  my($onlyGC) = $thiswin;
  $onlyGC =~ s/[^GC]//g;
  my($percent) = length($onlyGC)/$sitelen*100;
  if ($gcp) {
    next if ($percent < 30 || $percent > 70);
  }

  if ($dsseopt) {
    # include only sites with DSSE > -1 kcal/mol
    next if ($dsse[$i] <= -1);
  }

  my($prob) = 0;
  for (my $j=0; $j<$sitelen; $j++) {
    $prob += $uprob[$i+$j];
  }
  $prob /= $sitelen;
  next if ($prob < $threshold);

  push @spos, $i;
  $sites{$i} = $thiswin;
  $percent_gc{$i} = $percent;
  $avgprob{$i} = $prob;
  $site_dsse{$i} = $dsse[$i];
}


if ($combinesites) {
  # we try to eliminate "duplicate" sites here...

  # locate all peaks above the threshold probability
  my @pstart = ();
  my @pend = ();
  my($prev) = &max(0, $uprob[1]-$threshold);
  if ($prev > 0) {
    push @pstart, 1;
  }
  for (my $i=2; $i<=$seqlen; $i++) {
    my($current) = &max(0, $uprob[$i]-$threshold);

    if ($current > 0) {
      if ($i != $seqlen) {
        if ($prev <= 0) {
          push @pstart, $i;
        }
      } else {
        if ($prev > 0) {
          push @pend, $i;
        }
      }
    } else {
      if ($prev > 0) {
        push @pend, $i-1;
        if ($pend[$#pend]-$pstart[$#pstart] <= 0) {
          pop @pstart;
          pop @pend;
        }
      }
    }

    $prev = $current;
  }

  die "Error: inconsistent numbers of starting and ending positions of peaks"
    if ($#pstart != $#pend);

  # now we are going to sort the 2 arrays at the same time 
  # according to the difference between the values in the
  # arrays. we are sorting first by the difference (large 
  # to small), then by pstart (small to large).
  #
  my @sorted_ind = sort { (($pend[$b]-$pstart[$b]) <=> ($pend[$a]-$pstart[$a]))
                    || ($pstart[$a] <=> $pstart[$b]) } 0..$#pstart;
  @pstart = @pstart[@sorted_ind];
  @pend = @pend[@sorted_ind];

  #
  # at this point all peaks are stored in (pstart[i], pend[i])
  # and they are in descending order of peak width
  #

  %pincluded = ();
  foreach  $i (@pstart) {
    $pincluded{$i} = 0;
  }
  %sincluded = ();
  foreach  $i (@spos) {
    $sincluded{$i} = 0;
  }

  for (my $i=0; $i<=$#pstart; $i++) {
    next if ($pincluded{$pstart[$i]} == 1);

    my($maxwin) = sprintf("%d", ($pend[$i]-$pstart[$i]+1)/$sitelen);
#    $maxwin = 1 if ($maxwin == 0);
#
# talked to Ye about this. Ye suggested to pick at least 2 sites
# for each peak. code changed to reflect this update.
#
    $maxwin = &max($maxwin, 2);
    my(@allsites) = ();
    my(@allprob) = ();
    for (my $j=0; $j<=$#spos; $j++) {
  #
  # we are not going to check for the already-marked sites here,
  # which means a site might be checked more than once by different
  # peaks. that's ok, cuz a good site that works for more than
  # one peak is good!
  #
  #    next if ($sincluded{$spos[$j]} == 1);
      my($sposend) = $spos[$j] + $sitelen - 1;

      if (($spos[$j] >= $pstart[$i] && $spos[$j] < $pend[$i]) || 
          ($sposend > $pstart[$i] && $sposend <= $pend[$i]) ||
          ($spos[$j] <= $pstart[$i] && $sposend >= $pend[$i])) {
        # there is overlap between this peak and this site
        push @allsites, $spos[$j];
        push @allprob, $avgprob{$spos[$j]};
      }
    }

    if ($DEBUG) {
      print STDERR $pstart[$i], " ", $pend[$i], "----", "@allsites", "\n" if ($#allsites != -1);
    }

    # sort possible sites by probability
    my(@sorted_ind) = sort { $allprob[$b] <=> $allprob[$a] } 0..$#allprob;
    @allsites = @allsites[@sorted_ind];
    @allprob = @allprob[@sorted_ind];

    for (my $j=0; $j<$maxwin && $j<=$#allsites; $j++) {
      $sincluded{$allsites[$j]} = 1;
    }

    $pincluded{$pstart[$i]} = 1;
  }

#  for (my $i=0; $i<=$#pstart; $i++) {
#    print $pstart[$i], "    ", $pend[$i], "    ", $pend[$i]-$pstart[$i], "\n";
#  }

}

print OUTPUTFILE <<EndOfHeader;
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
    probability, all sites targeted to this same peak are
    ranked by their average unpaired probability (the higher the better) and
    at most n sites are selected for each peak, where n is determined by
    max([width of peak/site length], 2);
 F) Among sites satisfying criteria A-E, the top unique ones with
    the highest average unpaired probability are listed.

NOTE:
 i) The average unpaired probability is used in filter criteria D, E and F to
    cut down the number of reported sites in order to make the disruption
    energy calculation manageable on our web servers.
--------------------------------------------------------------------------------

EndOfHeader

foreach $i (@spos) {
  if (($combinesites && $sincluded{$i}==1) || $combinesites==0) {
    printf(OUTPUTFILE "%5d  %s  %.3f  %5.1f%%  %5.1f\n", $i, $sites{$i}, $avgprob{$i},
      $percent_gc{$i}, $site_dsse{$i});
  }
}

exit;


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


sub round() {
  my($num) = @_;

  return int($num + 0.5 * ($num <=> 0));
}


#
# Message about this program and how to use it
#
sub usage()
{
  my($gcp_def, $nx4_def, $dsseopt_def, $combinesites_def);
  if ($gcp) {
    $gcp_def = "on";
  } else {
    $gcp_def = "off";
  }
  if ($nx4) {
    $nx4_def = "on";
  } else {
    $nx4_def = "off";
  }
  if ($dsseopt) {
    $dsseopt_def = "on";
  } else {
    $dsseopt_def = "off";
  }
  if ($combinesites) {
    $combinesites_def = "on";
  } else {
    $combinesites_def = "off";
  }

  print STDERR << "EOF";

Usage: $0 [-cdgp] [-h] -i file [-l length] -s file [-t threshold]

 -c           : combine overlapping sites which share the same peak in the
                probability profile and report only the best one, i.e. the
                one with the highest unpaired probability [default = $combinesites_def]
 -d           : exclude sites that do not satisfy DSSE > -1 kcal/mol
                [default = $dsseopt_def]
 -g           : exclude sites with at least one of AAAA, UUUU, CCCC or GGGG
                [default = $nx4_def]
 -h           : print this help message
 -i file      : Sfold output file sirna.out
 -l length    : target site length [default = $sitelen]
 -p           : filter sites by GC content (30% <= GC% <= 70%) [default = $gcp_def]
 -s file      : file containing accessibility profile (sstrand.out from Sfold)
 -t threshold : all target sites with an average unpaired probability equal to 
                OR greater than this threshold will be printed as output 
                [default = $threshold]
 -o output-fille
              : file to which output will be sent

Example: $0 -c -g -i ./sirna.out -l 19 -p -s ./sstrand.out -t 0.5 -o myoutfile.txt

EOF

  exit;
}
