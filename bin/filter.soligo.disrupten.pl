#!/usr/bin/perl -w

#
# Globals
#
use vars qw/ %opt /;

$DEBUG=0;

$apfile = "";
$sitelen = 20;
$threshold = 0.000;
$gcp = 0;
$gx4 = 0;
$combinesites = 0;

&usage() if ($#ARGV == -1);
&init();

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
  $gx4 = 1;
}

if ($opt{c}) {
  $combinesites = 1;
}

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

$nwin = $seqlen - $sitelen + 1;
$seq =~ tr/atcgun/ATCGUN/;

if ($DEBUG) {
  print STDERR "\n";
  print STDERR " Sequence length = $seqlen\n";
  print STDERR " Number of $sitelen nt-long windows = $nwin\n";
}

@spos = ();
%sites = ();
%avgprob = ();
%percent_gc = ();
for (my $i=1; $i<=$nwin; $i++) {
  my($thiswin) = substr($seq, $i-1, $sitelen);

  if ($gx4) {
    # no GGGG in the target sequence
    next if ($thiswin =~ /GGGG/);
  }

  my($percent);
  if ($gcp) {
    my($onlyGC) = $thiswin;
    $onlyGC =~ s/[^GC]//g;

    $percent = length($onlyGC)/$sitelen*100;
    next if ($percent < 40 || $percent > 60);
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
}


if ($combinesites) {
  # we try to eliminate "duplicate" sites here...

  # locate all peaks above the threshold probability
  @pstart = ();
  @pend = ();
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
  @sorted_ind = sort { (($pend[$b]-$pstart[$b]) <=> ($pend[$a]-$pstart[$a]))
                    || ($pstart[$a] <=> $pstart[$b]) } 0..$#pstart;
  @pstart = @pstart[@sorted_ind];
  @pend = @pend[@sorted_ind];

  #
  # at this point all peaks are stored in (pstart[i], pend[i])
  # and they are in descending order of peak width
  #

  %pincluded = ();
  foreach my $i (@pstart) {
    $pincluded{$i} = 0;
  }
  %sincluded = ();
  foreach $i (@spos) {
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

foreach $i (@spos) {
  if (($combinesites && $sincluded{$i}==1) || $combinesites==0) {
    printf("%5d  %s  %.3f  %5.1f%%\n", $i, $sites{$i}, $avgprob{$i}, $percent_gc{$i});
  }
}

exit;


#
# Command line options processing
#
sub init()
{
    use Getopt::Std;
    my $opt_string = 'cghl:ps:t:';
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


sub round() {
  my($num) = @_;

  return int($num + 0.5 * ($num <=> 0));
}


#
# Message about this program and how to use it
#
sub usage()
{
  my($gcp_def, $gx4_def);
  if ($gcp) {
    $gcp_def = "on";
  } else {
    $gcp_def = "off";
  }
  if ($gx4) {
    $gx4_def = "on";
  } else {
    $gx4_def = "off";
  }
  if ($combinesites) {
    $combinesites_def = "on";
  } else {
    $combinesites_def = "off";
  }

  print STDERR << "EOF";

Usage: $0 [-h] [-l length] [-p] -s file [-t threshold]

 -c           : combine overlapping sites which share the same peak in the
                probability profile and report only the best one, i.e. the
                one with the highest unpaired probability [default = $combinesites_def]
 -g           : exclude sites with GGGG in target sequence [default = $gx4_def]
 -h           : print this help message
 -l length    : target site length [default = $sitelen]
 -p           : filter sites by GC content (40% <= GC% <= 60%) [default = $gcp_def]
 -s file      : file containing accessibility profile (sstrand.out from Sfold)
 -t threshold : all target sites with an average unpaired probability equal to 
                OR greater than this threshold will be printed as output 
                [default = $threshold]

Example: $0 -c -g -l 19 -p -s ./sstrand.out -t 0.5

EOF

  exit;
}
