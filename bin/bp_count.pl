#!/usr/bin/perl

#
# Calculate base pair frequencies from sampled ensemble
#

$infile = $ARGV[0];
if ($infile eq "") {
  print "Usage: $0 Sfold_bp_file\n";
  exit;
}

@a = ();
$max_x = 0;
$max_y = 0;
$nsamples = `grep -i Structure $infile | wc -l`;
$nsamples =~ s/^\s*(\d+)\D*$/$1/g;

open(INFILE, "$infile") || die "Unable to open file '$infile'!";
while (<INFILE>) {
  chomp;

  if (/^\s*\d+\s+\d+\s*$/) {
    my($x, $y) = split;

    if ($a[$x][$y] =~ /^\d+$/) {
      $a[$x][$y] += 1;
    } else {
      $a[$x][$y] = 1;
    }

    $max_x = $x if ($x > $max_x);
    $max_y = $y if ($y > $max_y);
  }
}
close(INFILE);

#$total = 0;
for ($i=1; $i<=$max_x; $i++) {
  for ($j=$i+1; $j<=$max_y; $j++) {
    print $i, " ", $j, " ", $a[$i][$j]," $nsamples\n" if ($a[$i][$j] =~ /^\d+$/);
#    $total += $a[$i][$j];
  }
}

#print "Total = $total\n";
