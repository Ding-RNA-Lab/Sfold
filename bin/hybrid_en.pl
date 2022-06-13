#!/usr/bin/perl -w

#
# This script determines the energies of a hybrid.
# it takes an alignment of the two strands from standard
# input and returns the energy. An example of the input
# alignment format is:
#
# miRNA strand
#   GACACUCCACCAUGAAU-CACUCCC
#    |||   ||||||     |||||| 
#   -UGUUU-GUGGUAACAGUGUGAGGU
# Target strand
#
# Update (20070710):
#  Corrected the way dangling end energies are included
#  in the final energy. To do so, we need to know which
#  strand is the miRNA and which is the target. Here we
#  will assume the top (positive) strand is the miRNA
#  and the bottom (negative) strand is the target. Since
#  there is uncertainty regarding the structure on the
#  dangling end of the target sequence, we will not add
#  any dangling term to single-stranded bases on the
#  target.
#
# Update (20070711):
#  Since there are several ways to handle the dangling
#  ends - no dangling at all; dangling just on the 
#  miRNAs; sum of dangling ends for both miRNAs and
#  target - the script has been updated to output three
#  energies at a time, one for each case specified above.
#

# specify location of the findfe program
my($findfe) = "/home/chanc/Sfold/bin/findfe";
die " Error: unable to locate program '$findfe'"
  if (!-e $findfe);
my($dangling) = "/home/chanc/Sfold/bin/dangling.I386.SUNOS";
die " Error: unable to locate program '$dangling'"
  if (!-e $dangling);
my($sfold_param) = "/home/chanc/Sfold/param/";
die " Error: unable to locate folder '$sfold_param'"
  if (!-e $sfold_param);

# terminal A-U penalty
my($terminal_AU_penalty) = 0.5;

# set auto flush
$| = 1;

use strict;
use Cwd qw/ cwd /;

my($basedir) = &cwd();
my($loopseq) = "GGGAAAAAAAAAACCC";
my(@lspos) = ();  # save the position ID of the loop seq.

# temporary file name
my($tmpfname) = join("", "/tmp/", &randname());

while (<>) {
  chomp;

  # skip blank lines
  next if (/^\s*$/);

  my($posstrand, $negstrand, $pairing);
  $posstrand = $_;
  $pairing = <>;
  chomp($pairing);
  $negstrand = <>;
  chomp($negstrand);

  $posstrand =~ tr/a-z/A-Z/;
  $posstrand =~ s/T/U/g;
  $negstrand =~ tr/a-z/A-Z/;
  $negstrand =~ s/T/U/g;

  # do some error checkings on the input alignment
  die " Error: invalid base found on first strand"
    if ($posstrand =~ /[^ATCGU-]/);
  die " Error: invalid base found on second strand"
    if ($negstrand =~ /[^ATCGU-]/);
  die " Error: lengths of strands do not match"
    if (length($posstrand) != length($negstrand));
  die " Error: invalid symbol found on pairing line"
    if ($pairing =~ /[^\s|]/);

  my($lastbp) = rindex($pairing, "|");
  die " Error: pairing line is longer than strand lines"
    if ($lastbp+1 > length($posstrand));
  die " Error: no base pairing specified"
    if ($lastbp < 0);
  my($firstbp) = index($pairing, "|");

  # add enough spaces to the end of the pairing line to 
  # make it the same length as the two strands
  $pairing .= " "x(length($posstrand) - length($pairing));

  # chop unnecessary white spaces at the end of the pairing line
  # unnecessary = anything longer than the two strands
  while (length($pairing) > length($posstrand)) {
    chop($pairing);
  }

  # position ID for the positive and negative strands
  my(@pspos, @nspos);
  my($seqlen) = 0;

  # set position ID for bases on the positive strand
  for (my $i=0; $i<length($posstrand); $i++) {
    if ($i <= $lastbp && substr($posstrand, $i, 1) ne "-") {
      $seqlen++;
      $pspos[$i] = $seqlen;
    } else {
      $pspos[$i] = 0;
    }
  }

  # set position ID for bases on the artificial loop seq.
  for (my $i=0; $i<length($loopseq); $i++) {
    if (substr($loopseq, $i, 1) ne "-") {
      $seqlen++;
      $lspos[$i] = $seqlen;
    } else {
      $lspos[$i] = 0;
    }
  }

  # set position ID for bases on the negative strand
  for (my $i=length($negstrand)-1; $i>=0; $i--) {
    if ($i <= $lastbp && substr($negstrand, $i, 1) ne "-") {
      $seqlen++;
      $nspos[$i] = $seqlen;
    } else {
      $nspos[$i] = 0;
    }
  }

# Debug messages:
#  print "Sequence length = $seqlen\n";
#  print "$pspos[0],$pspos[$#pspos-1]\n";
#  print "$nspos[1],$nspos[$#nspos-1]\n\n";
#
#  print $posstrand, "\n";
#  print $pairing, "\n";
#  print $negstrand, "\n";

  open(OUTFILE, ">$tmpfname")
    || die " Error: unable to create temp file";
  print OUTFILE "SFOLD of:  [initially 0.00]  Temp Check: 0 from: 1 to: $seqlen\n";
  print OUTFILE "Length: $seqlen Energy: 0.00 ..\n";

  # print GCG connect file, first strand
  for (my $i=0; $i<length($posstrand); $i++) {
    if ($pspos[$i] != 0) {
      my($bpid) = 0;
      if (substr($pairing, $i, 1) eq "|") {
        $bpid = $nspos[$i];

        # check for valid base pairing, i.e. only W-C and wobble pairs
        if (substr($posstrand, $i, 1) eq "A") {
          die " Error: invalid base pairing at position $pspos[$i]"
            if (substr($negstrand, $i, 1) ne "U");

        } elsif (substr($posstrand, $i, 1) eq "U") {
          die " Error: invalid base pairing at position $pspos[$i]"
            if (substr($negstrand, $i, 1) ne "A" && substr($negstrand, $i, 1) ne "G");

        } elsif (substr($posstrand, $i, 1) eq "C") {
          die " Error: invalid base pairing at position $pspos[$i]"
            if (substr($negstrand, $i, 1) ne "G");

        } elsif (substr($posstrand, $i, 1) eq "G") {
          die " Error: invalid base pairing at position $pspos[$i]"
            if (substr($negstrand, $i, 1) ne "C" && substr($negstrand, $i, 1) ne "U");
        }
      }

      print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $pspos[$i],
        substr($posstrand, $i, 1),
        $pspos[$i]-1, ($pspos[$i]+1)%($seqlen+1), $bpid, $pspos[$i]);
    }
  }

  # print GCG connect file, artificial loop, assuming 3 BPs plus 10 single-stranded
  for (my $i=0; $i<3; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $lspos[$i],
      substr($loopseq, $i, 1),
      $lspos[$i]-1, ($lspos[$i]+1)%($seqlen+1), $lspos[length($loopseq)-1-$i], $lspos[$i]);
  }
  for (my $i=3; $i<13; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $lspos[$i],
      substr($loopseq, $i, 1),
      $lspos[$i]-1, ($lspos[$i]+1)%($seqlen+1), 0, $lspos[$i]);
  }
  for (my $i=13; $i<16; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $lspos[$i],
      substr($loopseq, $i, 1),
      $lspos[$i]-1, ($lspos[$i]+1)%($seqlen+1), $lspos[length($loopseq)-1-$i], $lspos[$i]);
  }

  # print GCG connect file, second strand
  for (my $i=length($negstrand)-1; $i>=0; $i--) {
    if ($nspos[$i] != 0) {
      my($bpid) = 0;
      $bpid = $pspos[$i] if (substr($pairing, $i, 1) eq "|");

      print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $nspos[$i],
        substr($negstrand, $i, 1),
        $nspos[$i]-1, ($nspos[$i]+1)%($seqlen+1), $bpid, $nspos[$i]);
    }
  }

  close(OUTFILE);

  my($complete_dG) = `$findfe $tmpfname -p $sfold_param | grep "Free energy" | awk '{print \$4}'`;
  chomp($complete_dG);
  `rm -f $tmpfname`;


#
# Now construct the GCG connect file for the artificial hairpin loop
#

  open(OUTFILE, ">$tmpfname")
    || die " Error: unable to create temp file";
  print OUTFILE "SFOLD of:  [initially 0.00]  Temp Check: 0 from: 1 to: 18\n";
  print OUTFILE "Length: 18  Energy: 0.00 ..\n";

  # position on the positive strand for the last base pair
  print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", 1,
    substr($posstrand, $lastbp, 1), 0, 2, 18, 1);

  # print GCG connect file, artificial loop, assuming 3 BPs plus 10 single-stranded
  for (my $i=0; $i<3; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $i+2,
      substr($loopseq, $i, 1), $i+1, $i+3, 17-$i, $i+2);
  }
  for (my $i=3; $i<13; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $i+2,
      substr($loopseq, $i, 1), $i+1, $i+3, 0, $i+2);
  }
  for (my $i=13; $i<16; $i++) {
    print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $i+2,
      substr($loopseq, $i, 1), $i+1, $i+3, 17-$i, $i+2);
  }

  # position on the negative strand for the last base pair
  print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", 18,
    substr($negstrand, $lastbp, 1), 17, 0, 1, 18);

  close(OUTFILE);

  my($loop_dG) = `$findfe $tmpfname -p $sfold_param | grep "Free energy" | awk '{print \$4}'`;
  chomp($loop_dG);
  `rm -f $tmpfname`;


  # now we are ready to find the energy of the hybrid
  my($diff_dG) = $complete_dG - $loop_dG;

  # add terminal A-U penalty if necessary
  if ( (substr($posstrand, $lastbp, 1) eq "A" && 
        substr($negstrand, $lastbp, 1) eq "U") ||
       (substr($posstrand, $lastbp, 1) eq "U" && 
        substr($negstrand, $lastbp, 1) eq "A") ) {
    $diff_dG += $terminal_AU_penalty;
  }

  my($pos_dangle5) = 0.0;
  my($pos_dangle3) = 0.0;
  my($neg_dangle5) = 0.0;
  my($neg_dangle3) = 0.0;

  # 5' dangling end of the positive strand
  #       5'  --> 3'
  #           2
  #           13
  #       3' <--  5'
  if ($firstbp > 0) {
    # go backwards on the positive strand to look for a single-stranded base
    my($i) = $firstbp - 1;
    while ($i >= 0) {
      last if (substr($posstrand, $i, 1) ne "-");
      $i--;
    }

    if ($i >= 0) {
      # 5' dangling end exists on the positive strand
      my($cmd) = join("", "$dangling -p $sfold_param -e 5 -a ", 
        substr($posstrand, $firstbp, 1), " -b ",
        substr($negstrand, $firstbp, 1), " -c ", substr($posstrand, $i, 1));
      $pos_dangle5 = `$cmd`;
      chomp($pos_dangle5);
    }
  }

  # 3' dangling end of the positive strand
  #       5'  --> 3'
  #           23
  #           1
  #       3' <--  5'
  if (length($posstrand) > $lastbp+1) {

    # go forward on the positive strand to look for a single-stranded base
    my($i) = $lastbp + 1;
    while ($i < length($posstrand)) {
      last if (substr($posstrand, $i, 1) ne "-");
      $i++;
    }

    if ($i < length($posstrand)) {
      # 3' dangling end exists
      my($cmd) = join("", "$dangling -p $sfold_param -e 3 -a ", substr($negstrand, $lastbp, 1), " -b ",
        substr($posstrand, $lastbp, 1), " -c ", substr($posstrand, $i, 1));
      $pos_dangle3 = `$cmd`;
      chomp($pos_dangle3);
    }
  }

  # 5' dangling end of the negative strand
  #       5'  --> 3'
  #           2
  #           13
  #       3' <--  5'
  if (length($negstrand) > $lastbp+1) {

    # go forward on the negative strand to look for a single-stranded base
    my($i) = $lastbp + 1;
    while ($i < length($negstrand)) {
      last if (substr($negstrand, $i, 1) ne "-");
      $i++;
    }

    if ($i < length($negstrand)) {
      # 5' dangling end exists
      my($cmd) = join("", "$dangling -p $sfold_param -e 5 -a ", substr($negstrand, $lastbp, 1), " -b ",
        substr($posstrand, $lastbp, 1), " -c ", substr($negstrand, $i, 1));
      $neg_dangle5 = `$cmd`;
      chomp($neg_dangle5);
    }
  }

  # 3' dangling end on the negative strand
  #       5'  --> 3'
  #           23
  #           1
  #       3' <--  5'
  if ($firstbp > 0) {
    # go backwards on the negative strand to look for a single-stranded base
    my($i) = $firstbp - 1;
    while ($i >= 0) {
      last if (substr($negstrand, $i, 1) ne "-");
      $i--;
    }

    if ($i >= 0) {
      # 3' dangling end exists on the negative strand
      my($cmd) = join("", "$dangling -p $sfold_param -e 3 -a ", 
        substr($posstrand, $firstbp, 1), " -b ",
        substr($negstrand, $firstbp, 1), " -c ", substr($negstrand, $i, 1));
      $neg_dangle3 = `$cmd`;
      chomp($neg_dangle3);
    }
  }

  # no dangling at all
  my($hybrid_dG_v0) = $diff_dG - $pos_dangle5 - $neg_dangle3;

  # 5' and 3' dangling ends on the miRNA
  my($hybrid_dG_v1) = $diff_dG + $pos_dangle3 - $neg_dangle3;

  # 5' and 3' dangling ends on both miRNA and target
  my($hybrid_dG_v2) = $diff_dG + $pos_dangle3 + $neg_dangle5;

  print join("\t", $hybrid_dG_v0, $hybrid_dG_v1, $hybrid_dG_v2), "\n";
}

exit;


sub randname {
  my($lim) = 10000;
  my($fn) = "tmp";
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);

  return $fn;
}
