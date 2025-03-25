#!/usr/bin/perl -w

use Getopt::Std;
use Cwd qw/ cwd /;
use Cwd qw/abs_path/;

use strict;
use warnings;

my($bindir) = `dirname $0`;
chomp($bindir);
#my($sfold_root) = abs_path ("$bindir/../");
my($sfold_root) = abs_path ("$bindir/../");

my($sfold_param) = join("", $sfold_root, "/param/");
die " Error: unable to locate folder '$sfold_param'"
  if (!-e $sfold_param);
my($sfold_bin) = join("", $sfold_root, "/bin/");
die " Error: unable to locate folder '$sfold_bin'"
    if (!-e $sfold_bin);

my $opt_string = "a:b:m:s:";
my %opt = ();
my  $tg;
my  (@sites, @h3, @h1);

############################################################################
#   read target sequence
############################################################################
&getopts ($opt_string, \%opt) || &usage;
&usage unless defined $opt{m};

# read target sequence
&readFastaSeq ($opt{m}, \$tg);

############################################################################
#   specify target sites and rz arm lengthes
############################################################################
# starting position of the target sites
@sites  = split(',', $opt{s});

# arm1 length
@h3     = split(',', $opt{a});

# arm2 length
@h1     = split(',', $opt{b});

die " Error: numbers of entries to options -a,-b and -s do not match!"
  if ($#sites != $#h3 || $#h3 != $#h1);

############################################################################
#   read energy terms and set constants
############################################################################
my  (@dangle3, @dangle5);
&read_dangle_dat (\@dangle3, \@dangle5);

# constant term
my  $init   =   4.1;
my  $multi  =   5.2;
my  $term1  =   0.5;
my  $term2  =   0.5;

my  $dangle =   $dangle5[3][0][0] + $dangle3[1][2][3] + $dangle5[1][2][0];

my      $size   =   scalar(@sites);
my  $alph = "ACGU";
my  ($i, $s);

############################################################################
#   calculate dGhybrid for each rz
############################################################################
for ($i = 0; $i < $size; $i++) {
    my  ($stack3, $stack1);
    my  $s  =   $sites[$i];

    my  $pos    =   $s  -   1;
    $stack3 =   &HybrEnRNA ($tg, $pos, $h3[$i], 'r', 0, 0);
    $pos    +=  $h3[$i] + 1;
    $stack1 =   &HybrEnRNA ($tg, $pos, $h1[$i], 'r', 0, 0);

    my  $N2 =   rindex ($alph, substr($tg, $pos, 1));
    my  $N1 =   rindex ($alph, &revcompl(substr($tg, $pos, 1)));

    my  $d1 =   $dangle5[$N1][$N2][1] < $dangle3[3][0][1] ? $dangle5[$N1][$N2][1] : $dangle3[3][0][1];
    $dangle +=  $d1 + $dangle3[$N1][$N2][1];

    # dGhybrid
    my  $sum    =$init + $multi + $term1 + $term2 + $stack1 + $stack3 + $dangle;

    print join "\t", $s, $sum;
    print "\n";
}

sub usage
{
    print STDERR "\n";
    print STDERR "Usage: rz_dghybrid.pl -m target.fasta -s s1,s2,... -a a1,a2,... -b b1,b2,...\n";
    print STDERR "  where s1,s2,... are the starting positions of target sites\n";
    print STDERR "        a1,a2,... are the lengths of Helix III at target sites\n";
    print STDERR "        b1,b2,... are the lengths of Helix I at target sites\n";
    print STDERR "\n";
    exit;
}

sub revcompl {
    my ($str, $dna) = @_;
    my $i;
    for ($i=0; $i<length($str); $i++) {
        if (substr ($str, $i, 1) eq 'a') {
            if (!defined($dna)) {
                substr ($str, $i, 1) = 'u';
            } else {
                substr ($str, $i, 1) = 't';
            }
        } elsif (substr ($str, $i, 1) eq 'u') {
            substr ($str, $i, 1) = 'a';
        } elsif (substr ($str, $i, 1) eq 'g') {
            substr ($str, $i, 1) = 'c';
        } elsif (substr ($str, $i, 1) eq 'c') {
            substr ($str, $i, 1) = 'g';
        } elsif (substr ($str, $i, 1) eq 'T') {
            substr ($str, $i, 1) = 'A';
        }
    }
    my $newstr = "";
    for ($i = 0; $i < length($str); $i++) {
        $newstr .= substr ($str, length($str) - $i - 1, 1);
    }
    return $newstr;
}

# HybrEnRNA (02/24/2005):
#   calculate the hybridization energy an RNA/DNA oligo and an mRNA.
# (assume only wc cononical basepair). The binding site starts at position
# $pos of the mRNA and with length $len. $type specifies the oligo is 
# either DNA ("d") or RNA ("r"). And $overhang tells whether the program 
# penalize the 3' overhang or not (in case of siRNA with 3' TT overhang).
# NOTE: $len doesn't include overhang (if exists)!
sub HybrEnRNA {
  my ($rna, $pos, $len, $type, $overhang, $pel) = @_;
  # $type: "r" -- rna-rna, "d" -- rna-dna.
  # $overhang: 1 -- with overhang (penalize), 0 -- without overhang.
  # $pel: 1 -- apply terminal AU penalty, 0 -- not penalize

  if (($pos + $len - 1 > length($rna))) { 
    print STDERR "HybrEnRNA: invalid oligo position ($pos) or ",
      "length ($len) -- out of bound(." . length($rna) . ").\n"; 
    exit;
  }

  if (!defined($type) || ($type ne "r" && $type ne "d")) {
    print STDERR "HybrEnRNA: must specify the type: r for rna-rna ",
      "hybridiztion, d for rna-dna hybridization.\n";
    exit;
  }

  $overhang = 0 unless (defined($overhang));

  # OligoWalk doesn't penalizing terminal AU for RNA/RNA
  $pel = 0 unless (defined($pel));

  my ($i, $alph, $m, $n, $fe);

  # put reference here.
  my $rd_init = 3.1;
  my @rd_eij = ( [-1.0, -2.1, -1.8, -0.9],
                 [-0.9, -2.1, -1.7, -0.9],
                 [-1.3, -2.7, -2.9, -1.1],
                 [-0.6, -1.5, -1.6, -0.2]
               );
  my $rr_init = 4.09;
  my @rr_eij = ( [-0.9, -2.2, -2.1, -1.1],
                 [-2.1, -3.3, -2.4, -2.1],
                 [-2.4, -3.4, -3.3, -2.2],
                 [-1.3, -2.4, -2.1, -0.9]
               );
  my @rr_mismatch_eij
             = ( [-0.9, -2.2, -2.1, -1.1],
                 [-1.0, -1.5, -1.5, -1.0],
                 [-1.3, -2.5, -2.1, -1.4],
                 [-2.0, -2.5, -2.5, -2.0]
               );
  # 0.45 for INN model, and 0.4 for INN-HB model.
  my $AUpel = 0.45;

  $alph = "ACGU";

  if ($type eq "d") {
    $fe = $rd_init;
  } else {
    $fe = $rr_init;
  }
  for ($i = $pos; $i < ($pos + $len - 1); $i++) {
    $m = rindex ($alph, substr($rna, $i, 1));
    $n = rindex ($alph, substr($rna, $i + 1, 1));
    die "rna_dna_energy: invalid rna sequence (only A/C/G/U allowed)." 
        . substr($rna, $i, 2)
        unless ($m >= 0 && $n >= 0);

    if ($type eq "d") {
      $fe = $fe + $rd_eij[$m][$n];
    } else {
      $fe = $fe + $rr_eij[$m][$n];
    }

    # debug:
    # printf "%s: %f\n", substr ($rna, $i, 2), $rd_eij[$m][$n];
  }

  if ($overhang != 0) { 
    if ($pos <= 0) { 
      print STDERR "Warning[HybrEnRNA]: 3' overhang free energy adjustment not available at position 1!\n";
    } else {
      $m = rindex ($alph, substr($rna, $pos - 1, 1));
      $n = rindex ($alph, substr($rna, $pos, 1));
      $fe += $rr_mismatch_eij[$m][$n];
    }
  }

  # if no overhang, check both ends for termial AU penalty, otherwise, only
  # check the 3' (or the 5' of the siRNA).
  # ONLY penalize the rna-rna hybridization!
  # Updated: both ends should be checked -- the method is symmetrical.
  if ($pel != 0 && $type eq "r") {
    if (substr($rna, $pos + $len - 1, 1) eq "A" ||
        substr($rna, $pos + $len - 1, 1) eq "U") {
      $fe += $AUpel;
    }
    if (substr($rna, $pos, 1) eq "A" ||
        substr($rna, $pos, 1) eq "U") {
      $fe += $AUpel;
    }
# does 3' of the as siRNA matter in binding the target site?
#   if ($overhang == 0 && 
#       (substr($rna, $pos, 1) eq "A" || substr($rna, $pos, 1) eq "U")) {
#     $fe += $AUpel;
#   }
  }

  return $fe;
}
############################################################
## A set of subroutines that read the free energy tables
## from mfold package (in dat/ directory). They should also
## work for sfold versions of the energy tables.
## -- Yu Shao 04/14/2005
############################################################
# read the single stacking energy table (dangle.dat)
sub read_dangle_dat
{
  my ($dangle3, $dangle5, $datf) = @_;

  @$dangle3 = ();
  @$dangle5 = ();
  $datf = "$sfold_param/dangle.dat"
    unless defined($datf);

  my ($fh, $line, @fields, $i, $j, $k);

# 3' dangling
# 5' --> 3'
#   i k
#   j
# 3' <-- 5'

# 5' dangling
# 5' --> 3'
#   i
#   j k
# 3' <-- 5'
  $i = 0;
  open ($fh, $datf) || die "in read_dangle_dat: couldn't open $datf ($!)";
  while ($line = <$fh>) {
    $line =~ s/^\s+//;
    if ($line =~ /[A-Z'>]/ || $line =~ /--/ || $line =~ /^$/) { next; }
#   print $line;
    @fields = split (/\s+/, $line);
#   print @fields;
    for ($j = 0; $j < 4; $j++) { 
      for ($k = 0; $k < 4; $k++) {
        if ($i <= 3) {
          $$dangle3[$i][$j][$k] = $fields[$j * 4 + $k];
        } else {
          $$dangle5[$i - 4][$j][$k] = $fields[$j * 4 + $k];
        }
      }
    }
    $i++;
  }
  close ($fh);

# for ($i = 0; $i < 4; $i++) {
#   for ($j = 0; $j < 4; $j++) {
#     for ($k = 0; $k < 4; $k++) {
#       printf "%5s  ", $$dangle3[$i][$j][$k];
#     }
#   }
#   print "\n";
# }
# for ($i = 0; $i < 4; $i++) {
#   for ($j = 0; $j < 4; $j++) {
#     for ($k = 0; $k < 4; $k++) {
#       printf "%5s  ", $$dangle5[$i][$j][$k];
#     }
#   }
#   print "\n";
# }

}

## 1. read in RNA sequence from fasta format.
## 2. ignore header info.
## 3. conver all T's to U
## 4. $seq passed by reference
sub readFastaSeq {
      my ($fname, $seq) = @_;
        my ($infh, $line);

          $$seq = "";

            open ($infh, $fname) || die "$fname: file open failed.";
              while ($line = <$infh>) {
                      chomp ($line);
                          $line =~ s/^\s*//;
                              if ($line =~ /^>/ || $line =~ /^#/) {
                                        next;
                                            }
                                                $$seq .= $line;
                                                  }
                                                    close ($infh);

                                                      $$seq =~ s/\s*//g;
                                                        $$seq =~ s/\d*//g;
                                                          $$seq = uc ($$seq);

                                                            $$seq =~ s/T/U/g;
}

