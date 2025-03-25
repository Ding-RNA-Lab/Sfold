#!/usr/bin/perl -w

# moved the use statements to the top of the
# file as this is my custom
use strict;
use warnings;
use Cwd qw/ cwd /;
use Cwd 'abs_path';

my($basedir) = &cwd();
my($sfoldbin) = abs_path($0);
    
#
# This script determines the deltaG(switch) for a 
# hammerhead ribozyme. deltaG(switch) is the free
# energy cost for the ribozyme to switch from one
# conformation to the conformation that is most
# favorable for target binding. In this script,
# dG{switch} = dG{s} - dG{b}, where dG{s} is the
# average free energy of the 1000 structures predicted
# by Sfold for the ribozyme, and dG{b} is the free
# energy of the binding conformation.
#

# Again a lot of hard coding of directories here
# implys there was some problem with deducing locations

# specify location of the findfe program
    #my($findfe) = "/nfs/sfold/bin/findfe";
my($findfe) = "/$sfoldbin/findfe";    
die " Error: unable to locate program '$findfe'"
  if (!-e $findfe);
my($sfold_param) = "/home/williamrennie/Desktop./Sfoldrun1/Sfold-main/param";
die " Error: unable to locate folder '$sfold_param'"
  if (!-e $sfold_param);

# set auto flush
$| = 1;


&usage() if (scalar(@ARGV) != 3);

my($fe_out) = $ARGV[0];
my($arm_5p) = $ARGV[1];
my($arm_3p) = $ARGV[2];

$arm_5p =~ tr/a-z/A-Z/;
$arm_5p =~ s/T/U/g;
$arm_3p =~ tr/a-z/A-Z/;
$arm_3p =~ s/T/U/g;

die " Error: invalid base found in the 5' arm"
  if ($arm_5p =~ /[^AUCG]/);
die " Error: invalid base found in the 3' arm"
  if ($arm_3p =~ /[^AUCG]/);

my($Rz_core) = "CUGAUGAGUCCGUGAGGACGAA";
my($len_arm_5p) = length($arm_5p);
my($len_arm_3p) = length($arm_3p);
my($seqdata) = join("", $arm_5p, $Rz_core, $arm_3p);
my($seqlen) = $len_arm_5p + length($Rz_core) + $len_arm_3p;
my(@Rz_core_bp) = ();
for (my $i=1; $i<=$seqlen; $i++) {
  $Rz_core_bp[$i] = 0;
}
# set the 4 default base pairs in the binding conformation
$Rz_core_bp[$len_arm_5p+8] = $len_arm_5p+19;
$Rz_core_bp[$len_arm_5p+9] = $len_arm_5p+18;
$Rz_core_bp[$len_arm_5p+10] = $len_arm_5p+17;
$Rz_core_bp[$len_arm_5p+11] = $len_arm_5p+16;
$Rz_core_bp[$len_arm_5p+16] = $len_arm_5p+11;
$Rz_core_bp[$len_arm_5p+17] = $len_arm_5p+10;
$Rz_core_bp[$len_arm_5p+18] = $len_arm_5p+9;
$Rz_core_bp[$len_arm_5p+19] = $len_arm_5p+8;

# temporary file name
my($tmpfname) = &randname();

# build GCG connect file for the binding conformation
open(OUTFILE, ">$tmpfname")
  || die " Error: unable to create temp file";
print OUTFILE "SFOLD of:  [initially 0.00]  Temp Check: 0 from: 1 to: $seqlen\n";
print OUTFILE "Length: $seqlen Energy: 0.00 ..\n";

for (my $i=1; $i<=$seqlen; $i++) {
  print OUTFILE sprintf("%5d %s %5d %5d %5d %5d\n", $i, substr($seqdata, $i-1, 1),
    $i-1, ($i+1)%($seqlen+1), $Rz_core_bp[$i], $i);
}

close(OUTFILE);

# run findfe to determine free energy of binding conformation
my($binding_dG) = &run_findfe($tmpfname);
chomp($binding_dG);
#`rm -f $tmpfname`;

# determine average free energy from Sfold sample
my($avg_fe) = 0;
my($nstruct) = 0;
open(FE_OUT, "$fe_out") || die " Error: unable to read from file '$fe_out'";
while (<FE_OUT>) {
  chomp;
  my(@e) = split;

  $avg_fe += $e[1];
  $nstruct++;
}
close(FE_OUT);
$avg_fe = ($nstruct == 0) ? 0 : ($avg_fe/$nstruct);

printf("%.4f %.4f %.4f\n", $avg_fe, $binding_dG, ($avg_fe-$binding_dG));

exit;


sub randname {
  my($lim) = 10000;
  my($fn) = "tmp";
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);
  $fn .= sprintf "%05d", (rand($lim) + 1);

  return $fn;
}


sub run_findfe {
  my($fn) = @_;
  my($fe);

  # first check to see if the input file contains the energy
  # if so, we don't need to run findfe again.
  my($line) = `grep initially $fn`;
  $line =~ /\[initially\s+([\-\.\d]+)\]/;
  $fe = $1;

  if ($fe ne "" && $fe != 0.0) {
    return $fe;
  }

  # ok, we could not find the free energy inside the input file
  # so we need to run findfe to determine the energy...
  #
  my($cur_dir) = `sh -c 'echo \$PWD'`;
  chomp($cur_dir);
  $cur_dir .= "/" if ($cur_dir !~ /\/$/);
  my($findfe_dir) = $findfe;
  $findfe_dir =~ s/\/[^\/]+$/\//;

  my($path_to_fn);
  if ($fn =~ /^\s*\//) {
    $path_to_fn = $fn;
  } else {
    $path_to_fn = join("", $cur_dir, $fn);
  }

  chdir "$findfe_dir" || die "Error: unable to change to directory '$findfe_dir'";
  $fe = `./findfe $path_to_fn | grep "Free energy" | awk '{print \$4}'`;
  chomp($fe);
  # 2011-02-11 Adam. Repaired regex range [\d\-\.] (was [\d-\.] ).
  $fe = "0.00" if ($fe !~ /^[\-\d\.]+$/);
  chdir "$cur_dir" || die "Error: unable to change to directory '$cur_dir'";

  return $fe;
}


#
# display help information for this program and then exit
#
sub usage {
  print <<EndOfOutput;

Usage: $0 <fe.out> <Rz_5p_arm> <Rz_3p_arm>

EndOfOutput

  exit;
}
