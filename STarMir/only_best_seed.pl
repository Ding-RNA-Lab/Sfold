#!/usr/bin/perl -w
#This script is used by the STarMir module.

use strict;
use warnings;

#use Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

#only_best_seed.pl by Adam Wolenc for Wadsworth Center 2010-10-04
# 2022 -4- 14 (war)
# edited to remove requirement for the Carp::Assert package

#Filter rnahybrid output so there is at most one site for within +/- fudge factor (1)
# for every position on the mRNA target.
#When there is competition for two sites in the same pos on the mRNA (within fudge factor 1)
# (use ending position of last helix rather than starting position, for this purpose),
# use the following criteria to decide on the best one:
# 1. if site A has a better seed label than site B, as determined by previewing the seed labels
#  that would be applied to the given conformations, then chose site B
# 2. otherwise, if site A's label is worse than site B's seed label, choose site A.
# 3. otherwise, both sites have the same seed type. Choose the site with the lower (more 
#  favorable) free energy.
# 2010-12-16 Added %conf_distance to allow offset-6mer's to have a 2nt range of conflict
#  the other seed types have only a 1nt range of conflict.

my $InFname;
my $OutFname;
if ($#ARGV == 1) {
    $InFname = $ARGV[0] ;
    $OutFname= $ARGV[1] ;
}
else {
   die("usage: $0 <Input_File> <Output_File>\n"); 
}

open(IN,"<$InFname") || die "Unable to open the $InFname file to read" ;

my $fudge = 2;  #look +/- this number for potential conflicts.
my %data; #store all data from the input file here

#ranking of seed labels
my %conformation_code = ('8mer', 5, '7mer-m8', 4, '7mer-A1', 3, '6mer', 2, 'offset-6mer', 1, 'other', 0);
my %conf_code_rev  = (5, '8mer', 4, '7mer-m8', 3, '7mer-A1', 2, '6mer', 1, 'offset-6mer', 0, 'other');
my %conf_distance  = (5, 1,      4, 1,         3, 1,         2, 1,      1, 2,             0, 1);
#                                                                          ^ special offset-6mer rule

#keep rough track of why sites are rejected
my $better_conf = 0;
my $better_energy = 0;
my %better_than;

while (<IN>) {
  # 0  X66610-3pUTR
  # 1  317
  # 2  hsa-let-7a
  # 3  22
  # 4  -22.0
  # 5  0.488194
  # 6  36
  # 7  C    GUC        UAGA    C
  # 8   ACUA    ACCUACU    CUCA 
  # 9   UGAU    UGGAUGA    GAGU 
  #10  U    AUGU       UG       

  chomp ;
  my @rec = split(/:/) ;
  my $tar_name = $rec[0] ; 
  my $tar_len = $rec[1] ;
  my $mir_name = $rec[2] ; 
  my $h_ener = $rec[4] ;
  my $be_tar = $rec[6] ;

  my $tar_misses = lc(reverse($rec[7]));
  my $tar_matches = lc(reverse($rec[8]));
  my $mir_matches = lc(reverse($rec[9]));
  my $mir_misses = lc(reverse($rec[10]));

  my $tar_either = merge($tar_matches, $tar_misses);
  my $tar_either_no_space = $tar_either;
  $tar_either_no_space =~ s/ //g;

  #compute length of site on target
  my $end_tar = $be_tar+ length($tar_either_no_space)- 1 ;

  my $conformation="other"; #default

  #consider only positions 1-10 for seed-labeling
  my $mir_either = merge($mir_matches, $mir_misses);
  my $i=0;
  my $count=0;
  my @pos_map;
  for (my $i=1; $i<=10; $i++){
    $pos_map[$i] = -1;
  }
  while ($pos_map[10] == -1 && $i < length $mir_either) {
    my $c = substr($mir_either, $i, 1);
    if ($c ne " ") {
      ++$count;
      $pos_map[$count] = $i;
    }
    ++$i;
  }

  next if ($pos_map[1] == -1); #hopeless.

  #Look for A in position 1 of target. It can be either miss or match.
  my $tA1 = (substr($tar_either, $pos_map[1], 1) eq "a") ? 1 : 0;

  #gather all possible seed annotations.
  my %possible_annotation;
  #walk through the diagram to determine seed len and pos
  for (my $start_pos=2; $start_pos<=8; ++$start_pos) {

    my $mirlen = 0; #to start

    for (my $i=$pos_map[$start_pos]; $i<=$pos_map[8] && $i != -1 ; ++$i) {

      my $tar_miss_ch = substr($tar_matches, $i, 1);
      my $tar_match_ch = substr($tar_matches, $i, 1);
      my $mir_match_ch = substr($mir_matches, $i, 1);
      my $mir_miss_ch = substr($mir_misses, $i, 1);

      #there are 3 possibilities. 0 means space, 1 means base
      #tar_miss
      # mir_match
      #  mir_miss
      #010-WC = match, mir len increments
      #010-GU, some GU's left = match, mir len increments, gus left decrments
      #010-GU, no GU's left = essentially a bulge in both, exit loop
      #?0? = mismatch or bulge. exit loop.

      if ($mir_match_ch ne " ") {
        if (($tar_match_ch eq "g" && $mir_match_ch eq "u") ||
            ($tar_match_ch eq "u" && $mir_match_ch eq "g")) {
          #010-GU, no GU's allowed
          last;
        } else {
          #010-WC
          ++$mirlen;
        }
      } else {
        #?0?
        last;
      }

      my $end_pos = $start_pos + $mirlen - 1;
      my $label = "$mirlen-mer $start_pos-$end_pos";
      $possible_annotation{$label} = $mirlen;
    }
  }

  #look for primary conformations
  if ($tA1 == 1 && $possible_annotation{"7-mer 2-8"}) {
    $conformation = "8mer";
  } elsif ($tA1 == 0 && $possible_annotation{"7-mer 2-8"}) {
    $conformation = "7mer-m8";
  } elsif ($tA1 == 1 && $possible_annotation{"6-mer 2-7"}) {
    $conformation = "7mer-A1";
  } elsif ($tA1 == 0 && $possible_annotation{"6-mer 2-7"}) {
    $conformation = "6mer";
  } elsif ($possible_annotation{"6-mer 3-8"}) {
    $conformation = "offset-6mer";
  }

  my $seen_before = 0;
  for (my $o = -$fudge; $o <= $fudge; ++$o) {
    my $pos=$end_tar+$o;  #used for finding previous matches within this area +/- $fudge
    if (defined $data{$mir_name}->{$tar_name}->{$pos}) {
      #possible conflict.  Is it in range?
      # determine the conflict range for this pair. The rules are,
      # conflict range is always 1 unless one site is at least offset-6mer and
      # the other site is better than offset-6mer.
      my $conflict_range = 1;
      if ( ($conformation eq 'offset-6mer'
         && $data{$mir_name}->{$tar_name}->{$pos}->{'type'} > $conformation_code{'offset-6mer'})
        || ($conf_code_rev{$data{$mir_name}->{$tar_name}->{$pos}->{'type'}} eq 'offset-6mer'
         && $conformation_code{$conformation} > $conformation_code{'offset-6mer'}) ) {

        if ($conf_distance{$data{$mir_name}->{$tar_name}->{$pos}->{'type'}} > $conflict_range) {
          $conflict_range = $conf_distance{$data{$mir_name}->{$tar_name}->{$pos}->{'type'}};
        }

        #print "$_ versus\n";
        #print $data{$mir_name}->{$tar_name}->{$pos}->{'record'} . "\n";
        #print "" . <STDIN>;
      }
      next if (abs($o) > $conflict_range);  # sites are too far away.
      
      #confirmed conflict. Which one wins?
      my $better = 0;

      if (   $conformation_code{$conformation} >  $data{$mir_name}->{$tar_name}->{$pos}->{'type'}) {
        $better = 1;
        #current has better conformation than previous
        ++$better_than{$conformation_code{$conformation}}->{$data{$mir_name}->{$tar_name}->{$pos}->{'type'}};
        ++$better_conf;
      } else {
        if ($conformation_code{$conformation} < $data{$mir_name}->{$tar_name}->{$pos}->{'type'}) {
          #previous has better conformation than current
          ++$better_than{$data{$mir_name}->{$tar_name}->{$pos}->{'type'}}->{$conformation_code{$conformation}};
          ++$better_conf;
        } elsif ($h_ener <  $data{$mir_name}->{$tar_name}->{$pos}->{'deltaG'}) {  #lower is better
          $better = 1;
          #current has better energy than previous
          ++$better_energy;
        } else {
          #previous has better energy than current
          ++$better_energy;
        }
      }

      if ($better) {
        #replace.
        delete $data{$mir_name}->{$tar_name}->{$pos};
        $data{$mir_name}->{$tar_name}->{$end_tar}->{'type'} = $conformation_code{$conformation};
        $data{$mir_name}->{$tar_name}->{$end_tar}->{'deltaG'} = $h_ener;
        $data{$mir_name}->{$tar_name}->{$end_tar}->{'record'} = $_;
      }
      $seen_before = 1;
      last;
    }
  }
  if (!$seen_before) {
    #novel entry.
    $data{$mir_name}->{$tar_name}->{$end_tar}->{'type'} = $conformation_code{$conformation};
    $data{$mir_name}->{$tar_name}->{$end_tar}->{'deltaG'} = $h_ener;
    $data{$mir_name}->{$tar_name}->{$end_tar}->{'record'} = $_;
  }
}
close IN;

my $accepted=0;
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write" ;
for my $mir_name (keys %data) {
  for my $tar_name (keys %{$data{$mir_name}}) {
    for my $end_tar (keys %{$data{$mir_name}->{$tar_name}}) {
      print OUT $data{$mir_name}->{$tar_name}->{$end_tar}->{'record'} . "\n";
      ++$accepted;
    }
  }
}

my $rejected = ($better_conf + $better_energy);
print "$accepted Accepted\n";
print "$rejected Rejected\n";
print "The following breakdown is an  estimate and may vary depending on the sort order of the input data.\n";
print "($better_conf Rejected for inferior conformation)\n";
print "($better_energy Rejected for inferior energy)\n";
for my $a (sort keys %better_than) {
  for my $b (sort keys %{$better_than{$a}}) {
    print "\"" . $conf_code_rev{$a} . "\" beats \"" . $conf_code_rev{$b} . "\" "
      . $better_than{$a}->{$b} . " times.\n";
  }
}
close OUT ;

print STDERR "Exiting $0\n" if $debug;

exit 0;


sub merge {
  #fill in blank spaces in $line1 with characters from
  # $line2, resulting in a single, merged sequence.
  #line1 and line2 must be the same length for this to
  # make any sense.
  my $line1 = shift;
  my $line2 = shift;
  my $ret = "";
  for (my $i=0; $i<length $line1; ++$i) {
    if (substr($line1, $i, 1) ne " ") {
      $ret.=substr($line1, $i, 1);
    } else {
      $ret.=substr($line2, $i, 1);
    }
  }
  return $ret;
}
