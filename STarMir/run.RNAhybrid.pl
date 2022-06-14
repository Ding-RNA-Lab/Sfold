#!/usr/bin/perl -w

#This script is used by the STarMir module.
# Run RNAhybrid to identify candidate sites.
# 2011-02-01 Adam Synchronized with Prediction Steps.
# 2011-02-02 TODO (style) Convert to shell script.
# 2012-05-20 (war)
# adapted to conform to changes made to research code by ccliu
# 2022-04-14 (war)
# removed "Carp::Assert references

# reguired the BioPerl package

# usage > perl run.RNAhybrid.pl <DB_file_name> <miR_file_name> <Hyb_Ener_threshold> <Out_file_name> <RNAhyb_bindir> <training species>
# -----------------------------------------------------

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

# use Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# globals holding arguments
my $DBFname; 
my $MFname; 
my $Hyb_th; 
my $OutFname; 
my $RNAhyb_bindir; 
my $Hyb_species = 'human'; # this is the default

if ($#ARGV == 4 or $#ARGV == 5) {
  $DBFname =  $ARGV[0];
  $MFname = $ARGV[1];
  $Hyb_th = $ARGV[2];
  $OutFname = $ARGV[3];
  $RNAhyb_bindir = $ARGV[4];
  $Hyb_species = $ARGV[5] if defined $ARGV[5]; # we are assuming value is valid
} 
else {
  die("usage: > perl run.RNAhybrid.pl <DB_file_name> <miR_file_name> <Hyb_Ener_threshold> <Out_file_name> <RNAhyb_bindir> <Training species>\n");
}

### repaired, needs to be removed.
# the problem is related to the way the Starmir webserver processes input sequences.
# It is probably possible to eliminate the bioperl sequence below if you can be certain
# line lengths in the input fasta files are under 1000 nts long.\

# kludge 4/24/2012 (war)
# The file passed to us is in pseudo FASTA format.  The sequence is in one line which
# can be thousands of nts long.  RNAhybrid does not tolerate lines more than about 1,000 nts long.
# We read the file in and rewrite a local file to clean up the format.
# This should be corrected a lot sooner than here, but this is the first place the starmir protocol uses
# the sequence file.
# 2022-4-25 (war)
my $input_seqio = Bio::SeqIO->new(
    -file => $DBFname,
    -format => 'fasta');

my $output_seqio = Bio::SeqIO->new(
    -file => ">$DBFname.corr",
    -format => 'fasta');

my $seq = $input_seqio->next_seq();

$output_seqio = $output_seqio->write_seq($seq);

# run RNAhybrid 4 times to make sure we get all MIT-style seeds regardless of Ghyb

system("$RNAhyb_bindir/RNAhybrid -c -s 3utr_$Hyb_species -f 2,7 -e -1.0 -m 10000   -t $DBFname.corr -q $MFname > $OutFname-1");
system("$RNAhyb_bindir/RNAhybrid -c -s 3utr_$Hyb_species -f 2,8 -e -1.0 -m 10000  -t $DBFname.corr -q $MFname > $OutFname-2");
system("$RNAhyb_bindir/RNAhybrid -c -s 3utr_$Hyb_species -f 3,8 -e -1.0 -m 10000  -t $DBFname.corr -q $MFname > $OutFname-3");
system("$RNAhyb_bindir/RNAhybrid -c -s 3utr_$Hyb_species        -e $Hyb_th -m 10000 -t $DBFname.corr -q $MFname > $OutFname-5");

# removed because the merge step in the protocol removes duplicates, which this does not.
# the merge_123.sh script calls the collator executable, which in fact does not collate, but
# removes duplicates. 4/19/2012 (war)
# system("cat $OutFname-? | sort -u > $OutFname");
# system("rm $OutFname-?");

# merge files
print STDERR "running collator on RNAhybrid output\n" if $debug;
system( "./collator $OutFname-1 $OutFname-2 $OutFname-3 $OutFname-5 > $OutFname");
system("rm $OutFname-?");

#apply overlap filter
# not part of mouse/human protocol removed 4/19/2012 (war)
# system("/nfs/sfold/bin/post_filter_overlap.pl $OutFname");
# system("mv Passed-$OutFname $OutFname");
# system("rm *-$OutFname");

#apply only best seed filter

print STDERR "running only_best_seed.p. $OutFname Filtered-$OutFname-Filtered\n" if $debug;
system("./only_best_seed.pl $OutFname $OutFname-Filtered");
system("mv $OutFname-Filtered $OutFname");

print STDERR "Exiting $0\n" if $debug;

exit 0;
