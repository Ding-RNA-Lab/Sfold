#!/usr/bin/perl -w
#This script is used by the STarMir module.
# Calculate disruption Energy 
# usage > perl Calc_disruptEn.pl <target_name> <query_name> <SFold_outdir>

# 2012-5-20 (war)
# ported ccliu changes.
# -----------------------------------------------------
# 2022-4-14 (war)
# removed calls to Carp::Assert package

use Cwd qw/ abs_path /;
use warnings;
use strict;

#use Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# get the directory in which this script is located
my $base_dir  = `dirname $0`;
chomp $base_dir;
$base_dir = abs_path($base_dir);

# load parameters contained in option file
require "$base_dir/starmir_param.pl";
our $SFOLDBIN;
our $disruptEn_prog;

# Globals to hold arguments
my $folding_seq_fname;
my $bposf;
my $outfname;
my $SFold_outdir;

if ($#ARGV == 3) {
    $folding_seq_fname = $ARGV[0] ; 
    $bposf = $ARGV[1] ;  
    $outfname = $ARGV[2] ;
    $SFold_outdir  = $ARGV[3] ;  
} else {
    die("usage: > perl Calc_disruptEn_research.pl <folding_mRNA_fname> <bsite_file> <output_file> <SFold_outdir>\n") ;
}

#check the existence of input files for disruptEn_prog.
die " Error: unable to locate $bposf" if (!-e $bposf) ;

my $structfname = "bp.out" ; 
my $febfname= "fe.out" ;

my $structf = join("",$SFold_outdir,$structfname) ;
die " Error: unable to locate $structf" if (!-e $structf) ;

my $febf = $SFold_outdir.$febfname ;
die " Error: unable to locate $febf" if (!-e $febf) ;

# $folding_seq_fname = "$folding_seq_fname" ;
die " Error: unable to locate $folding_seq_fname" if (!-e $folding_seq_fname) ;

system("SFOLDBIN=$SFOLDBIN $disruptEn_prog -m $folding_seq_fname -p $bposf -s $structf -c 0 -f $febf -t 1 > $outfname") ;

print STDERR "Exiting $0\n" if $debug;

exit 0;

