#!/usr/bin/perl 

# starmir_research.pl
#
# used to predict the binding sites of one mRNA paired with a list of miRNAs and compute all features.
# ONLY used for internal research purpose, modified from starmir_2_call.pl written by Adam, War,
# which is used for StarMir v2.0 for webserver.

# Usage: starmir_research.pl <full_miRNA_file> <prediction mRNA_file> <folding mRNA filename> <Sfold_outdir> <localStore> <RNAhyb_species> <start coding region><end coding region>
# Arguments:
#     $ARGV[0] -- FASTA file containing miRNA sequence
#     $ARGV[1] -- FASTA file containing prediction/target mRNA sequence
#     $ARGV[2] -- Sfold output directory, top level output 
#     $ARGV[3] -- species for RNAhybrid (value: human, worm, fly);
# New parameters added on 11/2/2012 by ccliu
#     $ARGV[7] -- species for model: human or mouse (optional)
#     $ARGV[8] -- start coding region
#     $ARGV[9] -- end coding region
#     $ARGV[10] -- working directory     

# 5/16/2012 (ccliu)
# modified from starmir_2_call.pl for internal research purpose.
# Preconditions:
# 
# sfold must have been run.
# processing is dependent specifically
# on the existence of bp.out, fe.out  and sstrand.out
# which are assumed to be present in the sfold output directory sent as parameter
#

#
# 3/26/2012 (war)
# modified to reflect the change in the name of the files
# output by the regression R script.
#
# 5/20/2012 (war)
# added some of the changes in the research scripts.
#
# 2013-4-30 (war)
# modifications for the inclusion of worm model
#
# 2014-3-19 (war)
# contains modifications for multiple miRNAs
# mostly it changes the sense of a few variables.  The processing
# is almost identical.
#
# changes for distribution
# including working directories and
# removal of temporary files.
#
use warnings;
use strict;

#use Carp::Assert; # use to define $debug no for silent
use File::Basename;
use Cwd;
my $printCommand = 0;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

use Cwd qw/ abs_path /;

# get the directory in which this script is located
# that will be the base directory for all other scripts
$0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;
my $base_dir = abs_path($1);

# derive the starmir-params directory from the 
# execution directory $base_dir
# starmir-params contains files needed for execution of the
# logistic probability calculation

$base_dir =~ m/^([\/.]+)\/.+\/bin$/;
my $starmir_param_dir = $1;

print STDERR "starmir-param directory: $starmir_param_dir\n" if $debug; 

# load parameters contained in starmir parameters file
# confirm these before launch and edit if necessary
#   principly locations of binaries,
#   the hybirdization energy threshold
#   and the nucleation energy threshold.
require "$base_dir/starmir_param.pl";


# Arguments: 
die(" Usage: $0 <working directory <miRNA_file> <full_mRNA_file>  <sfold output directory> <target species> /n /t <model species> <CDS start> <CDS end> <delete temps(defaults to true)>\n") 
    unless scalar @ARGV >= 7;

(my $workDir, my $MirFname, my $prediction_seq_fname, my $SFold_outdir, my $targetspecies, my $trainspecies, my $cdsstart, my $cdsend) = @ARGV;

my $delTemps = 1;
#if (define $ARGV[7]) $delTemps = $ARGV[7];



# 1: run RNAhybrid and parsing the output
my $mRNA_name = $prediction_seq_fname;  #may contain directory.
   $mRNA_name = (split(/\//,$mRNA_name))[-1]; #remove directory
   $mRNA_name =~ s/\.\w+$//;  
# mRNA name will be the non path part of the filename.  $mRNA_name

our $RNAhyb_bindir;  # from starmir_param.pl
our $Hyb_threshold;  # from starmir_param.pl

# 1a: Run RNAhybrid
#     input -- prediction sequence filename
#                essentially the same as the total sequence in
#                this version
#              miRNA sequence filename
my $HybFname = "$workDir/Hyb-$mRNA_name.out"; # output name

# Note: run.RNAhybrid.pl runs the RNAHybrid binary and
# one  other script run to clean output
#     only_best_seed.pl
if ($trainspecies ne 'mouse') {
    if ($printCommand) {
	print STDERR "perl $base_dir/run.RNAhybrid.pl $prediction_seq_fname $MirFname $Hyb_threshold $HybFname $RNAhyb_bindir $trainspecies\n";
    }
    else {
	system("perl $base_dir/run.RNAhybrid.pl $prediction_seq_fname $MirFname $Hyb_threshold $HybFname $RNAhyb_bindir $trainspecies");
    } # end if else printCommand
} # end if not mouse
else {
    if ($printCommand) {
	print STDERR "perl $base_dir/run.RNAhybrid.pl $prediction_seq_fname $MirFname $Hyb_threshold $HybFname $RNAhyb_bindir\n";
    }
    else {
	system("perl $base_dir/run.RNAhybrid.pl $prediction_seq_fname $MirFname $Hyb_threshold $HybFname $RNAhyb_bindir");
    } # end if else printCommand
} # end else mouse

# at this point we have RNA output in compact form.  These will be our candidate sites.
# below this we are calculating features.

# 1b: run ParseHybFile.pl
#     input -- Hyb-$mRNA_name.out, output of run.RNAHybrid
my $ResHybFname = "$workDir/ResHyb-$mRNA_name.out"; # output name

if ($printCommand){
    print STDERR "perl $base_dir/ParseHybFile.pl $HybFname $ResHybFname\n";
}
else {
    system("perl $base_dir/ParseHybFile.pl $HybFname $ResHybFname");
}

# 2: run  ParseHybFile_2.pl 
#    input -- Hyp-$mir_name.out, from run.RNAHybrid 
my $SiteFname = "$workDir/Site-$mRNA_name.out"; # output name
if ($printCommand) {
    print STDERR "perl $base_dir/ParseHybFile_2.pl $HybFname $SiteFname\n";
}
else {
    system("perl $base_dir/ParseHybFile_2.pl $HybFname $SiteFname");
}

# 3: run SeedRegionFeaturesMIT.pl
#   input -- Hyp-$mir_name.out, from run.RNAHybrid output
my $SeedFname = "$workDir/Seed-$mRNA_name.out"; # output name
if($printCommand) {
    print STDERR "perl $base_dir/SeedRegionFeaturesMIT.pl $HybFname $SeedFname \n";
}
else {
    system("perl $base_dir/SeedRegionFeaturesMIT.pl $HybFname $SeedFname");
}

# 4: assemble seed and site data ***
#    input -- Site-$mir_name.out
#             Seed-$mir_name.out
#             sstrand.out, from SFold
my $SiteSeqFeatFname = "$workDir/SiteFeatures-$mRNA_name.out";
if($printCommand) {
    print STDERR "perl $base_dir/GetSiteSeqFeaturesMIT.pl $SiteFname $SeedFname $SFold_outdir/sstrand.out $SiteSeqFeatFname\n";
}
else {
    system("perl $base_dir/GetSiteSeqFeaturesMIT.pl $SiteFname $SeedFname $SFold_outdir/sstrand.out $SiteSeqFeatFname");
}


# 5: write out the possible binding sites
#    input -- ResHyb-$mir_name.out, ParseHybFile.pl output
#             Pos_Adj, is always 0
#    output $mir_name-bsites.txt
my $BsitesFname = "$workDir/$mRNA_name-bsites.txt"; # output file
#2011-02-03 Adam. Updated to Bsites.pl, which resolves off-by-one error in disen input ranges.
if ($printCommand){
    print STDERR "perl $base_dir/Bsites.pl $ResHybFname 0 $BsitesFname\n";
}
else {
    system("perl $base_dir/Bsites.pl $ResHybFname 0 $BsitesFname");
}

# 6: calculate disruption energy
#    input - fe.out, bp.out, from SFold
#    output - disruptEn_$mir_name.out
# Note: $prediction_seq_fname is the 
#  file containing the sequence we
#  gave to Sfold
my $tar_name = $prediction_seq_fname;
$tar_name =~ s/\.\w+$//;
my $disruptEnFname = "$workDir/disruptEn-$mRNA_name.out" ;  #output filename
if($printCommand) {
    print STDERR "perl $base_dir/Calc_disruptEn.pl $prediction_seq_fname $BsitesFname $disruptEnFname $SFold_outdir\n";
}
else {
    system("perl $base_dir/Calc_disruptEn.pl $prediction_seq_fname $BsitesFname $disruptEnFname $SFold_outdir");
}

# 7: calculate nucleation energies
#    input -- ResHyb-$mRNA_name.out, ParseHybFile.pl output
my $HybFilFname = "$workDir/Hyb-Fil-$mRNA_name.out"; # output name

#2011-02-02 Adam. Updated to CalcNE 4.2, which matches the publication.

if($printCommand){
    print STDERR "perl $base_dir/CalcNE.4.2.pl $ResHybFname 0 $HybFilFname $SFold_outdir";
}
else {
    system("perl $base_dir/CalcNE.4.2.pl $ResHybFname 0 $HybFilFname $SFold_outdir");
}

# 8: assemble the data to calculate the total interaction energy
#    input -- $mir_name-bsites.txt, from Bsites.pl
#             disruptEn_$mRNA_name.out, from Calc_disruptEn.pl
my $TotalEnFname = basename($HybFilFname); # output filename
$TotalEnFname = dirname($HybFilFname) . "/" . "TotalEn-" . $TotalEnFname;

# 2011-02-04 Adam. Updated to use essentially a functional combination of 
# Prediction Steps' Match_disEn and Match_totalEn
# starmir_2 doesn't care if there is not enough hyb to overcome init
if ($printCommand){
    print STDERR "perl $base_dir/Get_totalEn.pl $HybFilFname 0 $BsitesFname $disruptEnFname $TotalEnFname 1\n";
}
else {
    system("perl $base_dir/Get_totalEn.pl $HybFilFname 0 $BsitesFname $disruptEnFname $TotalEnFname 1");
}
# Now sort the output
if ($printCommand){
    print STDERR "perl $base_dir/Sort_TotalEn.pl $TotalEnFname $TotalEnFname\n";
}
else {
    system("perl $base_dir/Sort_TotalEn.pl $TotalEnFname $TotalEnFname");
}

# 9: Match total energy to sites  ***
#    input -- TotalEn-$mRNA_name, from Get_TotalEn.pl
#             SiteFeatures-$mir_name.out, from ParseHyp_2.pl
my $SiteFeaturesEnFname = "$workDir/SiteFeatures-En-$mRNA_name.out"; # output filename

if ($printCommand){
    print STDERR "perl $base_dir/Match_TotalEn.pl $TotalEnFname $SiteSeqFeatFname $SiteFeaturesEnFname\n";
}
else {
    system("perl $base_dir/Match_TotalEn.pl $TotalEnFname $SiteSeqFeatFname $SiteFeaturesEnFname");
}

# 10: Create final Feature files  ***
#    input SiteFeatures_En-$mRNA_name.out, from Match_TotalEn.pl
#          Seed-$mir_name.out, from SeedRegionFeatures.pl
my $SiteFeatEnSeedFname = "$workDir/SiteFeatures-En-Seed-$mRNA_name.out";

system("perl $base_dir/Match_SiteFeatures_Seed.pl $SiteFeaturesEnFname $SeedFname $SiteFeatEnSeedFname");


# 11: Split files
#   splits files by region and the presence of a site that
#   qualifies as a seed.
#   also adds additional columns related to accessibility.
#
if($printCommand) {
    print STDERR "perl $base_dir/SplitFilesWeb.pl $base_dir $base_dir/../starmir-param $SiteFeatEnSeedFname $cdsstart $cdsend $trainspecies\n";
}
else {
    system("perl $base_dir/SplitFilesWeb.pl $base_dir $base_dir/../starmir-param $SiteFeatEnSeedFname $cdsstart $cdsend $trainspecies");
}


# The logit regression step is the point of all.  It takes the features and generates a
# likelyhood for the binding site

# 12: run logistic regression
#    input - six files of the form
#      bsite-$type-$component-$trainspecies.txt
#        where $type == ("seed", "seedless")
#        and $component == ("3pUTR","CDs", 5pUTR)
#    output -- six files of the form
#      Output-bsite-$type-$component-$trainspecies.txt
#
my @args = (
    'perl',
    "$base_dir/starmir-param/run.logitWeb.pl", # location perl script that executes the model, also locates conservation info if present
    "$workDir/",             # working directory for this job
    "$base_dir/starmir-param",                   # root directory for starmir parameter files
    $targetspecies,
    $trainspecies,
    $cdsstart, $cdsend,
    $prediction_seq_fname,
    'Full',  # folded region, report hits from this region
    0        # did we lookup the sequence information. Always false in distribution version
    );

if ($printCommand) {
    print STDERR join(',', @args) . "\n" ;
}
else {
    system(@args) == 0 or die("logistic regression script did not run\n");
}

# 13 run scoring script
#    input six files of the form Output-$type-$component-$trainspecies.txt
#    output -- six files of the form
#     Final-bsite-$type-$component-$trainspecies.txt
#


system("perl $base_dir/run.score.pl $workDir  $trainspecies");


print STDERR "Exiting $0\n" if $debug;

exit 0;
