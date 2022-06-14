#This script is used by the STarMir module.
#
# Changeable options for the STarMir module
#

##
# Locations of programs and files
##

# directory of RNAhybrid executable
$RNAhyb_bindir = "/usr/local/bin/";
# directory of disruptEn executable
$disruptEn_bindir = "/home/bill/Desktop/finalFiles/Sfold_2010/bin/";
# directory of the findfe executable
$SFOLDBIN = "/home/bill/Desktop/finalFiles/Sfold_2010/bin/";

# path to the disruptEn executable
$disruptEn_prog = join("", $disruptEn_bindir, "disruptEn");


##
# Tunable values for the execution of the programs
##
# threshold for hybridization energy (RNAhybrid)
$Hyb_threshold = -15.0 ;
# nucleation energy threshold
$NE_threshold = 4.09 ;


1;
