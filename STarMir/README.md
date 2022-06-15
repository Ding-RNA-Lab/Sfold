# STarMir
 STarMir scripts and tables

Requirements

STarMir requires that RNAhybrid be installed to generate candidate
binding site.  You will also need the most recent version of R.

Sfold must be installed and tested.

The scripts are mostly in Perl, Perl 5 has to be installed and functioning.
The BioPerl library, must also be installed.

Configuring.

The file starmir_param.pl must be edited to give paths to the
sfold and RNAhybrid binaries.  Some tunable values for RNA hybrid must be
included here.

Running

Sfold must be run on the mRNA sequence first.  Make note of the directory the
Sfold output is in.


The top level script is starmir_research.pl.  It must be run in the directory
that the remaining scripts are located in.  You should not move or rename
any directory or file in that directory.

The instruction for running starmir_research.pl are in the header comment
to the file.

