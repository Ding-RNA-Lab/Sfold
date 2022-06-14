#!/usr/bin/perl -w

# MatchSiteFeatures_Seed.pl
#
# a script to match Site's Features values in the first file
# to the seed information in the second file
# usage > perl Match_SiteFeatures_Seed.pl <Feature_file> <Seed_file>  

# 2010-12-29 Adam. strict, cleanup, warnings, etc.
#
# 2011-8-2 (war) adapted for web
# Added an explicit argument for the output file.
# Added a header line to the output
#
# 2011-9-26 (war)
# Added 'SiteID' to cols created in the output file
# This supports using the site ID to maintain a relationship
# between the data and the conformation in web output.
#
# 2022-4-20 (war)
# edited to remove references to Carp
use warnings;
use strict;

#no Carp::Assert;
my $debug = 0;

print STDERR "\nEntering $0\n" if $debug;
print STDERR "Arguments:\n\t" if $debug;
print STDERR join("\n\t", @ARGV) . "\n" if $debug;

# Globals
my $InFname1;
my $InFname2;
my $OutFname;

# Fields in File 1
my $SITEID1 = 0;
my $SEQNAME1 = 1;
my $SEQLEN1 = 2;
my $GHYB1 = 4;
my $BEG1 = 5;
my $END1 = 6;

# Fields in File 2
my $SEQNAME2 = 0;
my $SEQLEN2 = 1;
my $GHYB2 = 3;
my $POS2 = 4;

if (scalar @ARGV == 3) {
    $InFname1 =  $ARGV[0];
    $InFname2 = $ARGV[1];
    $OutFname = $ARGV[2];
} else {
   die("usage:>perl Match_SiteFeatures_Seed.pl <Ener_file> <Seed_file> <Output_name>\n");
}

open(IN1,"<$InFname1") || die "Unable to open the $InFname1 file to read";
open(IN2,"<$InFname2") || die "Unable to open the $InFname2 file to read";
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";

output_header(*OUT);

my %Sites = ();
while (<IN1>) {
  #   1
  #   1  NM_175116
  #   2  706
  #   3  let-7a
  #   4  -20.4
  #   5  8
  #   6  37
  #   7  0.495
  #   8  23
  #   9  26
  #  10  0.791
  #  11  35
  #  12  37
  #  13  0.496
  # ...
  #  37  0.633
  #  38  -6.549
  #  39  -3.718
  chomp;
  my @rec = split(/:/);
  my $seqname = $rec[$SEQNAME1];

  # In a sequence name, if the entry is
  # in the format of a FASTA header
  # Take only the field after the first
  # "|".  for example.
  # gi|12345|gb|abcd gives "12345"

  if ($seqname =~ /\|/) {
    $seqname = (split(/\|/,$seqname))[1];
  }

  my $siteID = $rec[$SITEID1]; 
  my $seqlen = $rec[$SEQLEN1];
  my $Ghyb = $rec[$GHYB1];
  my $beg = $rec[$BEG1];
  my $end= $rec[$END1];

  #assert($beg > 0);
  #assert($end >= $beg);
  #assert($seqlen > 0);
  #assert($Ghyb <= 0.0);
  #assert($seqname ne "");

  my $key = join(":", $seqname, $seqlen, $Ghyb, $beg, $end);

  #assert(!defined $Sites{$key});  # enforce unique keys

  $Sites{$key}->{'SiteID'} = $siteID;

  splice(@rec, 0, 7); # originally 0,6 and now it is again . . . 
  my $value = join(":", @rec);
  $Sites{$key}->{'value'} = $value;
  $Sites{$key}->{'used'} = 0;
}
close IN1;

my $numline = 0;
<IN2>;
while (my $line = <IN2>) {
  #   0  NM_024451-3pUTR
  #   1  1288
  #   2  miR-125b-5p
  #   3  -26.4
  #   4  899-926
  #   5  923-926
  #   6  other
  #   7  -1
  #   8  1
  #   9  0
  #  10  0
  #  11  912-925
  #  12  910-911
  #  13  899-909
  chomp $line;
  my ($seqname, $seqlen, $Ghyb, $pos) = (split(/:/, $line))[$SEQNAME2,$SEQLEN2,$GHYB2,$POS2];

  # Sequence name in FASTA header format. see above
  if ($seqname =~ /\|/) {
    $seqname = (split(/\|/, $seqname))[1];
  }

  my ($beg, $end) = (split(/-/, $pos))[0,1];

  #assert($beg >= 1);  # non-positive position
  #assert($end >= $beg);
  #assert($seqlen > 0);
  #assert($Ghyb <= 0.0);
  #assert($seqname ne "");

  my $key = join(":", $seqname, $seqlen, $Ghyb, $beg, $end);

  #assert(defined $Sites{$key});
  #assert($Sites{$key}->{'used'} == 0);

  my @outrec = split(/:/, $line);
  splice(@outrec, 1, 1);
  my $outline = join(":", @outrec);
  print OUT $Sites{$key}->{'SiteID'} . ":" . $outline . ":" . $Sites{$key}->{'value'} . ":" .  $seqlen . "\n";
  $Sites{$key}->{'used'} = 1;
  ++$numline;
}

close IN2;
close OUT;

print "Matched $numline lines\n";

#warnings about unused
for my $key (keys %Sites) {
  if ($Sites{$key}->{'used'} == 0) {
    print "Warning: " . $Sites{$key}->{'value'} . " was not used.\n";
  }
}

# Output Header
# Outputs a header line.  This is needed by processing 
# That follows this script on the Website.
#    Arguments
#        -- output file handle
#    Returns
#        -- nothing
#
sub output_header 
{
    my $fh = shift;

    my @fields =
	(
	 "SiteID",       # Generated in Match_TotalEn.pl
	 "Transcript",	 # RNAhybrid
	 "miRNA",        # RNAhybrid
	 "Ghybrid",      # RNAhybrid
	 "Site_range",   # ResHyb- or TotalEn-
	 "Seed_range",   # Seed-
	 "Seed_kind",    # Seed-
	 "GU_pair",      # Seed-
	 "tA1",          # Seed-
	 "region12-17",  # Seed-
	 "region4-15",   # Seed-
	 "2-8_range",    # Seed-
	 "9-11_range",   # Seed-
	 "12+_range",    # Seed-
	 "SiteAcc",      # SiteFeatures-
	 "Old Seed_be",  # SiteFeatures-, not used
	 "Old Seed_en",  # SiteFeatures-, not used
	 "Old SeedAcc",  # SiteFeatures-, not used
	 "New Seed be",  # SiteFeatures-
	 "New Seed en",  # SiteFeatures-
	 "SeedAcc",      # SiteFeatures-
	 "UpsAcc_5nt",   # SiteFeatures-
	 "UpAU_5nt",     # SiteFeatures-
	 "DwsAcc_5nt",   # SiteFeatures-
	 "DwAU_5nt",     # SiteFeatures-
	 "UpsAcc_10nt",  # SiteFeatures-
	 "UpAU_10nt",    # SiteFeatures-
	 "DwsAcc_10nt",  # SiteFeatures-
	 "DwAU_10nt",    # SiteFeatures-
	 "UpsAcc_15nt",  # SiteFeatures-
	 "UpAU_15nt",    # SiteFeatures-
	 "DwsAcc_15nt",  # SiteFeatures-
	 "DwAU_15nt",    # SiteFeatures-
	 "UpsAcc_20nt",  # SiteFeatures-
	 "UpAU_20nt",    # SiteFeatures-
	 "DwsAcc_20nt",  # SiteFeatures-
	 "DwAU_20nt",    # SiteFeatures-
	 "UpsAcc_25nt",  # SiteFeatures-
	 "UpAU_25nt",    # SiteFeatures-
	 "DwsAcc_25nt",  # SiteFeatures-
	 "DwAU_25nt",    # SiteFeatures-
	 "UpsAcc_30nt",  # SiteFeatures-
	 "UpAU_30nt",    # SiteFeatures-
	 "DwsAcc_30nt",  # SiteFeatures-
	 "DwAU_30nt",    # SiteFeatures-
	 "Gtotal",       # TotalEn-
	 "Gnucl",        # TotalEn-
	 "UTR_length"    # RNAhybrid
	);

    my $tmpstr = join("\t", @fields);

    print $fh "$tmpstr\n";

} # end sub output_header

print STDERR "Exiting $0\n" if $debug;

exit 0;
