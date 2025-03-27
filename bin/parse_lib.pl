# Perl routines to manipulate different sequence formats
# Clarence Chan @ 2003
#

use warnings;

#
# Sample usage:
#
#  # undefine \n as the delimiting character
#  undef $/;
#  my($input) = <>;
#
#  $seq = &parsePlain($input);
#  $seq = &parseFASTA($input);
#  $seq = &parseGenBank($input);
#  $seq = &parseCT($input);
#
#  $seq = &formatSeq($seq, 60) if ($seq ne "");
#  print $seq, "\n";
#

#
# extract sequence from a plain sequence format
#
sub parsePlain {
  my($seq) = @_;
  my($data);

  # remove header line if found, even though we say this is a plain seq...
  $seq =~ /^\s*(>[^\n]*\n)*(.*)$/s;
  $data = $2;
  $data =~ tr/atcgun/ATCGUN/;
  $data =~ s/[^ATCGUN]//gi;

  return $data;
}

#
# extract sequence from the FASTA format
#
sub parseFASTA {
  my($in) = @_;
  my($seq) = "";

  my(@lines) = split(/\n/, $in);

  if ($lines[0] =~ /^>/) {
    for ($i=1; $i<=$#lines; $i++) {
      $lines[$i] =~ tr/atcgun/ATCGUN/;
      $lines[$i] =~ s/[^ATCGUN]//gi;
      $seq .= $lines[$i];
    }
  }

  return $seq;
}

#
# extract sequence from the FASTA format
#
sub parseMultipleFASTANames {
  my $in = pop(@_);
  my @names = ();
  my @lines = split(/\n/, $in);
	for (@lines) {
		if (/^>/) {
			tr/;|`><" \n\r\b\f\\//d;	#transobliterate
			push @names, $_;			
		}
  	}
  return @names;
}

sub parseMultipleFASTASeqs {
  my $in = pop(@_);
  my $seq = "";
  my @seqs = ();

  my @lines = split(/\n/, $in);

	for (@lines) {
		if (/^>/) {
			#start a new line
			if ($seq ne "") { push @seqs, $seq; }
			$seq = "";			
		} else {
			#add line to seq so far
			tr/atcgun/ATCGUN/;
			tr/ATCGUN//cd;	#transobliterate compliment
			$seq .= $_;
		}
	}
	#one last push
	if ($seq ne "") { push @seqs, $seq; }
	
	return @seqs;
}

#
# extract sequence from the GenBank format
#
sub parseGenBank {
  my($in) = @_;
  my($seq) = "";

  my(@lines) = split(/\n/, $in);

  my($i) = 0;
  while ($i <= $#lines) {
    if ($lines[$i] =~ /^\s*(\d+)\s+(.*)/) {
      if ($1 eq "1") {
        $seq = "";
      }
      my($l) = $2;
      $l =~ s/[^ATCGUN]//gi;
      $seq .= $l;
    }
    $i++;
  }

  return $seq;
}

#
# extract sequence from the GCG Connect format
#
sub parseCT {
  my($in) = @_;
  my($seq) = "";

  my(@lines) = split(/\n/, $in);

  my($i) = 0;
  while ($i <= $#lines) {
    if ($lines[$i] =~ /^\s*(\d+)\s+(.*)/) {
      if ($1 eq "1") {
        $seq = "";
      }
      $seq .= substr($2, 0, 1);
    }
    $i++;
  }

  $seq =~ s/[^ATCGUN]//gi;

  return $seq;
}

#
# format seq $s with $x bases in a line
#
sub formatSeq {
  my($s, $x) = @_;
  $s =~ s/(.{$x})/$1\n/g;
  $s =~ tr/atcgun/ATCGUN/;

  return $s;
}

#
# reverse and complement the input sequence
#
sub revcomp {
  my($seq) = @_;

  $seq =~ tr/atcgun/ATCGUN/;
  $seq =~ s/[^ATCGUN]//g;
  my($seq_size) = length($seq);

  my($new_seq) = $seq;
  $new_seq =~ tr/ATCGUN/UAGCAN/;

  $seq = "";
  for ($i=0; $i<$seq_size; $i++) {
    $seq = substr($new_seq, $i, 1) . $seq;
  }

  return $seq;
}

1; # return true
