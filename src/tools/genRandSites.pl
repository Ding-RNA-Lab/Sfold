#!/usr/bin/perl

#generate n fantasy binding sites each with start site p between a, inclusive and b, exclusive and end site between p+5 and p+25.
#genRandSites n a b

my $n = $ARGV[0];
my $a = $ARGV[1];
my $b = $ARGV[2];
srand;

for (my $i=0; $i<$n; $i++) {
	my $p = int(rand ($b - $a)) + $a;
	my $q = int(rand (20)) + 5 + $p;
	print "$p\t$q\n";
}
