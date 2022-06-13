#!/usr/bin/perl

#
# Take GCG connect file as input and convert structure
# to base-pair format
#

# skip connect file headers
while (<>) {
  last if (/^\s*1\s+/);
}

print " Structure 1001\n";
do {
  @fields = split;
  if ($fields[4] > $fields[0] && $fields[4] != 0) {
    print " $fields[0]\t$fields[4]\n";
  }
} while (<>);
