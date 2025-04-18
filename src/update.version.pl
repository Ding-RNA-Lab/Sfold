#!/usr/bin/perl

#
# This script is for updating the version number
# information contained in the header file
# ./version.h .
#

use Cwd qw(abs_path);

$scriptdir = $0;
$scriptdir =~ s/\/[^\/]*$/\//;
$CDIR = &abs_path($scriptdir);

$hdrfile = "$CDIR/version.h";

$cur_date = `date +\%Y\%m\%d`;
chomp($cur_date);

$line = `grep COMPILE_DATE $hdrfile`;
chomp($line);
$line =~ /\#define\s+COMPILE_DATE\s+(.*)/;
$last_date = $1;
$last_date =~ s/^\"//;
$last_date =~ s/\"$//;

if ($cur_date eq $last_date) {
  print "Version info is up-to-date. No change to version.h is necessary.\n";

} else {
  my($lines) = "";

  open(VERSION_FILE, $hdrfile) || die " Error: unable to locate file '$hdrfile'";
  while (<VERSION_FILE>) {
    if (/^\s*#define\s+COMPILE_DATE\s+.+$/) {
      $lines .= "#define COMPILE_DATE    \"$cur_date\"\n";
    } else {
      $lines .= $_;
    }
  }
  close(VERSION_FILE);

  open(VERSION_FILE, ">$hdrfile") || die " Error: unable to write to file '$hdrfile'";
  print VERSION_FILE $lines;
  close(VERSION_FILE);

  print "version.h has been updated. Please recompile the package.\n";
}
