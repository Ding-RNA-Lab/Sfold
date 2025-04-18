#!/usr/bin/perl

#
# This script is for updating the system information
# definitions inside ./sysinfo.h . This script
# should be called before every compilation to ensure
# the right platform parameters are used
#

#ctrl_c() {
#  echo "Compilation aborted by user."
#}
#trap  ctrl_c   2

use Cwd qw(abs_path);

$scriptdir = $0;
$scriptdir =~ s/\/[^\/]*$/\//;
$CDIR = &abs_path($scriptdir);
$showsysinfo = join("", $CDIR, "/../bin/showsysinfo.sh");

# check current operating system
$SYSINFO_OSTYPE = `$showsysinfo OSTYPE`;
chomp($SYSINFO_OSTYPE);

# check current hardware type
$SYSINFO_MACH = `$showsysinfo MACH`;
chomp($SYSINFO_MACH);

open(SYSINFO_FILE, "$CDIR/sysinfo.h") || die "Error: unable to locate file '$CDIR/sysinfo.h'";
while (<SYSINFO_FILE>) {
  if (/^\s*#define\s+OSTYPE\s+(\w+)\s*$/) {
    $SYSINFO_OSTYPE_ORI = $1;
  }
  if (/^\s*#define\s+MACH\s+(\w+)\s*$/) {
    $SYSINFO_MACH_ORI = $1;
  }
}
close(SYSINFO_FILE);

die "Error: '$CDIR/sysinfo.h' is not in the right format"
  if ($SYSINFO_OSTYPE_ORI eq "" || $SYSINFO_MACH_ORI eq "");


$UPDATED = 0;
if ($SYSINFO_OSTYPE_ORI ne $SYSINFO_OSTYPE) {
  my($lines) = "";

  open(SYSINFO_FILE, "$CDIR/sysinfo.h") || die "Error: unable to locate file '$CDIR/sysinfo.h'";
  while (<SYSINFO_FILE>) {
    if (/^\s*#define\s+OSTYPE\s+$SYSINFO_OSTYPE_ORI\s*$/) {
      $lines .= "#define OSTYPE\t\t$SYSINFO_OSTYPE\n";
    } else {
      $lines .= $_;
    }
  }
  close(SYSINFO_FILE);

  open(SYSINFO_FILE, ">$CDIR/sysinfo.h") || die "Error: unable to write to file '$CDIR/sysinfo.h'";
  print SYSINFO_FILE $lines;
  close(SYSINFO_FILE);

  $UPDATED = 1;
}

if ($SYSINFO_MACH_ORI ne $SYSINFO_MACH) {
  my($lines) = "";

  open(SYSINFO_FILE, "$CDIR/sysinfo.h") || die "Error: unable to locate file '$CDIR/sysinfo.h'";
  while (<SYSINFO_FILE>) {
    if (/^\s*#define\s+MACH\s+$SYSINFO_MACH_ORI\s*$/) {
      $lines .= "#define MACH\t\t$SYSINFO_MACH\n";
    } else {
      $lines .= $_;
    }
  }
  close(SYSINFO_FILE);

  open(SYSINFO_FILE, ">$CDIR/sysinfo.h") || die "Error: unable to write to file '$CDIR/sysinfo.h'";
  print SYSINFO_FILE $lines;
  close(SYSINFO_FILE);

  $UPDATED = 1;
}

if ($UPDATED == 1) {
  print "sysinfo.h has been updated. Please recompile the package.\n";
} else {
  print "System info is up-to-date. No change to sysinfo.h is necessary.\n";
}

