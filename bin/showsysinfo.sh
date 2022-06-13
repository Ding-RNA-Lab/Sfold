#!/bin/sh

VAR=$1

# function to convert all characters in a string to uppercase
toupper() {
  echo "$@" | tr '[a-z]' '[A-Z]'
}

# check operating system
os_str=`uname -s`
os_str=`toupper $os_str`
# (August 2008 Adam) removed SUNOS and OSF1 as possibilities
# (July 2009 Adam) restored SUNOS  as possibility
if test "$os_str" = "LINUX" || test "$os_str" = "DARWIN" || test "$os_str" = "SUNOS"; then
  OSTYPE=$os_str
else
  OSTYPE="UNKNOWN"
fi

# check hardware type
# (July 2008 Adam) changed this from uname -p to uname -m
# (July 2009 Adam) added s/I86PC/I386/
hw_str=`uname -m`
hw_str=`toupper $hw_str | sed 's/I.86/I386/' | sed 's/I86PC/I386/'`
# (August 2008 Adam) removed SPARC, ALPHA, and POWERPC as possibilities
if test "$hw_str" = "X86_64" || test "$hw_str" = "I386"; then
  MACH=$hw_str
else
  MACH="UNKNOWN"
fi

#(July 2009 Adam) bioqq4 only! force X86_64 despite uname
if test "$os_str" = "SUNOS"; then
  MACH="X86_64"
fi

if test "$OSTYPE" = "UNKNOWN" || test "$MACH" = "UNKNOWN"; then
  OSTYPE="UNKNOWN"
  MACH="UNKNOWN"
fi

if test "$VAR" = "OSTYPE"; then
  echo $OSTYPE
elif test "$VAR" = "MACH"; then
  echo $MACH
else
  echo $OSTYPE
  echo $MACH
fi
