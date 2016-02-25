#!/usr/bin/env bash

#
# Script: timecmd.sh
#
# Measure the time, and optionally the memory usage, of a specified command.
# Redirect the command's output to specified files, and save the time/memory usage information
# to other specified files.
#
# Usage:
#
#    timecmd.sh OPTIONS command arg1 args2 ...
#
# where OPTIONS are
#
#   -o, --outfile FILENAME            send command's stdout to this file (defaults to /dev/stdout)
#   -e, --errfile FILENAME            send command's stderr to this file (defaults to /dev/stderr)
#   -t, --timefile FILENAME           send timing information to this file, including the hostname of the host
#                                     on which timing was done (defaults to /dev/null)
#   -T, --timeout timeout             time out the command after this time
#   -k, --timeoutCode  val            in case of timeout return this code, rathern than the default timeout code of 124
#


################################
# Section: Parse arguments
################################

set -e -o pipefail -o nounset

ARGS=$(getopt -o "o:e:t:T:k:" -l "outfile:,errfile:,timefile:,timeout:,timeoutCode:" -n "timecmd.sh" -- "$@");

eval set -- "$ARGS";

OUTFILE=/dev/stdout
ERRFILE=/dev/stderr
TIMEFILE=/dev/null
TIMEOUT=
TIMEOUTCODE=

while true; do
  case "$1" in
    -o|--outfile)
      shift;
      if [ -n "$1" ]; then
        OUTFILE=$1
        shift;
      fi
      ;;
    -e|--errfile)
      shift;
      if [ -n "$1" ]; then
        ERRFILE=$1
        shift;
      fi
      ;;
    -t|--timefile)
      shift;
      if [ -n "$1" ]; then
        TIMEFILE=$1
        shift;
      fi
      ;;
    -T|--timeout)
      shift;
      if [ -n "$1" ]; then
        TIMEOUT=$1
        shift;
      fi
      ;;
    -k|--timeoutCode)
      shift;
      if [ -n "$1" ]; then
        TIMEOUTCODE=$1
        shift;
      fi
      ;;
    --)
      shift;
      break;
      ;;
  esac
done

REST=$@

echoerr() { echo "$@" 1>&2; }
echoerr TIMEFILE is $TIMEFILE
echoerr ERRFILE is $ERRFILE
echoerr OUTFILE is $OUTFILE

#
# End section: Parse arguments
#

rm -f $TIMEFILE
# Save the host on which we ran: on a faster host the same amount of time indicates more computation
#echo $HOST > $TIMEFILE
touch $TIMEFILE

# If a timeout is specified, prepend the timeout command
TIMEOUTCMD=
if [[ -n "$TIMEOUT" ]]; then
		TIMEOUTCMD="timeout $TIMEOUT"
fi

set +e
{ time $TIMEOUTCMD $REST >$OUTFILE 2>$ERRFILE; } 2>>$TIMEFILE

ret=$?
set -e
echoerr FINISHED writing to $TIMEFILE
ls -l $TIMEFILE 1>&2
cat $TIMEFILE 1>&2
echo retcode $ret >> $TIMEFILE
echo -n 'host ' >> $TIMEFILE
hostname >> $TIMEFILE
echo -n 'uname ' >> $TIMEFILE
uname -a >> $TIMEFILE

# If timeout was specified, exit code 124 indicates that the user's command timed out rather than
# exited with its own error.  In that case, if requested, adjust the error code we return from the script.

if [[ -n "$TIMEOUT" && -n "$TIMEOUTCODE" && $ret -eq 124 ]];
then		
		ret=$TIMEOUTCODE
fi

exit $ret


