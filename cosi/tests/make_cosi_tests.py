#!/usr/bin/env python

#
# Script: make_cosi_tests.py
#
# Generate the scripts to call cosi test cases. 
#

import sys, os, stat, logging

def SystemSucceed( cmd, dbg = False, exitCodesOk = ( 0, ) ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from ' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    exitCode = os.system( cmd )
    logging.info( 'Finished command ' + cmd + ' with exit code ' + str( exitCode ) )
    if exitCodesOk != 'any' and exitCode not in exitCodesOk:
        raise IOError( "Command %s failed with exit code %d" % ( cmd, exitCode ) )

logging.basicConfig( level = logging.DEBUG, format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

ntests = 11
include_xfail_tests = False
suffixes = ( '', )

xfail = [  ]
tests = [ i for i in range( 1, ntests+1 ) if i not in xfail ]
if not include_xfail_tests: xfail = []

for testNum in range( 1, ntests+1 ):
    for sfx in suffixes:
        scriptFN = 'tests/t%d/run_%d%s' % ( testNum, testNum, sfx )
        with open( scriptFN, 'w' ) as out:
            out.write( '''#!/bin/sh
set -e -v
TN=%d
TS=%s
if [ "$(expr substr $(uname -s) 1 6)" == "CYGWIN" -o "$(expr substr $(uname -s) 1 5)" == "MINGW" ]; then
  PS=win
else
  PS=
fi
TD=t${TN}test${TS}
rm -rf $TD
mkdir $TD
srcdirabs=$(cd $srcdir && pwd)
pf=$srcdirabs/tests/t${TN}/0_simple
cp $pf.cosiParams $pf.model $TD/
pushd $TD
$COSI_TEST_VALGRIND ../coalescent$TS -p 0_simple.cosiParams -o 0_simple_test
pwd
diff -q 0_simple_test.hap-1 ${pf}$TS$PS.hap-1 
diff -q 0_simple_test.hap-4 ${pf}$TS$PS.hap-4
diff -q 0_simple_test.hap-5 ${pf}$TS$PS.hap-5
diff -q 0_simple_test.pos-1 ${pf}$TS$PS.pos-1 
diff -q 0_simple_test.pos-4 ${pf}$TS$PS.pos-4
diff -q 0_simple_test.pos-5 ${pf}$TS$PS.pos-5
#if [ -x ../sample_stats_extra ]; then            
#   $COSI_TEST_VALGRIND ../coalescent$TS -p 0_simple.cosiParams -n 10 -m | ../sample_stats_extra -a 1,2,3-100,101-360 -l .200-.201 > sample_stats_out.txt
#   diff -q sample_stats_out.txt $srcdirabs/tests/t${TN}/sample_stats_out.txt         
#fi            
popd
rm -rf $TD
''' % ( testNum, sfx ) )
            os.fchmod( out.fileno(), stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO )

testData = []            
with open( 'cositests.am', 'w' ) as mf:
    is_first = [ True, True ]
    for sfx in suffixes:
        if sfx: mf.write( 'if ' + sfx[1:].upper() + '\n' )
        for in_xfail in ( False, True ):
            mf.write( ( '%sTESTS %s= ' % ( ( 'XFAIL_' if in_xfail else '' ),'' if is_first[ in_xfail ] else '+' ) ) + ' '.join([ 'tests/t%d/run_%d%s' % ( testNum, testNum,
                                                                                                                                                          sfx )
                                                                                                                                 for testNum in ( xfail if in_xfail else tests + xfail ) ]) + '\n' )
            is_first[ in_xfail ] = False
        if sfx: mf.write( 'endif\n' )

    for platformSuffix in ( '', 'win' ):
        testData += [ 'tests/t%d/0_simple%s.%s' % ( testNum, ( sfx + platformSuffix) if ext not in ( 'cosiParams', 'model' ) else '', ext )
                      for testNum in tests+xfail for sfx in suffixes for ext in ( 'cosiParams', 'model', 'hap-1', 'hap-4', 'hap-5', 'pos-1', 'pos-4', 'pos-5' ) ]

    testData = sorted( set( testData ) )
        
#    for f in testData: SystemSucceed( 'svn add --force ' + f )

    mf.write( 'COSI_TESTDATA = ' + ' '.join( testData ) + '\n' )
