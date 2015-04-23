#
# Script: update_cosi_tests.py
#
# Update cosi regtest cases to reflect the current output of the code.
#

import sys, os, stat, logging, platform

def SystemSucceed( cmd, dbg = False, exitCodesOk = ( 0, ) ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from ' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    exitCode = os.system( cmd )
    logging.info( 'Finished command ' + cmd + ' with exit code ' + str( exitCode ) )
    if exitCodesOk != 'any' and exitCode not in exitCodesOk:
        raise IOError( "Command %s failed with exit code %d" % ( cmd, exitCode ) )

logging.basicConfig( level = logging.DEBUG, format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

min_test = 11
max_test = 11
cosi_suffixes = ( '', )
platform_suffix = ''
if platform.platform().startswith( 'CYGWIN' ): platform_suffix = 'win'

for testNum in range( min_test, max_test+1 ):
    for cosiSfx in cosi_suffixes:
        SystemSucceed( 'pushd tests/t%d && ../../../builds/opt/cosi/coalescent%s -p 0_simple.cosiParams -o 0_simple%s && popd' % ( testNum, cosiSfx, platform_suffix ) )
        SystemSucceed( 'pushd tests/t%d && ( ../../../builds/opt/cosi/coalescent%s -p 0_simple.cosiParams -n 10 -m -t | ../../../builds/opt/cosi/sample_stats_extra -a 1,2,3-100,101-360 -l .200-.201 > sample_stats_out%s.txt )  && popd' % ( testNum, cosiSfx, platform_suffix ) )
