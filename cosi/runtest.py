#!/usr/bin/env python

"""Create or check a cosi test case."""

from __future__ import division
import os, sys, argparse, logging, subprocess, random, string, tempfile, inspect, re, stat, contextlib, fcntl, errno, time, \
    platform, shutil, socket, subprocess

try:
    import lsfutil
    have_lsfutil = True
except ImportError:
    have_lsfutil = False

logging.basicConfig( level = logging.DEBUG, format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

def parseArgs():    

    parser = argparse.ArgumentParser( description = 'Create or check a cosi test case',
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter )

    parser.add_argument( '--srcdir', default = os.environ.get( 'srcdir', '../../../cosi' ),
                         help = 'directory containing the source tree, under which the test cases are located' )
    parser.add_argument( '--builddir', default = os.environ.get( 'builddir', '.' ),
                         help = 'directory where we run the tests' )
    parser.add_argument( '--test-num', help = 'id of the test', type = int )
    parser.add_argument( '--test-name', help = 'test name (determines test dir)' )
    parser.add_argument( '--test-dir', help = 'the test dir' )
    parser.add_argument( '--update-exact', action = 'store_true',
                         default = bool( int( os.environ.get( 'COSI_TEST_UPDATE_EXACT', 0 ) ) ),
                         help = 'update the exact test case, rather than checking it' )
    parser.add_argument( '--update-stoch', action = 'store_true',
                         default = bool( int( os.environ.get( 'COSI_TEST_UPDATE_STOCH', 0 ) ) ),
                         help = 'update the stochastic test case, rather than checking it' )
    parser.add_argument( '--variant-exact',
                         help = 'variant string for the exact test case',
                         default = '_'.join((SysTypeString(),
                                             os.environ.get( 'CXX', 'g++' ))) )
    parser.add_argument( '--variant-stoch',
                         help = 'variant string for the stochastic test case',
                         default = 'dflt' )
    parser.add_argument( '--seed', help = 'random seed for cosi (0 to use time)', type = int, default = 0 )
    parser.add_argument( '--nsims-exact', help = 'number of simulations for the exact check', type = int,
                         default = int( os.environ.get( 'COSI_TEST_NSIMS_EXACT', 3 ) ) ) 
    parser.add_argument( '--nsims-stoch', help = 'number of simulations for the stochastic check', type = int,
                         default = int( os.environ.get( 'COSI_TEST_NSIMS_STOCH', 1000 ) ) )
    parser.add_argument( '--force-stoch', help = 'run the stochastic check even if exact check succeeds',
                         action = 'store_true',
                         default = bool( int( os.environ.get( 'COSI_TEST_FORCE_STOCH', 0 ) ) ) )
    parser.add_argument( '--cosi-sim-params', help = 'additional args to pass to cosi', default = '' )
    parser.add_argument( '--max-minutes', type = float, default = float( os.environ.get( 'COSI_TEST_MAX_MINUTES', 0 ) ),
                         help = 'if >0, stop sims after this many minutes' )
    parser.add_argument( '--max-slowdown',
                         type = float,
                         help = 'max slowdown factor compared to reference; test fails if slowdown is more than this',
                         default = float( os.environ.get( 'COSI_TEST_MAX_SLOWDOWN', 2 ) ) )
    parser.add_argument( '--use-orig-seed', action = 'store_true',
                         default = bool( int( os.environ.get( 'COSI_TEST_USE_ORIG_SEED', 0 ) ) ), 
                         help = 'for stochastic test running, use original random seed used to generate the test' )
    parser.add_argument( '--use-lsf', action = 'store_true', help = 'dispatch stochastic tests to LSF',
                         default = bool( int( os.environ.get( 'COSI_TEST_USE_LSF', 0 ) ) ) )

    miscArgs = parser.add_argument_group( 'misc', 'miscellaneous args' )
    miscArgs.add_argument( '--cosi-binary', help =  'path to cosi binary',
                           default = os.environ.get( 'COSI_TEST_COSI_BINARY', './coalescent' ) )
    miscArgs.add_argument( '--stats-binary', help =  'path to sample_stats_extra binary',
                           default = os.environ.get( 'COSI_TEST_STATS_BINARY', './sample_stats_extra' ) )

    args = parser.parse_args()
    logging.info( 'args=' + str( args ) + ' env=' + str( os.environ ) )

    return args

def SystemSucceed( cmd, dbg = False, useLSF = False, lsfJobName = 'runtest',
                   lsfPreExecCmd = '' ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from ' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    scriptFN = None
    outFN = None
    try:
        # scriptFN = ReserveTmpFileName( executable = True, prefix = 'runtestsysscrpt' )
        # with open( scriptFN, 'w' ) as out:
        #     out.write( '#!/usr/bin/env bash\n' )
        #     out.write( 'set -e -o pipefail -o nounset\n' )
        #     out.write( cmd + '\n' )
        #     os.fchmod( out.fileno(), stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO )

        if useLSF:
            outFN = ReserveTmpFileName( prefix = 'runtestlsfoutput' )
            grpId = '/ilya/cosi/runtest/' + socket.getfqdn() + ':' + str( os.getpid() )
            preExecOpts = '' if not lsfPreExecCmd else " -E '" + lsfPreExecCmd + "'"
            cmd = "bsub -K -q hour -W '04:00' -P sabeti_cosi -g '%(grpId)s' -J '%(lsfJobName)s' %(preExecOpts)s -R 'rusage[argon_io=1,mem=3]' -o %(outFN)s '%(cmd)s'" % locals()
            logging.info( 'outFN=' + outFN + ' scriptCmd=' + cmd )
        subprocess.check_call( cmd, shell = True )
        logging.info( 'Finished command ' + cmd + ' with exit code 0' )
    finally:
        if outFN:
            logging.info( 'printing output file ' + outFN )
            os.system( 'cat ' + outFN )
        if 'COSI_KEEP_TMP' not in os.environ:
            for f in ( scriptFN, outFN ):
                if f and os.path.isfile( f ) and os.access( f, os.W_OK ):
                    try:
                        logging.info( 'trying to delete tempfile ' + f )
                        os.remove( f )
                        logging.info( 'deleted tempfile ' + f )
                    except IOError as e:
                        logging.warning( "Error deleting temporary file; I/O error({0}): {1}".format(e.errno, e.strerror) )

    logging.info( 'exiting SystemSucceed for command ' + cmd )

def DumpFile( fname, contents, dbg = False, getio = None ):
    """Write out a file whose entire contents is the given string"""
    if getio: return dict( depends_on = (), creates = fname, attrs = dict( piperun_short = True ) )
    logging.info( 'Writing ' + fname  )
    with open( fname, 'wt' ) as f:
        f.write( str( contents ) )
    logging.info( 'Wrote ' + fname  )

def SlurpFile( fname ):
    """Read entire file into one string"""
    if not os.path.isfile( fname ): raise IOError( 'File not found: %s' % fname )
    with open( fname ) as f:
        return f.read()

def joinstr( *args ): return ' '.join( map( str, [_f for _f in args if _f] ) )

def subst( s, localEnv ): return string.Template( s ).substitute( localEnv )

def EnsureDirExists( d ):
    """Create the given dir (and any parent dirs) if it does not already exist"""
    if not os.path.isdir( d ):
        os.makedirs( d )
        assert os.path.isdir( d )

def readTime( timeFN ):
    with open( timeFN ) as inF:
        #host = inF.readline().strip()
        #cpuFactor = lsfutil.getHostCPUFactor( host ) if have_lsfutil else 1.0
        #inF.readline()
        cpuFactor = 1.0

        userTime = None
        sysTime = None
        
        for line in inF:
            if line.startswith( 'user' ) or line.startswith( 'sys' ):
                r = re.match( r'(?P<mins>\d+\.?\d*)m(?P<secs>\d+\.?\d*)s', line.strip().split()[1] )
                assert r
                t = float( r.group( 'mins' ) ) * 60.0 + float( r.group( 'secs' ) )
                if line.startswith( 'user' ):
                    userTime = t
                else:
                    sysTime = t
                if userTime is not None and sysTime is not None: break

        return userTime + sysTime
    
def ReserveTmpFileName( prefix = 'mytmp', suffix = '.tmp', text = True,
                        tmpDir = os.environ.get( '__LSF_JOB_TMPDIR__', '/tmp' ),
                        executable = False ):
    """Create a unique temp file, close it, and return its name.
    The caller then would typically overwrite this file,
    being assured that no other ReserveTmpFileName() call can return
    the same filename to another process or thread.
    The file is automatically deleted at program exit.
    """
    if tmpDir: EnsureDirExists( tmpDir )
    fileDescr, fileName = tempfile.mkstemp( prefix = prefix, suffix = suffix, dir = tmpDir, text = text )
    if executable:
        os.fchmod( fileDescr, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO )
        
    os.close( fileDescr )
    return fileName


class SimpleFlock(object):
   """Provides the simplest possible interface to flock-based file locking. Intended for use with the `with` syntax. It will create/truncate/delete the lock file as necessary.

   Adapted from https://github.com/derpston/python-simpleflock.git
   """

   def __init__(self, path, exclusive = True, timeout = None, minCheckInterval = 0.1, maxCheckInterval = 10, disable = False ):
      self._path = path
      self._exclusive = exclusive
      self._timeout = timeout
      self._fd = None
      self._minCheckInterval = minCheckInterval
      self._maxCheckInterval = maxCheckInterval
      self._rand = random.Random()
      self._rand.seed()
      self._disable = disable

   def __enter__(self):
      if self._disable: return 
      if not self._exclusive and not os.path.isfile( self._path ):
          raise IOError( "Lockfile for shared lock not found: " + self._path )
      self._fd = os.open(self._path, ( ( os.O_CREAT | os.O_EXCL | os.O_WRONLY ) if self._exclusive else os.O_RDONLY ) )
      start_lock_search = time.time()
      checkInterval = self._minCheckInterval
      while True:
         try:
            fcntl.lockf(self._fd, ( fcntl.LOCK_EX if self._exclusive else fcntl.LOCK_SH ) | fcntl.LOCK_NB)
            # Lock acquired!
            logging.info( 'acquired ' + ( 'exclusive' if self._exclusive else 'shared' ) + ' lock ' + self._path )
            return
         except IOError as ex:
            if ex.errno not in ( errno.EAGAIN, errno.EACCES ): # error other than "lock already held"
               raise
            elif self._timeout is not None and time.time() > (start_lock_search + self._timeout):
               # Exceeded the user-specified timeout.
               raise
         
         # TODO It would be nice to avoid an arbitrary sleep here, but spinning
         # without a delay is also undesirable.
         if checkInterval < self._maxCheckInterval: checkInterval *= 2
         time.sleep(checkInterval + self._rand.random() )

   # end: def __enter__() 

   def __exit__(self, *args):
      if self._disable: return 
      fcntl.lockf(self._fd, fcntl.LOCK_UN)
      os.close(self._fd)
      self._fd = None
      logging.info( 'released ' + ( 'exclusive' if self._exclusive else 'shared' ) +' lock ' + self._path )

      # Try to remove the lock file, but don't try too hard because it is
      # unnecessary. This is mostly to help the user see whether a lock
      # exists by examining the filesystem.
      # try:
      #    os.unlink(self._path)
      # except:
      #    pass

# end: class SimpleFlock

@contextlib.contextmanager
def EmptyContextMgr():
    """A context manager that does nothing"""
    yield

def SysTypeString( ):
    """Returns a string identifying the processor type and OS type.
    Useful for choosing among platform-specific binaries to execute.
    """

    if sys.platform == 'darwin':
        machType, osType = 'i386', 'darwin'
    else:    
        machType = platform.uname()[4]
        if machType.startswith('x86_64'): machType = 'x86_64'
        osType = platform.system().lower()
        if osType.startswith( 'darwin' ): osType = 'darwin'

    # linuxType = ''
    # if platform.system() == 'Linux':
    #   linuxType += '-' + MakeAlphaNum( platform.linux_distribution()[0].strip() )
        
    return machType + '-' + osType


def runTest( args ):
    """Create/update or run a cosi test case"""

    if not args.test_name: args.test_name = 't%03d' % args.test_num

    srcdir = args.srcdir
    testDir = args.test_dir or os.path.join( args.srcdir, 'tests', 'dist', args.test_name )
    exactDir = os.path.join( testDir, 'exact', args.variant_exact )
    stochDir = os.path.join( testDir, 'stoch', args.variant_stoch )

    updating = args.update_exact or args.update_stoch
    if updating:
        EnsureDirExists( testDir )
        if args.update_exact: EnsureDirExists( exactDir )
        if args.update_stoch: EnsureDirExists( stochDir )

    logging.info( 'getting lock' )
    with SimpleFlock( os.path.join( testDir, 'updating.lck' ), timeout = 0, exclusive = updating,
                      disable = not updating ):
        logging.info( 'got lock' )

        paramFN = os.path.join( testDir, '..', 'test.cosiParams' )
        genMapFN = os.path.join( testDir, '..', 'test.genmap' )

        totSampleSize = 0
        with open( paramFN ) as p:
            for line in p:
                if line.startswith( 'sample_size ' ):
                    totSampleSize += int( line.strip().split()[2] )

        if not args.seed:
            r = random.Random()
            # r.jumpahead( args.test_num )
            args.seed = r.randrange( 37, sys.maxsize )

        logging.info( 'runTest: seed is ' + str( args.seed ) )

        cosiBinary = args.cosi_binary
        cosiCmd = '$cosiBinary -p $paramFN -R $genMapFN -m ' + args.cosi_sim_params

        if args.nsims_exact > 0:
            cosiCmdFN = os.path.join( exactDir, 'exactcmd.txt' )
            if args.update_exact:
                # save the result of the exact check
                shutil.copy( cosiBinary, os.path.join( exactDir, 'coalescent' ) )
                shutil.copy( cosiBinary + '.flags.txt', os.path.join( exactDir, 'coalescent.flags.txt' ) )
                cosiCmdExact = joinstr( cosiCmd, '-n', args.nsims_exact, '--seed', args.seed )
                DumpFile( cosiCmdFN, cosiCmdExact )
            else:
                cosiCmdExact = SlurpFile( cosiCmdFN )

            sumFN = os.path.join( exactDir, 'exactsum.sha512' )
            SystemSucceed( joinstr( subst( cosiCmdExact, locals() ), ' | sha512sum', '>' if args.update_exact else '-c', sumFN ),
                           lsfPreExecCmd = 'ls -l %(srcdir)s/runtest.py' % locals() )

        # * Stochastic check
        if ( args.force_stoch or args.update_stoch ) and args.nsims_stoch > 0:
            afsMinBinSize = 4
            afsBinCount = 20
            afsBinSize = max( afsMinBinSize, int( totSampleSize / afsBinCount ) )
            statsBinary = args.stats_binary
            statsCmd = joinstr( '$statsBinary -a',
                                ','.join( map( str, [1,2,3,4,5] +
                                               [ '%d-%d' % ( f, min( f+afsBinSize-1, totSampleSize ) )
                                                 for f in range( 6, totSampleSize, afsBinSize )  ] ) ),
                                '--ld-seps 5,50,100,200,300,500,1000,2000,3000,5000,10000 -g 10' )

            cmpDistScript = os.path.join( args.srcdir, 'cmpdist.py' )
            stochSummaryFN = os.path.join( stochDir, 'stochsumm.tsv.bz2' )
            nsimsStoch = args.nsims_stoch

            stochCmdFN = os.path.join( stochDir, 'stochcmd.txt' )
            stochTimeFN = os.path.join( stochDir, 'stochtime.txt' )
            stochTimeFN_orig = None
            stochCountFN = os.path.join( stochDir, 'stochcount.txt' )
            timeScript = os.path.join( args.srcdir, 'timecmd.sh' )
            try:
                if args.update_stoch:
                    shutil.copy( cosiBinary, os.path.join( stochDir, 'coalescent' ) )
                    shutil.copy( cosiBinary + '.flags.txt', os.path.join( stochDir, 'coalescent.flags.txt' ) )
                    shutil.copy( statsBinary, os.path.join( stochDir, 'sample_stats_extra' ) )
                    shutil.copy( statsBinary, os.path.join( stochDir, 'sample_stats_extra.flags.txt' ) )

                    cosiCmdStoch = joinstr( '$timeScript -t $stochTimeFN', '--',
                                            cosiCmd, '--output-sim-times --output-end-gens -n $nsimsStoch', '--seed', args.seed,
                                            '|', statsCmd,
                                            '| $cmpDistScript --file1 stdin.tsv --file2 $stochSummaryFN --exclude-cols time' )
                    DumpFile( stochCmdFN, cosiCmdStoch )
                    DumpFile( stochCountFN, nsimsStoch )
                    cosiCmdStoch = joinstr( cosiCmdStoch, '--record' )
                else:
                    cosiCmdStoch = SlurpFile( stochCmdFN )
                    stochTimeFN_orig = stochTimeFN
                    stochTimeFN = ReserveTmpFileName( prefix = 'runteststochtime' )
                    print( 'stochTimeFN', stochTimeFN )
                    if args.max_minutes:
                        cosiCmdStoch = str.replace( cosiCmdStoch, '-m ', '-m --stop-after-minutes %f ' % args.max_minutes )
                    if not args.use_orig_seed:
                        cosiCmdStoch = re.sub( '--seed \d+', '--seed ' + str( args.seed ), cosiCmdStoch, count = 1 )

                SystemSucceed( subst( cosiCmdStoch, locals() ), useLSF = args.use_lsf, lsfJobName = testDir,
                               lsfPreExecCmd = 'ls -l %(srcdir)s/runtest.py' % locals() )

                if not args.update_stoch:
                    if os.path.isfile( os.path.join( stochDir, 'stochcount.long.txt' ) ):
                        stochCountFN = os.path.join( stochDir, 'stochcount.long.txt' )
                    origTime = readTime( stochTimeFN_orig ) / float( SlurpFile( stochCountFN ) )
                    curTime = readTime( stochTimeFN ) / float( nsimsStoch )
                    logging.info( 'origTime=' + str( origTime ) + ' curTime=' + str( curTime ) + ' slowdown=' +
                                  str( curTime / origTime ) )
                    assert  curTime <=  origTime * args.max_slowdown

                if updating:
                    checksumFN = '%(testDir)s/testchecksums.sha512' % locals()
                    if os.path.isfile( checksumFN ): os.remove( checksumFN )
                    SystemSucceed( "sha512sum `find %(testDir)s -type f -not -name testchecksums.sha512 -and -not -name '*~' | sort` > %(checksumFN)s" % locals(),
                                   lsfPreExecCmd = 'ls -l %(srcdir)s/runtest.py' % locals() )

            finally:
                pass
#                if not args.update_stoch and stochTimeFN_orig and os.path.isfile( stochTimeFN ) \
#                   and 'COSI_KEEP_TMP' not in os.environ:
#                    os.remove( stochTimeFN )

    # end: with lock
# end: def runTest

if __name__ == '__main__':
    runTest( parseArgs() )
