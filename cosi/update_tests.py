#!/usr/bin/env python

"""
 Script: update_tests.sh

 For each test directory specified on the command line, update the test.
"""

import argparse, os, sys, logging, stat, tempfile, string, types, collections, re

def updateCosiTests():

    parser = argparse.ArgumentParser( description = 'Update cosi test cases',
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter )

    parser.add_argument( '--min-model-num', type = int, help = 'minimum model id to check', required = True )
    parser.add_argument( '--max-model-num', type = int, help = 'maximum model id to check', required = True )
    parser.add_argument( '--cosi-sim-param-variants', action = DefaultAppend,
                         default = [ '', '-u .001', '-u .0001' ],
                         help = 'cosi simulation param variants' )
    parser.add_argument( '--update-stoch', action = 'store_true', help = 'update stochastic as well as exact test' )
    parser.add_argument( '--nsims-stoch', type = int, default = int( os.environ.get( 'COSI_TEST_NSIMS_STOCH', 10000 ) ) )
    
    args = parser.parse_args()
    logging.basicConfig( level = logging.DEBUG,
                         format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

    for modelNum in range( args.min_model_num, args.max_model_num+1 ):
        modelDir = os.path.join( 'tests', 'model%03d' % modelNum )
        for simParamVariant in args.cosi_sim_param_variants:
            testName = 'model%03d_t%s' % ( modelNum, Sfx( simParamVariant ) )
            testDir = os.path.join( modelDir, AddFileSfx( 't', simParamVariant ) )
            EnsureDirExists( testDir )

            nsims_stoch = args.nsims_stoch
            updateStochEnv = ''
            if args.update_stoch: updateStochEnv = 'COSI_TEST_UPDATE_STOCH=1'
            SystemSucceed( "bsub -q hour -sp 2 -W '04:00' -P sabeti_cosi -G sabeti_labfolk -E 'pushd /idi/sabeti-scratch/ilya/gsvn/Temp && python ../Operations/Ilya_Operations/PipeRun/python/testpython.py && popd' -R 'rusage[argon_io=1,mem=5]' -J 'cositest_%(testName)s' -o %(testDir)s/update.log1  'env COSI_TEST_UPDATE_EXACT=1 %(updateStochEnv)s COSI_TEST_NSIMS_STOCH=%(nsims_stoch)d srcdir=../../../cosi ../../../cosi/%(testDir)s/run.sh >& %(testDir)s/update.out1'" % locals() )

def EnsureDirExists( d ):
    """Create the given dir (and any parent dirs) if it does not already exist"""
    if not os.path.isdir( d ):
        os.makedirs( d )
        assert os.path.isdir( d )

def ReserveTmpFileName( prefix = 'mytmp', suffix = '.tmp', text = True, tmpDir = None, executable = False ):
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

def SystemSucceed( cmd, dbg = False, exitCodesOk = ( 0, ) ):
    """Run a shell command, and raise a fatal error if the command fails."""
    logging.info( 'Running command ' + cmd + ' ; called from ' + sys._getframe(1).f_code.co_filename + ':' +
                  str( sys._getframe(1).f_lineno ) )
    scriptFN = None
    try:
        scriptFN = ReserveTmpFileName( executable = True )
        with open( scriptFN, 'w' ) as out:
            out.write( '#!/usr/bin/env bash\n' )
            out.write( 'set -e -o pipefail -o nounset\n' )
            out.write( cmd + '\n' )
            os.fchmod( out.fileno(), stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO )
        
        exitCode = os.system( scriptFN )
        logging.info( 'Finished command ' + cmd + ' with exit code ' + str( exitCode ) )
        if exitCodesOk != 'any' and exitCode not in exitCodesOk:
            raise IOError( "Command %s failed with exit code %d" % ( cmd, exitCode ) )
    finally:
        if scriptFN and os.path.isfile( scriptFN ) and 'COSI_KEEP_TMP' not in os.environ:
            os.remove( scriptFN )

class DefaultAppend(argparse._AppendAction):
    def __call__(self, parser, namespace, values, option_string=None):
        items = argparse._copy.copy(argparse._ensure_value(namespace, 
                                                           self.dest, []))
        try:
            self._not_first
        except AttributeError:
            self._not_first = True
            del items[:]
        items.append(values)
        setattr(namespace, self.dest, items)

def MakeAlphaNum( str ):
    """Return a version of the argument string, in which all non-alphanumeric chars have been replaced
    by underscores.
    """
    return re.sub( '\W+', '_', str )

def AddFileSfx( fullFileName, *sfx ):
    """Add the specified suffix to a filename, inserting it before the file extension.
    So, if the file was named ../Data/myfile.tsv and the suffix is 'verA',
    this returns '../Data/myfile_verA.tsv' .
    """

    fileName, fileExt = os.path.splitext( fullFileName if not fullFileName.endswith( '/' ) else fullFileName[:-1] )
    return ( ( fileName + Sfx( sfx ) ) if fileExt not in ( '.gz', 'bz2' ) else \
               AddFileSfx( fileName, *sfx ) ) + fileExt + ( '' if not fullFileName.endswith( '/' ) else '/' )

def Sfx( *vals ):
    """Return a suffix string based on the value: if the value is non-empty, return '_' + val; but if it is
    already empty, then just return the empty string.
    """

    def underscorify( val ):
        """prepend undescrore if needed"""
        if not val and val != 0: return ''
        noPrefix = isinstance( val, types.StringTypes) and val.startswith('#')
        alphaVal = MakeAlphaNum( str( val ) )
        return alphaVal[1:] if noPrefix else ( alphaVal if alphaVal.startswith( '_' ) else '_' + alphaVal )
    
    return ''.join( map( underscorify, flatten( vals ) ) )

def flatten(*args):
    """flatten(sequence) -> tuple

    Returns a single, flat tuple which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, tuple((8,9,10))])
    (1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10)

    Taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    """

    x = args[0] if len( args ) == 1 else args

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if IsSeq( el ): result.extend(flatten(el))
        else: result.append(el)
    return tuple( result )

class AtomicForIsSeq(object):
    """Instances of classes derived from this class will be called _not_ sequences
    by IsSeq(), even if they otherwise look like a sequence."""
    pass

def IsSeq( val ):
    """Test if the value is a sequence but not a string or a dict.
    Derive your class from AtomicForIsSeq to force its instances to be called
    _not_ sequences by this function, regardless of anything else.
    """
    return ( isinstance( val, ( collections.Sequence, types.GeneratorType ) ) or
             ( hasattr( val, '__getitem__' ) and hasattr( val, '__len__' ) ) )  \
             and not isinstance( val, ( types.StringTypes, collections.Mapping, AtomicForIsSeq ) )


if __name__ == '__main__':
    updateCosiTests()
