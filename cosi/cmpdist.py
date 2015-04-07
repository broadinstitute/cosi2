#!/usr/bin/env python

"""Test whether the distribution of any of a set of summary statistics differs between two 
sets of simulations."""


import sys, argparse, logging, bz2, os, re
import numpy as np
from scipy.stats import ks_2samp

parser = argparse.ArgumentParser( description = 'Test for differences in distribution of summary statistics',
                                  formatter_class = argparse.ArgumentDefaultsHelpFormatter )

parser.add_argument( '--file1', help = 'first table of summary stats', required = True )
parser.add_argument( '--file2', help = 'second table of summary stats', required = True )
parser.add_argument( '-p', type = float, default = float( os.environ.get( 'COSI_TEST_P', .0001 ) ),
                     help = "maximum p level at which difference between distributions is tolerated"  )
parser.add_argument( '--record', action = 'store_true', help = 'save the stats to an .npz file and quit' )
parser.add_argument( '--exclude-cols', help = 'exclude columns whose names match this regexp', action = 'append' )
parser.add_argument( '-v', action =  'store_true', help = 'increase verbosity' )
parser.add_argument( '-c', action = 'store_true', help = 'check that this script runs' )

args = parser.parse_args()

if args.c:
    print('cmpdist.py script for comparing distributions starts up ok')
    sys.exit( 0 )

logging.basicConfig( level = logging.DEBUG, format='%(process)d %(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

def loadTable( f ):
    """Load a numpy recarray from file ``f``.  If f is a string that starts with ``stdin``, then stdin is used.
    f will be automatically uncompressed if it is compressed with .bz2 . tsv and npz files can be loaded.
    """
#    logging.info( 'type(f)=' + str( type( f ) ) )
#    logging.info( 'f=' + str( f ) )
    f_orig = f
    name, ext = os.path.splitext( f )
    try:
        if f_orig.startswith( 'stdin.' ): f = sys.stdin.buffer
        if ext == '.bz2':
            f = bz2.BZ2File( f )
            name, ext = os.path.splitext( name )
        if ext == '.tsv': return np.genfromtxt( f, names = True, delimiter = '\t', comments = '#', invalid_raise = True )
        elif ext == '.npz':
            result = np.load( f )
            assert len( result ) == 1
            return list( result.items() )[0][1]
        else:
            raise RuntimeError( 'unknown extension for input file %s: "%s"' % ( f_orig, ext ) )
    finally:
        if not f_orig.startswith( 'stdin.' ) and hasattr( f, 'close' ): f.close()
        logging.info( 'loaded stats table from ' + f_orig )

# end: def loadTable( f )

def saveTable( f, t ):
    """Save a numpy array ``t`` to file ``f``.  If f is a string that starts with ``stdout``, then stdout is used.
    f will be automatically compressed if it ends with .bz2 . tsv and npz files can be saved.
    """
    
    f_orig = f
    name, ext = os.path.splitext( f )
    try:
        if f_orig.startswith( 'stdout.' ): f = sys.stdout
        if ext == '.bz2':
            f = bz2.BZ2File( f, mode = 'w' )
            name, ext = os.path.splitext( name )
        if ext == '.tsv': np.savetxt( f, t, delimiter = '\t', comments = '', header = '\t'.join( t.dtype.names ) )
        elif ext == '.npz':
            result = np.savez_compressed( f, t )
            assert len( result ) == 1
            return list( result.items() )[0][1]
        else:
            raise RuntimeError( 'unknown extension for output file %s: "%s"' % ( f_orig, ext ) )
    finally:
        if not f_orig.startswith( 'stdout.' ) and hasattr( f, 'close' ): f.close()
        logging.info( 'saved stats table to ' + f_orig )

# end: def saveTable( f, t )        
        
if args.record:
    saveTable( args.file2, loadTable( args.file1 ) )
    logging.info( 'saved ' + args.file1 + ' to ' + args.file2 + ', exiting.' )
    sys.exit( 0 )

z1 = loadTable( args.file1 )
z2 = loadTable( args.file2 )
for n in z1.dtype.names:
    if n not in z2.dtype.names:
        print( 'column ' + str( n ) + ' is in z1 but not in z2' )
for n in z2.dtype.names:
    if n not in z1.dtype.names:
        print( 'column ' + str( n ) + ' is in z2 but not in z1' )
logging.info( 'z1.dtype.names=' + str( z1.dtype.names ) + ' z2.dtype.names=' + str( z2.dtype.names ) )
logging.info( 'checking two files, sizes %d and %d, for %d stats' % ( len( z1 ), len( z2 ), len( z1.dtype.names ) ) )
#assert z1.dtype.names == z2.dtype.names
pVals = []

excludeColsRegexps = list(map( re.compile, args.exclude_cols )) if args.exclude_cols else ()

for c in z1.dtype.names:
    if c not in z2.dtype.names or any([ ecr.match( c ) for ecr in excludeColsRegexps  ]): continue
    
    D, p = ks_2samp( z1[ c ], z2[ c ] )
    if p < args.p:
        raise RuntimeError( 'distribution of column %s does not match: p=%s, D=%s' % ( c, p, D ) )
    pVals.append( ( p, c ) )

sorted_pvals = sorted( pVals ) 
    
if not args.v: sorted_pvals = sorted_pvals[:5]
    
logging.info( '\n' + '\n'.join( map( str, sorted_pvals) ) )

logging.info( 'cmpdist completed successfully' )
        

                            
    

