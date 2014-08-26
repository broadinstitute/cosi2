#!/usr/bin/env python

"""For a given cosi model file, create a regression test for checking that the distribution of various
statistics of cosi output for this file has not changed."""

import argparse, os, sys, logging
import numpy as np
print 'importing stats'
from scipy import stats
print 'imported stats'
from Classes.DotData import DotData
from Operations.MiscUtil import SystemSucceed, EnsureDirExists, dbg
print 'imported rest'

logging.basicConfig( level = logging.DEBUG, format='%(asctime)s %(levelname)-8s %(filename)s:%(lineno)s %(message)s' )

print 'making parser'
parser = argparse.ArgumentParser( description = 'Create regtest for cosi based on given model',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( '--test-name', help = 'name of test; is name of directory under tests/ where test is saved' )
parser.add_argument( '--bootstrap-count', type = int, default = 1000, help = 'number of boostrap iters' )
parser.add_argument( '--cosi-binary', default = './coalescent', help = 'cosi binary to run' )

print 'calling parser'
args = parser.parse_args()
print 'parser done'
# do a reference run

print 'generating reference'
SystemSucceed( ' '.join(map( str, (args.cosi_binary, '-p', '1_simple.cosiParams', '-n', 100, '-m' ))) + ' | sample_stats_extra > ref.tsv' )
refData = DotData( SVPath = 'ref.tsv' )
min_p = np.ones( len( refData.dtype.names ) )
max_D = np.repeat( -np.inf, len( refData.dtype.names ) )
for i in range( 10 ):
    dbg( 'i' )
    refFN = 'reftest%d.tsv' % i
    SystemSucceed( ' '.join(map( str, ( args.cosi_binary, '-p', '0_simple.cosiParams', '-n', 100, '-m' ))) + ' | sample_stats_extra > ' + refFN )
    z = DotData( SVPath = refFN )
    for colNum, col in enumerate( z.dtype.names ):
        ks_D, ks_p = stats.ks_2samp( refData[ col ], z[ col ] )
        min_p[ colNum ] = np.min(( min_p[ colNum ], ks_p ))
        max_D[ colNum ] = np.max(( max_D[ colNum ], ks_D ))
    dbg( 'i min_p max_D' )

    

        
        
    
    








