/* $Id: mutate.c,v 1.4 2011/06/07 15:29:25 sfs Exp $ */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cosi/mutate.h>
#include <cosi/mutlist.h>
#include <cosi/demography.h>
#include <cosi/utils.h>
#include <cosi/seglist.h>

namespace cosi {

MutProcessor::~MutProcessor() {}
void MutProcessor::processMut(loc_t /*loc*/, leafset_p /*leaves*/, genid /*gen*/, popid /*popName*/) {}
void MutProcessor::postprocess() {}

Mutate::Mutate( RandGenP randGen_, prob_per_bp_per_gen_t mu_, len_bp_t chrom_len_bp_, bool_t /*treeSizeOnly_*/ ):
	HasRandGen( randGen_ ), mutate_mu( mu_ * chrom_len_bp_ ), mutate_chrom_len_bp( chrom_len_bp_ ),
	waitInitialized( False ) {
	setMutProcessor( boost::make_shared<MutProcessor>() );
}

// Method: setMutProcessor
// Sets the object to receive mutations as they are generated during the simulation.
void Mutate::setMutProcessor( MutProcessorP mutProcessor_ ) { mutProcessor = mutProcessor_; }

void Mutate::mutate_print_leafset (loc_t loc, leafset_p leaves, genid gen, popid popName) {
#ifndef COSI_DEV_BINOUTPUT
	mutProcessor->processMut( loc, leaves, gen, popName );
#else	
	fwrite( &loc.val, sizeof( loc.val ), 1, stdout );
#endif	
}


/*
 * Func: mutate_put_muts_on_seglist
 *
 * Place mutations on one edge of the <ARG>.
 *
 * Params:
 *
 *   seglist - the seglist that gives the segments inherited along the ARG edge
 *   edge_gens - the length of the ARG edge, in generations.
 *
 * Input/output global vars:
 *
 *   wait_left - the "wait time" left until the next mutation
 */
#ifndef COSI_DISABLE_MUTS
void Mutate::mutate_put_muts_on_seglist( const Seglist *seglist, gens_t edge_gens,
																				 genid edgeMinGen, popid popName ) {

	if ( !waitInitialized ) {
		/* Choose the "wait time" until the first mutation. */
		/* Wait times between mutations are exponentially distributed. */
		wait_left = mutime_t( expdev() / mutate_mu );
		waitInitialized = True;
	}
	
  if ( !seglist_is_empty( seglist ) && edge_gens > ZERO_GENS ) {

		/* Local var: len_left */
		/* The length left to walk along the segs of the seglist, until the next mutation. */
		len_t len_left = wait_left / edge_gens;

		/* Check whether the wait time until next mutation is so large that */
		/* we must skip this seglist entirely. */
		len_t new_len_left = len_left - seglist_tot_len( seglist );
		if ( is_positive( new_len_left ) ) {
			wait_left = new_len_left * edge_gens;
			return;
		}

		static seglist::Finger *finger;
		if ( !finger ) finger = seglist::seglist_finger_alloc();
	 
		seglist::seglist_init_finger( const_cast<Seglist *>(seglist), finger, /* init_lens= */ True );

		/* Local var: cur_loc */
		/* The current location as we move along the seglist. */
		loc_t cur_loc = MIN_LOC;

		/* Local var: cur_len */
		/* The length traveled so far along the seglist. */
		len_t cur_len = ZERO_LEN;
		leafset_p cur_leafset = LEAFSET_NULL;

		while( seglist::seglist_advance_by_len( seglist, finger, &cur_loc, &cur_len, &cur_leafset, &len_left ) ) {
			mutate_print_leafset( cur_loc, cur_leafset, edgeMinGen + factor_t( .5 ) * edge_gens, popName );
			/* Choose the length we must travel along the seglist's segs until the next mutation */
			len_left = expdev() / ( mutate_mu * edge_gens );
		}

		/* Record the wait time left until the next mutation */
		wait_left = len_left * edge_gens;
  }  /* if ( !seglist_is_empty( seglist ) && edge_gens > 0 ) */
}  /* Mutate::mutate_put_muts_on_seglist() */
#endif // #ifndef COSI_DISABLE_MUTS

}  // namespace cosi




