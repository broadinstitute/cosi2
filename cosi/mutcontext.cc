
#include <boost/foreach.hpp>
#include <cosi/general/utils.h>
#include <cosi/seglist.h>
#include <cosi/mutcontext.h>
#include <cosi/leafset.h>

#define ForEach BOOST_FOREACH

namespace cosi {
namespace mutcontext {

using std::vector;
using util::chkCond;

//
// Function: computeMutContexts
//
// For a mutation placed at location 'loc' of 'seglist', for leaves in the leafset of the seg
// of 'seglist' containing 'loc', return the boundaries of the mutcontext neighborhood
// around that mutation at that leaf.  In other words: the mutation was born on a specific
// ancestral chromosome; this shows the part of that ancestral chromosome 
// 
//
vector< BasicSeg_loc > computeMutContexts( const Seglist *seglist, loc_t loc ) {
	// Save mutcontext neighborhoods
			
	// at this node, find the seg containing the selpos.
	// then look at adjacent segs:
	// eg going left:
	// for each leaf in the core
	// so, when we get to the next left seg
	// intersect the set so far with the new set
	// and extend it to them.

	vector< BasicSeg_loc > mutContexts( leafset_get_max_leaf_id(), BasicSeg_loc( loc, loc ) );

	using seglist::Seg;
	std::vector< const Seg * > allSegs;
	ForEach( const Seg& s, *seglist )
		 allSegs.push_back( &s );

	size_t coreSegIdx = 0;
	while ( coreSegIdx < allSegs.size() ) {
		if ( allSegs[ coreSegIdx ]->contains( loc ) ) break;
		coreSegIdx++;
	}
	chkCond( coreSegIdx < allSegs.size() && allSegs[ coreSegIdx ]->contains( loc ) ,
					 "nowhere to put selected mut?" );

	for ( int dir = -1; dir <= +1; dir += 2 ) {

		size_t curSegIdx = coreSegIdx;
		//PRINT( coreSegIdx );

		leafset_p inters = allSegs[ curSegIdx ]->getLeafset();

		while( curSegIdx < allSegs.size() ) {
			//PRINT2( curSegIdx, allSegs[ curSegIdx ] );
			vector<leaf_id_t> leaves;
#ifdef COSI_FREQONLY
			assert(0);
#else			
			inters = leafset_intersection( inters, allSegs[ curSegIdx ]->getLeafset() );
			leafset_get_leaves( inters, std::back_inserter( leaves ) );
#endif			
			ForEach( leaf_id_t leaf, leaves ) {
				if ( dir == -1 )
					mutContexts[ leaf ].setBeg( allSegs[ curSegIdx ]->getBeg() );
				else
					mutContexts[ leaf ].setEnd( allSegs[ curSegIdx ]->getEnd() );
			}

			if ( ( dir == -1 &&  
						 ( curSegIdx == 0 || allSegs[ curSegIdx-1 ]->getEnd() < allSegs[ curSegIdx ]->getBeg() ) )
					 ||
					 ( dir == +1 &&  
						 ( curSegIdx+1 == allSegs.size() || allSegs[ curSegIdx+1 ]->getBeg() > allSegs[ curSegIdx ]->getEnd() ) ) )							 
								 
				 break;
					
			curSegIdx += dir;
		}

#if 0		
		if ( False ) {
			vector< leaf_id_t > leaves;
			leafset_get_leaves( sel_leaves, std::back_inserter( leaves ) );
			std::sort( leaves.begin(), leaves.end() );

			vector< len_bp_int_t > neighbLens;
			ForEach( leaf_id_t leaf, leaves ) {
				neighbLens.push_back( static_cast<len_bp_int_t>( 1000000.0 * ToDouble( dists[ leaf ] ) ) );
			}
			std::sort( neighbLens.begin(), neighbLens.end() );
			using namespace util;
				
			PRINT( neighbLens );
		}
#endif		
	}  // for each dir

	return mutContexts;
}  // computeMutContexts()
			

static mutContexts_t g_mutContexts;

void saveMutContexts( const Seglist *seglist, loc_t loc ) {
	vector< BasicSeg_loc > mutContexts( computeMutContexts( seglist, loc ) );
	g_mutContexts.insert( make_pair( loc, mutContexts ) );
}

const mutContexts_t& getSavedMutContexts() { return g_mutContexts; }

}  // namespace mutcontext
}  // namespace cosi



