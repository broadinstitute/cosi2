//
// notes
//
//   - if two begs at same place, the pair gets double-counted.

//   - if initial portion coalesces, have to subtract from number at zero.

//   - could unify the handling with beg0:
//       - at each pt, keep a count of how many beg here
//       - if >=2, count the pairs
//       - then, how many cross this pt?  (tot - endbefore - begin-here-or-later)
//       - if smth ends at or past 1, it doesn't endbefore.

//       - if two begs at same pt coalesce, decrement count of begs there
//
//    - what happens if a portion coalesces?
//
//         - only matters if two ends or two begs match?
//
//                         ==================
//                         =================
//
//               - then, both begs might disappear, and then possibly a new hullbeg be inserted (or not, if no new node gets created).
//               - for purpose of indexing count these as two nodes (or actually keep them as separate).
//
//               - if two ends match, count of ends there goes down.  also, possibly, a new hullend is inserted which was not there before.
//
//    - so, this deletes the old beg and end, correctly
//
//    - node can hold iterator to the beg/end records of its hull
//
//    - need ability to write the weight field, to adjust when inserting.
//    - if range being updated isn't large, simpler to just manually update?
//        - how common are large-range updates?  because if not common, may be more efficient to just update manually.
//        - delayed insertion: just store things till needed for query?
//
//     - we're storing node ids... ok better to just allow duplicates in the skiplist.
//       should never be more than two here.  just make sure it's counted correctly.
//         - when deleting, does it matter which?  it does if they have different inters-counts.
//
//      - adjust the new-level constant as needed.   can rebalance the list if searches get long.
//

#undef NDEBUG

#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <cosi/hullmgr.h>
#include <cosi/utils.h>


namespace cosi {

#ifdef COSI_SUPPORT_COALAPX

static const loc_t INF_LOC( std::numeric_limits<cosi_double>::infinity() );

HullMgr::HullMgr( ): margin( MAX_LOC - MIN_LOC ), ninters( 0 ) {
	std::cerr.precision(14);

	// begs.setIntrusive( true );
	// ends.setIntrusive( true );

	// begs.insert( HullBeg( INF_LOC, ninters_t( 0 ), (Node *)NULL ) );
	
	// ost_t t;

	// ost_iter_t v1 = t.insert( loc_t( .5 ) ).first;
	// ost_iter_t v2 = t.insert( loc_t( .5 ) ).first;
	// ost_iter_t v3 = t.insert( loc_t( .5 ) ).first;
	// PRINT3( v1.position(), v2.position(), v3.position() );
	// PRINT( t.lower_bound( loc_t( .5 ) ).position() );
	// PRINT( t.upper_bound( loc_t( .5 ) ).position() );
	// t.erase( v2 );
	// PRINT2( v1.position(), v3.position() );
}

HullMgr::~HullMgr() { }

// Function: last_among_equal
// Returns an iterator pointing to the rightmost item in a range of equal items.
template <typename ForwardIterator>
ForwardIterator last_among_equal( ForwardIterator i, ForwardIterator end ) {
	while ( i != end ) {
		ForwardIterator n = boost::next( i );
		if ( n != end && *n == *i )
			 i = n;
		else
			 break;
	}
	return i;
}

// Function: first_among_equal
// Returns an iterator pointing to the leftmost item in a range of equal items.
template <typename BidiIterator>
BidiIterator first_among_equal( BidiIterator i, BidiIterator begin ) {
	while ( i != begin ) {
		BidiIterator p = boost::prior( i );
		if ( *p == *i )
			 i = p;
		else
			 break;
	}
	return i;
}

const len_t smallestLen( std::numeric_limits<cosi_double>::min() );

HullMgr::Hull HullMgr::addHull( Node *node, loc_t beg, loc_t end ) {
	//PRINT6( "bef_add_hull", beg, end, begs.size(), ends.size(), getNumIntersections() );
	ninters += ( begs.lower_bound( loc_t( end + margin ) ).position() - ends.lower_bound( loc_t( beg - margin ) ).position() );
	return begs.insert( HullBeg( beg, node, ends.insert( end ).first ) ).first;
	
	//PRINT6( "aft_add_hull", beg, end, begs.size(), ends.size(), getNumIntersections() );
}

void HullMgr::removeHull( Hull hull ) {
	//PRINT6( "bef_rem_hull", beg, end, begs.size(), ends.size(), getNumIntersections() );

	loc_t beg = hull->beg, end = *hull->end_it;
	ends.erase( hull->end_it );
	begs.erase( hull );

	ninters -= ( begs.lower_bound( loc_t( end + margin ) ).position() - ends.lower_bound( loc_t( beg - margin ) ).position() );
	//PRINT6( "aft_rem_hull", beg, end, begs.size(), ends.size(), getNumIntersections() );
}

// Method: recomb
// Called after each recombination to update the hull information.
// Params:
//    beg1,end1 and beg2,end2 - the two hulls resulting from the recombination.
//    The original hull was [beg1,end2].  If the recomb point fell in the middle of a seg,
//    then end1 == beg2.  If the recomb point fell in a gap, then end1 < beg2.
//
// Complexity: lg n, where n is the current number of hulls.
//
void HullMgr::recomb( Hull hull, Node *newNode1, Node *newNode2,
											loc_t beg1, loc_t end1, loc_t beg2, loc_t end2, Hull *newHull1, Hull *newHull2 ) {
//	PRINT11( "bef_recomb", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2, loc, end1 == beg2 );
	removeHull( hull );
	*newHull1 = addHull( newNode1, beg1, end1 );
	*newHull2 = addHull( newNode2, beg2, end2 );
//	PRINT11( "aft_recomb", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2, loc, end1 == beg2 );	
	
#if 0	
	PRINT11( "bef_recomb", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2, loc, end1 == beg2 );
	
	assert( MIN_LOC <= beg1 && beg1 < end1 && end1 <= beg2 && beg2 < end2 && end2 <= MAX_LOC );
	
	using std::make_pair;
	
	//
	// What happens here:
	//
	//    - a new beg is added, beg2.  We must compute the number of intersections starting at beg2,
	//    i.e. the number of hulls in the middle of which beg2 falls.
	//    Beware that there may already be one hull beginning at beg2, if the seg immediately left of beg2 has
	//    fully coalesced; if there is a beg at beg2, its intersection
	//    with beg2 must be counted.
	//
	//    beg2 cannot be MIN_LOC, so numBeg0 is not affected.
	// 
	//    - a new end is added, end1.  We add it to <ends>.  Note that there may already be one end here,
	//    if the seg immediately right of end1 has fully coalesced.
	//
	//    - any existing begs in [beg2, end1+margin) previously intersected one hull, now they intersect two;
	//      so increment the count for them.
	//    
	//
	
	// Step 1: Count how many intersections that begin at beg2 are added.
	//
	// We subtract from the total:
	//     - any hulls that end before or at beg2
	//     - any hulls that begin after beg2.
	//
	// For hulls that begin at beg2, there may be one such existing hull (if seg left of beg2 has
	// fully coalesced).  We can say that beg2 landed slightly to the right of it,
	// so we do not subtract it from the total.  The alternative would have been to
	// say that we land slightly to the left; which can't, because it has coalesced.
	//
	assert( begs.find( loc ) == begs.end() );
	assert( ends.find( loc ) == ends.end() );
	
	assert(  // loc either broke an existing seg ...
		( end1 == beg2 && beg2 == loc )
		// ...or fell in a gap
		|| ( end1 < loc && loc < beg2 ) );
	
	begs.check();
	ends.check();
	
	// Determine the number of intersection that will start at beg2.
	// (Possibly including the intersection with the new hull [beg1,end1+margin]).
	
	// How many hulls begin left of loc?
	// How many end left of loc-margin?
	
	// (both of these can be answered on the same traversal that inserted loc;
	// initially, can just ask these from scratch)
	
	loc_t end1margin( end1 + margin );
	
	//if ( end1margin < MAX_LOC )
	ends.insert( end1 );
	
	PRINT3( "cntbef", begs.count( beg2 ), ends.count( end1 ) );
	
	PRINT( end1margin );
	
	ost_begs_t::const_iterator beg2_it = begs.lower_bound( beg2 );
	
#ifndef	NDEBUG
	ost_begs_t::const_iterator end1margin_it = begs.lower_bound( end1margin );
	assert( end1margin_it == begs.end() || end1margin_it->first > end1margin );
#endif	
	
	begs.addWeight( beg2, end1margin, 1 );
	begs.addWeight( end1margin, len_t( beg2 - smallestLen ), -1 ) ;
	
	if ( equal_eps( beg2, loc_t( 0.63281095310391 ) ) ) {
		if ( beg2_it != begs.end() ) {
			PRINT3( "nu", beg2_it->first, beg2_it->second );
			if ( beg2_it.position() > 0 ) PRINT2( boost::prior( beg2_it )->first, boost::prior( beg2_it )->second );
			if ( boost::next( beg2_it ) != begs.end() ) PRINT2( boost::next( beg2_it )->first, boost::next( beg2_it )->second );
		}
	}
	
	nchroms_t n_begs_bef_beg2 = nbegs0 + beg2_it.position();
	assert( ends.find( loc_t( beg2 - margin ) ) == ends.end() );
	nchroms_t n_ends_bef_beg2 = ends.upper_bound( loc_t( beg2 - margin ) ).position();
	
	ninters_t ninters_at_beg2 = n_begs_bef_beg2 - n_ends_bef_beg2;
	
	PRINT6( n_begs_bef_beg2, n_ends_bef_beg2, ninters_at_beg2, beg2_it == begs.end(), beg2_it == begs.end() ? 0 : beg2_it->second, beg2_it == begs.end() ? MIN_LOC : beg2_it->first );

	static size_t counter = 0;
	if ( counter == 2248 ) {
		
	}

	PRINT( counter++ );
	
	
	assert( beg2_it == begs.end() || beg2_it->first > beg2  ||  ( beg2_it->first == beg2 && beg2_it->second == ninters_at_beg2+1 ) );
	
	ost_begs_t::const_iterator beg2_new_it = begs.insert( make_pair( beg2, ninters_at_beg2 ) ).first;

	assert( beg2_new_it->first == beg2 );
	assert( beg2_new_it->second == ninters_at_beg2 );
	assert( beg2_it == begs.end() || beg2_it->first > beg2 || boost::next( beg2_it ) == beg2_new_it );
	
	//
	// 
	//


	PRINT3( "cntaft", begs.count( beg2 ), ends.count( end1 ) );
		
	// debugging / reporting
	
	nchroms_t tot_chroms = nbegs0 + begs.size();
	PRINT13( "aft_recomb", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2, n_begs_bef_beg2, n_ends_bef_beg2, tot_chroms, loc );
#endif	
}  // HullMgr::recomb()


void HullMgr::coalesce( Hull hull1, Hull hull2 ) {
//	PRINT9( "bef_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2 );
	removeHull( hull1 );
	removeHull( hull2 );
//	PRINT9( "aft_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2 );
}
	
//
// Method: coalesce
//
// Called after two nodes coalesce, to update the hull information and intersection count.
HullMgr::Hull HullMgr::coalesce( Hull hull1, Hull hull2, Node *newNode,
																 loc_t newbeg, loc_t newend ) {
	// PRINT11( "bef_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2,
	// 				newbeg, newend );
	removeHull( hull1 );
	removeHull( hull2 );
	return addHull( newNode, newbeg, newend );
	// PRINT11( "aft_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2,
	// 				newbeg, newend );
	
#if 0	
	PRINT9( "bef_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2 );

	loc_t overlap_beg = std::max( beg1, beg2 );
	loc_t overlap_end = std::min( end1, end2 );

	loc_t overlap_end_m( overlap_end + margin ); 
	assert( overlap_beg < overlap_end_m );

	assert( overlap_beg == MIN_LOC || begs.find( overlap_beg ) != begs.end() );
	assert( overlap_end_m > MAX_LOC || ends.find( overlap_end ) != ends.end() );
	
	if ( overlap_beg == MIN_LOC )
		nbegs0--;

	//
	// overlap_beg is being deleted from begs.
	// Let's count how many hull intersections there were starting at overlap_beg; we'll need to subtract them
	// from the inters count.  If there are multiple begs at overlap_beg, we say that the one we're deleting
	// was the rightmost.
	//

	PRINT2( ends.count( overlap_end ), begs.count( overlap_beg ) );

	ost_ends_iter_t i_end = ends.find( overlap_end );
	if ( i_end != ends.end() ) ends.erase( i_end );

	ost_begs_iter_t i_beg = begs.find( overlap_beg );
	if ( i_beg != begs.end() ) begs.erase( i_beg );

	PRINT2( ends.count( overlap_end ), begs.count( overlap_beg ) );
	

	begs.addWeight( overlap_beg, overlap_end_m, -1 );

	//if ( overlap_end_m < MAX_LOC ) ends.erase( overlap_end );
	
	
	PRINT11( "aft_coal", getNumIntersections(), nbegs0, begs.size(), ends.size(), beg1, end1, beg2, end2, overlap_beg, overlap_end  );
#endif
}  // HullMgr::coalesce()

std::pair< cosi::node::Node *, cosi::node::Node * > HullMgr::chooseRandomIntersection( RandGenP randGen ) {

	int nattempts = 0;


	
	while ( true ) {
		nattempts++;

//		if ( nattempts == 6 ) std::cerr << "\n";


		size_t hull1_idx = randGen->random_idx( begs.size() );
		
		Hull hull1 = begs[ hull1_idx ];
		size_t hull1_beg_succ_pos = hull1.position() + 1;
		size_t hull1_end_pos = begs.lower_bound( loc_t( hull1->getEnd() + margin ) ).position();
		ninters_t nbegs_inside_hull1 = hull1_end_pos - hull1_beg_succ_pos;

		assert( nbegs_inside_hull1 >= 0 );

		if ( nbegs_inside_hull1 == 0 ) continue;

		size_t hull2_idx;

		if ( ninters == 1 ) {
			assert( nbegs_inside_hull1 == 1 );
			hull2_idx = 0;
		} else {

			prob_t rejectProb = 1.0 - std::min( ((double)begs.size()), ((double)ninters) / begs.size() ) * ( ((double)nbegs_inside_hull1) / ((double)ninters ) );
			frac_t rejectCoinResult = randGen->random_double();

			if ( nattempts > 5 ) {
			 	 PRINT10( nattempts, nbegs_inside_hull1, ninters, begs.size(), ends.size(), hull1_idx,
			 						hull1_beg_succ_pos, hull1_end_pos, rejectProb, rejectCoinResult );
			}
			if ( rejectCoinResult < rejectProb ) {
				extern unsigned long nretry;
				nretry++;
				continue;
			}

			hull2_idx = randGen->random_idx( nbegs_inside_hull1 );
		}  // if ninters > 1

		Hull hull2 = begs[ hull1_beg_succ_pos + hull2_idx ];

		assert( hull1->node != hull2->node );
		
		return std::make_pair( hull1->node, hull2->node );
	}
}


#endif  // #ifdef COSI_SUPPORT_COALAPX

}  // namespace cosi
