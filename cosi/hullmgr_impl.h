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

//#undef NDEBUG

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility.hpp>
#include <boost/phoenix.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <cosi/general/utils.h>
#include <cosi/general/cosirand.h>
#include <cosi/general/datastruct/order_statistics.hpp>
//#include <cosi/seglist.h>
#include <cosi/hullmgr.h>

namespace cosi {

#ifdef COSI_SUPPORT_COALAPX

//static const loc_t INF_LOC( std::numeric_limits<cosi_double>::infinity() );


template<>
struct
pimpl<HullMgr>::implementation {

	 implementation(): margin( MAX_LOC - MIN_LOC ), ninters( 0 ), isRestrictingCoalescence( false ) { }

	 typedef util::order_statistics_tree< /* ValueType= */ loc_t > ost_ends_t;

	 // Struct: HullBeg
	 // Information about one hull.
	 struct HullBeg {
			// Field: beg
			// The location of the hull's beginning.  More than one hull may begin at a given loc,
			// but they're still ordered left-to-right by the time of addition.
			loc_t beg;

			// Field: end_it
			// Iterator pointing to this hull's endpoint in <HullMgr::ends>.
			ost_ends_t::iterator end_it;

			// Field: ninters
			// The number of intersections which begin at this hull; that is, the number of other hulls
			// with beg to the left and end to the right of <beg>.
			nchromPairs_t ninters;

			void *hullPtr;

			HullBeg(): ninters( 0 ), hullPtr( NULL ) { }
			HullBeg( loc_t beg_, ost_ends_t::iterator end_it_, void *hullPtr_ ):
				beg( beg_ ), end_it( end_it_ ), ninters( 0 ), hullPtr( hullPtr_ ) { }

			loc_t getBeg() const { return beg; }
			loc_t getEnd() const { return *end_it; }
	 };  // struct HullBeg

	 typedef util::order_statistics_tree< /* ValueType= */ HullBeg, /* Comparator= */ std::less<loc_t>,
																				/* KeyExtractor= */ cosi::util::FieldRef< HullBeg, loc_t,
																																									&HullBeg::beg >,
																				/* WeightExtractor= */ cosi::util::FieldRef< HullBeg, nchromPairs_t,
																																										 &HullBeg::ninters >								 
																				> ost_begs_t;
	 typedef ost_begs_t::iterator Hull;

	 
	 // Field: margin
	 //
	 // Specifies what near-intersections are considered to be intersections;
	 // specifically, for a seglist [beg,end] the corresponding hull is [beg, end+margin].
	 // So, two seglists [beg1,end1] and [beg2,end2] where end1 < beg2 have intersecting hulls
	 // if end1 is within 'margin' of beg2.
	 len_t margin;

	 // Field: begs
	 // The locations of all hull begs.
	 ost_begs_t begs;

	 // Field: ends
	 // The locations of hull ends.
	 ost_ends_t ends;

	 // Field: ninters
	 // Current number of intersections.
	 nchromPairs_t ninters;

	 // Field: isRestrictingCoalescence
	 // Whether we're restricting pairs that may coalesce.
	 bool isRestrictingCoalescence;
	 

};

	 // Constructor: HullMgr
	 // Construct a HullMgr.  It is not yet initialized (see <addLeaves()> for that).
	 // Params:
	 //    margin_ - specifies what near-intersections are considered to be intersections;
	 //      specifically, for a seglist [beg,end] the corresponding hull is [beg, end+margin].

	 // Method group: ARG operations
	 
	 // Method: recomb
	 //
	 // Called after each recombination to update the hull information.
	 //
	 // Params:
	 //
	 //    beg1,end1 and beg2,end2 - the two hulls resulting from the recombination.
	 //      The original hull was [beg1,end2].  If the recomb point fell in the middle of a seg,
	 //      then end1 == beg2.  If the recomb point fell in a gap, then end1 < beg2.
	 //
	 //    loc - the location of recombination
	 //
	 // Complexity: lg n, where n is the current number of hulls.
	 //
	 // void recomb( Hull hull, Node *newNode1, Node *newNode2,
	 // 							loc_t beg1, loc_t end1, loc_t beg2, loc_t end2, Hull *newHull1, Hull *newHull2 );
	 
	 // void coalesce( Hull hull1, Hull hull2 );
	 // Hull coalesce( Hull hull1, Hull hull2, 
	 // 								Node *newNode,
	 // 								loc_t newbeg, loc_t newend );

	 // End method group: ARG operations


	 // Method group: Debugging / testing

#endif  // #ifdef COSI_SUPPORT_COALAPX

}  // namespace cosi

