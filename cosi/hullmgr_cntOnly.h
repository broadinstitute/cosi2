#ifndef __COSI_INCLUDE_HULLMGR_H
#define __COSI_INCLUDE_HULLMGR_H

#include <cstdlib>
#include <utility>
#include <cosi/defs.h>
#include <cosi/general/utils.h>
#include <cosi/general/datastruct/order_statistics.hpp>
#include <cosi/general/math/cosirand.h>

namespace cosi {

namespace node { class Node; }

#ifdef COSI_SUPPORT_COALAPX

//
// Class: HullMgr
//
// Keeps track of the convex hulls of the seglists of all <Nodes>, and of the pairs of intersecting hulls.
// Supports choosing a pair of intersecting hulls uniformly at random from all such intersecting pairs, and
// updating the hull information through all standard operations (recombinations, coalescences,
// gene conversions, migrations).
//
// Complexity:
//
//    space complexity - for n hulls, O(n) storage.
//    time complexity - O(lg n) for update operations.
//
// Implementation notes:
//
// Because recombination points are chosen randomly on a uniform line, we assume that no two such points
// fall into exactly the same place.  So, at each loc except <MIN_LOC>, nbegs+nends <= 2 where nbegs is
// the number of hulls starting at that loc and nends is the number of hulls ending at that loc.
//
// Design rationale:
//
//    - why a simple weighted interval map is not enough?  we need to keep partial sums.  but, presumably
//      a map does that?   
//
// Implementation todo:
//
//   - handle additional cases
//
//       - when something fully coalesces
//       - migration in and out
//       - gene conversion (can/should we re-implement it as a pair of recombs?)
//
//   - keep track of node ids
//
//   - consider using knowledge of the genetic map to optimize access times to elements in hotspots.
//     note though, that for coalescence operations we might need different access probabilities.
//     but, we have some info on them if we keep the intersection counts.
//
//   - zero things out instead of deleting?
//
//   - remove seglist_summaries if no longer helpful
//
//   - avoid multiple tree traversals
//     - support finding several cose-by things at once
//
class HullMgr {
public:

	 typedef node::Node Node;

	 // Logical type: ninters_t
	 // A count of hull intersections; has units of nchroms * nchroms.
	 typedef int ninters_t;  // signed, for now, to allow negative deltas.

	 typedef util::order_statistics_tree< /* ValueType= */ loc_t > ost_ends_t;

	 struct HullBeg {
			loc_t beg;
			Node *node;
			ost_ends_t::iterator end_it;

			HullBeg(): node( NULL ) { }
			HullBeg( loc_t beg_, Node *node_, ost_ends_t::iterator end_it_ ):
				beg( beg_ ), node( node_ ), end_it( end_it_ ) { }

			loc_t getBeg() const { return beg; }
			loc_t getEnd() const { return *end_it; }
	 };

	 typedef util::order_statistics_tree< HullBeg, std::less<loc_t>,
																				cosi::util::FieldRef< HullBeg, loc_t, &HullBeg::beg >
																				> ost_begs_t;
	 typedef ost_begs_t::iterator Hull;

			
	 // Constructor: HullMgr
	 // Construct a HullMgr.  It is not yet initialized (see <addLeaves()> for that).
	 // Params:
	 //    margin_ - specifies what near-intersections are considered to be intersections;
	 //      specifically, for a seglist [beg,end] the corresponding hull is [beg, end+margin].
	 HullMgr();
	 virtual ~HullMgr();

	 void setMargin( len_t margin_ ) { this->margin = margin_; }
	 void setMargin( plen_t margin_ ) { this->setMargin( len_t( margin_ ) ); }

	 len_t getMargin() const { return this->margin; }

	 // Method group: ARG operations
	 
	 Hull addHull( Node *node, loc_t beg, loc_t end ); 
	 void removeHull( Hull hull ); 

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
	 void recomb( Hull hull, Node *newNode1, Node *newNode2,
								loc_t beg1, loc_t end1, loc_t beg2, loc_t end2, Hull *newHull1, Hull *newHull2 );
	 
	 void coalesce( Hull hull1, Hull hull2 );
	 Hull coalesce( Hull hull1, Hull hull2, 
									Node *newNode,
									loc_t newbeg, loc_t newend );

	 // End method group: ARG operations


	 // Method group: Debugging / testing
	 
private:
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
	 ninters_t ninters;

public:

	 // Method: getNumIntersections
	 // Returns the current number of intersecting hull pairs.
	 ninters_t getNumIntersections() const { return ninters; }

	 std::pair< Node *, Node * > chooseRandomIntersection( RandGenP randGen );
	 

};  // class HullMgr

#endif  // #ifdef COSI_SUPPORT_COALAPX

}  // namespace cosi

#endif  // #ifndef __COSI_INCLUDE_HULLMGR_H
