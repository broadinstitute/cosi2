
/*
 * Header: leafset.h
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 */

#ifndef __INCLUDE_COSI_LEAFSET_COUNTS_H
#define __INCLUDE_COSI_LEAFSET_COUNTS_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stack>
#include <valarray>
#include <vector>
#include <boost/intrusive_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>

namespace cosi {
namespace leafset_counts {

/* Type: leaf_id_t */
/* Identifier of one leaf node of the ARG, representing a present-day <chrom>. */
typedef nodeid leaf_id_t;

const leaf_id_t NULL_LEAF_ID = -1;

struct leafset_struct;

typedef leafset_struct leafset_t;
typedef boost::intrusive_ptr<leafset_t> leafset_p;
leafset_t * const LEAFSET_NULL = ((leafset_t * const)NULL);

extern const std::vector< popid > *leafset_leaf2popName;
extern const std::vector<pop_idx_t> *leafset_popname2idx;
extern int leafset_npops;


// Struct: leafset_t
// A set of leaves, i.e. present-day chroms.
//
// This implementation represents leafset as trees of nodes.
// Leaf nodes correspond to singleton leafsets.  Inner nodes
// represent unions of their children.
//
struct leafset_struct: public util::refcounted<leafset_struct> {
public:
	 // Field: chromCounts
	 // For each popidx, the number of chroms in that pop in this leafset.
	 std::valarray<nchroms_t> chromCounts;
	 
   void* operator new (size_t size);
   void operator delete (void *p);

	 leafset_struct( leaf_id_t leafId_ ):
		 chromCounts( 0, leafset_npops ) {
		 ::cosi::util::chk( leafset_leaf2popName, "did not set leafset info" );
		 popid popName = (*leafset_leaf2popName)[ leafId_ ];
		 pop_idx_t popidx = (*leafset_popname2idx)[ ToInt( popName ) ];
		 chromCounts[ popidx ] = 1;
	 }

	 leafset_struct( leafset_p childA_, leafset_p childB_ ):
		 chromCounts( childA_->chromCounts + childB_->chromCounts ) { }
	 
#ifdef COSI_R2
	 
	 // Func: compute_r2
	 // Compute the r^2 measure of linkage disequilibrium between the two leafsets.
	 static cosi_double compute_r2( leafset_p leafset1, leafset_p leafset2 );

private:
	 
	 // Field: leaves
	 // A sorted vector of the leafIds in this leafset, or NULL if not computed yet.
	 mutable boost::scoped_ptr< vector< leaf_id_t > > leaves;

	 void computeLeaves() const;
#endif	 
	 
};  // leafset_struct

/* FuncP: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
/* Must be called once, prior to the creation of any leafsets. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ );
leaf_id_t leafset_get_max_leaf_id(void);  

/* FuncP: make_singleton_leafset */
/* Returns the leafset containing only the specified leaf. */
inline leafset_p make_singleton_leafset( leaf_id_t leaf ) { return new leafset_t( leaf ); }

inline leafset_p make_empty_leafset() { return LEAFSET_NULL; }

inline bool_t leafset_is_empty( leafset_p leafset ) { return leafset == LEAFSET_NULL; }

/* FuncP: leafset_union */
/* Returns the leafset containing the union of the given leafsets. */
inline leafset_p leafset_union( leafset_p leafset1, leafset_p leafset2 ) {
	if ( leafset_is_empty( leafset1 ) ) return leafset2;
	if ( leafset_is_empty( leafset2 ) ) return leafset1;
	return new leafset_t( leafset1, leafset2 );
}

bool leafset_equal( leafset_p leafset1, leafset_p leafset2 );

leafset_p leafset_intersection( leafset_p leafset1, leafset_p leafset2 );
leafset_p leafset_difference( leafset_p leafset1, leafset_p leafset2 );

const char *leafset_str( leafset_p  );

//inline bool_t leafset_is_singleton( leafset_p leafset, leaf_id_t leaf ) { return leafset->leafId == leaf; }

void leafset_print( leafset_p leafset );

int leafset_num_alloced(void);

leafset_p leafset_from_str( const char *s );

inline bool_t leafset_may_differ( leafset_p leafset1, leafset_p leafset2) { return leafset1 != leafset2 ; }

inline
nchroms_t leafset_size( leafset_p leafset ) {
	return leafset_is_empty( leafset ) ? 0 : leafset->chromCounts.sum();
}

inline bool leafset_is_full( leafset_p leafset ) { return leafset_size( leafset ) == leafset_get_max_leaf_id(); }

// template <class OutputIter>
// void leafset_get_leaves( leafset_p leafset, OutputIter result ) {
// 	assert(0);
// }

#define COSI_FOR_LEAFSET(leafset,leaf_var,body) do {	assert(0); } while(0)

ostream& operator<<( std::ostream& s, leafset_p leafset );

} // namespace leafset_counts

using namespace leafset_counts;

} // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_LEAFSET_COUNTS_H
