
/*
 * Header: leafset.h
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 */

#ifndef __INCLUDE_COSI_LEAFSET_TREE_H
#define __INCLUDE_COSI_LEAFSET_TREE_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stack>
#include <boost/intrusive_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>

namespace cosi {
namespace leafset_tree {

/* Type: leaf_id_t */
/* Identifier of one leaf node of the ARG, representing a present-day <chrom>. */
typedef nodeid leaf_id_t;

const leaf_id_t NULL_LEAF_ID = -1;

struct leafset_struct;

typedef leafset_struct leafset_t;
typedef boost::intrusive_ptr<leafset_t> leafset_p;
leafset_t * const LEAFSET_NULL = ((leafset_t * const)NULL);

// Struct: leafset_t
// A set of leaves, i.e. present-day chroms.
//
// This implementation represents leafset as trees of nodes.
// Leaf nodes correspond to singleton leafsets.  Inner nodes
// represent unions of their children.
//
struct leafset_struct: public util::refcounted<leafset_struct> {
public:
	 // Field: size
	 // Number of leaves in the leafset.
	 nchroms_t size;
	 
	 // Field: leafId
	 // If this is a singleton leafset representing one leaf,
	 // the id of that leaf; otherwise, <NULL_LEAF_ID>.
	 leaf_id_t leafId;

	 // Fields: children
	 // If this is a singleton leafset representing one leaf,
	 // then LEAFSET_NULL, else the two leafsets whose disjoint
	 // union equals this leafset.
	 leafset_p childA,childB; 

   void* operator new (size_t size);
   void operator delete (void *p);

	 leafset_struct( leaf_id_t leafId_ ):
		 size( 1 ), leafId( leafId_ ) {}

	 leafset_struct( leafset_p childA_, leafset_p childB_ ):
		 size( childA_->size + childB_->size ), leafId( NULL_LEAF_ID ), childA( childA_ ), childB( childB_ ) { }

	 
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

leafset_p make_range_leafset( leaf_id_t fromLeaf, leaf_id_t toLeaf );


bool leafset_equal( leafset_p leafset1, leafset_p leafset2 );

leafset_p leafset_intersection( leafset_p leafset1, leafset_p leafset2 );
leafset_p leafset_difference( leafset_p leafset1, leafset_p leafset2 );

const char *leafset_str( leafset_p  );

inline bool_t leafset_is_singleton( leafset_p leafset, leaf_id_t leaf ) { return leafset->leafId == leaf; }

void leafset_print( leafset_p leafset );

int leafset_num_alloced(void);

leafset_p leafset_from_str( const char *s );

inline bool_t leafset_may_differ( leafset_p leafset1, leafset_p leafset2) { return leafset1 != leafset2 ; }

inline
nchroms_t leafset_size( leafset_p leafset ) {
	return leafset_is_empty( leafset ) ? 0 : leafset->size;
}

inline bool leafset_is_full( leafset_p leafset ) { return leafset_size( leafset ) == leafset_get_max_leaf_id(); }

template <class OutputIter>
void leafset_get_leaves( leafset_p leafset, OutputIter result ) {
	if ( !leafset_is_empty( leafset ) ) {
		std::stack<const leafset_t *> s;
		s.push( leafset.get() );
		while( !s.empty() ) {
			const leafset_t *p = s.top();
			s.pop();
			if ( p->leafId != NULL_LEAF_ID ) {
				*result++ = p->leafId;
			} else {
				s.push( p->childA.get() );
				s.push( p->childB.get() );
			}
		}
	}
}

#define COSI_FOR_LEAFSET(leafset,leaf_var,body) do {	  \
	  if ( !leafset_is_empty( leafset ) ) {               \
		  std::stack<const leafset_t *> s;									\
		  s.push( leafset.get() );													\
		  while( !s.empty() ) {															\
			  const leafset_t *p = s.top();										\
			  s.pop();																				\
			  if ( p->leafId != NULL_LEAF_ID ) {							\
		  		leaf_id_t leaf_var = p->leafId;								\
				 { body; }																			\
			  } else {																				\
				  s.push( p->childA.get() );										\
				  s.push( p->childB.get() );										\
			  }  																							\
		  }  /* while( !s.empty() ) */											\
	  }  /* if leafset nonempty */										    \
  } while(0)


ostream& operator<<( std::ostream& s, leafset_p leafset );

} // namespace leafset_tree

using namespace leafset_tree;

} // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_LEAFSET_TREE_H
