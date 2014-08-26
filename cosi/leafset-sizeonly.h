//
// Header: leafset-sizeonly.h
//
// An implementation of <leafsets> which keeps only the size of the leafset,
// and not the identity of the leaves in it.  Can be used when we're only interested
// in the frequency of each mutation (as when computing the allele frequency spectrum).
//

#ifndef __INCLUDE_COSI_LEAFSET_SIZEONLY_H
#define __INCLUDE_COSI_LEAFSET_SIZEONLY_H

#include <iostream>
#include <cosi/defs.h>

namespace cosi {

typedef nodeid leaf_id_t;

typedef int leafset_t;
typedef int leafset_p;

const leaf_id_t NULL_LEAF_ID = -1;

const leafset_p LEAFSET_NULL(-1);

/* FuncP: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
/* Must be called once, prior to the creation of any leafsets. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ );
leaf_id_t leafset_get_max_leaf_id(void);

/* FuncP: make_singleton_leafset */
/* Returns the leafset containing only the specified leaf. */
inline leafset_p make_singleton_leafset( leaf_id_t /*leaf*/ ) { return leafset_t( 1 ); }

inline leafset_p make_empty_leafset() { return 0; }

inline bool_t leafset_is_empty( leafset_p leafset ) { return leafset == LEAFSET_NULL; }

/* FuncP: leafset_union */
/* Returns the leafset containing the union of the given leafsets. */
inline leafset_p leafset_union( leafset_p leafset1, leafset_p leafset2 ) { return leafset1 + leafset2; }

inline const char *leafset_str( leafset_p  ) { return "leafset"; }

//inline bool_t leafset_is_singleton( leafset_p leafset, leaf_id_t leaf ) { return leafset->leafId == leaf; }

inline void leafset_print( leafset_p /*leafset*/ )  { }

//int leafset_num_alloced(void);

inline leafset_p leafset_from_str( const char * /*s*/ ) { return LEAFSET_NULL; }

inline bool_t leafset_may_differ( leafset_p leafset1, leafset_p leafset2) { return leafset1 != leafset2 ; }

inline
nchroms_t leafset_size( leafset_p leafset ) { return leafset; }

inline bool leafset_is_full( leafset_p leafset ) { return leafset_size( leafset ) == leafset_get_max_leaf_id(); }
//inline std::ostream& operator<<( std::ostream& s, leafset_p leafset ) { s << "(" << leafset << ")"; return s; }

}

#endif //  __INCLUDE_COSI_LEAFSET_SIZEONLY_H
