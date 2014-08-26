//
// Header: leafset-sizeonly.cc
//
// An implementation of <leafsets> which keeps only the size of the leafset,
// and not the identity of the leaves in it.  Can be used when we're only interested
// in the frequency of each mutation (as when computing the allele frequency spectrum).
//

#include <cosi/leafset-sizeonly.h>

namespace cosi {

namespace {
leaf_id_t max_leaf_id = NULL_LEAF_ID;
}

/* FuncP: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
/* Must be called once, prior to the creation of any leafsets. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ ) {
	max_leaf_id = max_leaf_id_;
}

leaf_id_t leafset_get_max_leaf_id(void) { return max_leaf_id; }

}

