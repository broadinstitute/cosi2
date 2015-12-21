/*
 * File: leafset.c
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 */

#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <boost/function_output_iterator.hpp>
#include <cosi/general/utils.h>
#include <cosi/general/mempool.h>
#include <cosi/leafset-counts.h>

namespace cosi {
namespace leafset_counts {

#define ForEach BOOST_FOREACH
using util::chkCond;


const std::vector< popid > *leafset_leaf2popName = NULL;
const std::vector<pop_idx_t> *leafset_popname2idx = NULL;
int leafset_npops = 0;
	
static Mempool mempool_leafsets( sizeof( leafset_struct ), 65536, 4096 );

/*static void leafset_clear( leafset_p  );*/

/* Static var: max_leaf_id */
/* The maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
static leaf_id_t max_leaf_id = NULL_LEAF_ID;

/* Func: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ ) {
  max_leaf_id = max_leaf_id_;
}  /* leafset_set_max_leaf_id() */

leaf_id_t leafset_get_max_leaf_id(void) { return max_leaf_id; }

void* leafset_struct::operator new (size_t /*size*/) {
  return mempool_leafsets.mempool_alloc_item();
}


void leafset_struct::operator delete (void *p) {
  mempool_leafsets.mempool_free_item( p );
}

ostream& operator<<( std::ostream& s, leafset_p leafset ) {
	if ( leafset_is_empty( leafset ) ) { s << "{}"; return s; }
	s << "{(" << leafset_size( leafset ) << ") ";
	vector<leaf_id_t> leaves;
//	leafset_get_leaves( leafset, back_inserter( leaves ) );
	std::sort( leaves.begin(), leaves.end() );
	std::copy( leaves.begin(), leaves.end(), std::ostream_iterator<leaf_id_t>( std::cout, " " ) );
	s  << "}";
	return s;
}

#ifdef COSI_R2
void leafset_struct::computeLeaves() const {
	if ( !leaves ) {
		leaves = boost::scoped_ptr< vector< leaf_id_t > >( new vector< leaf_id_t >() );
		leafset_get_leaves( leafset_p( (leafset_t *)this ), std::back_inserter( *leaves ) );
		std::sort( leaves->begin(), leaves->end() );
	}
}

namespace {
template <typename CountT, typename ValT>
struct counter_struct {
	 CountT& count;

	 counter_struct( CountT& count_ ): count( count_ ) { count = 0; }
	 void operator()( ValT ) { count++; } 
};
}

// Func: compute_r2
// Compute the r^2 measure of linkage disequilibrium between the two leafsets.
cosi_double leafset_struct::compute_r2( leafset_p leafset1, leafset_p leafset2 ) {
	leafset1->computeLeaves();
	leafset2->computeLeaves();

	nchroms_t num_ab;
	std::set_intersection( leafset1->leaves->begin(), leafset1->leaves->end(),
												 leafset2->leaves->begin(), leafset2->leaves->end(),
												 boost::make_function_output_iterator( counter_struct< nchroms_t, leaf_id_t >( num_ab ) ) );
//	PRINT3( leafset1->leaves->size(), leafset2->leaves->size(), num_ab );

	nchroms_t numIndivs = max_leaf_id;
	nchroms_t num_a = leafset1->leaves->size();
	nchroms_t num_b = leafset2->leaves->size();


	cosi_double p_a = ((cosi_double) num_a ) / ((cosi_double) numIndivs );
	cosi_double p_b = ((cosi_double) num_b ) / ((cosi_double) numIndivs );
	cosi_double p_ab = ((cosi_double) num_ab ) / ((cosi_double) numIndivs );

	cosi_double num_sqrt = p_ab - p_a * p_b;
	return num_sqrt * num_sqrt / ( p_a * ( 1.0 - p_a ) * p_b * ( 1.0 - p_b ) );
}  // leafset_struct::compute_r2
#endif // ifdef COSI_R2

}  // namespace leafset_counts
}  // namespace cosi

