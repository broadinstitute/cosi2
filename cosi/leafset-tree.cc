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
#include <cosi/leafset-tree.h>

namespace cosi {
namespace leafset_tree {

#define ForEach BOOST_FOREACH
using util::chkCond;
	
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

leafset_p leafset_intersection( leafset_p leafset1, leafset_p leafset2 ) {
	using namespace std;
	vector<leaf_id_t> leaves1, leaves2;
	leafset_get_leaves( leafset1, back_inserter( leaves1 ) );
	sort( leaves1.begin(), leaves1.end() );
	leafset_get_leaves( leafset2, back_inserter( leaves2 ) );
	sort( leaves2.begin(), leaves2.end() );
	vector<leaf_id_t> inters;
	set_intersection( leaves1.begin(), leaves1.end(), leaves2.begin(), leaves2.end(),
										back_inserter( inters ) );
	leafset_p result = make_empty_leafset();
	ForEach( leaf_id_t leaf, inters )
		 result = leafset_union( result, make_singleton_leafset( leaf ) );
	return result;
}


leafset_p leafset_difference( leafset_p leafset1, leafset_p leafset2 ) {
	using namespace std;
	vector<leaf_id_t> leaves1, leaves2;
	leafset_get_leaves( leafset1, back_inserter( leaves1 ) );
	sort( leaves1.begin(), leaves1.end() );
	leafset_get_leaves( leafset2, back_inserter( leaves2 ) );
	sort( leaves2.begin(), leaves2.end() );
	vector<leaf_id_t> inters;
	set_difference( leaves1.begin(), leaves1.end(), leaves2.begin(), leaves2.end(),
										back_inserter( inters ) );
	leafset_p result = make_empty_leafset();
	ForEach( leaf_id_t leaf, inters )
		 result = leafset_union( result, make_singleton_leafset( leaf ) );
	return result;
}




bool leafset_equal( leafset_p leafset1, leafset_p leafset2 ) {
	using namespace std;
	vector<leaf_id_t> leaves1, leaves2;
	leafset_get_leaves( leafset1, back_inserter( leaves1 ) );
	sort( leaves1.begin(), leaves1.end() );
	leafset_get_leaves( leafset2, back_inserter( leaves2 ) );
	sort( leaves2.begin(), leaves2.end() );
	return leaves1 == leaves2;
}



void leafset_struct::operator delete (void *p) {
  mempool_leafsets.mempool_free_item( p );
}
leafset_p leafset_from_str( const char *s ) {
	using util::cosi_strtok_r;
	
	PRINT2( "leafset_from_str", s );
  leafset_p leafset = LEAFSET_NULL;
  char *str = util::cosi_strdup( s );
  char *saveptr;
  for ( char *token = cosi_strtok_r( str, ",", &saveptr ); token; token = cosi_strtok_r( NULL, ",", &saveptr ) ) {
	 leaf_id_t leaf;
	 int num_scanned = sscanf( token, "%d", &leaf );
	 chkCond( num_scanned == 1, "leafset_from_str: bad set representation - %s", s );
	 leafset = leafset_union( leafset, make_singleton_leafset( leaf ) );
  }
  free( str );
	PRINT( leafset );
  /*printf( "leafset from str: %s is ", s ); leafset_print( leafset );*/
  return leafset;
}

const char *leafset_str( leafset_p /*leafset*/ ) {
  return "unimpl";
}

ostream& operator<<( std::ostream& s, leafset_p leafset ) {
	if ( leafset_is_empty( leafset ) ) { s << "{}"; return s; }
	s << "{(" << leafset_size( leafset ) << ") ";
	vector<leaf_id_t> leaves;
	leafset_get_leaves( leafset, back_inserter( leaves ) );
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

/* FuncP: make_range_leafset */
/* Returns the leafset containing the given range of leaves. */
leafset_p make_range_leafset( leaf_id_t fromLeaf, leaf_id_t toLeaf ) {
	leafset_p result = make_empty_leafset();
	for( leaf_id_t leaf = fromLeaf; leaf < toLeaf; ++leaf )
		 result = leafset_union( result, make_singleton_leafset( leaf ) );
	return result;
}

}  // namespace leafset_tree
}  // namespace cosi
