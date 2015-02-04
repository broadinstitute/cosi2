/* $Id: node.c,v 1.4 2011/06/07 15:29:25 sfs Exp $ */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <boost/utility/swap.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <cosi/node.h>
#include <cosi/nodelist.h>
#include <cosi/seglist.h>
#include <cosi/utils.h>
#include <cosi/mempool.h>
#include <cosi/recomb.h>
#include <cosi/genmap.h>
#include <cosi/pop.h>
#include <cosi/hullmgr.h>

unsigned long tot_recombs = 0;

namespace cosi {

namespace node {

using util::chkCond;
using util::isize;
	
NodePool::NodePool(void):
  node_index(0),
	mempool_node( sizeof( Node ), /* block_size_in_items= */ 8192,
								/* initial_freelist_size_in_items= */ 512 ),
	allNodes( 1, (Node *)NULL ),
	outputMuts( True ),
	enableGeneConv( True )
{
	allNodes.reserve( 16385 );
}

#ifdef COSI_STATS
unsigned long tot_coal = 0, tot_recomb = 0, tot_recomb_attempt = 0, tot_gc = 0, tot_gc_attempt = 0, tot_breakoff_attempt = 0,
tot_breakoff = 0, tot_alloced = 0;
#endif



Node *NodePool::alloc_node() {
  Node *n = (Node *)mempool_node.mempool_alloc_item();
  n->name = node_index++;
	n->pop = NULL;
	n->recombRate = ZERO_GLEN;
	n->gcRate = ZERO_GLEN;
	n->idx_in_allNodes = allNodes.size();
	allNodes.push_back( n );
  return n;
}

void NodePool::node_delete( Node **np ) {
  Node *n = *np; *np = NULL;

	assert( allNodes.size() >= 2 );

	setNodeRecombRate( n, ZERO_GLEN );
	setNodeGcRate( n, ZERO_GLEN );
	int idx = n->idx_in_allNodes;
	if ( idx < isize(allNodes)-1 ) {
		Node *last_node = allNodes.back();
		glen_t last_node_recomb_rate = last_node->recombRate;
		glen_t last_node_gc_rate = last_node->gcRate;
		setNodeRecombRate( last_node, ZERO_GLEN );
		setNodeGcRate( last_node, ZERO_GLEN );
		allNodes[ idx ] = last_node;
		last_node->idx_in_allNodes = idx;
		setNodeRecombRate( last_node, last_node_recomb_rate );
		setNodeGcRate( last_node, last_node_gc_rate );
	}
	allNodes.pop_back();

  seglist_delete( &( n->segs ) );
	
  mempool_node.mempool_free_item( n );
}

Node* 
NodePool::make_new_leaf ( ) 
{
  Node *newnodeptr = alloc_node();
#ifdef COSI_FREQONLY
  leafset_p leafset = make_singleton_leafset( newnodeptr->name );
#else
  leafset_p leafset = make_singleton_leafset( outputMuts ? newnodeptr->name : ((leaf_id_t)1) );
#endif	

  newnodeptr->segs = seglist::seglist_make_full( leafset, /* lastCoalGen= */ genid(0) );

  newnodeptr->gen = ZERO_GEN;

	getHooks()->fire_make_leaf( newnodeptr->name );

  node_chk( newnodeptr );

  return newnodeptr;
}

Node * 
NodePool::make_empty_node (genid gen) 
{
  Node *newnodeptr = alloc_node();
  newnodeptr->gen = gen;
  return newnodeptr;
}

/*
 * Func: node_coalesce
 *
 * Compute the coalescence of two nodes.  The two child nodes are destroyed.
 *
 * Input params:
 *
 *    node1p, node2p - the child nodes.   They're destroyed by this call.
 *    gen - the generation of the new coalesced node
 *
 * Returns:
 *
 *    The coalesced node.
 */
Node * 
NodePool::node_coalesce (Node ** node1p, Node ** node2p, genid gen) 
{
  Node *node1 = *node1p, *node2 = *node2p;
  *node1p = *node2p = NULL;

  node_chk( node1 );
  node_chk( node2 );

  IF_COSI_STATS( tot_coal++ ); 

  seglist_chk( node1->segs );
  seglist_chk( node2->segs );
	//PRINT2( "coal", gen );

	nodeid next_node_name = node_index;

	
	getHooks()->fire_add_edge( node1->name, next_node_name, node1->gen, gen, node1->segs, EDGE_COAL );
	getHooks()->fire_add_edge( node2->name, next_node_name, node2->gen, gen, node2->segs, EDGE_COAL );

	Node *newnodeptr = NULL;
														 
  Seglist *segs = seglist_union( &( node1->segs ), &( node2->segs ), gen
																 //IF_COSI_TRACK_LAST_COAL(, /* seg_union_callback= */ boost::bind( &Hooks::fire_tree_edge, hooks, _1, _2, _3, _4, _5 ) )
																 );

	if ( seglist_is_empty( segs ) ) {
		 seglist_delete( &segs );
		 node_index++;
	}
	else {
		newnodeptr = make_empty_node(gen);
		newnodeptr->segs = segs;
		setNodeRecombRate( newnodeptr, compute_node_recomb_rate( newnodeptr ) );
		setNodeGcRate( newnodeptr, compute_node_gc_rate( newnodeptr ) );
	}
	
  node_delete( &node1 );
  node_delete( &node2 );

  if ( newnodeptr ) { node_chk( newnodeptr ); }

  return newnodeptr;
}  // NodePool::node_coalesce()

/*
 * Func: node_recombine
 *
 * Computes the recombination of a node.  If the recombination point falls
 * inside the node's seglist, returns two parent nodes and destroys the 
 * child node.  If the recombination point falls completely to the left
 * or completely to the right of the node's seglist (i.e. if one of the
 * parent nodes resulting from the recombination would have an empty seglist),
 * the child node is left unchanged.
 *
 * Input params:
 *
 *    nodep - the node to be recombined
 *    gen - the generation of the two new parent nodes
 *    loc - the location of the recombination.
 *
 * Output params:
 *
 *    newnode1, newnode2 - the parent nodes resulting from the recombination.
 *
 * Returns:
 *
 *    The number of parent nodes of the recombination, 1 or 2.  If 2,
 *    the input node has been destroyed; if 1, it is left untouched.
 */
int
NodePool::node_recombine (Node **nodep, Node **newnode1, Node **newnode2, 
													genid gen, loc_t loc, Node **nodes_out) 
{
  Node *node = *nodep;

  node_chk( node );

  int node_idx_in_pop = node->get_idx_in_pop();
	Pop *popptr = node->getPop();

	if ( nodes_out ) { nodes_out[0] = node; nodes_out[1] = nodes_out[2] = (Node *)NULL; }

  IF_COSI_STATS( tot_recomb_attempt++ );

	assert( !seglist_is_empty( node->segs ) );
	//assert( seglist_beg( node->segs ) <= loc && loc <= seglist_end( node->segs ) );
	//PRINT2( "recomb", gen );

  if ( loc <= seglist_beg( node->segs ) ) {
		if ( nodes_out ) nodes_out[2] = node;
		*newnode1 = NULL;
		*newnode2 = node;
		return 1;
	} else if ( loc >= seglist_end( node->segs ) ) {
		if ( nodes_out ) nodes_out[1] = node;
		*newnode1 = node;
		*newnode2 = NULL;
		return 1;
	}
	
	tot_recombs++;

  IF_COSI_STATS( tot_recomb++ );

  *newnode1 = make_empty_node (gen);
  *newnode2 = make_empty_node (gen);

  seglist_split( &( node->segs ), loc, &( (*newnode1)->segs ), &( (*newnode2)->segs ), /* split_seg= */ True );

	setNodeRecombRate( *newnode1, compute_node_recomb_rate( *newnode1 ) );
	setNodeRecombRate( *newnode2, compute_node_recomb_rate( *newnode2 ) );

	setNodeGcRate( *newnode1, compute_node_gc_rate( *newnode1 ) );
	setNodeGcRate( *newnode2, compute_node_gc_rate( *newnode2 ) );

	getHooks()->fire_add_edge( node->name, (*newnode1)->name, node->gen, gen, (*newnode1)->segs,
														 EDGE_RECOMB );
	getHooks()->fire_add_edge( node->name, (*newnode2)->name, node->gen, gen, (*newnode2)->segs,
														 EDGE_RECOMB );
	
  node_delete( nodep );

  node_chk( (*newnode1) );
  node_chk( (*newnode2) );

	if ( nodes_out ) { nodes_out[1] = *newnode1; nodes_out[2] = *newnode2; }

	/* STEP 4 */
	popptr->pop_remove_node_by_idx (node_idx_in_pop);

	 
	/* STEP 5 */
	popptr->pop_add_node (*newnode1);
	popptr->pop_add_node (*newnode2);

	return 2;
	
}  // NodePool::node_recombine()



/* GENE CONVERSION */

/*
 * Func: node_gc
 *
 * Computes the gene conversion of a node.  If both resulting parents of the gene
 * conversion node have non-empty seglists, destroys the child node passed in;
 * otherwise, does not touch the child node.
 *
 * Input params:
 *
 *    nodep - the node to be recombined
 *    gen - the generation of the two new parent nodes
 *    loc, locend - the location of the "inside" segment of the gene conversion.
 *
 * Output params:
 *
 *    newnode1, newnode2 - the parent nodes resulting from the gene conversion
 *
 * Returns:
 *
 *    The number of parent nodes of the gene conversion, 1 or 2.  If 2,
 *    the input node has been destroyed; if 1, it is left untouched.
 */
int
NodePool::node_gc (Node **nodep, Node **newnode1, Node **newnode2, 
									 genid gen, loc_t loc, loc_t locend) 
{
  Node *node = *nodep;

  node_chk( node );

  IF_COSI_STATS( tot_gc_attempt++ );

  Seglist *inside = NULL, *outside = NULL;

  if ( !seglist_intersect_for_gc( &( node->segs ), loc, locend, &inside, &outside,
																	/* for_sentinels= */ False ) )
		 return 1;
  
  IF_COSI_STATS( tot_gc++ );
  *newnode1 = make_empty_node (gen);
  *newnode2 = make_empty_node (gen);
  
  (*newnode1)->segs = inside;
  (*newnode2)->segs = outside;

	setNodeRecombRate( *newnode1, compute_node_recomb_rate( *newnode1 ) );
	setNodeRecombRate( *newnode2, compute_node_recomb_rate( *newnode2 ) );

	setNodeGcRate( *newnode1, compute_node_gc_rate( *newnode1 ) );
	setNodeGcRate( *newnode2, compute_node_gc_rate( *newnode2 ) );

	getHooks()->fire_add_edge( node->name, (*newnode1)->name, node->gen, gen, (*newnode1)->segs, EDGE_GC );
	getHooks()->fire_add_edge( node->name, (*newnode2)->name, node->gen, gen, (*newnode2)->segs, EDGE_GC );

  node_delete( nodep );

  seglist_chk( (*newnode2)->segs );
  
  node_chk( (*newnode1) );
  node_chk( (*newnode2) );
  
  return 2;
}  /* node_gc() */

void NodePool::setGenMap( GenMapP genMap_ ) { genMap = genMap_; }

void NodePool::finishLeaves() {
	recombPartialSums.ensureCapacity( allNodes.size() );
	for ( idx_t leafId = 1; leafId < isize( allNodes ); leafId++ ) {
		setNodeRecombRate( allNodes[ leafId ], MAX_GLOC - MIN_GLOC );
		setNodeGcRate( allNodes[ leafId ], MAX_GLOC - MIN_GLOC );
	}
}

// Method: findRecomb
//
// Find the node to recombine and the loc at which to recombine it.
//
// Input params:
//
//    recombFrac - fraction of the total recomb rate
//
// Output params:
//
//    loc - location of the recombination
//
// Returns:
//
//    index in <allNodes> of the node to recombine.
//    loc is guaranteed to be strictly between the beg and end of
//    this node's seglist.
//   
Node *NodePool::findRecomb( frac_t recombFrac, loc_t *loc ) const {
#ifndef COSI_SLOW_RECOMB	
	glen_t recombAmt;
	idx_t nodeIdx = recombPartialSums.findCumulativeFraction( recombFrac, &recombAmt );
	Node *n = allNodes[ nodeIdx ];
	loc_t leftLoc = seglist_beg( n->segs );
	// if ( selPos )
	// 	 util::updateMin( leftLoc, *selPos );
	
	gloc_t leftGdPos = genMap->getGdPos( leftLoc );
	*loc = genMap->getLoc( leftGdPos + recombAmt );
	return n;
#else	// COSI_SLOW_RECOMB defined
	chkCond( recombFrac < 1.0 );
	cosi_double intPart;
	*loc = loc_t( modf( recombFrac * (allNodes.size()-1), &intPart ) );
	//PRINT4( recombFrac, (allNodes.size()-1), *loc, intPart );
	return allNodes.at( ((int)intPart) + 1 );
#endif // #ifdef COSI_SLOW_RECOMB
	
}  // NodePool::findRecomb()

//
// Method: findGcOrigin
//
// Find a node at which to do gene conversion, and the location of <gene conversion origin>.
//
// Called from <GeneConversion::gc_execute()>.  
//
// Params:
//
//    gcFrac - value indicating in what node, and at what location in that node, to
//      initiate the gene conversion.  More specifically: for each Node we keep the
//      probability of initiating a gene conversion within that node; frac is a
//      fraction of the total sum of these per-node probabilities.
//
//    tractLenLeft/tractLenRight - the distance from the initiation point of the gene conversion tract
//      to its left/right endpoint
//
Node *NodePool::findGcOrigin( frac_t gcFrac, gloc_t *gcOrigin ) const {
	//
	// Determine the origination point of the gene conversion.
	//

	// Var: gcAmt
	// Where within the chosen node will we initiate the gene conversion?
	glen_t gcAmt;
	idx_t nodeIdx = gcPartialSums.findCumulativeFraction( gcFrac, &gcAmt );

	Node *n = allNodes[ nodeIdx ];

	*gcOrigin = seglist_find_glen( n->segs, gcAmt );
	assert( MIN_GLOC <= *gcOrigin && *gcOrigin <= MAX_GLOC );
	return n;
}  // NodePool::findGcOrigin()


// Method: compute_node_recomb_rate
// Returns the prob of recomb between beg and end of the node's segs.
glen_t NodePool::compute_node_recomb_rate( Node *n ) const {
	loc_t segs_end = seglist_end( n->segs );
	loc_t segs_beg = seglist_beg( n->segs );
	
	// if ( selPos ) {
	// 	//PRINT4( "bef", *selPos, segs_beg, segs_end );
	// 	util::updateMin( segs_beg, *selPos );
	// 	util::updateMax( segs_end, *selPos );
	// 	//PRINT4( "aft", *selPos, segs_beg, segs_end );
	// }
	gloc_t gd_end = genMap->getGdPos( segs_end );
	gloc_t gd_beg = genMap->getGdPos( segs_beg );
	PRINT2( gd_beg, gd_end );
	return gd_end - gd_beg;
}

// Method: compute_node_gc_rate
// Returns the prob of a gene conversion within this node, such that the
// gene conversion tract and its complement would split this node's <Node::segs>
// into two non-empty seglists.
glen_t NodePool::compute_node_gc_rate( Node *n ) const {
	if ( !enableGeneConv ) return ZERO_GLEN;
	return seglist_tot_glen( n->segs );
}


/*
 * Section: Profiling and debugging
 *
 * Code for profiling and debugging; is inactive during normal execution.
 */

void node_print_stats() {
#ifdef COSI_STATS
  extern unsigned long leafsets_alloced, miss_inters, ok_inters;
  printf( "tot_alloced=%lu tot_coal=%lu tot_recomb=%lu tot_recomb_attempt=%lu tot_gc=%lu tot_gc_attempt=%lu tot_breakoff_attempt=%lu tot_breakoff=%lu leafsets_alloced=%lu\n",
					tot_alloced, tot_coal, tot_recomb, tot_recomb_attempt, tot_gc, tot_gc_attempt, tot_breakoff_attempt, tot_breakoff,
					leafsets_alloced );

  extern int n_seglist_empty, n_left_or_right, n_containing, n_in_gap, n_nonempty, n_segsumm, n_oneseg;

  printf( "n_seglist_empty=%d n_left_or_right=%d n_containing=%d n_in_gap=%d n_nonempty=%d n_segsumm=%d n_oneseg=%d\n",
					n_seglist_empty, n_left_or_right, n_containing, n_in_gap, n_nonempty, n_segsumm, n_oneseg );

  extern clock_t time_gc, time_coal_union;

  printf( "time_gc=%f\n", ((double)time_gc) / CLOCKS_PER_SEC );  
  printf( "time_coal_union=%f\n", ((double)time_coal_union) / CLOCKS_PER_SEC );
  
  extern unsigned long num_reused;
  printf( "num_reused=%lu\n", num_reused );

#if 0	
  extern Mempool mempool_seglist;
  extern Mempool *mempool_seg;
  extern Mempool mempool_header;

  printf( "max live seglists = %ld\n", mempool_seglist.max_live_items );
  for ( level_t level = seglist_min_level; level <= seglist_max_level; level++ )
		 printf( "level=%d max live items=%ld\n", level, mempool_seg[ level ].max_live_items );
  printf( "max live nodes=%ld\n", mempool_node.max_live_items );
  printf( "max live headers=%ld\n", mempool_header.max_live_items );
#endif	

  printf( "tot nodes alloced=%d\n", node_index );

  extern int level_hist[1024];

  printf( "\n\nlevel histogram:\n" );
  for ( level_t level = seglist_min_level; level <= seglist_max_level; level++ )
		 printf( "level=%d count=%d\n", level, level_hist[ level ] );
  
#endif  
}

void NodePool::node_chk_helper( const Node * /*n*/, const char * /*fname*/, int /*line*/  ) {
  /*seglist_chk_all();*/
  /*printf( "checking node %d\n", n->name );*/

#if 0	
  seglist_chk( n->segs );
  seglist_chk( n->sentinel_for );
#endif	
}

ostream& operator<<( ostream& s, const Node *n ) {
	s << "[node: " << n->getName() << " gen: " << n->getGen() << " segs: " << n->getSegs() << "]";
	return s;
}


#ifdef COSI_DEV_NODE_DEBUG
void node_chk_all( void ) {
#if 0	
  for ( Node *node = all_nodes; node; node = node->next ) {
		assert( !node->next || node->next->prev == node );
		assert( !node->prev || node->prev->next == node );
		assert( node->segs );
		seglist_set_node( node->segs, NULL );
		assert( node->sentinel_for );
		seglist_set_node( node->sentinel_for, NULL );
  }
  for ( Node *node = all_nodes; node; node = node->next ) {
		assert( !seglist_get_node( node->segs ) );
		seglist_set_node( node->segs, node );
		assert( !seglist_get_node( node->sentinel_for ) );
		seglist_set_node( node->sentinel_for, node );
  }
#endif	
}
#endif

int
nodelist_add(NodeList *nlptr, Node *nodeptr) {
  nlptr->push_back( nodeptr );
  return nlptr->size()-1;
}

void nodelist_remove_idx( NodeList *nlptr, int idx ) {
  if ( idx < ((int)nlptr->size())-1 ) {
		Node *last_node = nlptr->back();
		(*nlptr)[ idx ] = last_node;
		last_node->idx_in_pop = idx;
  }
  nlptr->pop_back();
}  

} // namespace node
}  // namespace cosi
