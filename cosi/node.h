/* $Id: node.h,v 1.4 2011/06/07 15:29:25 sfs Exp $ */

/**
 * Header: node.h
 *
 * Defines the <Node> struct, a node of the <ARG>, and related operations.
 */

#ifndef __INCLUDE_COSI_NODE_H
#define __INCLUDE_COSI_NODE_H

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/intrusive/parent_from_member.hpp>
#include <boost/optional.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/utils.h>
#include <cosi/leafset.h>
#include <cosi/mempool.h>
#include <cosi/cosirand.h>
#include <cosi/hooks.h>
#include <cosi/nodefwd.h>
#include <cosi/hullmgr.h>

namespace cosi {

using std::ostream;

class Pop;

namespace node {

/**
 * Struct: Node
 *
 * A node of the ARG.  Represents a chromosome instance that existed at a particular point in time.
 *
 * Nodes are assigned to <populations>.  During the backwards-in-time simulation in <Simulator::sim_execute()>, nodes can be moved between
 * populations.
 *
 * Conceptually, each node of the ARG has children (node(s) that inherited genetic material from it) and parents (node(s) from which
 * this it inherited genetic material; because of recombination and gene conversion a node may have two parents).
 * Since during simulation we only keep the "frontier" of the ARG (nodes that do not yet have a parent), the Node data structure
 * does not keep pointers to child or parent nodes.
 *
 * *Invariants*
 *
 *   - all genetic material in a Node gets inherited by its children.
 *   - everything a node has, it gets from its parent(s)
 *   - right now, parents of a recombined node have only that node as a child.
 *     so you can get the material inherited along an edge.
 *
 */
class Node {
public:
	 typedef seglist::Seglist Seglist;
	 
	 nodeid getName() const { return name; }
	 genid getGen() const { return gen; }
	 Seglist *getSegs() const { return segs; }
	 Pop *getPop() const { return pop; }
	 
	 int get_idx_in_pop() const { return idx_in_pop; }
	 void set_idx_in_pop( int idx_in_pop_ ) { idx_in_pop = idx_in_pop_; }

private:
	 /* Field: name */
	 /* A unique id of the node. */
	 nodeid name;

	 /* Field: gen */
	 /* The time point (generation) at which this chromosome existed. */
	 genid gen;

	 /* Field: idx_in_pop */
	 /* Index of this node in its population.  Used when removing the node */
	 /* from its population's nodelist. */
	 int idx_in_pop;

	 /* Field: segs */
	 /* The parts of this chromosome inherited by _some_ present-day (sample) chromosome.  */
	 Seglist *segs;

	 // Field: pop
	 // The population to which the node belongs.
	 Pop *pop;

	 // Field: idx_in_allNodes
	 // Index of this node in <NodePool::allNodes>
	 int idx_in_allNodes;

	 // Field: recombRate
	 // Total recomb rate within this node.  Equal to the genetic distance between
	 // seglist_beg(segs) and seglist_end(segs).
	 glen_t recombRate;

	 // Field: gcRate
	 // Total gc rate within this node.  Equal to the sum of the genetic lengths of
	 // the segments of this node's <Seglist>, <segs>.
	 glen_t gcRate;

#ifdef COSI_SUPPORT_COALAPX	 
	 HullMgr::Hull hull;
#endif	 

	 friend class NodePool;
	 friend void nodelist_remove_idx( NodeList *nlptr, int idx );

public:

	 // Class: PopAccess
	 // Lets methods in class Pop set <Node::pop>
	 class PopAccess {
			friend class ::cosi::Pop;
			static void SetNodePop( Node *n, Pop *p ) { n->pop = p; }
#ifdef COSI_SUPPORT_COALAPX			
			static HullMgr::Hull* GetNodeHullPtr( Node *n ) { return &( n->hull ); }

			static Node *GetNodeFromHullPtr( const HullMgr::Hull *hull ) {
				return boost::intrusive::get_parent_from_member< Node >( const_cast< HullMgr::Hull * >( hull ),
																																 &Node::hull );
			}
					 
#endif			
	 };
	 
	 friend class Node::PopAccess;
};  // class Node

#ifndef COSI_DEV_NODE_DEBUG
#define node_chk_all()
#else
void node_chk_all( void );
#endif


#ifndef COSI_DEV_NODE_DEBUG
#define node_chk(n)
#else
#define node_chk(n) node_chk_helper(n,__FILE__,__LINE__)
#endif

void node_print_stats(void);

void node_chk_helper( const Node *, const char *fname, int line );

//
// Class: NodePool
//
// Keeps track of all existing <Nodes>, in all <Pops>, during a simulation.
// Deals with node memory management.  Implements the basic operations on
// the nodes: coalescence, recombination, gene conversion.  (Of course,
// the most complex part of these operations is performing the corresponding
// operation on the nodes' <Seglists>; see the module <seglist.h> for that.)
//
class NodePool: public HasRandGen, public Hookable {
public:

	 //
	 // Method group: Initialization
	 //
	
	 NodePool( );

	 void setGenMap( GenMapP );
	 void setEnableGeneConv( bool_t enableGeneConv_ ) { enableGeneConv = enableGeneConv_; }
	 void setOutputMuts( bool_t outputMuts_ ) { outputMuts = outputMuts_; }
	 void setMaxCoalDist( plen_t maxCoalDist_ );

	 //
	 // Method group: Manipulating <Nodes>
	 //
	 
	 Node * make_new_leaf (void);
	 void finishLeaves();
	 
	 Node * node_coalesce (Node **, Node **, genid  gen);
	 
	 /* returns the number of recombinant ancestors */
	 int node_recombine( Node**, Node**, Node**, genid  gen, loc_t loc, Node **nodes_out);
	 int
	 node_gc (Node **  node , Node **  newnode1 , Node **  newnode2 , 
						genid  gen , loc_t  loc , loc_t locend);

	 //
	 // Method group: Keeping track of per-node recomb and gc rates
	 //
	 // This information is used to for incrementally updating the region-wide probability of a recomb or
	 // gc event, and -- when a recomb or gc event is chosen to happen -- to quickly find the node
	 // undergoing the event and the location within the node.
	 //

	 // MethodP: getAllNodesRecombRate
	 // Returns the total rate of recombination over all nodes, conditioned on the recombination
	 // falling between start and end of a node's seglist.
	 // Complexity: constant-time.
	 glen_t getAllNodesRecombRate() const {
#ifndef COSI_SLOW_RECOMB		 
		 return recombPartialSums.getTotalSum();
#else
		 return glen_t( allNodes.size()-1 );
#endif		 
	 }

	 glen_t getAllNodesGeneConvRate() const { return gcPartialSums.getTotalSum(); }
	 
	 nchroms_t getNumNodes() const { return allNodes.size()-1; }

	 Node *pickRandomNode() const { return allNodes[ 1 + ((int)(random_double() * (allNodes.size()-1))) ]; }

	 // MethodP: findRecomb
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
	 Node *findRecomb( frac_t recombFrac, loc_t *loc ) const;

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
	 //
	 Node *findGcOrigin( frac_t gcFrac, gloc_t *gcOrigin ) const;

//	 void setSelPos( boost::optional<loc_t> selPos_ ) { this->selPos = selPos_; }
	 GenMapP getGenMap() const { return genMap; }

private:
	 
	 /* Field: node_index */
	 /* Unique index (name) to be given to the next <Node> we allocate. */
	 int node_index;
	 
	 /* Field: mempool_node */
	 /* <Mempool> from which <Node>s are allocated. */
	 Mempool mempool_node;

	 // Field: allNodes
	 // All currently active nodes. NOTE: 0th element of this vector is unused!
	 vector<Node *> allNodes;

	 // Field: recombPartialSums
	 // Partial sum tree of recomb rates.  Keeps the
	 // recomb rate for each node in <allNodes> (identified by
	 // <Node::idx_in_allNodes>), as well as partial sums for subsets of nodes.
	 util::PartialSumTree<glen_t> recombPartialSums;

	 // Field: gcPartialSums
	 // Partial sum tree of recomb rates.  Keeps the
	 // gc rate for each node in <allNodes> (identified by
	 // <Node::idx_in_allNodes>), as well as partial sums for subsets of nodes.
	 util::PartialSumTree<glen_t> gcPartialSums;

	 // Field: genMap
	 // The genetic map
	 GenMapP genMap;

	 // Field: outputMuts
	 // Will we be outputting the full haplotypes?  (If not, can save computation.)
	 bool_t outputMuts;

	 // Field: enableGeneConv
	 // If true, we keep track of gene conversion rate for each node;
	 // if false, we don't.
	 bool_t enableGeneConv;

//	 boost::optional<loc_t> selPos;

	 //
	 // Private methods
	 //

	 Node *alloc_node();
	 Node * make_empty_node (genid  gen);
	 void node_delete( Node ** );

	 // MethodP: node_recomb_rate
	 // Returns the prob of recomb between beg and end of the node's segs.
	 glen_t compute_node_recomb_rate( Node * ) const;

	 
	 // Method: compute_node_gc_rate
	 // Returns the prob of a gene conversion within this node, such that the
	 // gene conversion tract and its complement would split this node's <Node::segs>
	 // into two non-empty seglists.
	 glen_t compute_node_gc_rate( Node * ) const;

	 void setNodeRecombRate( Node *n, glen_t recombRate ) {
		 recombPartialSums.add( n->idx_in_allNodes, recombRate - n->recombRate );
		 n->recombRate = recombRate;
	 }
	 void setNodeGcRate( Node *n, glen_t gcRate ) {
		 if ( enableGeneConv ) {
			 gcPartialSums.add( n->idx_in_allNodes, gcRate - n->gcRate );
			 n->gcRate = gcRate;
		 }
	 }

	 void node_chk_helper( const Node *n, const char *fname, int line  );

};  // class NodePool

ostream& operator<<( ostream&, const Node *n );

} // namespace node

}  // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_NODE_H
