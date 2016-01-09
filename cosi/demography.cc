/* $Id: demography.c,v 1.8 2011/06/08 11:23:54 sfs Exp $ */

/*
 * File: demography.c
 *
 * Data structures and methods for representing the recombination sites, and the <populations>.
 */

#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <cassert>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <cosi/general/utils.h>
#include <cosi/defs.h>
#include <cosi/pop.h>
#include <cosi/node.h>
#include <cosi/seglist.h>
#include <cosi/demography.h>
#include <cosi/leafset.h>
#include <cosi/seglist.h>
#include <cosi/mutate.h>
#include <cosi/hooks.h>

namespace cosi {

#define ForEach BOOST_FOREACH

using util::chkCond;
using util::isize;

Demography::Demography():
  LOGGING( False ), logfp( NULL ), verbose( True ), maxCoalDist( plen_t( 1.0 ) ), maxCoalDistCvxHull( False ) {

}

void Demography::dg_setMutate( MutateP mutate_ ) { dg_mutate = mutate_; mutate_->setDemography( this ) ; }
void Demography::dg_setNodePool( NodePoolP nodePool_ ) { dg_nodePool = nodePool_; }

/* ERROR HANDLING */
static void dg_exit (const char *funct_name, const char *error_string);
static void dg_error_nonfatal (const char *funct_name, const char *error_string);

/* POP FUNCTIONS 
 * "index" refers to the location in the population list.
 * "name" refers to the numerical name of the population, as
 *        defined in the param file.
 * 
 * How to create a population:
 * 1. call dg_create_pop to create an empty population
 * 2. call dg_set_pop.. to set the population size
 * 3. (opt) call dg_populate_by_name to fill in the nodes
 */
void
Demography::dg_create_pop (popid popname, const std::string& label, genid gen) {
  Pop *temppop;
	
  /* check to make sure the popname is not taken and
		 is non negative
		 if it is, throw a fatal exception.  */
  if (is_null( popname )) {
		dg_exit("dg_create_pop","popname must not be null");
  }
  temppop = dg_get_pop_by_name(popname);

  if (temppop != NULL) {
		dg_exit("dg_create_pop","duplicate popname used.");
  }

	nchroms_t zero_chroms = 0;
	PopP pop = boost::make_shared<Pop>(popname, zero_chroms, label);
  dg_add_pop (pop, gen);
}


void
Demography::dg_populate_by_name (popid popname, int members, genid gen) {
  populate_request_t req;
  req.popname = popname;
  req.members = members;
  req.gen = gen;
  populate_requests.push_back( req );
}

void
Demography::dg_populate_by_name_do (popid popname, int members, genid gen) {
  Pop *popptr = dg_get_pop_by_name(popname);
  Node *tempnode;
  int i;	

  for (i = 0; i < members; i++) {
		tempnode = dg_nodePool->make_new_leaf();
		popptr->pop_add_node (tempnode);

		dg_log (ADD_NODE, gen, tempnode, popptr);
  }	
//  totmembers += members;
}

void Demography::dg_complete_initialization() {
  int totmembers = 0;

  ForVec( populate_request, req, populate_requests )
		 totmembers += req->members;

  leafset_set_max_leaf_id( totmembers );
	this->leaf2popName.resize( totmembers + 1 );

#ifdef COSI_FREQONLY
	// extern const std::vector< popid > *leafset_leaf2popName;
	// extern const std::vector<pop_idx_t> *leafset_popname2idx;

	leafset_counts::leafset_leaf2popName = &leaf2popName;
	leafset_counts::leafset_popname2idx = &popname2idx;
	leafset_counts::leafset_npops = populate_requests.size();
	
#endif
	

  leaf_id_t leafFrom = 0;
  ForVec( populate_request, req, populate_requests ) {
		popNames.push_back( req->popname );
		sampleSizes.push_back( req->members );

		pop2leaves.push_back( make_range_leafset( leafFrom, leafFrom + req->members ) );

		for ( leaf_id_t leaf = leafFrom; leaf < leafFrom + req->members; leaf++ )
			 this->leaf2popName[ leaf ] = req->popname;
		leafFrom += req->members;
		dg_populate_by_name_do( req->popname, req->members, req->gen );
		
  }


	dg_nodePool->finishLeaves();
}

void Demography::setMaxCoalDist( plen_t maxCoalDist_ ) {
	this->maxCoalDist = maxCoalDist_;
#ifdef COSI_SUPPORT_COALAPX	
	BOOST_FOREACH( PopP p, this->pops )
		 p->setCoalMargin( len_t( maxCoalDist_ ) );
#endif	
}

const char * 
Demography::dg_get_pop_label_by_name (popid popname) const
{
  return dg_get_pop_by_name(popname) -> get_label().c_str();
}

/* returns NULL_POPID if something goes wrong.  */
popid
Demography::dg_get_pop_name_by_label (const char *label) const
{
  for (size_t ipop = 0; ipop < pops.size(); ipop++) {
    if (strcmp(pops[ipop]->get_label().c_str(), label) == 0) {
      return pops[ipop]->pop_get_name();
    }
  } 
  return NULL_POPID;
}

/* POP_SIZE FUNCTIONS */
/*
 * Returns zero if the population specified does not exist.
 */
int 
Demography::dg_set_pop_size_by_name (genid gen, popid popname, nchroms_t newsize) 
{
  Pop *popptr =  dg_get_pop_by_name(popname);
  if (popptr == NULL) {
		dg_error_nonfatal("dg_set_pop_size_by_name", 
											"pop does not exist");
		return 0;
  }
	//PRINT6( "change_size", gen, popname, popptr->pop_get_size(), popptr->pop_get_num_nodes(), newsize );
	//chkCond( 2*newsize >= popptr->pop_get_num_nodes(), "setting invalid pop size: %d;%d", 2*newsize, popptr->pop_get_num_nodes() );
  popptr->pop_set_size(newsize);
  dg_log(CHANGE_SIZE, gen, popptr);
  return 1;
}

/* COALESCE */

// Method: dg_coalesce_by_index
//
// Choose a pair of nodes from the specified population, and coalesce them.
//
// Params:
//
//    popindex - the <pop index> of the population in which to do coalescence
//    gen - the <generation> in which the coalescence happens
void
Demography::dg_coalesce_by_index (int popindex, genid gen) 
{
  Pop *popptr = dg_get_pop_by_index (popindex);
  assert(popindex >= 0);
  dg_coalesce_by_pop (popptr, gen);
}

unsigned long nreject = 0, nrecombs = 0, ncoals = 0, nretry = 0;

double coalApxRejSamp_u( 1.0 );

class MyRejPrint {
public:
	 ~MyRejPrint() { /*PRINT4( nreject, ncoals, nrecombs, nretry );*/ }
} myitem;

// Method: dg_coalesce_by_pop
//
// Choose a pair of nodes from the specified population, and coalesce them.
//
// Params:
//
//    pop - the population in which to do coalescence
//    gen - the <generation> in which the coalescence happens
Node *
Demography::dg_coalesce_by_pop(Pop* popptr, genid gen, bool forceCoalescence ) 
{
	// int node1index, node2index;
  Node *node1, *node2, *newnode;
  
  /* Choose two unique nodes to coalesce */
	
#ifdef COSI_SUPPORT_COALAPX

	if ( !forceCoalescence ) {
		std::pair< Node *, Node * > pickedNodes = popptr->chooseRandomIntersection( getRandGen() );
		node1 = pickedNodes.first;
		node2 = pickedNodes.second;
	} else {
		nchroms_t node1index = (int) (random_double() * popptr->pop_get_num_nodes());
		nchroms_t node2index = (int) (random_double() * (popptr->pop_get_num_nodes() - 1));
		
		if (node2index >= node1index) node2index++;
		
		node1 = popptr->pop_get_node (node1index);
		node2 = popptr->pop_get_node (node2index);
	}

#else  // if not approximating the coalescent

	nchroms_t node1index = (int) (random_double() * popptr->pop_get_num_nodes());
	nchroms_t node2index = (int) (random_double() * (popptr->pop_get_num_nodes() - 1));
		
	if (node2index >= node1index) node2index++;
		
	node1 = popptr->pop_get_node (node1index);
	node2 = popptr->pop_get_node (node2index);
	
#endif  // #ifdef COSI_SUPPORT_COALAPX

	if ( coalApxRejSamp_u < 1.0 && !forceCoalescence ) {
		double b1 = ToDouble( get_ploc( seglist_beg( node1->getSegs() ) ) );
		double e1 = ToDouble( get_ploc( seglist_end( node1->getSegs() ) ) );
		double b2 = ToDouble( get_ploc( seglist_beg( node2->getSegs() ) ) );
		double e2 = ToDouble( get_ploc( seglist_end( node2->getSegs() ) ) );

		if ( ( e1 + coalApxRejSamp_u < b2 ) || ( e2 + coalApxRejSamp_u < b1 ) )
			 return NULL;
	}

  /* Put mutations on the edges going from these nodes to their (about-to-be-created */
  /* parent node */
  dg_mutate->mutate_put_muts_on_seglist( node1->getSegs(), gen - node1->getGen(), node1->getGen(), popptr->pop_get_name() );
  dg_mutate->mutate_put_muts_on_seglist( node2->getSegs(), gen - node2->getGen(), node2->getGen(), popptr->pop_get_name() );

  /* Remove the two child nodes from their respective pops */
  popptr->pop_remove_node (node1);
  popptr->pop_remove_node (node2);

  /* Compute the coalesced node.  This destroys the child nodes. */
  newnode = dg_nodePool->node_coalesce (&node1, &node2, gen);

	ncoals++;

	if ( newnode ) {

		/* Add the new coalesced node into the population, in place of its two child nodes */
		popptr->pop_add_node (newnode);
		
		hooks->fire_coal( newnode );
		
		dg_log (COALESCE, gen, node1, node2, newnode, popptr);
	}

	return newnode;
}  /* dg_coalesce_by_pop() */

/* RECOMBINE */

namespace node {
extern bool NodePool_selPosGiven;
extern loc_t NodePool_selPos;
}

/*
 * Func: dg_recombine
 *
 * Choose a node to recombine, and perform the recombination.
 * If the recombination results in two parent nodes with nonempty segs,
 * the chosen node is replaced in its population by the two parent nodes;
 * if one of the parent nodes would have an empty segs, nothing is done.
 *
 * Params:
 *
 *    popindex - index of the population in which the recombination takes place
 *    gen - the generation of the two parent nodes resulting from the recombination
 *    loc - the location of the recombination event
 *    nodes_out - if non-NULL, must be an array of three Nodes; if the recombination
 *      yields two non-empty parent nodes, nodes_out is filled with
 *      pointers to the original node and the two parent nodes, otherwise its first
 *      element is set to NULL.
 */
void Demography::dg_recombine (Node *node, genid gen, loc_t loc, Node** nodes_out)
{
  Node    *newnode1 = NULL, 
		 *newnode2 = NULL;
	int nr = -1;

	Pop *popptr = node->getPop();

  /*
   * 1. Pick a random node to recombine.
   * 2. Execute recombination.
   * 3. If recombination produces two nodes (i.e.
   *    recombination occurs in one of these locations:
   *        a. in the middle of an "active segment", or
   *        b. between two "active segments",
   *    then do following steps, otherwise exit.
   * 4. Remove old node.
   * 5. Add two new nodes.
   * 7. Log it.
   */
  
  /* STEP 1 */
  
  /* STEP 2 */
	genid node_gen = node->getGen();
  gens_t edge_len = gen - node->getGen();

	nr = dg_nodePool->node_recombine (&node, &newnode1, &newnode2, gen, loc, nodes_out);
	if ( nr == 2 ) {
  
		/* STEP 3 */
		nrecombs++;

		dg_mutate->mutate_put_muts_on_seglist( newnode1->getSegs(), edge_len, node_gen, popptr->pop_get_name() );
		dg_mutate->mutate_put_muts_on_seglist( newnode2->getSegs(), edge_len, node_gen, popptr->pop_get_name() );

		// seglist_add_straight_branch_length( newnode1->getSegs(), edge_len );
		// seglist_add_straight_branch_length( newnode2->getSegs(), edge_len );

		dg_log (RECOMBINE, gen, node, newnode1, newnode2, popptr, get_loc( loc ) );
	 
	}
	hooks->fire_recomb( newnode1, newnode2, loc, gen );
}  // Demography::dg_recombine()

/* GENE CONVERSION */

void Demography::dg_gc(Node *node, genid gen, loc_t loc1, loc_t loc2, Node** nodes_out)
{
  Node    *newnode1 = NULL, *newnode2 = NULL;
  Pop     *popptr = node->getPop();
  int nr;

  /*
	 * 1. Pick a random node to geneconvert.
	 * 2. Execute recombination.
	 * 3. If recombination produces two nodes (i.e.
	 *    the location of gene conversion overlaps an
	 *    active region)
	 *    then do following steps, otherwise exit.
	 * 4. Remove old node.
	 * 5. Add two new nodes.
	 * 7. Log it.
	 */

  /* STEP 1 */
  int node_idx_in_pop = node->get_idx_in_pop();
	genid node_gen = node->getGen();
  gens_t edge_len = gen - node->getGen();
  Node *node_address = node;
  
  /* STEP 2 */
  nr = dg_nodePool->node_gc (&node, &newnode1, &newnode2, gen, loc1, loc2);

  /* STEP 3 */
  if (nr == 2) {

		dg_mutate->mutate_put_muts_on_seglist( newnode1->getSegs(), edge_len, node_gen, popptr->pop_get_name() );
		dg_mutate->mutate_put_muts_on_seglist( newnode2->getSegs(), edge_len, node_gen, popptr->pop_get_name() );

		// seglist_add_straight_branch_length( newnode1->getSegs(), edge_len );
		// seglist_add_straight_branch_length( newnode2->getSegs(), edge_len );

		if ( nodes_out ) {
			nodes_out[0] = node_address;
			nodes_out[1] = newnode1;
			nodes_out[2] = newnode2;
		}

		dg_log (GENE_CONVERSION, gen, node, newnode1, newnode2, popptr, get_loc( loc1 ), get_loc( loc2 ) );
	 
		/* STEP 4 */
		popptr->pop_remove_node_by_idx (node_idx_in_pop);

		/* STEP 5 */
		popptr->pop_add_node (newnode1);
		popptr->pop_add_node (newnode2);
		
		/* STEP 7 */
		hooks->fire_gc( newnode1, newnode2, loc1, loc2, gen );
  }
  else {
		if ( nodes_out )
			 nodes_out[0] = NULL;
  }
}

/* MIGRATE */

void 
Demography::dg_migrate_one_chrom (Pop* from_popptr, Pop* to_popptr, genid gen) 
{
  Node *tempnode;
  int node_index;

	
  node_index = (int) (random_double() * from_popptr->pop_get_num_nodes());
  tempnode = from_popptr->pop_get_node (node_index);
  from_popptr->pop_remove_node (tempnode);
  to_popptr->pop_add_node (tempnode);

  dg_log(MOVE, gen, tempnode, from_popptr, to_popptr);
}


/* MOVING NODES */
/* dg_move_nodes_by_name is called from functions that implement 
 * admixing and splitting, to move a bunch of nodes from one 
 * population to another. The number moved is binomially distributed.
 */
void 
Demography::dg_move_nodes_by_name (popid frompop, popid topop, frac_t fractionToMove, genid gen,
																	 bool_t exactFraction ) 
{
  Pop     *from_popptr, 
		 *to_popptr;
  Node *tempnode;
  nchroms_t num_to_move;
  int node_index, i;

  from_popptr = dg_get_pop_by_name (frompop);
  to_popptr = dg_get_pop_by_name(topop);

	util::chk( from_popptr, "move_nodes: unknown frompop" );
	util::chk( to_popptr, "move_nodes: unknown topop" );

	if ( to_popptr->isInactive() ) {
		std::cerr << "WARNING: moving nodes to inactive pop " << topop << " at gen " << gen << " - ignored\n";
		return;
	}

  num_to_move = exactFraction ? nchroms_t( fractionToMove * from_popptr->pop_get_num_nodes() ) : ranbinom(from_popptr->pop_get_num_nodes(), fractionToMove);
	PRINT11( "moving_nodes", gen, frompop, topop, fractionToMove, bool_t(exactFraction), num_to_move,
					 from_popptr->pop_get_size(), from_popptr->pop_get_num_nodes(), to_popptr->pop_get_size(), to_popptr->pop_get_num_nodes() );

  for (i = 0; i < num_to_move; i++) {
		node_index = (int) (random_double() 
												* from_popptr->pop_get_num_nodes());
		tempnode = from_popptr->pop_get_node (node_index);
		from_popptr->pop_remove_node (tempnode);
		to_popptr->pop_add_node (tempnode);

		dg_log(MOVE, gen, tempnode, from_popptr, to_popptr);
  }
	PRINT3( "aft_moving_nodes", from_popptr->pop_get_num_nodes(), to_popptr->pop_get_num_nodes() );
}

/* NODE FUNCTIONS */

/* GET_NUM_NODES */
nchroms_t
Demography::dg_get_num_nodes (void) const
{
  nchroms_t     total = 0;

  for (size_t i = 0; i < pops.size(); i++)
		 total += dg_get_num_nodes_in_pop_by_index (i);
  return total;
}

bool_t
Demography::dg_done_coalescent() const
{
  /*
   * so i don't want to delete the pops right away, but i don't
   * want to count empty pops either. so we count non-empty pops
   * instead of deleting them. If there is more than one non-empty
   * pop, we know we're not done.
   */
#if 1
  ForEach( PopP p, pops ) if ( !p->empty() ) return False;
#else
	bool_t seenNonemptyPop = False;
	For( p, pops ) {
		if ( !(*p)->empty() ) {
			if ( seenNonemptyPop ) return False;
			else seenNonemptyPop = True;
		}
	}
#endif	
	return True;
}

#ifdef COSI_DEV_DO_LOGGING
/* LOGGING */
/* The dg_log function can take a variable number of arguments,
 * depending on what we are logging. The only external call of this
 * function occurs in historical.c, where it is needed to log the
 * historical events that occur. All other calls are from within
 * demography.c.
 */
void 
Demography::dg_log (event_kind_t type, genid gen, ...) {
  FILE   *outputptr = logfp;
  Pop     *popptr, 
		 *popptr2;
  Node    *nodeptr1, 
		 *nodeptr2, 
		 *nodeptr3;
  double  double1;
  loc_t  loc, loc2;
  char*   string1;
  va_list ap; /* points to unnamed args */

  if (LOGGING) {

		va_start(ap, gen);

		switch (type) {
		case ADD_NODE: 
			/* ADD_NODE newnode pop */
			nodeptr1 = va_arg(ap, Node *);
			popptr = va_arg(ap, Pop *);
			if (outputptr != NULL) {fprintf(outputptr, 
																			"%f\tADD\tnode: %d pop: %d\n", 
																			gen, nodeptr1->name, ToInt( popptr->pop_get_name( ) ));}
			break;
      
		case CHANGE_SIZE: 
			/* CHANGE_SIZE pop */
			popptr = va_arg(ap, Pop*);
			if (outputptr != NULL) {fprintf(outputptr, 
																			"%f\tchange_size\tpop: %d size: %d\n", 
																			gen, ToInt( popptr->pop_get_name() ), popptr->popsize);}
			break;

		case COALESCE:
			/* COALESCE oldnode1 oldnode2 newnode pop */
			nodeptr1 = va_arg(ap, Node *);
			nodeptr2 = va_arg(ap, Node *);
			nodeptr3 = va_arg(ap, Node *);
			popptr = va_arg(ap, Pop *);
			if (outputptr != NULL) {fprintf(outputptr, "%f\tC\t%d %d -> %d pop: %d\n", 
																			gen, nodeptr1->name, nodeptr2->name, 
																			nodeptr3->name, ToInt( popptr->pop_get_name() ));}
			break;
      
		case CREATE_POP:
			/* CREATE_POP pop */
			popptr = va_arg(ap, Pop*);
			if (outputptr != NULL) {fprintf(outputptr, 
																			"%f\tcreate_pop\tpop: %d size: %d\n", 
																			gen, ToInt( popptr->pop_get_name() ), popptr->popsize);}
			break;

		case DONE:
			/* DONE oldnode donenode [newnode|NULL] */
			nodeptr1 = va_arg(ap, Node *);
			nodeptr2 = va_arg(ap, Node *);
			nodeptr3 = va_arg(ap, Node *);

			if (nodeptr3 != NULL) {
				if (outputptr != NULL) {fprintf(outputptr, 
																				"%f\tD\t%d -> %d | %d %f %f\n",
																				gen, nodeptr1->name, 
																				nodeptr3->name, 
																				nodeptr2->name, 
																				seglist_first_seg_beg( nodeptr2->getSegs() ), 
																				seglist_first_seg_end( nodeptr2->getSegs() ) );}
			}
			else {
				if (outputptr != NULL) {fprintf(outputptr, "%f\tD\t%d\n", gen,
																				nodeptr1->name);}
			}
			break;

		case GENE_CONVERSION:
			/* RECOMBINE oldnode1 oldnode2 newnode pop */
			nodeptr1 = va_arg(ap, Node*);
			nodeptr2 = va_arg(ap, Node *);
			nodeptr3 = va_arg(ap, Node *);
			popptr = va_arg(ap, Pop*);
			loc = va_arg(ap, loc_t);
			loc2 = va_arg(ap, loc_t);
			
			if (nodeptr3 != NULL) {
				if (outputptr != NULL) {fprintf(outputptr, 
																				"%f\tG\t%d -> %d %d %d %f %f\n", 
																				gen, nodeptr1->name, 
																				nodeptr2->name, 
																				nodeptr3->name,
																				ToInt( popptr->pop_get_name() ), loc, loc2);}
			}
			else {
				if (outputptr != NULL) {fprintf(outputptr, 
																				"%f\tG\t%d -> %d  - %d %f %f\n", 
																				gen, nodeptr1->name, 
																				nodeptr2->name, 
																				ToInt( popptr->pop_get_name() ), loc, loc2);}
			}
			break;
			


		case HISTORICAL:
			/* HISTORICAL string_description */
			string1 = va_arg(ap, char*);
			if (outputptr != NULL) {fprintf(outputptr, "%f\tH\t%s\n",
																			gen, string1);}
			break;

		case MIG_RATE:
			/* MIG_RATE from-pop to-pop new_rate */
			/* note that "from" and "to" are in real time.
			 * chromosomes will move in the opposite direction */
			popptr = va_arg(ap, Pop *);
			popptr2 = va_arg(ap, Pop *);
			double1 = va_arg(ap, double);
			if (outputptr != NULL) {fprintf(outputptr, 
																			"%f\tmig_rate\t%d %d %f\n", 
																			gen, ToInt( popptr->pop_get_name() ), 
																			ToInt( popptr2->pop_get_name() ), double1);}
			break;

		case MOVE:
			/* MOVE node from-pop to-pop */
			nodeptr1 = va_arg(ap, Node *);
			popptr = va_arg(ap, Pop *);
			popptr2 = va_arg(ap, Pop *);			
			if (outputptr != NULL) {fprintf(outputptr, "%f\tM\t%d %d %d\n",
																			gen, nodeptr1->name, ToInt( popptr->pop_get_name() ), 
																			ToInt( popptr2->pop_get_name() ) );}
			break;

		case RECOMBINE:
			/* RECOMBINE oldnode1 oldnode2 newnode pop */
			nodeptr1 = va_arg(ap, Node*);
			nodeptr2 = va_arg(ap, Node *);
			nodeptr3 = va_arg(ap, Node *);
			popptr = va_arg(ap, Pop*);
			loc = va_arg(ap, loc_t);

			if (nodeptr3 != NULL) {
				if (outputptr != NULL) {fprintf(outputptr, 
																				"%f\tR\t%d -> %d %d %d %f\n", 
																				gen, nodeptr1->name, 
																				nodeptr2->name, 
																				nodeptr3->name,
																				ToInt( popptr->pop_get_name() ), loc);}
			}
			else {
				if (outputptr != NULL) {fprintf(outputptr, 
																				"%f\tR\t%d -> %d  - %d %f\n", 
																				gen, nodeptr1->name, 
																				nodeptr2->name, 
																				ToInt( popptr->pop_get_name() ), loc);}
			}
			break;
		       
		}
  }	
}
#endif

void 
Demography::dg_set_logfile (FILE *fp) 
{
  logfp = fp;
}

void 
Demography::dg_logging_on () 
{
  LOGGING = 1;
}

void 
Demography::dg_logging_off () 
{
  LOGGING = 0;
}

/***********************************************************/
/* UTILITIES FOR INTERNAL USE ONLY */

/* DG_EXIT
 */
void 
dg_exit(const char* funct_name, const char* error_string) 
{
  fprintf(stderr, "demography.c | %s: %s\n",
					funct_name, error_string);
	chkCond( False, "demography error" );
}

void 
dg_error_nonfatal(const char* funct_name, const char* error_string) 
{
  fprintf(stderr, "demography.c | %s: %s\n",
					funct_name, error_string);
}

/* POP FUNCTIONS */
Pop *
Demography::dg_add_pop (PopP newpop, genid gen) {
	int popid_int = ToInt( newpop->pop_get_name() );
	if ( popid_int >= isize( popname2idx ) )
		 popname2idx.resize( popid_int+1, NULL_POP_IDX );
	popname2idx[ popid_int ] = pops.size();

#ifdef COSI_SUPPORT_COALAPX
	newpop->setCoalMargin( len_t( maxCoalDist ) );
#endif
	
  pops.push_back( newpop );
  dg_log (CREATE_POP, gen, newpop.get());
	return newpop.get();
}

popid Demography::find_unused_popname() const {
	int max_id = 1;
	ForEach( PopP p, pops )
		 max_id = std::max( max_id, ToInt( p->pop_get_name() ) );
	return popid( max_id + 1 );
}

/* GET POPS */
/* both functions return NULL if the specified pop is not to
 * be found.
 */

Pop * 
Demography::dg_get_pop_by_name(popid popname) const {
	pop_idx_t idx = dg_get_pop_index_by_name( popname );
	return idx == NULL_POP_IDX ? ((Pop *)NULL) : pops[ idx ].get();
}

popid 
Demography::dg_get_pop_name_by_index (int popindex) const {
  return dg_get_pop_by_index(popindex)->pop_get_name();
}

pop_idx_t Demography::dg_get_pop_index_by_name (popid popname) const {
	int popname_int = ToInt( popname );
	return popname_int < isize( popname2idx ) ? popname2idx[ popname_int ] : NULL_POP_IDX ;
}


Pop * 
Demography::dg_get_pop_by_index (int index1) const { 
  return dg_get_pop_from_list(index1);
}

Pop *
Demography::dg_get_pop_from_list(int index1) const {
  return pops[index1].get();
}

}  // namespace cosi
