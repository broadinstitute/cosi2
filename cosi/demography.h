/* $Id: demography.h,v 1.6 2011/06/01 21:33:16 sfs Exp $ */

/*
 * Header: demography.h
 *
 * Data structures and methods for representing the recombination sites, and the <populations>.
 */

#ifndef __INCLUDE_COSI_DEMOGRAPHY_H
#define __INCLUDE_COSI_DEMOGRAPHY_H
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/cosirand.h>
#include <cosi/hooks.h>
#include <cosi/leafset.h>
#include <cosi/nodefwd.h>
#include <cosi/pop.h>

namespace cosi {

using std::vector;
using node::Node;
using node::NodePoolP;

//
// Struct: Demography
//
// The demographic model: the list of populations (and their current state during
// simulation), the sample sizes, etc.
//
class Demography: public Hookable, public HasRandGen {
	 
public:

/* definitions for dg_log */
	 enum event_kind_t { ADD_NODE = 1,
											 COALESCE = 2,
											 CREATE_POP = 3,
											 RECOMBINE = 4,
											 DONE = 5,
											 MOVE = 7,
											 CHANGE_SIZE = 8,
											 MIG_RATE = 9,
											 HISTORICAL = 10,
											 GENE_CONVERSION = 12 };

	 Demography();


	 // Method: getTotSamples
	 // Returns the total number of samples, from all populations.
	 nchroms_t getTotSamples() const { return std::accumulate( sampleSizes.begin(), sampleSizes.end(), nchroms_t(0) ); }

	 void dg_setMutate( MutateP mutate_ );
	 void dg_setNodePool( NodePoolP nodePool_ );

	 void setMaxCoalDist( plen_t maxCoalDist_ );
	 void setMaxCoalDistCvxHull( bool_t maxCoalDistCvxHull_ ) { this->maxCoalDistCvxHull = maxCoalDistCvxHull_; }

	 NodePoolP getNodePool() const { return dg_nodePool; }
	 
	 void dg_complete_initialization(void);

	 /* POP FUNCTIONS */

	 void dg_create_pop (popid popname, const std::string& label, 
											 genid gen);
	 void dg_populate_by_name (popid popname, 
														 int members,  genid gen = ZERO_GEN );

	 int dg_get_num_pops(void) const { return pops.size(); }

	 popid find_unused_popname() const;

	 const char* dg_get_pop_label_by_name (popid popname) const;
	 popid dg_get_pop_name_by_label (const char* label) const;
	 popid dg_get_pop_name_by_index (int popindex) const;
	 int dg_get_pop_index_by_name (popid popname) const;
	 Pop* dg_get_pop_by_index (int index1) const;
	 Pop* dg_get_pop_by_name(popid popname) const;

	 /* POP_SIZE FUNCTIONS */

	 nchroms_t dg_get_num_nodes (void) const;
	 
	 nchroms_t dg_get_pop_size_by_index(int popindex) const { return pops[popindex]->pop_get_size(); }

	 int dg_set_pop_size_by_name (genid gen, popid popname, 
																nchroms_t newsize);

	 /* POISSON EVENTS */

	 /* coalesce */

	 // Method: dg_coalesce_by_index
	 //
	 // Choose a pair of nodes from the specified population, and coalesce them.
	 //
	 // Params:
	 //
	 //    popindex - the <pop index> of the population in which to do coalescence
	 //    gen - the <generation> in which the coalescence happens
	 void dg_coalesce_by_index (int popindex, genid gen);
	 
	 // Method: dg_coalesce_by_pop
	 //
	 // Choose a pair of nodes from the specified population, and coalesce them.
	 //
	 // Params:
	 //
	 //    pop - the population in which to do coalescence
	 //    gen - the <generation> in which the coalescence happens
	 //    forceCoalescence - if True, force a coalescence even if in approximate-coalescence mode.
	 Node * dg_coalesce_by_pop (Pop* pop, genid gen, bool_t forceCoalescence = False );

	 /* recombine */  /* array[0]=old node  [1,2]=new nodes.  Caller must free ptr.  NULL => no visible recomb. */
	 void dg_recombine(Node *, genid, loc_t, Node** );

	 void dg_gc(Node *node, genid gen, loc_t loc1, loc_t loc2, Node** nodes_out);
	 
	 /* migrate */
	 void dg_migrate_one_chrom (Pop* frompop, Pop* topop, genid);

	 /* FOR HISTORICAL EVENTS */
	 void dg_move_nodes_by_name (popid frompop, popid topop, frac_t fractionToMove, genid gen, bool_t exactFraction = False);

	 /* NODE FUNCTIONS */

	 nchroms_t dg_get_num_nodes_in_pop_by_index(int popindex) const {
		 return pops[popindex]->pop_get_num_nodes();
	 }

	 /* DONE? */
	 bool_t dg_done_coalescent(void) const;

	 /* logging function */
	 /* since logging is slow, we have the option to turn it off. */

#ifdef COSI_DEV_DO_LOGGING
	 void dg_log (event_kind_t type, genid gen, ...);
#else
	 void dg_log (event_kind_t /*type*/, genid /*gen*/, ...) { }
#endif

	 void dg_logging_on (void);
	 void dg_logging_off (void);
	 void dg_set_logfile (FILE *);
	 
	 void setVerbose( bool_t verbose_ ) { verbose = verbose_; }
	 bool_t getVerbose() const { return verbose; }

	 Pop *dg_add_pop (PopP, genid gen);

	 MutateP getMutate() const { return dg_mutate; }

	 const vector< PopP >& getPops() const { return pops; }
	 const vector< popid >& getPopNames() const { return popNames; }
	 const vector< nchroms_t >& getSampleSizes() const { return sampleSizes; }
	 const vector< popid >& getLeaf2PopName() const { return leaf2popName; }

#ifndef COSI_FREQONLY	 
	 std::map< popid, nchroms_t > getLeafsetPopCounts( leafset_p leafset ) const {
		 std::map< popid, nchroms_t > result;
		 COSI_FOR_LEAFSET( leafset, leaf, {
				 result[ leaf2popName[ leaf ] ]++;
			 });
		 return result;
	 }
#endif	 

	 template <typename OutputIterator>
	 void getLeafsetPopCounts( leafset_p leafset, OutputIterator it ) const {
		 std::vector<nchroms_t> counts( sampleSizes.size(), 0 );
		 COSI_FOR_LEAFSET( leafset, leaf, {
				 counts[ popname2idx[ ToInt( leaf2popName[ leaf ] ) ] ]++;
			 });
		 boost::copy( counts, it );
	 }

	 nchroms_t get_N0_tot() const { return N0_tot; }
	 
	 
private:
	 // Field: pops
	 // All populations currently in the model.  Note that some of them may be empty of <Nodes>.
	 vector< PopP > pops;
	 //	 nchroms_t totmembers;
	 
	 // Field: popNames
	 // For each pop, the popid of that pop.
	 vector< popid > popNames;
	 
	 // Field: sampleSizes
	 // For each pop, the size of the present-day sample from that pop.
	 vector< nchroms_t > sampleSizes;

#ifdef COSI_EHH	 
	 // Field: pop2leaves
	 // For each pop, a leafset containing the leaves in that pop.
	 vector< leafset_p > pop2leaves;
#endif

	 // Field: leaf2popName
	 // For each leaf, the name of the pop from which that leaf comes
	 vector< popid > leaf2popName;
	 
	 MutateP dg_mutate;
	 NodePoolP dg_nodePool;

	 int LOGGING;
	 FILE* logfp;

	 Pop *dg_get_pop_from_list(int index) const;

	 typedef struct populate_request {
			popid popname;
			int members;
			genid gen;
	 } populate_request_t;
	 
	 vector<populate_request_t> populate_requests;

	 vector<pop_idx_t> popname2idx;
	 
	 // Field: verbose
	 // Print progress during sim
	 bool_t verbose;

	 plen_t maxCoalDist;
	 bool_t maxCoalDistCvxHull;
	 
	 void
	 dg_populate_by_name_do (popid popname, int members, genid gen);

	 nchroms_t N0_tot;

};  // class Demography

typedef boost::shared_ptr<Demography> DemographyP;

} // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_DEMOGRAPHY_H

