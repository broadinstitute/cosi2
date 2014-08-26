//
// Header: pop.h
//
// Representation of populations.
//

/* $Id: pop.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_POP_H 
#define __INCLUDE_COSI_POP_H

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/utility/declval.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/nodelist.h>
#include <cosi/generalmath.h>
#include <cosi/coalrate.h>

namespace cosi {

using std::string;
using node::Node;
using node::NodeList;

//
// Class: Pop
//
// The current state of one population (deme) during the simulation:
// the size of the population, as well as the list of nodes (lineages)
// from this population we're currently tracking.  The size of the population
// should be much larger than the number of tracked nodes, for the coalescent
// approximation to hold.
//
// <Demography> keeps the list of all pops in the simulation.
//
// Each population is assigned an integer id (also called pop name) in the parameter
// file; these ids, of type popid, are used throughout the code to refer to
// populations.  <Demography> keeps the mapping of pop names to Pop objects.
//
// A Pop object may also represent a sub-population of an actual population.
// For example, when simulating selective sweeps, for each possible allele at
// the causal mutation location, a separate Pop object is used
// to represent the subpopulation of chromosomes carrying that allele; see
// 
// Each <Node> belongs to one pop at any given time.
//
// Invariant: members.size() <= 2 * popsize
//
// In practice, the above invariant may be temporarily violated in several cases: eg a recombination or
// gene conversion may increase the number of nodes beyond the popsize.  But in such cases, the
// probability of coalescence within such a population will be so high that coalescences will be forsed
// and the situation will quickly rectify itself.
//
class Pop {
	 
public:

	 Pop(popid  popname_ , nchroms_t  popsize_ , 
			 const string& label_ );
	 ~Pop();

	 class PopListener {
	 public:
			virtual void nodeAdded( Pop *pop, nodeid nodeName, idx_t idx_in_pop, loc_t beg, loc_t end ) = 0;
			virtual void nodeRemoved( Pop *, idx_t idx_in_pop ) = 0;
			virtual ~PopListener() { }
	 };

	 void setPopListener( boost::shared_ptr<PopListener> popListener_ ) { this->popListener = popListener_; }
	 void clearPopListener() { this->popListener.reset(); }

	 void pop_remove_node ( Node *);
	 void pop_remove_node_by_idx ( int idx_in_pop);
	 void pop_add_node ( Node *);
	 nchroms_t pop_get_num_nodes () const { return members.size(); }
	 bool_t empty() const { return members.empty(); }
	 void pop_set_size (nchroms_t  popsize_  ) { popsize = popsize_; }
	 popid pop_get_name() const { return name; }

	 Node *pop_get_node( int nodeIdx ) const;

	 nchroms_t pop_get_size() const { return popsize; }
	 const std::string& get_label() const { return label; }

	 // Method: getMembers
	 // Returns the current list of <Nodes> in the population.
	 const NodeList& getMembers() const { return members; }

	 typedef math::ArrivalProcess<genid, math::Any< RandGen > > coal_arrival_process_type_ptr;

	 void setCoalArrivalProcess( coal_arrival_process_type_ptr coalArrivalProcess_ ) {
		 this->coalArrivalProcess = coalArrivalProcess_;
//				boost::make_shared< coal_arrival_process_type >( math::makeNonHomogeneousPoissonProcess( math::coalRateFunction( *sizeTraj_ ) ) );
	 }
	 coal_arrival_process_type_ptr getCoalArrivalProcess() const { return this->coalArrivalProcess; }
	 void clearCoalArrivalProcess() { this->coalArrivalProcess.reset(); }

	 template <typename TSpec>
	 void setSizeTraj( const math::Function< genid, popsize_float_t, TSpec >& traj, genid gen ) {
		 setCoalArrivalProcess(
			 math::ArrivalProcess< genid, math::Any< RandGen > >(
				 math::makeNonHomogeneousPoissonProcess (
					 math::coalRateFunction( traj ), gen, std::string( "coal in pop " + label )
					 )));
	 }

#ifdef COSI_DEV	 
	 prob_t getCoalRate() const { return coalRate; }
	 void setCoalRate( prob_t coalRate_ ) { coalRate = coalRate_; }
#endif

#ifdef COSI_SUPPORT_COALAPX
	 void setCoalMargin( len_t margin_ );

	 bool restrictingCoalescence() const;
	 nchromPairs_t getNumCoalesceableChromPairs() const;

	 std::pair< Node *, Node * > chooseRandomIntersection( RandGenP randGen );
	 
#endif	 // #ifdef COSI_SUPPORT_COALAPX
	 
private:
	 // Private field: name
	 // The id assigned to the population in the parameter file.
	 popid name;
	 
	 // Private field: members
	 // The current <Node>s belonging to this population.  
	 NodeList members;

	 // Private field: popsize
	 // The current (effective) size of this population.
	 // This is actually the number of diploid individuals in the population;
	 // the number of haploid chromosomes is twice this.
	 nchroms_t popsize;

	 // Private field: label
	 // A human-readable name for a population.  Does not affect
	 // any computations, used only for debugging and reporting.
	 std::string label;

	 // Field: sizeTraj
	 // Trajectory for the population size.
//	 size_traj_type_ptr sizeTraj;

	 coal_arrival_process_type_ptr coalArrivalProcess;

	 bool isRestrictingCoalescence;

#ifdef COSI_SUPPORT_COALAPX	 
	 // Field: hullMgr
	 // Keeps track of the convex hulls of the nodes in this pop.
	 HullMgrP hullMgr;

	 void chkHullMgr() const;
#endif	 

#ifdef COSI_DEV	 
	 prob_t coalRate;
#endif

	 bool useCoalApx() const;

	 boost::shared_ptr<PopListener> popListener;
	 
};  // class Pop


inline Node *Pop::pop_get_node( int nodeIdx ) const { return node::nodelist_get_node( nodeIdx, &members ); }

}  // namespace cosi

#endif /* POP_H */
