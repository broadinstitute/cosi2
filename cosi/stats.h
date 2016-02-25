/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright (2005) by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/

/*
	Header: stats.h

	Code for collecting and reporting various statistics about simulations.
*/

#ifndef __INCLUDE_COSI_STATS_H
#define __INCLUDE_COSI_STATS_H

#include <stdint.h>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <cosi/general/utils.h>
#include <cosi/defs.h>
#include <cosi/hooks.h>

namespace cosi {

using std::map;
using std::ostringstream;
using std::string;

//
// Class: TreeStatsHook
//
// Gathers statistics about the ARG during the ARG's construction.
//
class TreeStatsHook: public Hook {
public:
	 TreeStatsHook();
	 virtual ~TreeStatsHook();

	 virtual void handle_recomb( Node *, Node *, loc_t, genid);
	 virtual void handle_coal( Node * );

	 // VMethodP: handle_make_leaf
	 //
	 // Called whenever a leaf node is created.
	 virtual void handle_make_leaf( nodeid leafNodeId );

	 // VMethodP: handle_add_edge
	 //
	 // Called after adding an edge to the ARG.  The new edge may be the result of a coalescence,
	 // a recombination or a gene conversion.
	 //
	 // Params:
	 //
	 //   nodeId_moreRecent - <nodeid> of the more recent (lower gen) node of the edge
	 //   nodeId_lessRecent - <nodeid> of the less recent (higher gen) node of the edge.
	 //   genId_moreRecent - <genid> of the more recent (lower gen) node of the edge
	 //   genId_lessRecent - <genid> of the less recent (higher gen) node of the edge.
	 //   seglist - the segments inherited along the edge.  NOTE: this seglist may be destroyed
	 //       after the call, so make a copy of you need to save it.
	 //       
	 virtual void handle_add_edge( nodeid nodeId_moreRecent,
																 nodeid nodeId_lessRecent,
																 genid genId_moreRecent,
																 genid genId_lessRecent, const Seglist *seglist,
																 edge_kind_t edgeKind );

	 // VMethodP: handle_fully_coalesced
	 //
	 // Called when, after a coalescence where part of the resulting node's segs has fully coalesced,
	 // a new node is created to carry the not-yet-coalesced parts of the old node's segs.
	 virtual void handle_fully_coalesced( nodeid nodeId_old, nodeid nodeId_new, const Seglist *seglist_new );

	 // MethodP: printStats
	 // Write, to standard output, any stats gathered during the simulation, after one simulation completes.
	 void printStats();

private:
	 // Field: totMutimeAtJoins
	 // Total mutime computed at seg join time.
	 mutime_t totMutimeAtJoins;

	 // Field: totMutime
	 // Total mutime computed directly.
	 mutime_t totMutime;

	 // Field: n_joins
	 // Number of segment joins
	 uintmax_t n_joins;

	 // Field: n_edges
	 // Number of (short) ARG edges.
	 uintmax_t n_edges;

	 // Field: n_coal
	 // Number of coalescent events
	 uintmax_t n_coal;

	 // Field: n_recomb
	 // Number of recomb events
	 uintmax_t n_recomb;
	 
};  // class TreeStatsHook

typedef boost::shared_ptr<TreeStatsHook> TreeStatsHookP;

}  // namespace cosi


#endif  // #ifndef __INCLUDE_COSI_STATS_H

