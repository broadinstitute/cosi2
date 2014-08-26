
#include <utility>
#include <iostream>
#include <cosi/stats.h>
#include <cosi/seglist.h>

namespace cosi {

using std::make_pair;
using std::map;
using std::max;
using std::cout;

//
// Class impl: TreeStatsHook
//

TreeStatsHook::TreeStatsHook(): totMutimeAtJoins( mutime_t(0) ), totMutime( mutime_t( 0 ) ), n_joins(0), n_edges(0),
																n_coal(0), n_recomb(0)
{
}

TreeStatsHook::~TreeStatsHook() { }

// Virtual method: handle_make_leaf
//
// Called whenever a leaf node is created.
void TreeStatsHook::handle_make_leaf( nodeid /*leafNodeId*/ ) {

}

void TreeStatsHook::handle_recomb( Node *, Node *, loc_t, genid) { n_recomb++; }
void TreeStatsHook::handle_coal( Node * ) { n_coal++; }

// Virtual method: handle_tree_edge
//
// Called when a tree edge is added.  Only coalescence nodes of the ARG are taken into account; recombination and gc nodes
// are ignored.
//
// Params:
//
//     beg, end - the <Seg> inherited along the edge
//     genMoreRecent, genLessRecent - generations at edge ednpoints
//     leafset - the leafset inherited along the edge
// void TreeStatsHook::handle_tree_edge( loc_t beg, loc_t end, genid genMoreRecent, genid genLessRecent, leafset_p leafset ) {
// 	cout << "ARG " << get_ploc( beg ) << " " << get_ploc( end ) << " " << genMoreRecent << " " << genLessRecent << " " << leafset_size( leafset ) << "\n";
// }


// Virtual method: handle_add_edge
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
void TreeStatsHook::handle_add_edge( nodeid /*nodeId_moreRecent*/,
																		 nodeid /*nodeId_lessRecent*/,
																		 genid /*genId_moreRecent*/,
																		 genid /*genId_lessRecent*/, const Seglist * /*seglist*/,
																		 edge_kind_t /*edgeKind*/ ) {
//	totMutime += seglist_tot_len( seglist ) * ( genId_lessRecent - genId_moreRecent );
	//PRINT9( "add_edge", n_edges, nodeId_moreRecent, genId_moreRecent, nodeId_lessRecent, genId_lessRecent, seglist, seglist_tot_len( seglist ), totMutime );
//	n_edges++;
}  // Hook::handle_add_edge()

// Virtual method: handle_fully_coalesced
//
// Called when, after a coalescence where part of the resulting node's segs has fully coalesced,
// a new node is created to carry the not-yet-coalesced parts of the old node's segs.
void TreeStatsHook::handle_fully_coalesced( nodeid /*nodeId_old*/, nodeid /*nodeId_new*/, const Seglist * /*seglist_new*/ ) {
}

// Method: printStats
// Write, to standard output, any stats gathered during the simulation, after one simulation completes.
void TreeStatsHook::printStats() {
	cout << "stat nrecomb " << n_recomb << "\n";
	cout << "stat ncoal " << n_coal << "\n";
	// if ( !equal_eps( totMutimeAtJoins, totMutime, 1e-10 ) )
	// 	 PRINT5( totMutimeAtJoins, totMutime, totMutimeAtJoins - totMutime, n_joins, n_edges );
}

//
// End class impl: TreeStatsHook
//


}
