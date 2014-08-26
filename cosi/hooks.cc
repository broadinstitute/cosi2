
#include <cosi/hooks.h>

namespace cosi {

//
// Class impl: Hook
//

Hook::Hook() { }
Hook::~Hook() { }
void Hook::handle_recomb( Node *, Node *, loc_t, genid ) { }
void Hook::handle_gc( Node *, Node *, loc_t , loc_t, genid ) { }
void Hook::handle_coal( Node * ) {}
void Hook::handle_set_migrate_rate( popid /*from*/, popid /*to*/, prob_per_chrom_per_gen_t /*rate*/ ) { }
void Hook::handle_sweep_end( leafset_p /* sel_leaves */ ) { };

// Virtual method: handle_make_leaf
//
// Called whenever a leaf node is created.
void Hook::handle_make_leaf( nodeid /*leafNodeId*/ ) { }


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
void Hook::handle_add_edge( nodeid /*nodeId_moreRecent*/,
														nodeid /*nodeId_lessRecent*/,
														genid /*genId_moreRecent*/,
														genid /*genId_lessRecent*/, const Seglist * /*seglist*/,
														edge_kind_t ) { }

// Virtual method: handle_fully_coalesced
//
// Called when, after a coalescence where part of the resulting node's segs has fully coalesced,
// a new node is created to carry the not-yet-coalesced parts of the old node's segs.
void Hook::handle_fully_coalesced( nodeid /*nodeId_old*/, nodeid /*nodeId_new*/, const Seglist * /*seglist_new*/ ) { }

}  // namespace cosi
