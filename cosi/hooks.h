//
// Header: hooks.h
//
// A mechanism for establishing callbacks to be called when certain events
// happen during simulation.
//
// To use this mechanism, you subclass the class <Hook>, overriding some of the handle_* methods for the events
// you want to handle.  All these methods have empty default implementations, so you only need to override the
// methods you care about.
//

#ifndef __INCLUDE_COSI_HOOKS
#define __INCLUDE_COSI_HOOKS

#include <list>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>
#include <cosi/nodefwd.h>
#include <cosi/seglistfwd.h>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

namespace cosi {

using std::list;
using node::Node;
using seglist::Seglist;

enum edge_kind_t { EDGE_RECOMB, EDGE_GC, EDGE_COAL };
 

//
// Class: Hook
//
// Base class from which user hooks should be derived.
// Users can override any subset of handle_* methods to handle specific events;
// other handle_* methods will default to no-ops.
class Hook {
protected:
	 Hook();
	 virtual ~Hook();

private:	 
	 // Virtual method: handle_recomb
	 // Called after each recombination.
	 virtual void handle_recomb( Node *, Node *, loc_t, genid);

	 // Virtual method: handle_gc
	 // Called after each gene conversion.
	 virtual void handle_gc( Node *, Node *, loc_t, loc_t, genid);

	 virtual void handle_coal( Node * );

	 // Virtual method: handle_set_migrate_rate
	 // Called after each setting of migrate rate.
	 virtual void handle_set_migrate_rate( popid from, popid to, prob_per_chrom_per_gen_t rate );


	 //
	 // Method group: NodePool events
	 //

	 // Virtual method: handle_make_leaf
	 //
	 // Called whenever a leaf node is created.
	 virtual void handle_make_leaf( nodeid leafNodeId );

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
	 //       after the call, so make a copy if you need to save it.
	 //       
	 virtual void handle_add_edge( nodeid nodeId_moreRecent,
																 nodeid nodeId_lessRecent,
																 genid genId_moreRecent,
																 genid genId_lessRecent, const Seglist *seglist,
																 edge_kind_t edgeKind );

	 // Virtual method: handle_fully_coalesced
	 //
	 // Called when, after a coalescence where part of the resulting node's segs has fully coalesced,
	 // a new node is created to carry the not-yet-coalesced parts of the old node's segs.
	 virtual void handle_fully_coalesced( nodeid nodeId_old, nodeid nodeId_new, const Seglist *seglist_new );

	 // Method: handle_sweep_end
	 //
	 // Called when the backwards simulation reaches the generation where the sweep ends.
	 virtual void handle_sweep_end( leafset_p sel_leaves );

	 friend class Hooks;
};  // class Hook

typedef boost::shared_ptr<Hook> HookP;

//
// Class: Hooks
//
// Manages a list of hooks.  Provides methods to signal a specified
// event to all hooks in the list.
//
class Hooks {
public:
	 void addHook( HookP hook ) { hooks.push_back( hook ); }
	 void removeHook( HookP hook ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				if ( *it == hook ) {
					hooks.erase( it );
					break;
				}
	 }

	 void fire_recomb( Node *node1, Node *node2, loc_t loc, genid gen ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_recomb( node1, node2, loc, gen );
	 }
	 void fire_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2, genid gen ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_gc( node1, node2, loc1, loc2, gen );
	 }
	 void fire_coal( Node *node ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_coal( node );
	 }

	 void fire_set_migrate_rate( popid from, popid to, prob_per_chrom_per_gen_t rate ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_set_migrate_rate( from, to, rate );
	 }

	 void fire_make_leaf( nodeid leafNodeId ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_make_leaf( leafNodeId );
	 }
	 
	 void fire_add_edge( nodeid nodeId_moreRecent, nodeid nodeId_lessRecent,
											 genid genId_moreRecent, genid genId_lessRecent, const Seglist *seglist, edge_kind_t edgeKind ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_add_edge( nodeId_moreRecent, nodeId_lessRecent,
																genId_moreRecent, genId_lessRecent,
																seglist, edgeKind );
	 }

	 void fire_fully_coalesced( nodeid nodeId_old, nodeid nodeId_new, const Seglist *seglist_new ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_fully_coalesced( nodeId_old, nodeId_new, seglist_new );
	 }

	 void fire_sweep_end( leafset_p sel_leaves ) {
		 for ( list<HookP>::iterator it = hooks.begin(); it != hooks.end(); it++ )
				(*it)->handle_sweep_end( sel_leaves );
	 }

	 
private:
	 list<HookP> hooks;
};  // class Hooks

typedef boost::shared_ptr<Hooks> HooksP;

//
// Class: Hookable
//
// Convenience base class to help classes include a list of hooks.
//
struct Hookable {
	 Hookable(): hooks( boost::make_shared<Hooks>() ) { }
	 void setHooks( HooksP hooks_ ) { hooks = hooks_; }
	 HooksP getHooks() const { return hooks; }

	 void addHook( HookP hook ) { hooks->addHook( hook ); }
	 void removeHook( HookP hook ) { hooks->removeHook( hook ); }

protected:
	 HooksP hooks;
};

}


#endif // #ifndef __INCLUDE_COSI_HOOKS
