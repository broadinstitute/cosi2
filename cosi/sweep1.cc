#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <ios>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <cosi/general/utils.h>
#include <cosi/decls.h>
#include <cosi/node.h>
#include <cosi/pop.h>
#include <cosi/demography.h>
#include <cosi/mutlist.h>
#include <cosi/mutate.h>
#include <cosi/historical.h>
#include <cosi/traj.h>
#include <cosi/hooks.h>
#include <cosi/seglist.h>
#include <cosi/mutcontext.h>
#include <cosi/sweep1.h>

namespace cosi {

namespace sweep1 {

#define ForEach BOOST_FOREACH  
  
using std::vector;
using std::map;
using std::set;
using std::ifstream;
using std::ofstream;
using util::map_get;
using util::STLContains;
using node::Node;
using node::NodeList;

/////////////////

//
// Class: Event_SweepOnePop
//
// A selective sweep in one population.
//
// The sweep is implemented as follows.  We take a joint allele frequency trajectory
// for the causal allele, specifying the allele's frequency in each population.
// (Or, we take parameters of the sweep and from them construct the frequency trajectory --
// how, is another subject).  The trajectory ends at <gen>, this sweep's ending generation,
// and starts at a generation at which the total frequency of the causal allele drops to zero
// (viewing time as going backwards.)
//
// We try to reuse as much of the neutral simulation machinery as possible when
// simulating sweeps.  At the end of the sweep, when we first encounter it during backwards
// simulation, we split each population into two: one will contain nodes from the original population
// that carry the derived allele at the causal mutation SNP, and the other nodes that carry the
// ancestral allele.  More specifically, for each population we create one new "companion" population
// and move the derived-allele nodes into it.  The original population, with the original <popid>,
// then contains nodes with the ancestral allele.  Keeping the original popid lets us catch events
// executed on the original population and forward them to the ancestral & derived subpops as appropriate.
// The number of nodes moved from the original pop to the companion derived-allele pop is determined by
// the frequency of the derived (causal) allele in the original pop at the end of the sweep,
// as determined by the frequency trajectory.
//
// (Let's call each original pop which we split into two 'orig-pop', and the two resulting pops
// 'anc-pop' and 'der-pop', with the understanding that 'anc-pop' is the same <Pop> object as 'orig-pop' --
// we just moved some of its <Nodes> into the companion 'der-pop'.)
//
// Then, we let the simulation run as usual, with the following adjustments:
//   - at each point of the frequency trajectory, for each orig-pop we adjust the relative sizes of the
//     anc-pop and der-pop to reflect the derived allele frequency in the orig-pop at that generation
//   - after each recombination or gene conversion, we take the resulting node that did NOT receive
//     the causal mutation location, and decide whether it has the derived or the ancestral allele
//     at the causal location by flipping a biased coin based on the frequency of the derived allele
//     (aka pop size of the der-pop); then, if the node is in der-pop and it was assigned the ancestral
//     allele, we move it to the anc-pop, and vice versa.
class Event_SweepOnePop: virtual public HistEvents::Event {
public:
	 Event_SweepOnePop( HistEvents *histEvents_, const string& label_, genid gen_, popid sweepPop_,
											gensInv_t selCoeff_, loc_t selPos_,
											freq_t final_sel_freq_ ):
		 Event( histEvents_, label_, gen_ ), sweepPop( sweepPop_ ),
		 selCoeff( selCoeff_ ), selPos( selPos_ ), final_sel_freq( final_sel_freq_ ),
		 ancPop( NULL ), derPop( NULL ), sel_leaves( LEAFSET_NULL ), sweepFinished( false )
			{ }
	 Event_SweepOnePop( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ),
																															ancPop( NULL ), derPop( NULL ),
																															sel_leaves( LEAFSET_NULL ), sweepFinished( false ) {
		 is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
	 }
	 static const char *typeStr() { return "sweep1"; }
	 virtual eventKind_t getEventKind() const { return E_SWEEP; }	 
				
	 virtual ~Event_SweepOnePop();

	 virtual genid execute();

	 virtual void processSimEnd( genid gen );


private:
	 // Field: sweepPop
	 // The sole population in which the sweep happens.
	 // The frequency of the derived allele in other populations always remains zero.
	 popid sweepPop;

	 // Field: selCoeff
	 // The selection coefficient.
	 gensInv_t selCoeff;

	 // Field: selPos
	 // The location of the causal mutation.
	 loc_t selPos;
	 
	 // Field: final_sel_freq
	 // The frequency of the derived allele at the end of the sweep.
	 freq_t final_sel_freq;

	 // Field: ancPop
	 // The Pop holding the nodes that carry the ancestral allele at the selected site;
	 // same object as the original Pop in which the sweep happens.
	 Pop *ancPop;

	 // Field: derPop
	 // The Pop holding the nodes that carry the derived allele at the selected site.
	 // Created when the sweep starts (pastward).
	 Pop *derPop;

	 // Field: freqTraj
	 // The frequency trajectory of the causal allele
	 FreqTrajP freqTraj;

	 // Field: sel_leaves
	 // Leaves which inherit the causal allele.
	 leafset_p sel_leaves;

	 // Field: totPopSize
	 // Total pop size of the population undergoing sweep.
	nchroms_t totPopSize;

	////////////////////////////////

	 // Class: SweepOnePopHook
	 // Callbacks invoked when certain events happen during a sweep.  We're using the standard
	 // neutral simulation machinery for most of the sweep simulation, but a few things are
	 // different -- this class implements these differences.
	 class SweepOnePopHook: public Hook {
	 public:
			SweepOnePopHook( Event_SweepOnePop *evt_): evt( evt_ ), inHook( False ) {}
			virtual ~SweepOnePopHook();
			virtual void handle_recomb( Node *, Node *, loc_t, genid);
			virtual void handle_gc( Node *, Node *, loc_t, loc_t, genid);
//			virtual void handle_coal( Node * );
			
	 private:
			// Private field: evt
			// The sweep event to which these hooks relate.
			Event_SweepOnePop *evt;

			// Private field: inHook
			// Whether we're already executing a handle_* method of this hook.
			// Used to prevent recursively invoking a handle_* method.
			bool_t inHook;

	 };  // class SweepOnePopHook

	 typedef boost::shared_ptr<SweepOnePopHook> SweepOnePopHookP;

	 // Field: hook
	 // Hook that makes adjustments to neutral coalescent machinery needed to simulate sweeps.
	 SweepOnePopHookP sweepHook;

	 bool sweepFinished;

	 // Private methodP: init_
	 // Called when the backward simulation first encounters the sweep
	 // (at the generation corresponding to the sweep's end).
	 // Splits each population existing at that time into anc-pop and der-pop
	 // containing nodes carrying the ancestral/derived allele at the causal site,
	 // and establishes hooks (callbacks) to be called during the simulation
	 // to handle events according to the sweep.
	 void init_();

	 bool_t initialized_() const { return ancPop != NULL; } 

	 // Private methodP: determineAlleleAtSelPos_
	 // Given a node (resulting from a recombination or gene conversion) for which
	 // <Node::segs> does not contain <Event_SweepOnePop::selPos>, choose the allele at selPos
	 // and if needed move node to the appropriate partial pop (anc-pop or der-pop).
	 void determineAlleleAtSelPos_( Node *node );
};  // class Event_SweepOnePop

//unsigned nrecomb_sel = 0, nrecomb_unsel = 0, nnew_sel = 0, nnew_unsel = 0, ncoal_sel = 0, ncoal_unsel = 0;

//
// Implementation of class: Event_SweepOnePop
//
const char *Event_SweepOnePop_typeStr() { return Event_SweepOnePop::typeStr(); }
HistEvents::EventP make_shared_Event_SweepOnePop( HistEvents *histEvents, istream& is );
HistEvents::EventP make_shared_Event_SweepOnePop( HistEvents *histEvents, istream& is ) {
  HistEvents::EventP ep( new Event_SweepOnePop( histEvents, is ) );
  return ep;
}

Event_SweepOnePop::~Event_SweepOnePop() {}

static filename_t sweep1_trajFN;

void sweep1_setTrajFN( filename_t fname ) {
	sweep1_trajFN = fname;
}

static double deltaTfactor(1.0);

void sweep1_set_deltaTfactor( double deltaTfactor_ ) { deltaTfactor = deltaTfactor_; }



// Private method: init_
// Called when the backward simulation first encounters the sweep
// (at the generation corresponding to the sweep's end).
// Splits each population existing at that time into anc-pop and der-pop
// containing nodes carrying the ancestral/derived allele at the causal site,
// and establishes hooks (callbacks) to be called during the simulation
// to handle events according to the sweep.
void Event_SweepOnePop::init_() {
	if ( sweep1_trajFN.empty() ) {
		nchroms_t popSize = getDemography()->dg_get_pop_by_name( sweepPop )->pop_get_size();
		freqTraj = boost::make_shared<DeterministicSweepTraj>( sweepPop, gen, selCoeff, final_sel_freq,
																													 popSize, deltaTfactor );
	} else
		freqTraj = boost::make_shared<TrajFromFile>( sweep1_trajFN, gen, final_sel_freq );

//	PRINT( freqTraj.get() );
	
	sel_leaves = make_empty_leafset();
	
	//
	// Split pop into two pops, anc-pop and der-pop.
	//

	Pop *pop = getDemography()->dg_get_pop_by_name( sweepPop );

	chkCond( pop != NULL, "pop is null" );
		
	nchroms_t final_pop_size = pop->pop_get_size();
	totPopSize = final_pop_size;
	freq_t cfreq = freqTraj->getCurFreq( pop->pop_get_name() );
	nchroms_t derPopSize( static_cast<double>(final_pop_size) * cfreq );
	std::string derPopLabel( pop->get_label() + "_sel" );

	ancPop = pop;

	// For nodes carrying the derived allele at selPos, create a new pop.
	// For nodes carrying the ancestral allele, reuse the original Pop.
	derPop =
		 getDemography()->dg_add_pop( boost::make_shared<Pop>( getDemography()->find_unused_popname(),
																													 derPopSize,
																													 derPopLabel ),
																	gen );
	
	
	getDemography()->dg_move_nodes_by_name( pop->pop_get_name(), derPop->pop_get_name(),
																					freqTraj->getCurFreq( pop->pop_get_name() ),
																					gen
																					IF_COSI_DEV( , /* exactFraction= */ True ) );
		
	ForEach (Node *it, derPop->getMembers() ) {
		ForEach( const seglist::Seg& seg, *it->getSegs() ) {
			if ( seg.contains( selPos ) )
				 sel_leaves = leafset_union( sel_leaves, seg.getLeafset() );
		}
	}
		
	// Set up a hook, so that we get called when certain events happen during the sweep
	// (specifically recombinations, gene conversions, and setting of migration rate.)
	Event_SweepOnePop *thisPtr = this;
	getDemography()->getHooks()->addHook( sweepHook = boost::make_shared<SweepOnePopHook>( thisPtr ) );

}  // Event_SweepOnePop::init()

//
// Method: execute
//
// Execute a selective sweep.  This method is invoked the first time the sweep is encountered during pastward simulation,
// and then each time the frequency of the selected allele in the population changes according to the frequency trajectory.
//
// When first called: splits each pop into der-pop and anc-pop; adds hooks to be called on recomb and gc events.
//
// When called for frequency changes:
// Updates the sizes of the der-pop and anc-pop making up each orig-pop
// (see <Event_SweepOnePop> for definitions of terms used here), based on the frequency of the selected allele in the orig-pop
// given by <freqTraj>.
//
genid Event_SweepOnePop::execute() {
	// If this is the first time this sweep is encountered,
	// initialize the processing of this sweep by splitting each orig-pop currently existing in the demographic model into
	// anc-pop and der-pop, etc.
	if ( !initialized_() ) init_();

	nchroms_t derPopSize( static_cast<double>( totPopSize ) * freqTraj->getCurFreq( ancPop->pop_get_name() ) );
		//nchroms_t ancSize = totSize - derSize;

	getDemography()->dg_set_pop_size_by_name( gen, derPop->pop_get_name(),
																						derPopSize );
		
	getDemography()->dg_set_pop_size_by_name( gen, ancPop->pop_get_name(),
																						totPopSize - derPopSize );

	genid curGen = this->gen;

	if ( derPopSize >= 1 ) {
		freqTraj->next();
		this->gen = freqTraj->getCurGen();
		addEvent( shared_from_this() );
	} else {

		assert( !this->sweepFinished );

		//
		// The sweep has completed.
		//

		// Force coalescence of any nodes carrying the derived allele, until there is at most one

		while( derPop->pop_get_num_nodes() > 1 ) {
			getDemography()->dg_coalesce_by_pop( derPop, curGen );
			curGen += gens_t( 1e-8 );
		}

		
		// Merge the split population back into a single pop.
		// Remove the hooks that got called after each recomb or gc event.
		// Put the selected mutation on the leaves.  (We could have done it at the start -- in the pastwards sense --
		// of the sweep, but there we did not know the generation at which the mutation originated, and
		// we want to be able to output the generation of each mutation.)
		//

#ifdef COSI_DEV_MUTCONTEXT
		mutcontext::saveMutContexts( derPop->pop_get_node( 0 )->getSegs(), selPos );
#endif		
		
		getDemography()->dg_move_nodes_by_name( derPop->pop_get_name(), ancPop->pop_get_name(), /* fractionToMove= */ 1.00, gen, /* exactFraction= */ True );
		chkCond( derPop->pop_get_num_nodes() == 0, "must clear derPop" );
		
		getDemography()->getHooks()->removeHook( sweepHook );

		getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
		//PRINT2( "sweep finished", gen );

		this->sweepFinished = true;

	}  // if sweep finished
	
	return curGen;
}  // Event_SweepOnePop::execute()

// Private method: determineAlleleAtSelPos_
// Given a node (resulting from a recombination or gene conversion) for which
// <Node::segs> does not contain <Event_SweepOnePop::selPos>, choose the allele at selPos
// and if needed move node to the appropriate partial pop (anc-pop or der-pop).
void Event_SweepOnePop::determineAlleleAtSelPos_( Node *node ) {
	Pop *curPop = node->getPop();
	if ( curPop == ancPop || curPop == derPop ) {
		//Pop *companionPop = ( curPop == ancPop ? derPop : ancPop );

		freq_t testFreq = double( derPop->pop_get_size() ) / double( totPopSize );
		{
			prob_t pval = random_double();
			Pop * shouldBeInPop = ( pval < testFreq ) ? derPop : ancPop;
			//( pval < testFreq ? nnew_sel : nnew_unsel )++;
			if ( curPop != shouldBeInPop ) {
//			PRINT5( "moving", curPop->pop_get_name(), shouldBeInPop->pop_get_name(), testFreq, pval );
				curPop->pop_remove_node( node );
				shouldBeInPop->pop_add_node( node );
			}
		}
	}  // if ( curPop == ancPop || curPop == derPop ) 
} // void Event_SweepOnePop::SweepOnePopHook::determineAlleleAtSelPos_( Node *node )

void Event_SweepOnePop::processSimEnd( genid gen ) {
	if ( !this->sweepFinished ) {
		getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
		getDemography()->getHooks()->removeHook( sweepHook );
		this->sweepFinished = true;
	}
}

	
//
// Implementation of class: Event_SweepOnePop::SweepOnePopHook
//

Event_SweepOnePop::SweepOnePopHook::~SweepOnePopHook() {
}


void Event_SweepOnePop::SweepOnePopHook::handle_recomb( Node *node1, Node *node2, loc_t loc, genid ) {
	evt->determineAlleleAtSelPos_( evt->selPos < loc ? node2 : node1 );
}
void Event_SweepOnePop::SweepOnePopHook::handle_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2,
																										genid ) {
	evt->determineAlleleAtSelPos_( loc1 <= evt->selPos && evt->selPos <= loc2 ? node2 : node1 );
}


// void Event_SweepOnePop::SweepOnePopHook::handle_coal( Node *n ) {
// //	( STLContains( evt->origPops, n->getPop()->pop_get_name() ) ? ncoal_unsel : ncoal_sel )++;
// }

}  // namespace sweep1



}  // namespace cosi
