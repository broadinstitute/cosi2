// Preamble
//         :PROPERTIES:
//         :ID:       f1b4c935-a3d5-48e2-9575-61d54ff5fca2
//         :END:


// [[file:sweep2.org::*Preamble][Preamble:1]]
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
#include <limits>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <cosi/decls.h>
#include <cosi/node.h>
#include <cosi/pop.h>
#include <cosi/demography.h>
#include <cosi/utils.h>
#include <cosi/mutlist.h>
#include <cosi/mutate.h>
#include <cosi/historical.h>
#include <cosi/traj.h>
#include <cosi/hooks.h>
#include <cosi/seglist.h>
#include <cosi/mutcontext.h>
#include <cosi/sweep2.h>
#include <cosi/coalrate.h>

namespace cosi {

namespace sweep2 {

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
// Preamble:1 ends here

// Class Event_SweepOnePop
//         :PROPERTIES:
//         :ID:       0ed47a21-1dd9-4e2a-9131-1ba029e7435a
//         :END:

//   A selective sweep in one population.

//         :DETAILS:
        
//         The sweep is implemented as follows.  We take a joint allele frequency trajectory
//         for the causal allele, specifying the allele's frequency in each population.
//         (Or, we take parameters of the sweep and from them construct the frequency trajectory --
//         how, is another subject).  The trajectory ends at <gen>, this sweep's ending generation,
//         and starts at a generation at which the total frequency of the causal allele drops to zero
//         (viewing time as going backwards.)
        
//         We try to reuse as much of the neutral simulation machinery as possible when
//         simulating sweeps.  At the end of the sweep, when we first encounter it during backwards
//         simulation, we split each population into two: one will contain nodes from the original population
//         that carry the derived allele at the causal mutation SNP, and the other nodes that carry the
//         ancestral allele.  More specifically, for each population we create one new "companion" population
//         and move the derived-allele nodes into it.  The original population, with the original <popid>,
//         then contains nodes with the ancestral allele.  Keeping the original popid lets us catch events
//         executed on the original population and forward them to the ancestral & derived subpops as appropriate.
//         The number of nodes moved from the original pop to the companion derived-allele pop is determined by
//         the frequency of the derived (causal) allele in the original pop at the end of the sweep,
//         as determined by the frequency trajectory.
        
//         (Let's call each original pop which we split into two 'orig-pop', and the two resulting pops
//         'anc-pop' and 'der-pop', with the understanding that 'anc-pop' is the same <Pop> object as 'orig-pop' --
//         we just moved some of its <Nodes> into the companion 'der-pop'.)
        
//         Then, we let the simulation run as usual, with the following adjustments:
        
//     - at each point of the frequency trajectory, for each orig-pop we adjust the relative sizes of the
//       anc-pop and der-pop to reflect the derived allele frequency in the orig-pop at that generation
                        
//     - after each recombination or gene conversion, we take the resulting node that did NOT receive
//       the causal mutation location, and decide whether it has the derived or the ancestral allele
//       at the causal location by flipping a biased coin based on the frequency of the derived allele
//       (aka pop size of the der-pop); then, if the node is in der-pop and it was assigned the ancestral
//       allele, we move it to the anc-pop, and vice versa.
                        
//     :END:                       
        

// [[file:sweep2.org::*Event_SweepOnePop][Event_SweepOnePop:1]]
class Event_SweepOnePop: virtual public HistEvents::Event {
// Event_SweepOnePop:1 ends here

// Constructor, destructor, housekeeping code

// [[file:sweep2.org::*Constructor,%20destructor,%20housekeeping%20code][Constructor\,\ destructor\,\ housekeeping\ code:1]]
public:
   Event_SweepOnePop( HistEvents *histEvents_, const string& label_, genid gen_, popid sweepPop_,
                      gensInv_t selCoeff_, loc_t selPos_,
                      freq_t final_sel_freq_ ):
     Event( histEvents_, label_, gen_ ), sweepPop( sweepPop_ ),
     selCoeff( selCoeff_ ), selPos( selPos_ ), final_sel_freq( final_sel_freq_ ),
     ancPop( NULL ), derPop( NULL ), sel_leaves( LEAFSET_NULL ),
     sweepStartTime( NULL_GEN ), sweepFinished( false )
      { }
   Event_SweepOnePop( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ),
                                                              ancPop( NULL ), derPop( NULL ),
                                                              sel_leaves( LEAFSET_NULL ),
                                                              sweepStartTime( NULL_GEN ),
         sweepFinished( false ) {
     is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
   }
   static const char *typeStr() { return "sweep2"; }
        
   virtual ~Event_SweepOnePop();
// Constructor\,\ destructor\,\ housekeeping\ code:1 ends here

// Method execute

//                         Execute a selective sweep.  This method is invoked the first time the sweep is encountered during pastward simulation,
//                         and then each time the frequency of the selected allele in the population changes according to the frequency trajectory.

//                         :DETAILS:
//                         When first called: splits each pop into der-pop and anc-pop; adds hooks to be called on recomb and gc events.
                        
//                         When called for frequency changes:
//                         Updates the sizes of the der-pop and anc-pop making up each orig-pop
//                         (see <Event_SweepOnePop> for definitions of terms used here), based on the frequency of the selected allele in the orig-pop
//                         given by <freqTraj>.
//                         :END:


// [[file:sweep2.org::*execute][execute:1]]
virtual genid execute()
// execute:1 ends here

// impl


// [[file:sweep2.org::*impl][impl:1]]
{
    // If this is the first time this sweep is encountered,
    // initialize the processing of this sweep by splitting each orig-pop currently existing in the demographic model into
    // anc-pop and der-pop, etc.
    //PRINT2( "in sweep2::execute", gen );
    genid curGen = gen;
    if ( !initialized_() ) {
      init_();
      
      // compute when sweep started
  
      // double epsilon = 1 / (2. * final_pop_size);
      // gens_t end_shift =  (final_sel_freq > 1-epsilon) ? ZERO_GENS :
      //   ( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / selCoeff );
      // genid tend = gen - ( 2 * log(epsilon) / selCoeff );
      
      // genid sweepStartGen = tend - end_shift;
      this->gen = this->sweepStartTime;
      addEvent( shared_from_this() );
      
      // this->gen = 
      // addEvent( 
      return curGen;
    }
    
    //
    // The sweep has completed.
    //
    
    // Force coalescence of any nodes carrying the derived allele, until there is at most one
    
    while( derPop->pop_get_num_nodes() > 1 ) {
      getDemography()->dg_coalesce_by_pop( derPop, curGen );
      curGen += gens_t( 1e-12 );
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
    ancPop->clearCoalArrivalProcess();
    derPop->clearCoalArrivalProcess();
    
    getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
    //PRINT2( "sweep finished", gen );
    
    this->sweepFinished = true;
    
    return curGen;
  }  // Event_SweepOnePop::execute()
// impl:1 ends here

// Method processSimEnd


// [[file:sweep2.org::*processSimEnd][processSimEnd:1]]
virtual void processSimEnd( genid gen )
// processSimEnd:1 ends here

// impl

// [[file:sweep2.org::*impl][impl:1]]
{
  if ( !this->sweepFinished ) {
    getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
    getDemography()->getHooks()->removeHook( sweepHook );
    this->sweepFinished = true;
  }
}
// impl:1 ends here

// private                                                                                                                                                                                                                                 :private:
//           :PROPERTIES:
//                 :ID:       2ccb5043-97f5-40b8-adb1-d189a762669c
//           :END:

// [[file:sweep2.org::*private][private:1]]
private:
// private:1 ends here

// Field sweepPop - The sole population in which the sweep happens.
//                                 :PROPERTIES:
//                                 :ID:       3c4bd8b9-f4c8-4a04-a10f-4573d4a33bc6
//                                 :END:

// [[file:sweep2.org::*sweepPop%20-%20The%20sole%20population%20in%20which%20the%20sweep%20happens.][sweepPop\ -\ The\ sole\ population\ in\ which\ the\ sweep\ happens\.:1]]
popid sweepPop;
// sweepPop\ -\ The\ sole\ population\ in\ which\ the\ sweep\ happens\.:1 ends here

// Field selCoeff - the selection coefficient
//                                 :PROPERTIES:
//                                 :ID:       25c7b4d7-d0c7-4a18-bdcf-5567f9b31d1b
//                                 :END:

// [[file:sweep2.org::*selCoeff%20-%20the%20selection%20coefficient][selCoeff\ -\ the\ selection\ coefficient:1]]
gensInv_t selCoeff;
// selCoeff\ -\ the\ selection\ coefficient:1 ends here

// Other fields
//                 :PROPERTIES:
//                                         :ID:       f9832672-559a-4d17-ae32-c07060b31a02
//                       :END:
                

// [[file:sweep2.org::*Other%20fields][Other\ fields:1]]
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
popsize_float_t totPopSize;

 // Field: sweepStartTime
 // Time when sweep started (going forward)
 genid sweepStartTime;

 typedef math::Function< genid, popsize_float_t, math::SweepPopSizeTraj > sweepTraj_t;
 typedef boost::shared_ptr<sweepTraj_t> sweepTraj_p;
 sweepTraj_p sweepTraj;
 bool sweepFinished;
// Other\ fields:1 ends here

// Class SweepOnePopHook
//                                 :PROPERTIES:
//                                 :ID:       de739866-b523-4774-8a6f-a422923663ee
//                                 :END:
                        
//                                 Callbacks invoked when certain events happen during a sweep.  We're using the standard
//                                 neutral simulation machinery for most of the sweep simulation, but a few things are
//                                 different -- this class implements these differences.


// [[file:sweep2.org::*SweepOnePopHook][SweepOnePopHook:1]]
class SweepOnePopHook: public Hook {
public:
   SweepOnePopHook( Event_SweepOnePop *evt_): evt( evt_ ), inHook( False ) {}
   virtual ~SweepOnePopHook() {}
   virtual void handle_recomb( Node *node1, Node *node2, loc_t loc, genid curGen ) {
      evt->determineAlleleAtSelPos_( evt->selPos < loc ? node2 : node1, curGen );
   }
   virtual void handle_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2,
                                                 genid curGen ) {
      evt->determineAlleleAtSelPos_( loc1 <= evt->selPos && evt->selPos <= loc2 ? node2 : node1,
                                     curGen );
   }
  
   
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
// SweepOnePopHook:1 ends here

// Method init_
//                                         :PROPERTIES:
//                                         :ID:       a972a681-98c2-46fa-9b92-2fbc8ca1d6d4
//                                         :END:


// [[file:sweep2.org::*init_][init_:1]]
void init_();
// init_:1 ends here

// Method initialized_
//                                         :PROPERTIES:
//                                         :ID:       0910734d-3108-4ef7-b4e2-825bee6d372a
//                                         :END:
                                

// [[file:sweep2.org::*initialized_][initialized_:1]]
bool_t initialized_() const { return ancPop != NULL; }
// initialized_:1 ends here

// Method determineAlleleAtSelPos_

//                                         Given a node (resulting from a recombination or gene conversion) for which
//                                         <Node::segs> does not contain <Event_SweepOnePop::selPos>, choose the allele at selPos
//                                         and if needed move node to the appropriate partial pop (anc-pop or der-pop).

// [[file:sweep2.org::*determineAlleleAtSelPos_][determineAlleleAtSelPos_:1]]
void determineAlleleAtSelPos_( Node *node, genid );
// determineAlleleAtSelPos_:1 ends here

// End of class Event_SweepOnePop


// [[file:sweep2.org::*End%20of%20class%20Event_SweepOnePop][End\ of\ class\ Event_SweepOnePop:1]]
};  // class Event_SweepOnePop
// End\ of\ class\ Event_SweepOnePop:1 ends here

// Misc routines
                

// [[file:sweep2.org::*Misc%20routines][Misc\ routines:1]]
static bool sweepFracSample = False;

void sweep2_set_sweepFracSample( bool sweepFracSample_ ) {
        sweepFracSample = sweepFracSample_;
}

//unsigned nrecomb_sel = 0, nrecomb_unsel = 0, nnew_sel = 0, nnew_unsel = 0, ncoal_sel = 0, ncoal_unsel = 0;

const char *Event_SweepOnePop2_typeStr() { return Event_SweepOnePop::typeStr(); }
HistEvents::EventP make_shared_Event_SweepOnePop2( HistEvents *histEvents, istream& is );
HistEvents::EventP make_shared_Event_SweepOnePop2( HistEvents *histEvents, istream& is ) {
  HistEvents::EventP ep( new Event_SweepOnePop( histEvents, is ) );
  return ep;
}

Event_SweepOnePop::~Event_SweepOnePop() {}

static filename_t sweep2_trajFN;

void sweep2_setTrajFN( filename_t fname ) {
        sweep2_trajFN = fname;
}

void Event_SweepOnePop::init_() {


if ( sweep2_trajFN.empty() ) {
        nchroms_t popSize = getDemography()->dg_get_pop_by_name( sweepPop )->pop_get_size();
        freqTraj = boost::make_shared<DeterministicSweepTraj>( sweepPop, gen, selCoeff, final_sel_freq,
                                                                                                                                                                                                                                 popSize );
} else
        freqTraj = boost::make_shared<TrajFromFile>( sweep2_trajFN, gen, final_sel_freq );

nchroms_t popSizeNow = getDemography()->dg_get_pop_by_name( sweepPop )->pop_get_size();
BOOST_AUTO( freqTraj2,
                                                ( boost::make_shared<DeterministicSweepTraj>( sweepPop, gen, selCoeff, final_sel_freq,
                                                                                                                                                                                                                                        popSizeNow ) ) );

this->sweepStartTime = gen;
while ( !freqTraj2->done() ) {
        this->sweepStartTime = freqTraj2->getCurGen();
        freqTraj2->next();
}


//PRINT( this->sweepStartTime );

sel_leaves = make_empty_leafset();

                                
        Pop *pop = getDemography()->dg_get_pop_by_name( sweepPop );

        chkCond( pop != NULL, "pop is null" );
                
        nchroms_t final_pop_size = pop->pop_get_size();
        totPopSize = static_cast<popsize_float_t>( final_pop_size );
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
                                                                                                                                                                        gen, /* exactFraction= */ sweepFracSample
#if 0                                                                                                                                                                   
                                                                                                                                                                        IF_COSI_DEV( , /* exactFraction= */ True )
#endif                                                                                                                                                                  
                );
        
                
        ForEach (Node *it, derPop->getMembers() ) {
                ForEach( const seglist::Seg& seg, *it->getSegs() ) {
                        if ( seg.contains( selPos ) )
                                 sel_leaves = leafset_union( sel_leaves, seg.getLeafset() );
                }
        }
        //PRINT2( derPop->getMembers().size(), leafset_size( sel_leaves ) );
        getDemography()->getHooks()->fire_sweep_end( sel_leaves );
                
        // Set up a hook, so that we get called when certain events happen during the sweep
        // (specifically recombinations, gene conversions, and setting of migration rate.)
        Event_SweepOnePop *thisPtr = this;
        getDemography()->getHooks()->addHook( sweepHook = boost::make_shared<SweepOnePopHook>( thisPtr ) );

        double epsilon = 1. / (2. * final_pop_size);
        gens_t end_shift =  (final_sel_freq > 1-epsilon) ? ZERO_GENS :
                 ( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / selCoeff );
        genid tend = gen - ( 2 * log(epsilon) / selCoeff );

        sweepTraj = boost::make_shared<sweepTraj_t>(
                                                                        popsize_float_t( final_pop_size ),
                                                                        selCoeff,
                                                                        epsilon,
                                                                        end_shift,
                                                                        tend
                );

//      genid tstart = 
        

        derPop->setCoalArrivalProcess(
                math::ArrivalProcess< genid, math::Any< RandGen > >(
                        math::makeNonHomogeneousPoissonProcess
                        (
                                math::coalRateFunction(
                                        *sweepTraj

																	), gen, "coal in derPop"
                                )));

        BOOST_AUTO( sweepPopSizeComplementTraj,
                                                        (
                                        math::Function< genid, popsize_float_t, math::SweepPopSizeComplementTraj >(
                                                popsize_float_t( final_pop_size ),
                                                selCoeff,
                                                epsilon,
                                                end_shift,
                                                tend
                                                )
                                                                ));

        ancPop->setCoalArrivalProcess(
                math::ArrivalProcess< genid, math::Any< RandGen > >(
                        math::makeNonHomogeneousPoissonProcess
                        (
                                math::coalRateFunction(
                                        sweepPopSizeComplementTraj
                                        )
                                , gen, "coal in ancPop"
                                )));





#if !defined(NDEBUG) && defined(COSI_DEV_PRINT)
                BOOST_AUTO( coalRate, math::coalRateFunction( *sweepTraj ) );
                BOOST_AUTO( coalRateComplement, math::coalRateFunction( sweepPopSizeComplementTraj ) );

                PRINT( coalRate( genid( 370.0 ) ) ); 
                PRINT( coalRateComplement( genid( 370.0 ) ) ); 

                BOOST_AUTO( coalRateIntegral, math::integralFromPoint( coalRate, genid(100.0) ) );
                BOOST_AUTO( coalRateComplementIntegral, math::integralFromPoint( coalRateComplement, genid(100.0) ) );

                PRINT( coalRateIntegral( genid( 370.0 ) ) );
                PRINT(( math::integrateNumerically( coalRate, genid(100.0), genid( 370.0 ), 10 ) ));
                // PRINT(( math::integrateNumerically( coalRate, genid(0), genid( 1.0 ), 10 ) ));
                // PRINT(( math::integrateNumerically( coalRate, genid(1.0), genid( 370.0 ), 10 ) ));
                PRINT( coalRateComplementIntegral( genid( 370.0 ) ) );
                PRINT(( math::integrateNumerically( coalRateComplement, genid(100.0), genid( 370.0 ), 10 ) ));
                // PRINT(( math::integrateNumerically( coalRateComplement, genid(0), genid( 1.0 ), 10 ) ));
                // PRINT(( math::integrateNumerically( coalRateComplement, genid(1.0), genid( 370.0 ), 10 ) ));
#endif // #if !defined(NDEBUG) && defined(COSI_DEV_PRINT)               




#if 0   
        if( 0 ) {
                BOOST_AUTO( coalRate, math::coalRateFunction( sweepPopSizeTraj ) );
                BOOST_AUTO( coalRateComplement, math::coalRateFunction( sweepPopSizeComplementTraj ) );

                BOOST_AUTO( coalRateIntegral, math::indefiniteIntegral( coalRate ) );
                BOOST_AUTO( coalRateComplementIntegral, math::indefiniteIntegral( coalRateComplement ) );

                PRINT2( coalRateIntegral( genid(0.0) ), coalRateComplementIntegral( genid(0.0) ) );

                std::cerr.precision(10);
                genid a(127.3737), b( 537.321 );

                PRINT(( math::integrateNumerically( coalRateComplement, a, b, 3 ) ));
                PRINT(( math::integrateNumerically( coalRateComplement, a, b, 4 ) ));
                PRINT(( math::integrateNumerically( coalRateComplement, a, b, 5 ) ));
                PRINT(( math::integrateNumerically( coalRateComplement, a, b, 10 ) ));
                
                PRINT(( coalRateComplementIntegral( genid( b ) ) - coalRateComplementIntegral( genid( a ) ) ) );
                
                
                //assert(false);


                gens_t delta(0.01);
                typedef BOOST_TYPEOF( coalRateIntegral( boost::declval<genid>() ) ) integral_t;
                integral_t coalRateManualIntegral(0.0);
                integral_t coalRateComplementManualIntegral(0.0);
                
                for ( genid g(0.0); g < genid(1000.0); g += delta ) {
                        PRINT10( g, sweepPopSizeTraj(g), sweepPopSizeComplementTraj(g),
                                                         sweepPopSizeTraj(g) + sweepPopSizeComplementTraj(g),
                                                         coalRate(g), coalRateComplement(g),
                                                         coalRateIntegral(g), coalRateManualIntegral,
                                                         coalRateComplementIntegral(g), coalRateComplementManualIntegral );

                        coalRateManualIntegral += delta * coalRate(g + 0.5*delta);
                        PRINT( coalRateManualIntegral - coalRateIntegral(g+delta) );

                        coalRateComplementManualIntegral += delta * coalRateComplement(g + 0.5*delta);
                        PRINT( coalRateComplementManualIntegral - coalRateComplementIntegral(g+delta) );
                        

                        assert( coalRate(g) > ZERO );
                        assert( coalRateComplement( g ) > ZERO );

                        assert( coalRateIntegral( g ) > ZERO );

                        assert(( equal_eps( math::integrateNumerically( coalRateComplement, genid(0.0), g, 10 ),
                                                                                                        coalRateComplementIntegral(g) - coalRateComplementIntegral(genid(0.0)) ) ));
                        

                        if ( g > ZERO ) {
                                assert( math::integrateNumerically( coalRateComplement, genid(0.0), g, 10 ) > ZERO );
                        }
                        
                        //assert( coalRateComplementIntegral( g ) > ZERO );
                        if ( g > ZERO_GEN ) {
                                assert( coalRateIntegral(g) > coalRateIntegral(g-gens_t(1)) );
                                assert( coalRateComplementIntegral(g) > coalRateComplementIntegral(g-gens_t(1)) );
                        }
                }
        }
#endif // #if 0 

}  // Event_SweepOnePop::init()

void Event_SweepOnePop::determineAlleleAtSelPos_( Node *node, genid curGen ) {
	if ( node ) {
		Pop *curPop = node->getPop();
		if ( curPop == ancPop || curPop == derPop ) {
			//Pop *companionPop = ( curPop == ancPop ? derPop : ancPop );

			freq_t testFreq = eval( *sweepTraj, curGen )  / totPopSize;
			{
				prob_t pval = random_double();
				Pop * shouldBeInPop = ( pval < testFreq ) ? derPop : ancPop;
				//( pval < testFreq ? nnew_sel : nnew_unsel )++;
				if ( curPop != shouldBeInPop ) {
//                      PRINT5( "moving", curPop->pop_get_name(), shouldBeInPop->pop_get_name(), testFreq, pval );
					curPop->pop_remove_node( node );
					shouldBeInPop->pop_add_node( node );
				}
			}
		}  // if ( curPop == ancPop || curPop == derPop ) 
	}
} // void Event_SweepOnePop::SweepOnePopHook::determineAlleleAtSelPos_( Node *node )
// Misc\ routines:1 ends here

// Postamble
//         :PROPERTIES:
//         :ID:       eee793ed-1ebf-44fd-a409-b55cb2ce5b30
//         :END:


// [[file:sweep2.org::*Postamble][Postamble:1]]
}  // namespace sweep2



}  // namespace cosi
// Postamble:1 ends here
