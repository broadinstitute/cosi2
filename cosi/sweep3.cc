
// #+TITLE: sweep3 - implementation of selective sweep in one population 

// * Preamble
// ** Includes

//#define COSI_DEV_PRINT

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <ios>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <limits>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/bind.hpp>
#include <boost/none.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/filesystem/fstream.hpp>
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
#include <cosi/sweep3.h>
#include <cosi/coalrate.h>
#include <cosi/genmap.h>

namespace cosi {

bool sweep3_no_oneSidedRecombs = False;

extern double poisPrec;
extern unsigned int poisMaxSteps;

namespace sweep3 {

#define ForEach BOOST_FOREACH  

class Event_SweepOnePop;
boost::shared_ptr<Event_SweepOnePop> theSweep;

using std::vector;
using std::map;
using std::set;
using std::ifstream;
using std::ofstream;
using util::map_get;
using util::STLContains;
using util::isize;
using node::Node;
using node::NodeList;

template <typename T> 
class PartialSumTree {
public:
	 typedef T value_type;
	 
	 PartialSumTree(): partialSums( 9, T(0.0) )  { partialSums.reserve( 16385 ); }

	 // Method: add
	 // To the weight of the object at index 'itemIdx', add the value 'delta'
	 // (which can be positive or negative).  Note that the weight of the
	 // object with index 0 is always 0 and cannot be changed; thus, object
	 // indexing effectively starts at 1.
   void add( idx_t itemIdx, T delta ) {
		 assert( itemIdx > 0 );
		 ensureCapacity( itemIdx+1 );
		 assert( itemIdx < static_cast< idx_t >( partialSums.size() ) );

		 while (itemIdx < isize( partialSums ) ){
			 partialSums[itemIdx] += delta;
			 itemIdx += (itemIdx & -itemIdx);
		 }

#ifndef NDEBUG
#if 0		 
		 T chkSum(0.0);
		 for ( size_t i = 1; i < partialSums.size(); i++ ) {
			 chkSum += read( i );
			 assert( chkSum == readCumulative( i ) );
		 }
#endif		 
#endif		 
		 
#if 0		 
	   for (; itemIdx < util::isize( partialSums ); itemIdx += ( itemIdx & -itemIdx ) ) {
				partialSums[ itemIdx ] += delta;
				assert( partialSums[ itemIdx ] >= T(0.0) );
		 }
#endif		 
   }

	 void set( idx_t itemIdx, T val ) {
		 add( itemIdx, val - read( itemIdx ) );
		 assert( equal_eps( read( itemIdx ), val, 1e-10 ) );
	 }

	 T readCumulative(idx_t idx) const {
		 T sum( 0.0 ) ;
		 while (idx > 0){
			 sum += partialSums[idx];
			 idx -= (idx & -idx);
		 }
		 return sum;
	 }

	 T read( idx_t idx ) const {
		 assert( idx > 0 );
		 const_cast< PartialSumTree * >(this)->ensureCapacity( idx+1 );
		 return readCumulative( idx ) - readCumulative( idx-1 );

		 // copied from TopCoder http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=binaryIndexedTrees
		 T sum = partialSums[idx]; // sum will be decreased
		 if (idx > 0){ // special case
			 idx_t z = idx - (idx & -idx); // make z first
			 idx--; // idx is no important any more, so instead y, you can use idx
			 while (idx != z){ // at some iteration idx (y) will become z
				 sum -= partialSums[idx]; 
// substruct tree frequency which is between y and "the same path"
				 idx -= (idx & -idx);
			 }
		 }
		 if ( !( equal_eps( sum, readCumulative( idx ) - readCumulative( idx-1 ), 1e-10 ) ) ) {
			 PRINT10( idx, sum, read( idx-1 ), read( idx ), read( idx+1 ), readCumulative( idx+1 ), readCumulative(idx), readCumulative(idx-1), readCumulative(idx)-readCumulative(idx-1),
							readCumulative(idx) - readCumulative( idx-1 ) - sum );
			 std::cerr << "psums: ";
			 std::ostream_iterator<glen_t> out_it (std::cerr,", ");
			 std::copy( partialSums.begin(), partialSums.end(), out_it );
			 std::cerr << "\n";
		 }
		 assert( equal_eps( sum, readCumulative( idx ) - readCumulative( idx-1 ), 1e-10 ) );
		 return sum;
	 }

	 // Method: getSotalSum
	 // Returns the sum of weights of all objects.  Recall that the
	 // weights of all objects are initially zero until changed
	 // by the <add()> method, so objects whose weights were never changed
	 // count as weight zero.
	 T getTotalSum() const { return partialSums.back(); }

	 // Method: ensureCapacity
	 // Prepares this PartialSumTree for storing weights of objects with indices
	 // of at least 'size' by pre-allocating internal storage.  Like <std::vector::ensureCapacity>,
	 // affects only performance and not functionality.
	 void ensureCapacity( size_t size ) {
		 while ( partialSums.size() <= size+3 ) {
			 size_t oldSize = partialSums.size()-1;
			 if ( size+2 >= oldSize ) {
				 size_t newSize = oldSize * 2;
				 partialSums.resize( newSize+1, T( 0.0 ) );
				 partialSums[ newSize ] = partialSums[ oldSize ];
			 }
		 }
	 }

	 //
	 // Method: findCumulativeFraction
	 //
	 // Given a fraction of the sum of all weights, find the largest index i
	 // such that the sum of weights through i divided by <getTotalSum()>
	 // is less than that fraction.  Also returns the residue, such that
	 // adding the residue to the sum of weights through i and dividing by
	 // <getTotalSum()> equals 'cumulativeFraction'.
	 //
	 // Input params:
	 //
	 //     cumulativeFraction - specifies a fraction of <getTotalSum()>
	 //
	 // Returns:
	 //
	 //     the highest index i such that the sum of weights of elements 1..i
	 //     is less than cumulativeFraction*getTotalSum().   Also, in 'residue'
	 //     stores the value r such that the sum of weights of elements 1..i
	 //     plus r equals cumulativeFraction*getTotalSum().
	 //
	 idx_t findCumulativeFraction( frac_t cumulativeFraction, T *residue ) const {
	
		 idx_t searchRegionBase = 0;
		 int searchRegionHalfSize = partialSums.size() >> 1;
		 T amt = cumulativeFraction * getTotalSum();
		 //T amtOrig = amt;
		 while( searchRegionHalfSize > 0 ) {
			 idx_t searchRegionMid = searchRegionBase + searchRegionHalfSize;
			 T amtInFirstHalf = partialSums[ searchRegionMid ];
			 if ( amt > amtInFirstHalf ) {
				 amt -= amtInFirstHalf;
				 searchRegionBase = searchRegionMid;
			 }
			 searchRegionHalfSize >>= 1;
		 }
		 idx_t itemIdx = searchRegionBase+1;
		 *residue = amt;
		 return itemIdx;
	 }  // findCumulativeFraction

	 const std::vector<T>& getPartialSums() const { return partialSums; }
	 
private:
	 // Field: partialSums
	 // Partial sums of sub-ranges of elements of the underlying vector of weights.
	 // (The vector of weights itself is not stored, though can be recovered from <partialSums>).
	 // See the references in the <PartialSumTree> class comment for details.
	 std::vector<T> partialSums;
};  // class PartialSumTree



// * Class Event_SweepOnePop

//   A selective sweep in one population.

// ** details
//   
  
//   The sweep is implemented as follows.  We take a joint allele frequency trajectory
//   for the causal allele, specifying the allele's frequency in each population.
//   (Or, we take parameters of the sweep and from them construct the frequency trajectory --
//   how, is another subject).  The trajectory ends at <gen>, this sweep's ending generation,
//   and starts at a generation at which the total frequency of the causal allele drops to zero
//   (viewing time as going backwards.)
  
//   We try to reuse as much of the neutral simulation machinery as possible when simulating
//   sweeps.  At the end of the sweep, when we first encounter it during backwards simulation,
//   we split each population into two: one will contain nodes from the original population
//   that carry the derived allele at the causal mutation SNP, and the other nodes that carry
//   the ancestral allele.  More specifically, for each population we create one new
//   "companion" population and move the derived-allele nodes into it.  The original
//   population, with the original <popid>, then contains nodes with the ancestral allele.
//   Keeping the original popid lets us catch events executed on the original population and
//   forward them to the ancestral & derived subpops as appropriate.  The number of nodes
//   moved from the original pop to the companion derived-allele pop is determined by the
//   frequency of the derived (causal) allele in the original pop at the end of the sweep, as
//   determined by the frequency trajectory.
  
//   (Let's call each original pop which we split into two 'orig-pop', and the two resulting
//   pops 'anc-pop' and 'der-pop', with the understanding that 'anc-pop' is the same <Pop>
//   object as 'orig-pop' -- we just moved some of its <Nodes> into the companion 'der-pop'.)
  
//   Then, we let the simulation run as usual, with the following adjustments:
  
//     - at each point of the frequency trajectory, for each orig-pop we adjust the relative sizes of the
//       anc-pop and der-pop to reflect the derived allele frequency in the orig-pop at that generation
      
//     - after each recombination or gene conversion, we take the resulting node that did NOT receive
//       the causal mutation location, and decide whether it has the derived or the ancestral allele
//       at the causal location by flipping a biased coin based on the frequency of the derived allele
//       (aka pop size of the der-pop); then, if the node is in der-pop and it was assigned the ancestral
//       allele, we move it to the anc-pop, and vice versa.
      
class Event_SweepOnePop: virtual public HistEvents::Event {

// *** Constructor, destructor, housekeeping code

public:
   Event_SweepOnePop( HistEvents *histEvents_, const string& label_, genid gen_, popid sweepPop_,
                      gensInv_t selCoeff_, loc_t selPos_,
                      freq_t final_sel_freq_ ):
     Event( histEvents_, label_, gen_ ), sweepPop( sweepPop_ ),
     selCoeff( selCoeff_ ), selPos( selPos_ ), final_sel_freq( final_sel_freq_ ),
     ancPop( NULL ), derPop( NULL ), sel_leaves( LEAFSET_NULL ),
     sweepStartTime( NULL_GEN ), sweepFinished( false ), sideRecombPop( NULL )
      { }
   Event_SweepOnePop( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ),
                                                              ancPop( NULL ), derPop( NULL ),
                                                              sel_leaves( LEAFSET_NULL ),
                                                              sweepStartTime( NULL_GEN ),
																															sweepFinished( false ) {
     is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
   }
   static const char *typeStr() { return "sweep"; }
	 virtual eventKind_t getEventKind() const { return E_SWEEP; }	 

	 virtual ~Event_SweepOnePop();

// **** Method execute

//      Execute a selective sweep.  This method is invoked the first time the sweep is encountered during pastward simulation,
//      and then each time the frequency of the selected allele in the population changes according to the frequency trajectory.
		 
//      When first called: splits each pop into der-pop and anc-pop; adds hooks to be called on recomb and gc events.
     
//      When called for frequency changes:
//      Updates the sizes of the der-pop and anc-pop making up each orig-pop
//      (see <Event_SweepOnePop> for definitions of terms used here), based on the frequency of the selected allele in the orig-pop
//      given by <freqTraj>.

	 virtual genid execute()
	 // ***** impl
			{
				// If this is the first time this sweep is encountered, initialize the processing of this sweep
				// by splitting each orig-pop currently existing in the demographic model into anc-pop and
				// der-pop, etc.  PRINT2( "in sweep3::execute", gen );
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
					PRINT4( "end_of_sweep", gen, derPop->pop_get_num_nodes(), ancPop->pop_get_num_nodes() );
					getDemography()->dg_coalesce_by_pop( derPop, curGen, /* forceCoalescence= */ True );
					curGen += gens_t( 1e-12 );
				}
    
				// Merge the split population back into a single pop.  Remove the hooks that got called after each
				// recomb or gc event.  Put the selected mutation on the leaves.  (We could have done it at the
				// start -- in the pastwards sense -- of the sweep, but there we did not know the generation at
				// which the mutation originated, and we want to be able to output the generation of each
				// mutation.)
				//
    
#ifdef COSI_DEV_MUTCONTEXT
				mutcontext::saveMutContexts( derPop->pop_get_node( 0 )->getSegs(), selPos );
#endif    
    
				getDemography()->dg_move_nodes_by_name( derPop->pop_get_name(),
																								ancPop->pop_get_name(),
																								/* fractionToMove= */ 1.00, gen, /* exactFraction= */ True );
				chkCond( derPop->pop_get_num_nodes() == 0, "must clear derPop" );
    
				getDemography()->getHooks()->removeHook( sweepHook );
				ancPop->clearCoalArrivalProcess();
				derPop->clearCoalArrivalProcess();
				ancPop->clearPopListener();
				derPop->clearPopListener();
    
				getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
				PRINT2( "sweep finished", gen );

				//getDemography()->getNodePool()->setSelPos( boost::none );
    
				this->sweepFinished = true;
				theSweep.reset();
    
				return curGen;

			}  // Event_SweepOnePop::execute()

// **** Method processSimEnd

	 virtual void processSimEnd( genid gen )

	 // ***** impl

			{
				if ( !this->sweepFinished ) {
					getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);
					getDemography()->getHooks()->removeHook( sweepHook );
					this->sweepFinished = true;
				}

			}

// *** private

private:

// **** Field sweepPop - The sole population in which the sweep happens.
	 popid sweepPop;

// **** Field selCoeff - the selection coefficient
	 gensInv_t selCoeff;

// **** Other fields

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

	 typedef math::Function< genid, popsize_float_t, math::Piecewise< math::Const<> > > sweepTraj_t;
	 typedef boost::shared_ptr<sweepTraj_t> sweepTraj_p;
	 sweepTraj_p sweepTraj;

	 bool sweepFinished;

	 PartialSumTree<glen_t> sideRecombRates_der, sideRecombRates_anc;
	 std::vector< boost::optional<glen_t> > sideRecombVals_der, sideRecombVals_anc;

	 typedef math::ArrivalProcess<genid, math::Any< RandGen > > side_recomb_process_type;
	 side_recomb_process_type sideRecombProcess_der, sideRecombProcess_anc;

// ***** Class SweepOnePopHook

//       Callbacks invoked when certain events happen during a sweep.  We're using the standard
//       neutral simulation machinery for most of the sweep simulation, but a few things are
//       different -- this class implements these differences.
			

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

	 class SweepOnePopListener: public Pop::PopListener {
	 public:
			SweepOnePopListener( Event_SweepOnePop *evt_ ): evt( evt_ ) { }
			virtual ~SweepOnePopListener() { }
			virtual void nodeAdded( Pop *pop, nodeid nodeName, idx_t idx, loc_t beg, loc_t end ) {
				(void)nodeName;
				PartialSumTree<glen_t> *pst = ( pop == evt->derPop ? &(evt->sideRecombRates_der) :
																							( pop == evt->ancPop ? &(evt->sideRecombRates_anc) : NULL ) );
				PRINT7( "nodeAdded", pop->pop_get_name(), pop->getMembers().size(), nodeName, idx, beg, end );
				if ( pst ) {

					std::vector< boost::optional<glen_t> > *srv =
						 ( pop == evt->derPop ? &(evt->sideRecombVals_der) : &(evt->sideRecombVals_anc) );

					
					idx_t nodeidx = idx + 1;
					assert( nodeidx > 0 );
					while ( int( srv->size() ) <= nodeidx+5 )
						 srv->resize( 2 * std::max( srv->size(), static_cast<size_t>(1024) ) );

					assert( !boost::get_pointer( (*srv)[ nodeidx ] ) );
					if ( !equal_eps( pst->read( idx+1 ), glen_t( 0.0 ), 1e-10 ) ) {
						PRINT3( idx+1, nodeidx, pst->read( idx+1 ) );
					}
					assert( equal_eps( pst->read( idx+1 ), glen_t( 0.0 ), 1e-10 ) );
					
					PRINT6( beg, end, evt->selPos, evt->selPos_gloc, get_gloc( beg ), get_gloc( end ) );
					
					if ( evt->selPos < beg ) {
						(*srv)[ nodeidx ] = get_gloc( beg ) - evt->selPos_gloc;
						assert( get( (*srv)[ nodeidx ] ) > glen_t( 0.0 ) );
						pst->set( nodeidx, boost::get( (*srv)[ nodeidx ] ) );
					} else if ( evt->selPos > end ) {
						(*srv)[ nodeidx ] = evt->selPos_gloc - get_gloc( end ); 
						assert( get( (*srv)[ nodeidx ] ) > glen_t( 0.0 ) );
						pst->set( nodeidx, boost::get( (*srv)[ nodeidx ] ) );
					} else
						 (*srv)[ nodeidx ] = glen_t( 0.0 );

					assert( boost::get_pointer( (*srv)[ nodeidx ] ) );

					if ( !( equal_eps( pst->read( idx+1 ), get( (*srv)[ idx+1 ] ), 1e-10 )  ) ) {
						PRINT3( idx, pst->read( idx+1 ), get( (*srv)[ idx+1 ] ) );
					}
					assert( equal_eps( pst->read( idx+1 ) , get( (*srv)[ idx+1 ] ), 1e-10 ) );
					
					
					PRINT3( idx+1, srv->size(), boost::get( (*srv)[ idx+1 ] ) );
				} else {
					PRINT( "added to unknownPop!" );
				}
			}
			virtual void nodeRemoved( Pop *pop, idx_t idx ) {
				PRINT3( "nodeRemoved", pop->pop_get_name(), idx );
				PartialSumTree<glen_t> *pst = ( pop == evt->derPop ? &(evt->sideRecombRates_der) :
																							( pop == evt->ancPop ? &(evt->sideRecombRates_anc) : NULL ) );
				assert( pst );
				if ( pst ) {
					std::vector< boost::optional<glen_t> > *srv = ( pop == evt->derPop ? &(evt->sideRecombVals_der) : &(evt->sideRecombVals_anc) );
					assert( boost::get_pointer( (*srv)[ idx+1 ] ) );
					assert( equal_eps( pst->read( idx+1 ), get( (*srv)[ idx+1 ] ), 1e-10 ) );
					PRINT3( idx+1, srv->size(), get( (*srv)[ idx+1 ] ) );
					assert( boost::get( (*srv)[ idx+1 ] ) >= glen_t(0.0) );
					if ( boost::get( (*srv)[ idx+1 ] ) > glen_t(0.0) )
						 pst->set( idx+1, glen_t( 0.0 ) );
					(*srv)[ idx+1 ] = boost::optional<glen_t>();
					assert( equal_eps( pst->read( idx+1 ), glen_t( 0.0 ), 1e-10 ) );
					
				} else {
					PRINT( "removed from unknownPop!" );
				}
			}

	 private:
			// Private field: evt
			// The sweep event to which these hooks relate.
			Event_SweepOnePop *evt;
	 };
  
	 typedef boost::shared_ptr<SweepOnePopHook> SweepOnePopHookP;
  
// ****** Field: hook
// Hook that makes adjustments to neutral coalescent machinery needed to simulate sweeps.
	 SweepOnePopHookP sweepHook;

	 typedef boost::shared_ptr<SweepOnePopListener> SweepOnePopListenerP;
	 SweepOnePopListenerP popListener;

// ****** Method init_
	 void init_();

// ****** Method initialized_
	 bool_t initialized_() const { return ancPop != NULL; }

// ****** Method determineAlleleAtSelPos_

//        Given a node (resulting from a recombination or gene conversion) for which
//        <Node::segs> does not contain <Event_SweepOnePop::selPos>, choose the allele at selPos
//        and if needed move node to the appropriate partial pop (anc-pop or der-pop).
			 
	 void determineAlleleAtSelPos_( Node *node, genid );

public:
	 gens_t getSideRecombWaitTime( genid gen_, gens_t maxWaitTime  ) {
		 PRINT5( "getSideRecombWaitTime", gen_, maxWaitTime, initialized_(), sweepFinished );
		 gens_t gens_inf( std::numeric_limits<double>::infinity() );
		 if ( !initialized_() || sweepFinished ) return gens_inf;
		 genid genid_inf( std::numeric_limits<double>::infinity() );
		 PRINT( sideRecombRates_der.getPartialSums().size() ); 
		 PRINT( sideRecombRates_anc.getPartialSums().size() );
#ifdef COSI_DEV_PRINT
		 if ( !cosi::util::noDbgPrint ) {
			 std::ostream_iterator<glen_t> out_it (std::cerr,", ");
			 std::cerr << "ratesDer: ";
			 std::copy( sideRecombRates_der.getPartialSums().begin(),
									sideRecombRates_der.getPartialSums().end(),
									out_it );
			 std::cerr << "\n";
			 std::cerr << "ratesAnc: ";
			 std::copy( sideRecombRates_anc.getPartialSums().begin(),
									sideRecombRates_anc.getPartialSums().end(),
									out_it );
			 std::cerr << "\n";
		 }
#endif		 
		 
								
		 glen_t derTot = sideRecombRates_der.getTotalSum();
		 PRINT( derTot );
		 glen_t ancTot = sideRecombRates_anc.getTotalSum();
		 PRINT( ancTot );
		 glen_t glen_eps( 1e-10 );
		 PRINT( glen_eps );
		 if ( derTot < glen_eps && ancTot < glen_eps ) return gens_inf;

		 sideRecombPop = NULL;
		 PRINT4( gen_, maxWaitTime, ancTot, ancTot < glen_eps );
		 genid nextCoalTime = ancTot < glen_eps ?  genid_inf : sideRecombProcess_der.nextArrivalTime( gen_, gen_ + maxWaitTime, ToDouble( ancTot ), *getRandGen(), /* eps= */ poisPrec, /* maxSteps= */ poisMaxSteps );
		 if ( nextCoalTime < gen_ + maxWaitTime ) {
			 maxWaitTime = nextCoalTime - gen_;
			 sideRecombPop = ancPop;
		 }
		 PRINT4( gen_, maxWaitTime, derTot, derTot < glen_eps );
		 nextCoalTime = derTot < glen_eps ? genid_inf : sideRecombProcess_anc.nextArrivalTime( gen_, gen_ + maxWaitTime, ToDouble( derTot ), *getRandGen(), /* eps= */ poisPrec, /* maxSteps= */ poisMaxSteps );
		 if ( nextCoalTime < gen_ + maxWaitTime ) {
			 maxWaitTime = nextCoalTime - gen_;
			 sideRecombPop = derPop;
		 }
		 return sideRecombPop ? maxWaitTime : gens_inf;
	 }

	 void sideRecombExecute() {
		 if ( sideRecombPop == derPop ) {
			 glen_t residue;
			 idx_t nodeIdx = sideRecombRates_der.findCumulativeFraction( random_double(), &residue );
			 Node *n = derPop->pop_get_node( nodeIdx-1 );
			 assert( n->get_idx_in_pop() == nodeIdx-1 );
			 assert( residue > glen_t( 0.0 ) );
			 assert( boost::get_pointer( sideRecombVals_der[ nodeIdx ] ) );
			 if ( !( boost::get( sideRecombVals_der[ nodeIdx ] ) > residue ) ) {
				 std::cerr << "nodeIdx=" << nodeIdx << " residue=" << residue << " sideRecombVals_der[nodeIdx]=" <<
						boost::get( sideRecombVals_der[ nodeIdx ] ) << "\n";
			 }
			 assert( boost::get( sideRecombVals_der[ nodeIdx ] ) > residue );
			 assert( ( ( selPos < seglist_beg( n->getSegs() ) )
								 && ( sideRecombVals_der[ nodeIdx ] == get_gloc( seglist_beg( n->getSegs() ) ) - selPos_gloc ) ) 
							 ||
							 ( ( selPos > seglist_end( n->getSegs() ) ) &&
								 ( sideRecombVals_der[ nodeIdx ] == selPos_gloc - get_gloc( seglist_end( n->getSegs() ) ) ) )
							 ||
							 ( sideRecombVals_der[ nodeIdx ] == glen_t( 0.0 ) ) );
																																																															
			 PRINT4( "sideRecombExecute - moving node ", nodeIdx, sideRecombVals_der[ nodeIdx ], residue );
			 derPop->pop_remove_node( n );
			 ancPop->pop_add_node( n );
		 } else {
			 glen_t residue;
			 idx_t nodeIdx = sideRecombRates_anc.findCumulativeFraction( random_double(), &residue );
			 assert( residue > glen_t( 0 ) );
			 assert( boost::get_pointer( sideRecombVals_anc[ nodeIdx ] ) );
			 
			 Node *n = ancPop->pop_get_node( nodeIdx-1 );
			 assert( n->get_idx_in_pop() == nodeIdx-1 );
			 assert( boost::get( sideRecombVals_anc[ nodeIdx ] ) > residue );

			 PRINT4( "sideRecombExecute - moving node ", nodeIdx, sideRecombVals_anc[ nodeIdx ], residue );
			 assert( ( ( selPos < seglist_beg( n->getSegs() ) )
								 && ( sideRecombVals_anc[ nodeIdx ] == get_gloc( seglist_beg( n->getSegs() ) ) - selPos_gloc ) ) 
							 ||
							 ( ( selPos > seglist_end( n->getSegs()  ) ) &&
								 ( sideRecombVals_anc[ nodeIdx ] == selPos_gloc - get_gloc( seglist_end( n->getSegs() ) ) ) )
							 ||
							 ( sideRecombVals_anc[ nodeIdx ] == glen_t( 0.0 ) ) );

			 
			 ancPop->pop_remove_node( n );
			 derPop->pop_add_node( n );
		 }
		 sideRecombPop = NULL;
	 }

	 Pop *sideRecombPop;
	 gloc_t selPos_gloc;
	 

};  // class Event_SweepOnePop

// *** Misc routines
                
static bool sweepFracSample = False;

void sweep3_set_sweepFracSample( bool sweepFracSample_ ) {
	sweepFracSample = sweepFracSample_;
}

// unsigned nrecomb_sel = 0, nrecomb_unsel = 0, nnew_sel = 0, nnew_unsel = 0, ncoal_sel = 0, ncoal_unsel = 0;



gens_t getSideRecombWaitTime( genid gen_, gens_t maxWaitTime  ) {
	return theSweep ? theSweep->getSideRecombWaitTime( gen_, maxWaitTime ) : static_cast<gens_t>( std::numeric_limits<double>::infinity() );
}
void sideRecombExecute() { if ( theSweep ) theSweep->sideRecombExecute(); } 


const char *Event_SweepOnePop3_typeStr() { return Event_SweepOnePop::typeStr(); }
HistEvents::EventP make_shared_Event_SweepOnePop3( HistEvents *histEvents, istream& is );
HistEvents::EventP make_shared_Event_SweepOnePop3( HistEvents *histEvents, istream& is ) {
	theSweep.reset( new Event_SweepOnePop( histEvents, is ) );
  HistEvents::EventP ep( theSweep );
  return ep;
}

Event_SweepOnePop::~Event_SweepOnePop() {}

static filename_t sweep3_trajFN;

void sweep3_setTrajFN( filename_t fname ) {
	sweep3_trajFN = fname;
}

void Event_SweepOnePop::init_() {

	PRINT2( "setting_selpos", selPos );
	if ( !sweep3_no_oneSidedRecombs ) {
		//getDemography()->getNodePool()->setSelPos( selPos );
		
	}
	selPos_gloc = getDemography()->getNodePool()->getGenMap()->getGdPos( selPos );
	PRINT2( selPos, selPos_gloc );
	
	sel_leaves = make_empty_leafset();
  
	Pop *pop = getDemography()->dg_get_pop_by_name( sweepPop );
	
	chkCond( pop != NULL, "pop is null" );
  
	nchroms_t final_pop_size = pop->pop_get_size();
	totPopSize = static_cast<popsize_float_t>( final_pop_size );
	nchroms_t derPopSize( static_cast<double>(final_pop_size) * final_sel_freq );
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
																					final_sel_freq,
																					gen, /* exactFraction= */ sweepFracSample );
  
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

	{
		popListener = boost::make_shared<SweepOnePopListener>( thisPtr );
		
		BOOST_FOREACH( Node *it, derPop->getMembers() )
			 popListener->nodeAdded( derPop, it->getName(), it->get_idx_in_pop(), seglist_beg( it->getSegs() ), seglist_end( it->getSegs() ) );
		BOOST_FOREACH( Node *it, ancPop->getMembers() )
			 popListener->nodeAdded( ancPop, it->getName(), it->get_idx_in_pop(), seglist_beg( it->getSegs() ), seglist_end( it->getSegs() ) );

		derPop->setPopListener( popListener );
		ancPop->setPopListener( popListener );
	}


	sweepTraj = boost::make_shared<sweepTraj_t>();
	if ( sweep3_trajFN.empty() ) {

		typedef math::Function< genid, popsize_float_t, math::SweepPopSizeTraj > sweepTraj_det_t;
		typedef boost::shared_ptr<sweepTraj_det_t> sweepTraj_det_p;

		double epsilon = 1. / (2. * final_pop_size);
		gens_t end_shift =  (final_sel_freq > 1-epsilon) ? ZERO_GENS :
			 ( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / selCoeff );
		genid tend = gen - ( 2 * log(epsilon) / selCoeff );

		sweepTraj_det_p sweepTraj_det = boost::make_shared<sweepTraj_det_t>(
			popsize_float_t( final_pop_size ),
			selCoeff,
			epsilon,
			end_shift,
			tend
			);
	
	

		genid g( 0.0 );
		gens_t g_step( 1.0 );
		while ( true ) {
			popsize_float_t psize = eval( *sweepTraj_det, g );
			if ( psize < static_cast< popsize_float_t >( .5 ) ) break;
			sweepTraj->getPieces().insert( std::make_pair( g, psize ) );
			//PRINT2( g, psize );
			g += g_step;
		}
		
	} else {
		// traj specified
	  boost::filesystem::ifstream trajFile( sweep3_trajFN );
		typedef math::Function< genid, freq_t, math::Piecewise< math::Const<> > > freq_traj_type;
		freq_traj_type genid2derFreq;
		loadFrom( trajFile, genid2derFreq );

		using boost::for_each;
		using boost::bind;
		using boost::adaptors::map_values;
		using boost::adaptors::transformed;
		using math::Function;
		using math::Const;
		using boost::algorithm::clamp;
			 

		freq_t freqDelta = static_cast<popsize_float_t>(1.0) / ( 2 * totPopSize );
		// for_each( genid2derFreq.getPieces() | map_values |
		// 					transformed( bind( &Function<genid, freq_t, Const<> >::getValRef, _1 ) ), 
		// 					_1 = bind( clamp<popsize_float_t>, _1, freqDelta, 1.0 - freqDelta ) );

		
		BOOST_FOREACH( freq_t& derFreq, genid2derFreq.getPieces() | map_values |
									 transformed( bind( &Function<genid, freq_t, Const<> >::getValRef, _1 ) ) )
			 derFreq = clamp( derFreq, freqDelta, 1.0 - freqDelta );
		
		*sweepTraj = math::Function< genid, popsize_float_t, math::Const<> >( totPopSize ) * genid2derFreq;
	}
	PRINT( *sweepTraj );

	this->sweepStartTime = sweepTraj->getPieces().begin()->first - gens_t( 1 );
	PRINT( this->sweepStartTime );
	
	//std::cerr << "sweepStartTime=" << this->sweepStartTime << "\n";

	PRINT( "creating derPop arrival process" );
	derPop->setCoalArrivalProcess(
		math::ArrivalProcess< genid, math::Any< RandGen > >(
			math::makeNonHomogeneousPoissonProcess
			(
				math::coalRateFunction(
					*sweepTraj
					
					), gen, "coal in derPop"
				)));
	
	BOOST_AUTO( sweepPopSizeComplementTraj,
							(( math::Function< genid, popsize_float_t,
								 math::Const<> >( totPopSize ) - *sweepTraj )) );

	PRINT( "creating ancPop arrival process" );
	ancPop->setCoalArrivalProcess(
		math::ArrivalProcess< genid, math::Any< RandGen > >(
			math::makeNonHomogeneousPoissonProcess
			(
				math::coalRateFunction(
					sweepPopSizeComplementTraj
					)
				, gen, "coal in ancPop"
				)));

	PRINT( derPop->getCoalArrivalProcess() );
	PRINT( ancPop->getCoalArrivalProcess() );

	BOOST_AUTO( totPopSizeInv, (math::Function< genid, popsizeInv_float_t, math::Const<> >( factor_t(1.0) / totPopSize ) ) );

	BOOST_AUTO( freqTraj, totPopSizeInv * ( *sweepTraj ) );

	math::Function< genid, factor_t, math::Const<> > totRecombRate( ToDouble( getDemography()->getNodePool()->getGenMap()->getRegionRecombRateAbs() ) );

	PRINT( totRecombRate );
	PRINT( freqTraj );
	
	PRINT( "creating sideRecombProcess_der" );
	sideRecombProcess_der = side_recomb_process_type(
		math::makeNonHomogeneousPoissonProcess(
			totRecombRate * freqTraj, gen, "sideRecomb anc->der" ) );
	PRINT( "creating sideRecombProcess_anc" );

	sideRecombProcess_anc = side_recomb_process_type(
		math::makeNonHomogeneousPoissonProcess(
			totRecombRate * ( math::Function< genid, factor_t, math::Const<> >( 1.0 ) - freqTraj ), gen,
			"sideRecomb der->anc" ) );

	PRINT( sideRecombProcess_der );
	PRINT( sideRecombProcess_anc );
	
	
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
			freq_t testFreq = eval( *sweepTraj, curGen )  / totPopSize;
			prob_t pval = random_double();
			Pop * shouldBeInPop = ( pval < testFreq ) ? derPop : ancPop;
			if ( shouldBeInPop != curPop ) {
				curPop->pop_remove_node( node );
				shouldBeInPop->pop_add_node( node );
			}
		}  // if ( curPop == ancPop || curPop == derPop ) 
	}  // if ( node ) 
} // void Event_SweepOnePop::SweepOnePopHook::determineAlleleAtSelPos_( Node *node )

// * Postamble
}  // namespace sweep3
}  // namespace cosi
