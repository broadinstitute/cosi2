/* $Id: historical.c,v 1.4 2011/05/26 21:36:22 sfs Exp $ */

/*
 * Module: historical.c
 *
 * 1. Processes historical events from datafile.
 * 2. Returns the time to the next historical event when called
 *    by simulator.c.
 * 3. Executes the historical event, by making calls to demography.c.
 *
 */

#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <ios>
#include <sstream>
#include <algorithm>
#include <boost/make_shared.hpp>
#include <boost/next_prior.hpp>
#include <cosi/utils.h>
#include <cosi/pop.h>
#include <cosi/historical.h>
#include <cosi/demography.h>
#include <cosi/migrate.h>
#include <cosi/sweep.h>
#include <cosi/sweep1.h>
#include <cosi/sweep2.h>
#include <cosi/sweep3.h>
#include <cosi/generalmath.h>

namespace cosi {

HistEvents::HistEvents( DemographyP demography_ ):
	demography( demography_ ) { }

HistEvents::Event::~Event() {}


double HistEvents::Event::to_4N0( genid gen ) const { return ToDouble( gen / ToDouble( 4 * getDemography()->get_N0_tot() ) ); }

void HistEvents::write_ms_flags( ostream& s ) const {
	BOOST_FOREACH( events_type::value_type e, events ) {
		e.second->write_ms_flags( s );
	}
}

void HistEvents::processSimEnd( genid simEndTime ) {
	BOOST_FOREACH( events_type::value_type e, events )
		 e.second->processSimEnd( simEndTime );
}

namespace histevents {

// Class: Event_PopSize
// Historical event that sets the size of a given pop at a given generation.  The size remains at this value going
// pastward, unless/until changed by another historical event.  Note that the present-day size of each population
// is set in the <parameter file>.
// Note also that this event affects only population size, not the sample size (number of <Nodes> belonging to that pop).
class Event_PopSize: public HistEvents::Event {
public:
	 Event_PopSize( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, nchroms_t popsize_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ), popsize( popsize_ ) { }
	 Event_PopSize( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> pop >> gen >> popsize;
	 }
	 virtual ~Event_PopSize();

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "change_size"; }
			
	 virtual genid execute();

private:
	 // Field: pop
	 // The population in which we're setting the new population size.
	 popid pop;
			
	 // Field: popsize
	 // The new population size
	 nchroms_t popsize;
};  // class Event_PopSize

// Class: Event_PopSizeExp
// Exponential expansion of population size over a range of generations.
class Event_PopSizeExp: public HistEvents::Event {
public:
	 Event_PopSizeExp( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, nchroms_t popsize_,
										 genid genBeg_, nchroms_t popSizeBeg_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ) , popsize( popsize_ ), genBeg( genBeg_ ), popSizeBeg( popSizeBeg_ ) {
		 calcExpansionRate();
	 }
	 Event_PopSizeExp( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> pop >> gen >> genBeg >> popsize >> popSizeBeg;
		 calcExpansionRate();
		 PRINT6( pop, gen, genBeg, popsize, popSizeBeg, expansionRate );
	 }
	 virtual ~Event_PopSizeExp();
			
	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "exp_change_size"; }

	 virtual genid execute();

	 // Const: STEP
	 // By how many generations to step each time while we implement the exponential expansion?
	 static const gens_t STEP;
			
private:
	 // Field: pop
	 // The population in which we're setting the new population size.
	 popid pop;
			
	 // Field: popsize
	 // Population size at <gen>.
	 popsize_float_t /* nchroms_t */  popsize;
			
	 // Field: genBeg
	 // Time when the increase started (going forward).  Since time numerically increases pastward,
	 // genBeg > gen.
	 genid genBeg;
			
	 // Field: popSizeAtBeg
	 // Population size at <genBeg>.
	 popsize_float_t /* nchroms_t */ popSizeBeg;

	 // Field: expansionRate
	 // The expansion rate at each step
	 gensInv_t expansionRate;

	 void calcExpansionRate() { expansionRate = log( popSizeBeg / popsize ) / ( genBeg - gen ); }

};  // class Event_PopSizeExp 


// Class: Event_PopSizeExp2
// Exponential expansion of population size over a range of generations.
class Event_PopSizeExp2: public HistEvents::Event {

	 static int procCount;
	 
public:
	 Event_PopSizeExp2( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, nchroms_t popsize_,
										 genid genBeg_, nchroms_t popSizeBeg_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ) , popsize( popsize_ ), genBeg( genBeg_ ), popSizeBeg( popSizeBeg_ ) , procName( -1 ) {
		 calcExpansionRate();
	 }
	 Event_PopSizeExp2( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ), procName( -1 ) {
		 is >> pop >> gen >> genBeg >> popsize >> popSizeBeg;
		 calcExpansionRate();
		 PRINT6( pop, gen, genBeg, popsize, popSizeBeg, expansionRate );
	 }
	 virtual ~Event_PopSizeExp2();
			
	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "exp_change_size2"; }

	 virtual genid execute();

private:
	 // Field: pop
	 // The population in which we're setting the new population size.
	 popid pop;
			
	 // Field: popsize
	 // Population size at <gen>.
	 popsize_float_t /* nchroms_t */  popsize;
			
	 // Field: genBeg
	 // Time when the increase started (going forward).  Since time numerically increases pastward,
	 // genBeg > gen.
	 genid genBeg;
			
	 // Field: popSizeAtBeg
	 // Population size at <genBeg>.
	 popsize_float_t /* nchroms_t */ popSizeBeg;

	 // Field: expansionRate
	 // The expansion rate at each step
	 gensInv_t expansionRate;

	 int procName;
	 std::string procNameStr;

	 void calcExpansionRate() { expansionRate = log( popSizeBeg / popsize ) / ( genBeg - gen ); }

};  // class Event_PopSizeExp2 


// Class: Event_Split
// Split a pop going forward (join two pops going backward).
class Event_Split: public HistEvents::Event {
public:
	 Event_Split( HistEvents *histEvents_, const string& label_, genid gen_, popid fromPop_, popid newPop_ ):
		 Event( histEvents_, label_, gen_ ), fromPop( fromPop_ ), newPop( newPop_ ) { }
	 Event_Split( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) { is >> fromPop >> newPop >> gen; }
	 virtual ~Event_Split();

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "split"; }

	 virtual genid execute();

	 virtual void write_ms_flags( std::ostream& s ) const {
		 s << " -ej " << to_4N0( gen ) << " " << ( getDemography()->dg_get_pop_index_by_name( newPop ) + 1 )
			 << " " << ( getDemography()->dg_get_pop_index_by_name( fromPop ) + 1 );
   }
	 
			
private:
	 // Field: fromPop
	 // Population that is split, going forward.
	 popid fromPop;

	 // Field: newPop
	 // New pop that is created from part of <fromPop>, going forward.
	 popid newPop;

};  // class Event_Split

// Class: Event_MigrationRate
// Set the migration rate between a pair of pops.
class Event_MigrationRate: public HistEvents::Event {
public:
	 Event_MigrationRate( HistEvents *histEvents_, const string& label_, genid gen_, popid fromPop_, popid toPop_, prob_t rate_ ):
		 Event( histEvents_, label_, gen_ ), fromPop( fromPop_ ), toPop( toPop_ ), rate( rate_ ) { }
	 Event_MigrationRate( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 // is >> fromPop >> toPop >> gen >> rate;
		 is >> toPop >> fromPop >> gen >> rate;  // deliberate reversal if to and from BUG FIXME
	 }
	 virtual ~Event_MigrationRate();

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "migration_rate"; }

	 virtual genid execute();
			
private:
	 // Field: fromPop
	 // Migration from which pop?
	 popid fromPop;

	 // Field: toPop
	 // Migration to which pop?
	 popid toPop;

	 // Field: rate
	 // Probability of migration, per chromosome per generation
	 prob_per_chrom_per_gen_t rate;
			
};  // class Event_MigrationRate

// Class: Event_Bottleneck
// An (instantaneous) population bottleneck.
class Event_Bottleneck: public HistEvents::Event {
public:
	 Event_Bottleneck( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, double inbreedingCoefficient_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ), inbreedingCoefficient( inbreedingCoefficient_ ) { }
	 Event_Bottleneck( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) { is >> pop >> gen >> inbreedingCoefficient; }
	 virtual ~Event_Bottleneck();

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "bottleneck"; }

	 virtual genid execute();
			
private:
	 // Field: pop
	 // Population in which the bottleneck happens.
	 popid pop;

	 // Field: inbreedingCoefficient
	 // Strength of bottleneck
	 double inbreedingCoefficient;

};  // class Event_Bottleneck

// Class: Event_Admix
// Admixture of two populations.
class Event_Admix: public HistEvents::Event {
public:
	 Event_Admix( HistEvents *histEvents_, const string& label_, genid gen_, popid admixedPop_, popid sourcePop_, frac_t admixFrac_ ):
		 Event( histEvents_, label_, gen_ ), admixedPop( admixedPop_ ), sourcePop( sourcePop_ ), admixFrac( admixFrac_ ) { }
	 Event_Admix( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) { is >> admixedPop >> sourcePop >> gen >> admixFrac; }

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "admix"; }

	 virtual ~Event_Admix();

	 virtual genid execute();
			
private:
	 // When does the 
	 // Field: admixedPop
	 // The pop that receives chroms from the source pop
	 popid admixedPop;
			
	 // Field: sourcePop
	 // The pop that sends chroms to the admixed pop
	 popid sourcePop;

	 // Field: admixFrac
	 // Fraction of admixed chroms from source
	 frac_t admixFrac;
}; // class Event_Admix

	 // Class: Event_Sweep
	 // A selective sweep in a single population.
class Event_Sweep: public HistEvents::Event {
public:
	 Event_Sweep( HistEvents *histEvents_, const string& label_, genid gen_, popid sweepPop_, gensInv_t selCoeff_, loc_t selPos_,
								freq_t final_sel_freq_ ):
		 Event( histEvents_, label_, gen_ ), sweepPop( sweepPop_ ), selCoeff( selCoeff_ ), selPos( selPos_ ), final_sel_freq( final_sel_freq_ ) { }
	 Event_Sweep( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
	 }

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "sweep_orig"; }
				
	 virtual ~Event_Sweep();

	 virtual genid execute();
			
private:
	 // Field: sweepPop
	 // The population in which the sweep happens.
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
};  // class Event_Sweep



Event_PopSize::~Event_PopSize() {}
Event_PopSizeExp::~Event_PopSizeExp() {}
Event_PopSizeExp2::~Event_PopSizeExp2() {}
Event_Bottleneck::~Event_Bottleneck() {}
Event_MigrationRate::~Event_MigrationRate() {}
Event_Admix::~Event_Admix() {}
Event_Split::~Event_Split() {}
Event_Sweep::~Event_Sweep() {}

genid Event_PopSize::execute() {
	getDemography()->dg_set_pop_size_by_name( gen, pop, popsize );
	return gen;
}

// Const: STEP
// By how many generations to step each time while we implement the exponential expansion?
const gens_t Event_PopSizeExp::STEP( getenv( "COSI_POPSIZEEXP_STEP" ) ? atof( getenv( "COSI_POPSIZEEXP_STEP" ) ) : 10.0 );

genid Event_PopSizeExp::execute() {
	PRINT( STEP );
	getDemography()->dg_set_pop_size_by_name( gen, pop, nchroms_t( ToDouble( popsize + static_cast< popsize_float_t >( .5 ) ) ) );
	genid oldgen = gen;
	if ( gen < genBeg ) {
		genid newgen = std::min( oldgen + STEP, genBeg );
		//		double popsize_bef = popsize;
		popsize *= exp( expansionRate * ( newgen - oldgen ) );
		gen = newgen;
		addEvent( shared_from_this() );
	}
	return oldgen;
}

int Event_PopSizeExp2::procCount = 0;

genid Event_PopSizeExp2::execute() {

	genid oldgen = gen;
	Pop *popPtr = getDemography()->dg_get_pop_by_name( pop );
	if ( procName == -1 ) {
		procName = ++procCount;
		std::ostringstream procNameStrm;
		procNameStrm << "exp" << ++procCount;
		procNameStr = procNameStrm.str();

		using namespace math;

		BOOST_AUTO( exponent,
								( Function< genid, gensInv_t, Const<> >( expansionRate ) *
									(
										Function< genid, gens_t, Const<> >( genBeg - static_cast<genid>( 0.0 ) )
										 -
										 Function< genid, gens_t, X_To<1> >()
										) ) );
		BOOST_AUTO( with_exp, exp_( exponent ) );
		BOOST_AUTO( f, ( Function<genid,popsizeInv_float_t,Const<> >( 1. / popSizeBeg ) * with_exp ) );

		PRINT4( popsize, popSizeBeg, 1.0 / eval( f, gen ), 1.0 / eval( f, genBeg ) );

		BOOST_AUTO( coalRateFn, ( Function< genid, double, Const<> >( .5 ) * ( f ) ) );

#ifndef NDEBUG		
		{
			std::cerr.precision( 15 );
			const genid g_beg( 10.0 ), g_end( 150.0 );
			const gens_t g_step( 1.0 );
			genid g( g_beg );
			while( g < g_end ) {
				PRINT2( g, eval( coalRateFn, g ) );
				g += g_step;
			}

			for( int n = 1; n < 25; ++n ) {
				PRINT2( n, ( integrateNumerically( coalRateFn, g_beg, g_end, n ) ) );
			}
			PRINT( definiteIntegral( coalRateFn, g_beg, g_end ) );
		}
#endif		
		
		popPtr->
			 setCoalArrivalProcess(
				 ArrivalProcess< genid, Any< RandGen > >(
					 makeNonHomogeneousPoissonProcess
					 ( coalRateFn, gen, procNameStr ) ) );

		PRINT( popPtr->getCoalArrivalProcess() );

		gen = genBeg;
		addEvent( shared_from_this() );
	} else {
		if ( popPtr->getCoalArrivalProcess() && popPtr->getCoalArrivalProcess().getLabel() == procNameStr ) {
			 popPtr->clearCoalArrivalProcess();
			 
			 PRINT( "clearedArrivalProcess" );
			 getDemography()->dg_set_pop_size_by_name( gen, pop, nchroms_t( ToDouble( popSizeBeg + static_cast< popsize_float_t >( .5 ) ) ) );			 
		}
	}
	return oldgen;
}


genid Event_Bottleneck::execute() {
  Pop* the_pop = getDemography()->dg_get_pop_by_name(pop);
  int num_nodes = the_pop->pop_get_num_nodes();
  gens_t t = ZERO_GENS;
	rate_t rate;
  
  if (num_nodes < 2) return gen;
  double effective_N( - 1.0 / ( 2.0 * log (1.0 - inbreedingCoefficient)) );
  rate = (4 * effective_N) / (num_nodes * (num_nodes - 1));
  t += gens_t( poisson_get_next (1 / rate) ) ;
	const gens_t GEN_ONE(1.0);
  while (t <= GEN_ONE) {
    getDemography()->dg_coalesce_by_pop (the_pop, gen + t, /* forceCoalescence= */ True);
		num_nodes --;
		if (num_nodes > 1) {
			rate = (double) (4 * effective_N)
				 / (num_nodes * (num_nodes - 1));
			t += gens_t( poisson_get_next (1/rate) );
		}
		else t += GEN_ONE;      /* escape loop */
	}

	return gen;  // FIXME - should return t?
}  // Event_Bottleneck::execute()

genid Event_Sweep::execute() {
	return getSweep()->sweep_execute( sweepPop, selCoeff, gen, selPos, final_sel_freq );
}

genid Event_MigrationRate::execute() {
	getMigrate()->migrate_set_rate (fromPop, toPop, rate);
	return gen;
}

genid Event_Admix::execute() {

	/* 27-sep-04.  Modified to simply shift the required # chroms to another pop. */

	/* join two pops forward = split pops backwards */
	/* change pop size to what is specified. */
	/*		  		  demography->dg_set_pop_size_by_name (currentevent->gen,
								currentevent->popindex[1],
								currentevent->params[0]);  */
		  
	getDemography()->dg_move_nodes_by_name (admixedPop, sourcePop, admixFrac, gen );
	return gen;
}

genid Event_Split::execute() {
	/* split two pops forward = join pops backwards */
	getDemography()->dg_move_nodes_by_name ( /* going backwards, move all nodes from */ newPop,
																					 /* into */ fromPop,
																					 /* fractionToMove= */ 1.0,
																					 gen, /* exactFraction= */ true );
	getMigrate()->migrate_delete_all_for_pop( newPop );
	/*		  demography->dg_end_pop_by_name (currentevent->popindex[1]); */
	return gen;
}

}  // namespace histevents

///////////////////////////////////////////////////////////////////////////////////////////////////

// MethodP: historical_event_execute
// Execute the next historical event, applying its effects to the demographic model.
//
// Input params:
//
//   historical_event_time - time of the next historical event
//
// Returns:
//
//   the time at the end of this historical event.  Most events are instantaneous,
//   but a few are not.
genid HistEvents::historical_event_execute
(genid COSI_IF_DEBUG(historical_event_time) ) 
{
	EventP curEvent = getCurEvent();
	COSI_IF_DEBUG( genid curEvent_gen = curEvent->getGen() );
	assert( equal_eps( curEvent_gen, historical_event_time ) );
	genid gen_after_event = curEvent->execute();
	assert( gen_after_event >= curEvent_gen );
	events.erase( events.begin() );
	assert( gen_after_event >= curEvent_gen );
	assert( events.empty() || gen_after_event <= getCurEvent()->getGen() );
	return gen_after_event;
}

// Method: historical_get_time_till_next_event
// Returns the amount of time until the next historical event, from the current time (passed in).
gens_t HistEvents::historical_get_time_till_next_event (genid gen) const {
	return events.empty() ? NULL_GENS : getCurEvent()->getGen() - gen;
}

namespace sweep1 {
 HistEvents::EventP make_shared_Event_SweepOnePop( HistEvents *histEvents, istream& is );
}
namespace sweep2 {
 HistEvents::EventP make_shared_Event_SweepOnePop2( HistEvents *histEvents, istream& is );
}
namespace sweep3 {
 HistEvents::EventP make_shared_Event_SweepOnePop3( HistEvents *histEvents, istream& is );
}

HistEvents::EventP HistEvents::parseEvent( const char *buffer ) {
	using std::istringstream;
	using std::string;
	using boost::make_shared;
	using std::ios;

	const char *Event_SweepNew_typeStr();
	EventP make_shared_Event_SweepNew( HistEvents *histEvents, istream& is );

	istringstream is( buffer );
	is.exceptions( ios::failbit | ios::badbit | ios::eofbit );
	EventP event;
	try {
		using namespace histevents;
		string typestr;
		is >> typestr;
		if ( typestr == Event_PopSize::typeStr() ) event.reset( new Event_PopSize( this, is ) );
		else if ( typestr == Event_PopSizeExp::typeStr() ) event.reset( new Event_PopSizeExp( this, is ) );
		else if ( typestr == Event_PopSizeExp2::typeStr() ) event.reset( new Event_PopSizeExp2( this, is ) );
		else if ( typestr == Event_Split::typeStr() ) event.reset( new Event_Split( this, is ) );
		else if ( typestr == Event_MigrationRate::typeStr() ) event.reset( new Event_MigrationRate( this, is ) );
		else if ( typestr == Event_Bottleneck::typeStr() ) event.reset( new Event_Bottleneck( this, is ) );
		else if ( typestr == Event_Admix::typeStr() ) event.reset( new Event_Admix( this, is ) );
		else if ( typestr == Event_Sweep::typeStr() ) event.reset( new Event_Sweep( this, is ) );
		else if ( typestr == Event_SweepNew_typeStr() ) event = make_shared_Event_SweepNew( this, is );
		else if ( typestr == sweep1::Event_SweepOnePop_typeStr() ) event = sweep1::make_shared_Event_SweepOnePop( this, is );
		else if ( typestr == sweep2::Event_SweepOnePop2_typeStr() ) event = sweep2::make_shared_Event_SweepOnePop2( this, is );
		else if ( typestr == sweep3::Event_SweepOnePop3_typeStr() ) event = sweep3::make_shared_Event_SweepOnePop3( this, is );
		else chkCond( False, "could not parse event %s", buffer );
	} catch( ios::failure e ) {
		chkCond( False, "could not parse event %s", buffer );
	}
	return event;
}

void HistEvents::addEvent( EventP event ) {
	events.insert( make_pair( event->getGen(), event ) );
}

}  // namespace cosi
