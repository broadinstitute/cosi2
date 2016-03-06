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

//#define COSI_DEV_PRINT

#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <ios>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/make_shared.hpp>
#include <boost/next_prior.hpp>
#include <boost/exception/all.hpp>
#include <cosi/general/utils.h>
#include <cosi/pop.h>
#include <cosi/historical.h>
#include <cosi/demography.h>
#include <cosi/migrate.h>
#include <cosi/sweep.h>
#include <cosi/sweep1.h>
#include <cosi/sweep2.h>
#include <cosi/sweep3.h>
//#include <cosi/generalmath.h>
#include <cosi/basemodel.h>
#include <cosi/msweep.h>

namespace cosi {

struct cosi_hist_event_error: virtual cosi_io_error {};
typedef boost::error_info<struct errinfo_hist_err_detail_,std::string> errinfo_hist_err_detail;

HistEvents::HistEvents( DemographyP demography_ ):
	demography( demography_ ) { }

HistEvents::Event::~Event() {}

namespace histevents {

// Class: Event_PopSize
// Historical event that sets the size of a given pop at a given generation.  The size remains at this value going
// pastward, unless/until changed by another historical event.  Note that the present-day size of each population
// is set in the <parameter file>.
// Note also that this event affects only population size, not the sample size (number of <Nodes> belonging to that pop).
class Event_PopSize: public HistEvents::Event {
public:
	 Event_PopSize( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, nchroms_t popsize_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ), popsize( popsize_ ) {
		 chkRep();
	 }
	 Event_PopSize( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> pop >> gen >> popsize;
		 chkRep();
	 }
	 virtual ~Event_PopSize();

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "change_size"; }
			
	 virtual genid execute();
	 virtual void addToBaseModel( BaseModel& ) const;

	 virtual eventKind_t getEventKind() const { return E_POPSIZE; }

private:
	 // Field: pop
	 // The population in which we're setting the new population size.
	 popid pop;
			
	 // Field: popsize
	 // The new population size
	 nchroms_t popsize;

	 void chkRep() const {
		 if ( popsize < nchroms_t(0) )
				BOOST_THROW_EXCEPTION( cosi_hist_event_error()
															 << error_msg( "pop size negative" ) );
	 }
};  // class Event_PopSize

// Class: Event_PopSizeExp
// Exponential expansion of population size over a range of generations.
class Event_PopSizeExp: public HistEvents::Event {
public:
	 Event_PopSizeExp( HistEvents *histEvents_, const string& label_, genid gen_, popid pop_, nchroms_t popsize_,
										 genid genBeg_, nchroms_t popSizeBeg_ ):
		 Event( histEvents_, label_, gen_ ), pop( pop_ ) , popsize( popsize_ ), genBeg( genBeg_ ), popSizeBeg( popSizeBeg_ ) {
		 PRINT6( pop, gen, genBeg, popsize, popSizeBeg, expansionRate );
		 calcExpansionRate();
	 }
	 Event_PopSizeExp( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> pop >> gen >> genBeg >> popsize >> popSizeBeg;
		 PRINT6( pop, gen, genBeg, popsize, popSizeBeg, expansionRate );
		 calcExpansionRate();
		 PRINT6( pop, gen, genBeg, popsize, popSizeBeg, expansionRate );
	 }
	 virtual ~Event_PopSizeExp();
			
	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "exp_change_size_orig"; }

	 virtual genid execute();
	 virtual void addToBaseModel( BaseModel& ) const;

	 virtual eventKind_t getEventKind() const { return E_POPSIZEEXP; }	 

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

	 void calcExpansionRate() {
		 if ( genBeg <= gen )
				BOOST_THROW_EXCEPTION( cosi_hist_event_error()
															 << error_msg( "pop size expansion ends before it starts" ) );
		 // if ( popSizeBeg >= popsize ) 
		 // 		BOOST_THROW_EXCEPTION( cosi_hist_event_error()
		 // 													 << error_msg( "pop size expansion: size does not increase" ) );
		 expansionRate = log( popSizeBeg / popsize ) / ( genBeg - gen );
	 }

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
	 static const char *typeStr2() { return "exp_change_size"; }

	 virtual genid execute();
	 virtual void addToBaseModel( BaseModel& ) const;
	 virtual eventKind_t getEventKind() const { return E_POPSIZEEXP; }	 

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

	 void calcExpansionRate() {
		 if ( genBeg <= gen )
				BOOST_THROW_EXCEPTION( cosi_hist_event_error()
															 << error_msg( "pop size expansion ends before it starts" ) );
		 // if ( popSizeBeg >= popsize ) 
		 // 		BOOST_THROW_EXCEPTION( cosi_hist_event_error()
		 // 													 << error_msg( "pop size expansion: size does not increase" ) );
		 
		 expansionRate = log( popSizeBeg / popsize ) / ( genBeg - gen );
	 }

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
	 virtual eventKind_t getEventKind() const { return E_SPLIT; }	 
	 virtual void addToBaseModel( BaseModel& ) const;
			
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
	 virtual eventKind_t getEventKind() const { return E_MIGRATIONRATE; }	 

	 virtual genid execute();
	 virtual void addToBaseModel( BaseModel& baseModel ) const;
			
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
	 virtual eventKind_t getEventKind() const { return E_BOTTLENECK; }	 

	 virtual genid execute();
	 virtual void addToBaseModel( BaseModel& ) const;
			
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
	 virtual eventKind_t getEventKind() const { return E_ADMIX; }	 

	 virtual ~Event_Admix();

	 virtual genid execute();
	 void addToBaseModel( BaseModel& baseModel ) const;
			
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
	 virtual eventKind_t getEventKind() const { return E_SWEEP; }	 
				
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


// Class: Event_MSweep
// A selective sweep .
class Event_MSweep: public HistEvents::Event {
public:
	 Event_MSweep( HistEvents *histEvents_, istream& is ): Event( histEvents_, is ) {
		 is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
		 if ( !( final_sel_freq.getMin() && final_sel_freq.getMax() &&
						 0. <= final_sel_freq.getMin() && final_sel_freq.getMin() <= 1. ) )
				BOOST_THROW_EXCEPTION( cosi_hist_event_error() << error_msg( "invalid final freq range" ) );
	 }

	 // Method: typeStr
	 // Return the string denoting this type of historical event in the <parameter file>.
	 // Historical events are specified in the parameter file by lines of the form
	 // pop_event <eventType> <eventParams>
	 // This method specifies the eventType.
	 static const char *typeStr() { return "sweep_mult"; }
	 virtual eventKind_t getEventKind() const { return E_SWEEP; }

	 virtual void addToBaseModel( BaseModel& m ) const {
		 setSweepInfo( m, gen, selCoeff, selPos, sweepPop, final_sel_freq  );
	 }
				
	 virtual ~Event_MSweep();

	 virtual genid execute();
			
private:
	 // Field: sweepPop
	 // The population in which the sweep happens.
	 popid sweepPop;

	 // Field: selCoeff
	 // The selection coefficient.
	 double selCoeff;

	 // Field: selPos
	 // The location of the causal mutation.
	 loc_t selPos;

	 // Field: final_sel_freq
	 // The frequency of the derived allele at the end of the sweep.
	 util::ValRange<freq_t> final_sel_freq;
};  // class Event_MSweep


Event_PopSize::~Event_PopSize() {}
Event_PopSizeExp::~Event_PopSizeExp() {}
Event_PopSizeExp2::~Event_PopSizeExp2() {}
Event_Bottleneck::~Event_Bottleneck() {}
Event_MigrationRate::~Event_MigrationRate() {}
Event_Admix::~Event_Admix() {}
Event_Split::~Event_Split() {}
Event_Sweep::~Event_Sweep() {}
Event_MSweep::~Event_MSweep() {}

genid Event_PopSize::execute() {
	getDemography()->dg_set_pop_size_by_name( gen, pop, popsize );
	return gen;
}

void Event_PopSize::addToBaseModel( BaseModel& baseModel ) const {
	baseModel.popInfos[ pop ].setSizeFrom( gen, popsize_float_t( popsize ) );
}


// Const: STEP
// By how many generations to step each time while we implement the exponential expansion?
const gens_t Event_PopSizeExp::STEP( getenv( "COSI_POPSIZEEXP_STEP" ) ? atof( getenv( "COSI_POPSIZEEXP_STEP" ) ) : 10.0 );

genid Event_PopSizeExp::execute() {
	//	PRINT( STEP );
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

void Event_PopSizeExp::addToBaseModel( BaseModel& baseModel ) const {
	using namespace math;

	BOOST_AUTO( &popInfo, baseModel.popInfos[ pop ] );
	BOOST_AUTO( &pieces, popInfo.popSizeFn.getPieces() );

	if ( cosi_fabs( popsize - popSizeBeg ) < popsize_float_t( 1. ) )
		 popInfo.setSizeFrom( gen, fn_const<genid>( popsize ) );
	else
		 popInfo.setSizeFrom( gen, 
													cval( popSizeBeg ) *
													exp_(
														cval( log( popsize / popSizeBeg ) / ( genBeg - gen ) ) *
														( ( cval( genBeg ) -
																fn_x<genid,genid>() ) ) ) );
	
	if ( pieces.find( genBeg ) == pieces.end() )
		 popInfo.setSizeFrom( genBeg, math::fn_const<genid>( popSizeBeg ) );
}

int Event_PopSizeExp2::procCount = 0;

genid Event_PopSizeExp2::execute() {

	genid oldgen = gen;
	Pop *popPtr = getDemography()->dg_get_pop_by_name( pop );
	util::chk( popPtr, "exp2: unknown pop!" );
	//std::cerr << "exp2 exec: gen=" << gen << " procName=" << procName << " pop=" << pop << " nodes=" << popPtr->pop_get_num_nodes() << "\n";
	if ( procName == -1 ) {
		procName = ++procCount;

		std::ostringstream procNameStrm;
		procNameStrm << "exp" << ++procCount;
		procNameStr = procNameStrm.str();

		getDemography()->dg_set_pop_size_by_name( gen, pop, nchroms_t( ToDouble( popsize + static_cast< popsize_float_t >( .5 ) ) ) );
		if ( cosi_fabs( popsize - popSizeBeg ) > popsize_float_t( .9 ) ) {
			using namespace math;
			//std::cerr << "SETTING EXP PROCESS\n";

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

#if 0
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
		}

		gen = genBeg;
		addEvent( shared_from_this() );
	} else {
		if ( popPtr->getCoalArrivalProcess() && popPtr->getCoalArrivalProcess().getLabel() == procNameStr ) {
			 popPtr->clearCoalArrivalProcess();
			 //std::cerr << "gen " << gen << ": CLEARED ARRIVAL PROCESS\n";
			 
			 PRINT( "clearedArrivalProcess" );
		}
		getDemography()->dg_set_pop_size_by_name( gen, pop, nchroms_t( ToDouble( popSizeBeg + static_cast< popsize_float_t >( .5 ) ) ) );
	}
	return oldgen;
}

void Event_PopSizeExp2::addToBaseModel( BaseModel& baseModel ) const {
	using namespace math;

	BOOST_AUTO( &popInfo, baseModel.popInfos[ pop ] );
	BOOST_AUTO( &pieces, popInfo.popSizeFn.getPieces() );

	if ( cosi_fabs( popsize - popSizeBeg ) < popsize_float_t( 1. ) )
		 popInfo.setSizeFrom( gen, fn_const<genid>( popsize ) );
	else
		 popInfo.setSizeFrom( gen, 
													cval( popSizeBeg ) *
													exp_(
														cval( log( popsize / popSizeBeg ) / ( genBeg - gen ) ) *
														( ( cval( genBeg ) -
																fn_x<genid,genid>() ) ) ) );
	
	if ( pieces.find( genBeg ) == pieces.end() )
		 popInfo.setSizeFrom( genBeg, math::fn_const<genid>( popSizeBeg ) );
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

void Event_Bottleneck::addToBaseModel( BaseModel& baseModel ) const {

	BOOST_AUTO( &popInfo, baseModel.popInfos[ pop ] );
	BOOST_AUTO( &pieces, popInfo.popSizeFn.getPieces() );
	BOOST_AUTO( &pieces_coalRate, popInfo.coalRateFn.getPieces() );

  popsize_float_t effective_N( - 1.0 / ( 2.0 * log (1.0 - inbreedingCoefficient)) );
	pieces[ gen + gens_t( 1 ) ] = pieces.lower_bound( gen )->second;
	pieces_coalRate[ gen + gens_t( 1 ) ] = pieces_coalRate.lower_bound( gen )->second;
	popInfo.setSizeFrom( gen, math::fn_const< genid >( effective_N ) );

	// baseModel.popInfos[ pop ].setSizeFrom( gen + gens_t(1),
	// 																			 baseModel.popInfos[ pop ].popSizeFn
	//baseModel.popInfos[ pop ].setSizeFrom( gen, effective_N );
}

genid Event_Sweep::execute() {
	return getSweep()->sweep_execute( sweepPop, selCoeff, gen, selPos, final_sel_freq );
}

genid Event_MSweep::execute() {
	return gen;
}

genid Event_MigrationRate::execute() {
	Pop *fromPopPtr = getDemography()->dg_get_pop_by_name( fromPop );
	Pop *toPopPtr = getDemography()->dg_get_pop_by_name( toPop );
	util::chk( fromPopPtr, "migrate::execute - unknown pop" );
	util::chk( toPopPtr, "migrate::execute - unknown pop" );
	//std::cerr << "migration rate setting: gen=" << gen << " fromPop=" << fromPop << " toPop=" << toPop << "\n";
	if ( !fromPopPtr->isInactive() && !toPopPtr->isInactive() )
		 getMigrate()->migrate_set_rate (fromPop, toPop, rate);
	return gen;
}

void Event_MigrationRate::addToBaseModel( BaseModel& baseModel ) const {
	baseModel.popInfos[ fromPop ].setMigrRate( toPop, gen, rate );
}

genid Event_Admix::execute() {

	/* 27-sep-04.  Modified to simply shift the required # chroms to another pop. */

	/* join two pops forward = split pops backwards */
	/* change pop size to what is specified. */
	/*		  		  demography->dg_set_pop_size_by_name (currentevent->gen,
								currentevent->popindex[1],
								currentevent->params[0]);  */
		  
	Pop *admixedPopPtr = getDemography()->dg_get_pop_by_name( admixedPop );
	Pop *sourcePopPtr = getDemography()->dg_get_pop_by_name( sourcePop );
	util::chk( admixedPopPtr, "admix::execute - unknown admixed pop" );
	util::chk( sourcePopPtr, "admix::execute - unknown source" );

	if ( !admixedPopPtr->isInactive() && !sourcePopPtr->isInactive() )
		 getDemography()->dg_move_nodes_by_name (admixedPop, sourcePop, admixFrac, gen );
	return gen;
}

void Event_Admix::addToBaseModel( BaseModel& baseModel ) const {
	BOOST_AUTO( &popInfo, baseModel.popInfos[ admixedPop ] );
	BOOST_AUTO( &pieces, popInfo.migrRateTo[ sourcePop ].getPieces() );

	if ( pieces.empty() )
		 pieces[ ZERO_GEN ] =
				math::fn_const< genid >( prob_per_chrom_per_gen_t( 0. ) );
	
	pieces.insert( std::make_pair( gen + gens_t( 1 ), pieces.lower_bound( gen )->second ) );
	pieces.insert( std::make_pair( gen, math::fn_const< genid >( prob_per_chrom_per_gen_t( admixFrac ) ) ) );
}

genid Event_Split::execute() {
	Pop *fromPopPtr = getDemography()->dg_get_pop_by_name( fromPop );
	Pop *newPopPtr = getDemography()->dg_get_pop_by_name( newPop );
	util::chk( fromPopPtr, "split: unknown from pop!" );
	util::chk( newPopPtr, "split: unknown new pop!" );
	//std::cerr << "bef split: gen=" << gen << " from=" << fromPop << " new=" << newPop << " fromPopNodes=" << fromPopPtr->pop_get_num_nodes() << " newPopNodes=" << newPopPtr->pop_get_num_nodes() << "\n";
	/* split two pops forward = join pops backwards */
	getDemography()->dg_move_nodes_by_name ( /* going backwards, move all nodes from */ newPop,
																					 /* into */ fromPop,
																					 /* fractionToMove= */ 1.0,
																					 gen, /* exactFraction= */ true );
	getMigrate()->migrate_delete_all_for_pop( newPop );
	/*		  demography->dg_end_pop_by_name (currentevent->popindex[1]); */
	newPopPtr->makeInactive();
	//std::cerr << "aft split: gen=" << gen << " from=" << fromPop << " new=" << newPop << " fromPopNodes=" << fromPopPtr->pop_get_num_nodes() << " newPopNodes=" << newPopPtr->pop_get_num_nodes() << "\n";
	return gen;
}

void Event_Split::addToBaseModel( BaseModel& baseModel ) const {
	baseModel.popInfos[ newPop ].setMigrRate( fromPop, gen, prob_per_chrom_per_gen_t( 1.0 ) );
}

////////////////////////

}  // namespace histevents

void HistEvents::constructBaseModel( BaseModelP baseModel ) {
	typedef std::pair<genid,EventP> pair_t;
	BOOST_FOREACH( pair_t e, events )
		 if ( e.second->getEventKind() != Event::E_BOTTLENECK &&
					e.second->getEventKind() != Event::E_ADMIX )
				e.second->addToBaseModel( *baseModel );

	//std::cerr << "\nb4 bnecks:" << *baseModel << "\n";

	BOOST_FOREACH( pair_t e, events )
		 if ( e.second->getEventKind() == Event::E_BOTTLENECK ||
					e.second->getEventKind() == Event::E_ADMIX )
				e.second->addToBaseModel( *baseModel );

	events.clear();
#if 0
	size_t prevSize;
	do {
		prevSize = events.size();

		for( BOOST_AUTO( e, events.begin() ); e != events.end();  ++e )  {
			if ( true || e->second->getEventKind() == Event::E_BOTTLENECK ||
					 e->second->getEventKind() == Event::E_POPSIZEEXP ) {
				events.erase( e );
				break;
			}
		}
		
	} while ( events.size() < prevSize );
#endif
	
	
}

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
		else if ( ( typestr == Event_PopSizeExp2::typeStr() ) ||
							( typestr == Event_PopSizeExp2::typeStr2() ) ) event.reset( new Event_PopSizeExp2( this, is ) );
		else if ( typestr == Event_Split::typeStr() ) event.reset( new Event_Split( this, is ) );
		else if ( typestr == Event_MigrationRate::typeStr() ) event.reset( new Event_MigrationRate( this, is ) );
		else if ( typestr == Event_Bottleneck::typeStr() ) event.reset( new Event_Bottleneck( this, is ) );
		else if ( typestr == Event_Admix::typeStr() ) event.reset( new Event_Admix( this, is ) );
		else if ( typestr == Event_Sweep::typeStr() ) event.reset( new Event_Sweep( this, is ) );
		else if ( typestr == Event_SweepNew_typeStr() ) event = make_shared_Event_SweepNew( this, is );
		else if ( typestr == sweep1::Event_SweepOnePop_typeStr() ) event = sweep1::make_shared_Event_SweepOnePop( this, is );
		else if ( typestr == sweep2::Event_SweepOnePop2_typeStr() ) event = sweep2::make_shared_Event_SweepOnePop2( this, is );
		else if ( typestr == sweep3::Event_SweepOnePop3_typeStr() ) event = sweep3::make_shared_Event_SweepOnePop3( this, is );
		else if (  typestr == Event_MSweep::typeStr() ) event.reset( new Event_MSweep( this, is ) );
		else chkCond( False, "could not parse event %s", buffer );
	} catch( ios::failure e ) {
		chkCond( False, "could not parse event %s", buffer );
	}
	return event;
}

void HistEvents::addEvent( EventP event ) {
	events.insert( std::make_pair( event->getGen(), event ) );
}

}  // namespace cosi
