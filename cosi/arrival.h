//
// Header: poissonprocess.h
//
// Generic handling of poisson processes.
//


#ifndef __INCLUDE_COSI_POISSONPROCESS
#define __INCLUDE_COSI_POISSONPROCESS

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/range.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <cosi/defs.h>
#include <cosi/cosirand.h>
#include <cosi/traj.h>
#include <cosi/utils.h>
#include <cosi/generalmath.h>

namespace cosi {

typedef cosi_double rate_t;

class InhomogeneousPoissonProcess {
	 typedef cosi_double integral_t;
public:

	 template <typename Range>
	 InhomogeneousPoissonProcess( const Range traj ) {
		 typedef typename boost::range_iterator<Range>::type iter_t;
		 typedef typename boost::iterator_value<iter_t>::type pair_t;

		 integral_t integralSoFar = 0;
		 integralToGen.addPt( integral_t(0.0), genid(0.0) );
		 genid lastGen(0.0);
		 BOOST_FOREACH( pair_t p, traj ) {
			 genid gen = p.first;
			 rate_t rate = p.second;
			 assert( ( lastGen == ZERO_GEN && gen >= lastGen ) || ( lastGen > ZERO_GEN && gen > lastGen ) );
			 if ( gen > lastGen && rate > 0 ) {
				 integralSoFar += rate * ToDouble( gen - lastGen );
				 //PRINT2( integralSoFar, lastGen );
				 integralToGen.addPt( integralSoFar, gen );
			 }
			 lastGen = gen;
		 }
	 }

	 genid getTimeOfNextEvent( RandGenP randGen, genid fromGen, factor_t factor ) {
		 double s = - log( randGen->random_double() );
		 //PRINT( s );

		 integral_t throughFromGen = integralToGen.getInverse()( fromGen );
		 //PRINT( throughFromGen );

		 //PRINT( integralToGen.evalInv( genid( 30 ) ) );
		 
		 return integralToGen( ( throughFromGen + s ) / ToDouble( factor ) );
	 }
	 
private:
	 
	 math::InterpBiFun<integral_t,genid> integralToGen;
	 
};  // class InhomogeneousPoissonProcess

void testArrival();

//
// Abstract Class: ArrivalProcess
//
class ArrivalProcess {
public:
	 virtual ~ArrivalProcess() { }

	 class Event {
	 public:
			virtual void execute() = 0;
			
			genid getGen() const { return gen; }
	 protected:
			Event( genid gen_ ): gen( gen_ ) { }
	 private:
			genid gen;
	 };

	 typedef boost::shared_ptr<Event> EventP;

	 virtual EventP getNextEvent( genid fromGen, genid maxGen ) = 0;

};  // class ArrivalProcess


//
// Logical type: rate_t
//
// The rate of a Poisson process at a given time point.  Multiplying the rate
// by the length of a time interval gives the expected number of events in that
// interval, so the rate has dimensions "eventCount / time".
//
//COSI_DEFINE_TYPEDVAL(rate_t);

#if 0
//
// Abstract Class: PoissonProcess
//
// A 
//
class PoissonProcess: public ArrivalProcess {

public:
	 virtual ~PoissonProcess() { }
	 
	 // Method: getRateAt
	 // Returns the instantaneous rate of this Poisson process at a given time
	 virtual rate_t getRateAt( genid gen ) = 0;
};  // class PoissonProcess

typedef boost::shared_ptr<PoissonProcess> PoissonProcessP;

class HomogeneousPoissonProcess: public PoissonProcess {
public:
	 HomogeneousPoissonProcess

	 virtual rate_t getRateAt( genid gen ) { return rate; }
private:
	 rate_t rate;
};

class CompoundPoissonProcess: public PoissonProcess {
public:
	 gens_t getTimeUntilNextEvent( genid gen, RandGenP randGen ) {
		 gens_t waitTime( 0.0 );
		 BOOST_FOREACH( PoissonProcessP p, componentProcesses ) 
				
	 }
private:
	 vector< PoissonProcessP > componentProcesses;
};
#endif
}  // namespace cosi


#endif // #ifndef __INCLUDE_COSI_POISSONPROCESS
