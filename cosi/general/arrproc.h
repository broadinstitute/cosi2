#ifndef COSI_INCLUDE_ARRPROC_H
#define COSI_INCLUDE_ARRPROC_H

//
// Header: arrproc.h
//
// Generic arrival processes.
//

#include <vector>
#include <cstdlib>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/value_init.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <cosi/general/typeutil.h>
#include <cosi/general/math/generalmath.h>

namespace cosi {
namespace arrival {

// Class: ArrivalProcess
// An arrival process of events.
template <typename TSpec> class ArrivalProcess;


// Metafunction: TimeType
// The time type used by a given ArrivalProcess
template <typename TSpec> struct TimeType { typedef double type; };

template <typename TSpec> struct TimeDiffType: public DiffType< typename TimeType<TSpec>::type > { };

// Function: waitTime
// Returns the waiting time until the next event of this process.
template <typename TSpec>
typename TimeDiffType<TSpec>::type waitTime( typename TimeType<TSpec>::type fromTime, const ArrivalProcess<TSpec>& );


// Func: executeNextEvent
// Execute the next event of the process.  It is guaranteed that waitTime() has been called first,
// and that the system time has not advanced by more than this time.
template <typename TSpec> void executeNextEvent( ArrivalProcess<TSpec>& p ) {
	p.nextEvent();
	updateAfterExecutingEvent( p );
}

template <typename TSpec>
boost::signals2::signal<void()>* getUpdateSignal( const ArrivalProcess<TSpec>& );

struct Basic;

template <typename TComponentSpec> struct Compound;

template <typename TSpec> class ArrivalProcessImpl;

template <>
class ArrivalProcessImpl< Basic > {
public:
	 std::string name;
	 boost::signals2::signal<void()> signal_update;
	 boost::function< void() > nextEvent;
};

template <typename TComponentSpec>
class ArrivalProcessImpl< Compound<TComponentSpec> >: public ArrivalProcessImpl< Basic > {
public:
	 typedef ArrivalProcess<TComponentSpec> component_process_type;
	 typedef boost::shared_ptr< component_process_type > component_process_p;
	 std::vector< component_process_p > _components;
	 size_t nextEventComponent;
};

template <typename TComponentSpec>
class ArrivalProcess< Compound<TComponentSpec> >: public ArrivalProcessImpl< Compound<TComponentSpec> > {
};

template <typename TComponentSpec, typename TTime>
TTime nextEventTime( TTime fromTime, ArrivalProcess< Compound<TComponentSpec> >& p ) {
	TTime t_next = -1;

	std::cerr << "computing nextEventTime for " << p.name << "\n";
	for ( size_t j = 0; j < p._components.size(); j++ ) {
		TTime t_here = nextEventTime( fromTime, *p._components[ j ] );
		std::cerr << "j=" << j << " p_here=" << p._components[ j ]->name << " t_here=" << t_here << "\n";
		if ( t_here > 0 && ( t_next < 0|| t_here < t_next ) ) {
			std::cerr << "updating t_next at j=" << j << ": t_next was " << t_next << " t_here is " << t_here << "\n";
			t_next = t_here;
			p.nextEventComponent = j;
			p.nextEvent = p._components[ j ]->nextEvent;
			std::cerr << "t_next is now " << t_next << "\n";
		} else {
			std::cerr << "skipping event at " << j << "\n";
		}
	}
	std::cerr << "returning t_next as "<< t_next << "\n";
	return t_next;
}

template <typename TComponentSpec>
void updateAfterExecutingEvent( ArrivalProcess< Compound< TComponentSpec > >& p ) {
	std::cerr << "updating after event: " << p.name << "\n";
	updateAfterExecutingEvent( *p._components[ p.nextEventComponent ] );
	p.nextEvent.clear();
	std::cerr << "updated after event: " << p.name << "\n";
}

template <typename TPoissonSpec=void> struct Poisson;

template <typename TPoissonSpec=void> struct RateFunctionType;

template <>
class ArrivalProcess< Compound< Poisson<> > >: public ArrivalProcessImpl< Compound< Poisson<> > >  {
	 
};

template <typename TSpec = void> struct OneEvent;

template <typename TSpec>
class ArrivalProcessImpl< OneEvent< TSpec > >: public ArrivalProcessImpl< Basic > {
public:
	 typedef typename TimeType< ArrivalProcess< OneEvent< TSpec > > >::type time_type;
	 time_type eventTime;
};

template <typename TSpec>
class ArrivalProcess< OneEvent< TSpec > >: public ArrivalProcessImpl< OneEvent< TSpec > > { };

template <typename TSpec, typename TTime>
TTime nextEventTime( TTime fromTime, const ArrivalProcess< OneEvent<TSpec> >& p ) {
	return p.eventTime >= fromTime ? p.eventTime : static_cast<TTime>( -1 );
}

template <typename TSpec>
void updateAfterExecutingEvent( ArrivalProcess< OneEvent<TSpec> >& p ) {
	std::cerr << "updating after event: " << p.name << "\n";
	p.eventTime = -1;
	p.nextEvent.clear();
	std::cerr << "updated after event: " << p.name << "\n";
}

}  // namespace arrival
}  // namespace cosi

#endif // #ifndef COSI_INCLUDE_ARRPROC_H
