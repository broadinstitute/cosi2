#ifndef COSI_INCLUDE_ARRPROC_H
#define COSI_INCLUDE_ARRPROC_H

//
// Header: arrproc.h
//
// Generic arrival processes.
//

////////////////
// interface 
////////////////

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {

namespace math {
template <typename TDomain, typename TRange, typename TSpec> struct Function;
}

template <typename T> inline
T& setLabel( T& x, std::string lbl ) { x.label = lbl; return x; }

namespace arrival2 {

using math::Function;

// * Class: ArrivalProcess
// An arrival process of events.  The process can give the time to its next event, and execute that
// event.
template <typename TTime, typename TSpec> class ArrivalProcess;

template <typename URNG, typename TSpec = void> struct Stoch;

// ** Func: nextEventTime
// Return the time of the next event of this process, starting with time 'fromTime'.  If the next
// event is past 'maxTime', return maxTime.  Possibly store details of the event in the process
// object, for use by a subsequent call to 'executeNextEvent'.
template <typename TTime, typename URNG, typename TSpec>
TTime nextEventTime( ArrivalProcess<TTime, Stoch< URNG, TSpec> >& p, TTime fromTime, TTime maxTime, URNG& );

// ** Func: executeNextEvent
// Execute the next event of the process.  It is guaranteed that generateNextEvent() has been called first,
// and that no other events (from this process or another) have been executed since then.
template <typename TTime, typename URNG, typename TSpec>
void executeNextEvent( ArrivalProcess<TTime, Stoch< URNG, TSpec > >&, TTime, URNG& );

//
// ** Poisson processes
//

struct CinlarsMethod;
template <typename TRateFnSpec, typename TRateFactor, typename TSpec=CinlarsMethod> class Poisson;

template <typename TTime, typename URNG, typename TRateFnSpec, typename TRateFactor> struct
ArrivalProcess<TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > >;

// template <typename TTime, typename TSpec> struct EventRunner;
// template <typename TTime, typename URNG, typename TSpec> struct EventRunner<TTime, Stoch<URNG, TSpec> > {
// 	 virtual void executeNextEvent( TTime eventTime, URNG& urng  ) = 0;
// };

// template <typename TTime, typename URNG>
// struct EventRunner<TTime, Stoch<URNG, Poisson<> > > {
// 	 virtual void executeNextEvent( TTime eventTime, URNG& urng  ) = 0;
// 	 virtual double getRateFactor() const = 0;
// };

//template <typename TTime, typename URNG, typename TRateFnSpec> struct result_of_makePoissonProcess;

// template <typename TTime, typename URNG, typename TRateFnSpec>
// typename result_of_makePoissonProcess<TTime,URNG,TRateFnSpec>::type
// makePoissonProcess( Function< TTime, typename RateType<TTime>::type, TRateFnSpec> const& rateFn,
// 										boost::function< double (TTime) > rateFactorFn,
// 										boost::function< void ( TTime, URNG ) > eventAction );

struct AnyProc;

template <typename TTime, typename URNG> struct ArrivalProcess<TTime, Stoch<URNG, AnyProc> >;

template <typename TTime, typename URNG, typename TSpec>
ArrivalProcess<TTime, Stoch<URNG, AnyProc> >
makeAnyProc( ArrivalProcess< TTime, Stoch< URNG, TSpec > > const& );

template <typename TComponentSpec, typename TSpec = void> struct Compound;

template <typename TTime, typename URNG, typename TComponentSpec>
struct ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > >;

template <typename TTime, typename URNG, typename TComponentSpec>
ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > >
makeCompoundProc( );

template< typename TTime, typename URNG, typename TComponentSpec >
void add( ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > >&,
					boost::shared_ptr< ArrivalProcess<TTime, Stoch<URNG, TComponentSpec> > > );

}  // namespace arrrival
}  // namespace cosi


/////////////////////////
// implementation
/////////////////////////

#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/value_init.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <cosi/general/typeutil.h>
#include <cosi/general/math/generalmath.h>

namespace cosi {
namespace arrival2 {

using math::IsFunctionSpec;
using math::SpecType;
using math::integralFromPoint;
using math::DiffType;

template <typename TTime, typename TSpec> struct ArrivalProcessImpl;


template <typename TTime, typename URNG>
struct ArrivalProcessImpl< TTime, Stoch< URNG > > {
	BOOST_CONCEPT_ASSERT(( boost::LessThanComparable<TTime> ));

};

template <typename TTime, typename TRateFnSpec, typename TRateFactor>
struct get_rate_fn_integral {
	 BOOST_MPL_ASSERT(( IsFunctionSpec<TRateFnSpec> ));
	 //BOOST_MPL_ASSERT(( boost::is_convertible(

   typedef BOOST_TYPEOF_TPL( boost::declval<TTime>() - boost::declval<TTime>()) time_diff_type;
	 typedef BOOST_TYPEOF_TPL(( boost::declval<double>() / boost::declval<time_diff_type>()
															/ boost::declval<TRateFactor>() )) rate_type;
	 
   typedef Function< TTime, rate_type, TRateFnSpec > rate_fn_type;
   typedef BOOST_TYPEOF_TPL(( integralFromPoint( boost::declval<rate_fn_type>(),
                                                 boost::declval<TTime>() ) ))
	 rate_fn_integral_type;
	 typedef rate_fn_integral_type type;

	 typedef typename SpecType<rate_fn_integral_type>::type rate_fn_integral_spec;
	 

   typedef typename rate_fn_integral_type::result_type integral_val_type;
   typedef typename DiffType< integral_val_type >::type
   integral_diff_type;
	 
};

template <typename TTime, typename URNG, typename TRateFnSpec, typename TRateFactor>
struct ArrivalProcessImpl< TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > >:
		 public ArrivalProcessImpl< TTime, Stoch< URNG > > {
};

template <typename TTime, typename URNG, typename TRateFactor>
struct ArrivalProcessDef {
	 virtual TRateFactor getRateFactor() const = 0;
	 virtual void executeEvent( TTime, URNG& ) = 0;
};

template <typename TTime, typename URNG, typename TRateFnSpec, typename TRateFactor>
class ArrivalProcess<TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > >:
		 public ArrivalProcessImpl< TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TRateFnSpec> ));
	 typedef ArrivalProcessImpl< TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > > PARENT;
public:

	 typedef typename get_rate_fn_integral<TTime, TRateFnSpec, TRateFactor>::rate_fn_type rate_fn_type;
	 typedef typename get_rate_fn_integral<TTime, TRateFnSpec, TRateFactor>::type rate_fn_integral_type;

	 ArrivalProcess( rate_fn_type const& rateFn_, TTime startTime_,
									 boost::shared_ptr< ArrivalProcessDef<TTime, URNG, TRateFactor> > processDef_,
									 std::string label_="" ):
		 rateFnIntegral( integralFromPoint( rateFn_, startTime_ ) ),
		 processDef( processDef_ ),
		 label( label_ ) { }

   TTime nextArrivalTime( TTime fromTime, TTime maxTime, TRateFactor rateFactor_, URNG& randGen,
													double eps = 1e-5, int maxSteps = 100000  ) const {
		 double rateFactor = ToDouble( rateFactor_ );
		 assert( rateFactor >= 0 );
     assert( fromTime >= TTime(0.0) );
     assert( maxTime >= fromTime );

		 typedef typename get_rate_fn_integral<TTime, TRateFnSpec, TRateFactor>::integral_diff_type integral_diff_type;

		 if ( rateFactor == 0  ||  fromTime == maxTime ) return maxTime;
     BOOST_AUTO_TPL(integralAtFromTime, ( eval( rateFnIntegral, fromTime )  ));
     BOOST_AUTO_TPL(integralAtMaxTime, ( eval( rateFnIntegral, maxTime ) ));
     assert( integralAtMaxTime >= integralAtFromTime );
		 //PRINT7( label, fromTime, maxTime, rateFactor, eps, integralAtFromTime, integralAtMaxTime );
		 if ( integralAtMaxTime == integralAtFromTime ) return maxTime;
#ifndef NDEBUG
#ifdef COSI_DEV_PRINT_GENERALMATH
		 PRINT( rateFnIntegral );
#endif		 
     BOOST_TYPEOF_TPL( integralAtFromTime ) zeroIntegral(0.0);
     assert( integralAtFromTime >= zeroIntegral );
     assert( integralAtMaxTime >= zeroIntegral );
#endif
		 
     typename boost::random::uniform_01<> uniformDistr;
     BOOST_AUTO_TPL( u,-std::log( uniformDistr( randGen ) ) );
     BOOST_AUTO_TPL(targetIntegralValue, integralAtFromTime +
                    static_cast< integral_diff_type >( u / rateFactor ) );
     if ( integralAtMaxTime < targetIntegralValue ) return maxTime;
     TTime  targetTime =  evalInverse( rateFnIntegral, fromTime, maxTime, targetIntegralValue,
                                       static_cast<integral_diff_type>( eps / ( 10 * rateFactor ) ),
																			 maxSteps );

     // if ( !(std::abs( targetTime - targetIntegralValue ) < eps ) ) {
     //     PRINT2( std::abs( targetTime - targetIntegralValue ), eps );
     // }
     assert(( equal_eps( rateFnIntegral( targetTime ), targetIntegralValue ), eps ));
#ifndef NDEBUG
		 if ( !( equal_eps( ( rateFnIntegral( targetTime ) - rateFnIntegral( fromTime ) )
												* rateFactor, static_cast< integral_diff_type >(u),  eps ) ) ) {
			 PRINT12( fromTime, targetTime, targetTime - fromTime, u, rateFnIntegral( targetTime ), rateFnIntegral( fromTime ),
								targetIntegralValue, rateFnIntegral( targetTime ) - targetIntegralValue,
								rateFactor, rateFnIntegral( targetTime ) - rateFnIntegral( fromTime ),
								( rateFnIntegral( targetTime ) - rateFnIntegral( fromTime ) ) * rateFactor,
								ToDouble( ( rateFnIntegral( targetTime ) - rateFnIntegral( fromTime ) ) * rateFactor ) - u );
							
		 }
#endif		 
     assert(( equal_eps( ( rateFnIntegral( targetTime ) - rateFnIntegral( fromTime ) )
                         * rateFactor, static_cast< integral_diff_type >(u), eps ) ));
     return targetTime;
//     return std::min( maxTime, findIntegralInverse( rateFnIntegral( fromTime ) ;
   }
   const rate_fn_integral_type& getRateFnIntegral() const { return rateFnIntegral; }
   
   friend std::ostream& operator<<( std::ostream& s, const ArrivalProcess& f ) {
		 s << "ArrivalProcess[label=" << f.label << ", rateFnIntegral=" << f.rateFnIntegral << "]";
		 return s;
	 }

	 std::string getLabel() const { return label; }

   rate_fn_integral_type rateFnIntegral;
	 boost::shared_ptr< ArrivalProcessDef<TTime, URNG, TRateFactor> > processDef;
	 std::string label;
   //integral_val_type integralAtZeroTime;
};  // class ArrivalProcess<TTime, Poisson< TRate, NonHomogeneous< TRateFnSpec > > >

template <typename TTime, typename URNG, typename TRateFnSpec, typename TRateFactor>
TTime nextEventTime( ArrivalProcess<TTime, Stoch< URNG, Poisson< TRateFnSpec, TRateFactor > > >& p,
										 TTime fromTime, TTime maxTime, URNG& urng ) {
	return p.nextArrivalTime( fromTime, maxTime, p.processDef->getRateFactor(), urng );
}

template <typename TTime, typename URNG, typename TSpec>
void executeNextEvent( ArrivalProcess<TTime, Stoch<URNG, TSpec> >& p, TTime t, URNG& urng ) {
	p.processDef->executeEvent( t, urng );
}


template <typename TTime, typename URNG, typename TRateFactor>
class ArrivalProcess<TTime, Stoch< URNG, Poisson< math::Const<>, TRateFactor > > > {
public:

   typedef BOOST_TYPEOF_TPL( boost::declval<TTime>() - boost::declval<TTime>()) time_diff_type;
	 typedef BOOST_TYPEOF_TPL(( boost::declval<double>() / boost::declval<time_diff_type>()
															/ boost::declval<TRateFactor>() )) rate_type;
	 
   typedef Function< TTime, rate_type, math::Const<> > rate_fn_type;

	 ArrivalProcess( rate_fn_type const& rateFn_,
									 boost::shared_ptr< ArrivalProcessDef<TTime, URNG, TRateFactor> > processDef_,
									 std::string label_="" ):
		 rateFn( rateFn_ ),
		 processDef( processDef_ ),
		 label( label_ ) { }

   
   friend std::ostream& operator<<( std::ostream& s, const ArrivalProcess& f ) {
		 s << "ArrivalProcess[label=" << f.label << ", rateFn=" << f.rateFn << "]";
		 return s;
	 }

	 std::string getLabel() const { return label; }

   rate_fn_type rateFn;
	 boost::shared_ptr< ArrivalProcessDef<TTime, URNG, TRateFactor> > processDef;
	 std::string label;
   //integral_val_type integralAtZeroTime;
};  // class ArrivalProcess<TTime, Poisson< TRate, NonHomogeneous< TRateFnSpec > > >


template <typename TTime, typename URNG, typename TRateFactor>
TTime nextEventTime( ArrivalProcess<TTime, Stoch< URNG, Poisson< math::Const<>, TRateFactor > > >& p,
										 TTime fromTime, TTime maxTime, URNG& urng ) {
	typedef BOOST_TYPEOF_TPL( boost::declval<TTime>() - boost::declval<TTime>()) time_diff_type;
	double poisRate = ToDouble( p.processDef->getRateFactor() * evalConst( p.rateFn ) );

	if ( poisRate < 1e-30 ) return maxTime;
															
	return fromTime +
		 time_diff_type( urng.poisson_get_next( poisRate ) );
}


// template <typename TTime, typename URNG, typename TRateFnSpec >
// struct result_of_makePoissonProcess {
//    BOOST_MPL_ASSERT(( IsFunctionSpec<TRateFnSpec> ));
//    typedef ArrivalProcess<TTime, Stoch< URNG, Poisson< TRateFnSpec > > > type;
// };

// template <typename TTime, typename URNG, typename TRateFnSpec >
// typename result_of_makePoissonProcess<TTime, URNG, TRateFnSpec>::type
// makePoissonProcess( const Function< TTime, typename RateType<TTime>::type, TRateFnSpec >& rateFn,
// 										TTime fromTime,
// 										boost::function< double( TTime ) > rateFactorFn_,
// 										boost::function< void( TTime, URNG& ) > eventAction_,
// 										std::string label ) {
//   typedef typename result_of_makePoissonProcess<TTime, URNG, TRateFnSpec>::type
//      result_type;
//   return result_type( rateFn, fromTime, rateFactorFn_, eventAction_, label );
// }

template <typename TTime, typename URNG, typename TComponentSpec>
struct ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > > {
	 typedef ArrivalProcess< TTime, Stoch< URNG, TComponentSpec > > component_type;
	 
	 std::vector< component_type > procs;
	 component_type *nextEvtProc;

	 ArrivalProcess(): nextEvtProc( NULL ), label("compound") { }

	 friend std::ostream& operator<<( std::ostream& s, ArrivalProcess const& p ) {
		 s << "[CompoundArrProc: ";
		 std::copy( p.procs.begin(), p.procs.end(), std::ostream_iterator< component_type >( s, "," ) );
		 s   << "]";
		 return s;
	 }

	 std::string getLabel() const { return label; }
	 std::string label;
	 
};  // ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > >

template< typename TTime, typename URNG, typename TComponentSpec >
void add( ArrivalProcess<TTime, Stoch<URNG, Compound< TComponentSpec > > >& compoundProcess,
					ArrivalProcess<TTime, Stoch<URNG, TComponentSpec> > const& newComponent ) {
	compoundProcess.procs.push_back( newComponent );
}


template <typename TTime, typename URNG, typename TComponentSpec>
TTime nextEventTime( ArrivalProcess< TTime, Stoch< URNG, Compound< TComponentSpec > > >& p, TTime fromTime,
										 TTime maxTime, URNG& urng ) {
	TTime curMaxTime = maxTime;
	p.nextEvtProc = NULL;
	typedef typename ArrivalProcess< TTime, Stoch< URNG, Compound< TComponentSpec > > >::component_type component_type;
	BOOST_FOREACH( component_type& c, p.procs ) {
		TTime nextTimeHere = nextEventTime( c, fromTime, curMaxTime, urng );
		//std::cerr << "process: " << c << " nextevt: " << nextTimeHere << " curMaxTime=" << curMaxTime << "\n";
		if ( nextTimeHere < curMaxTime ) {
			curMaxTime = nextTimeHere;
			p.nextEvtProc = &c;
			//std::cerr << "   -> new curMaxTime: " << curMaxTime << " from process " << c.getLabel() << "\n"; 
		}
	}
	return curMaxTime;
}

template <typename TTime, typename URNG, typename TComponentSpec>
void executeNextEvent( ArrivalProcess< TTime, Stoch< URNG, Compound< TComponentSpec > > >& p, TTime t, URNG& urng ) {
	if ( p.nextEvtProc ) executeNextEvent( *p.nextEvtProc, t, urng );
}

template <typename TTime, typename URNG>
class ArrivalProcess<TTime, Stoch< URNG, AnyProc > > {
   struct ArrivalProcessConcept {
      virtual ~ArrivalProcessConcept() {}
      virtual TTime do_nextEventTime( TTime fromTime, TTime maxTime, URNG& ) = 0;
      virtual void do_executeNextEvent( TTime t, URNG& ) = 0;
			virtual std::ostream& print( std::ostream& s ) const = 0;
			virtual std::string doGetLabel() const = 0;
   };
   template <typename TSpec>
   class ArrivalProcessModel: public ArrivalProcessConcept {
   public:
      ArrivalProcessModel( const ArrivalProcess<TTime, Stoch< URNG, TSpec> >& proc_ ): proc( proc_ ) {}
      virtual ~ArrivalProcessModel() {}

			
      virtual TTime do_nextEventTime( TTime fromTime, TTime maxTime, URNG& urng ) {
				return nextEventTime( proc, fromTime, maxTime, urng );
			}
      virtual void do_executeNextEvent( TTime t, URNG& urng ) {
				executeNextEvent( proc, t, urng );
			}
			
			virtual std::ostream& print( std::ostream& s ) const { s << proc; return s; }
			virtual std::string doGetLabel() const { return proc.getLabel(); }
   private:
      ArrivalProcess<TTime, Stoch< URNG, TSpec> > proc;
   };

   boost::shared_ptr<ArrivalProcessConcept> object;
public:
   ArrivalProcess() { }
   
   template< typename TSpec >
   ArrivalProcess( const ArrivalProcess< TTime, Stoch< URNG, TSpec > >& proc ):
     object( new ArrivalProcessModel<TSpec>( proc ) ) { }
   
   template< typename TSpec >
   void reset( const ArrivalProcess< TTime, Stoch< URNG, TSpec > >& proc ) {
     object.reset( new ArrivalProcessModel<TSpec>( proc ) );
   }
   bool empty() const { return !object.get(); }
   operator bool() const { return !empty(); }
   void reset() { object.reset(); }
   
   TTime nextArrivalTime( TTime fromTime, TTime maxTime, URNG& urng ) {
     return object->do_nextEventTime( fromTime, maxTime, urng );
   }

	 void execNextEvent( TTime t, URNG& urng ) { object->do_executeNextEvent( t, urng ); }
	 
	 std::ostream& print( std::ostream& s ) const { return object->print( s ); }

	 std::string getLabel() const { return object->doGetLabel(); }
   
   friend std::ostream& operator<<( std::ostream& s, const ArrivalProcess& f ) {
		 return f.print( s );
	 }
//	 friend std::string getLabel<>( const ArrivalProcess& p ); // { return p.object->doGetLabel(); }
}; // class ArrivalProcess<TTime, Any< TRand > >

template <typename TTime, typename URNG>
TTime nextEventTime( ArrivalProcess<TTime, Stoch< URNG, AnyProc > >& p,
										 TTime fromTime, TTime maxTime, URNG& urng ) {
	return p.nextArrivalTime( fromTime, maxTime, urng );
}

template <typename TTime, typename URNG>
void executeNextEvent( ArrivalProcess<TTime, Stoch<URNG, AnyProc> >& p, TTime t, URNG& urng ) {
	p.execNextEvent( t, urng );
}

#if 0


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

template <typename TPoissonSpec=void> struct RateFnType;

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

#endif  // #if 0

}  // namespace arrival
}  // namespace cosi

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::ArrivalProcess,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::Stoch,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::ArrivalProcess,2)
BOOST_TYPEOF_REGISTER_TYPE(cosi::arrival2::CinlarsMethod)
BOOST_TYPEOF_REGISTER_TYPE(cosi::arrival2::AnyProc)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::result_of_makePoissonProcess,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::RateType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::Compound,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::ArrivalProcessImpl,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::get_rate_fn_integral,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::arrival2::result_of_makePoissonProcess,3)

#endif // #ifndef COSI_INCLUDE_ARRPROC_H
