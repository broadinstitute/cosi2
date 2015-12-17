/* $Id: historical.h,v 1.2 2011/05/03 18:50:54 sfs Exp $ */

/*
 * Header: historical.h
 *
 * Representation of <historical events>: population size changes, bottlenecks, migrations, sweeps, etc.
 */


#ifndef __INCLUDE_COSI_HISTORICAL_H
#define __INCLUDE_COSI_HISTORICAL_H

//#define COSI_DEV_PRINT

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/foreach.hpp>
#include <cosi/utils.h>
#include <cosi/cosirand.h>
#include <cosi/defs.h>
#include <cosi/decls.h>

namespace cosi {

using std::string;
using std::istream;
using std::map;
using util::chkCond;

#define HE_LABEL_TOKENS "\""

//
// Class: HistEvents
//
// Keeps track of historical events during a simulation.  Events can affect parameters of the demographic model
// (pop sizes, migration rates), join or split pops, or execute entire segments of the simulation such as
// sweeps or bottlenecks.
//
class HistEvents: public HasRandGen {
public:

	 /////////////////////////////////////////////////////////////////////////////////

	 //
	 // Method group: Setup
	 //
	 // Methods called when initializing the simulation and parsing the parameter file

	 HistEvents( DemographyP demography_ );

	 class Event;
	 typedef boost::shared_ptr<Event> EventP;

	 EventP parseEvent( const char *buffer );
	 void addEvent( EventP event );

	 void constructBaseModel( BaseModelP );

	 void historical_setMigrate( MigrateP migrate_ ) { migrate = migrate_; }
	 void historical_setSweep( SweepP sweep_ ) { sweep = sweep_; }

	 // Method group: During sim
	 //
	 // Methods called during the simulation.
	 
	 // MethodP: historical_get_time_till_next_event
	 // Returns the amount of time until the next historical event, from the current time (passed in).
	 gens_t historical_get_time_till_next_event (genid gen) const;

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
	 genid historical_event_execute (genid historical_event_time);

	 void processSimEnd( genid simEndTime ) {
		 typedef std::pair<genid,EventP> pair_t;
		 BOOST_FOREACH( pair_t p, events )
			 p.second->processSimEnd( simEndTime );
	 }

	 // Abstract class: Event
	 // Base class for historical events.
	 class Event: virtual public boost::enable_shared_from_this<Event> {
	 public:

			enum eventKind_t { E_POPSIZE, E_POPSIZEEXP, E_SPLIT, E_MIGRATIONRATE,
												 E_BOTTLENECK, E_ADMIX, E_SWEEP };

			// Method: getGen
			// Returns the time when this events occurs (for instantenous events),
			// or -- for events taking place over a time interval -- the time when the event
			// starts (going pastward ) / ends (going forward).
			genid getGen() const { return gen; }

			// Method: getLabel
			// Returns the human-readable label of the event.  This is a comment used for messages,
			// does not affect the simulation.
			const string& getLabel() const { return label; }

			// Pure virtual method: execute
			// Execute the event (going pastward), modifying the state of the simulation (demography, migration rates etc) according
			// to the particular historical event.
			// Returns: the generation after the execution of the event (moving pastward), i.e. the generation at which the
			// backwards simulation should continue after this event has completed.
			virtual genid execute() = 0;

			virtual eventKind_t getEventKind() const = 0;

			virtual void processSimEnd( genid /*gen*/ ) { }

			virtual void addToBaseModel( BaseModel& ) const { }
			
	 private:
			// Field: histEventsP
			// The HistEvents object to which this event belongs.  Provides access
			// to <HistEvents::demography>, <HistEvents::migrate> etc, which can be used
			// in implementations of <Event::execute()>.
			HistEvents *histEvents;
			
			// Field: label
			// A human-readable description of this event; does not affect the simulation.
			string label;

	 protected:
			// Field: gen
			// When does this event occur (for instantaneous events), or -- for events taking place over a time interval --
			// when does it start (going pastward ) / end (going forward)?
			genid gen;
			
			Event( HistEvents *histEvents_, const string& label_, genid gen_ ): histEvents( histEvents_ ), label( label_ ), gen( gen_ ) {}
			Event( HistEvents *histEvents_, istream& is ): histEvents( histEvents_ ), gen( NULL_GEN ) {
				char quote;
				is >> quote;
				std::getline( is, label, '\"' );
			}
			virtual ~Event();

			// Method group: Forwarding methods
			// Methods used by concrete implementations of Event to affect the state of the simulation.
			DemographyP getDemography() const { return histEvents->getDemography(); }
			MigrateP getMigrate() const { return histEvents->getMigrate(); }
			SweepP getSweep() const { return histEvents->getSweep(); }
			
			void addEvent( EventP event ) { histEvents->addEvent( event ); }
			
			double poisson_get_next( double rate ) { return histEvents->poisson_get_next( rate ); }
			double random_double() { return histEvents->random_double(); }
			RandGenP getRandGen() { return histEvents->getRandGen(); }
			
	 };  // class Event
	 
private:
	 DemographyP demography;
	 MigrateP migrate;
	 SweepP sweep;

	 friend class Event;
	 friend class Event_SweepNew;

	 // Field: events
	 // All historical events, chronologically ordered (going pastward).  Map from starting point of an event
	 // (in the moving-pastward sense) to the <Event> object.  Note that several historical events may in principle
	 // happen at the same time point, so this must be a multimap; although, what happens in that case has not been
	 // fully checked out, so ideally should be avoided.  FIXME
	 std::multimap<genid, EventP> events;

	 EventP getCurEvent() const { chkCond( !events.empty(), "events empty!" ); return events.begin()->second; }

	 DemographyP getDemography() const { return demography; }
	 MigrateP getMigrate() const { return migrate; }
	 SweepP getSweep() const { return sweep; }
	 
};  // class HistEvents

}  // namespace cosi

#endif 
// #ifndef __INCLUDE_COSI_HISTORICAL_H

