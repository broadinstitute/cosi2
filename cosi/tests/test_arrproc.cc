
#define BOOST_TEST_MODULE arrproc_funcs test

#include <utility>
#include <iostream>
#include <string>
#include <boost/test/included/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <cosi/arrproc.h>

void doPrint() { std::cout << "hello1\n"; }
void doPrint2() { std::cout << "hello2\n"; }


BOOST_AUTO_TEST_CASE( arrproc_tests ) {
	using namespace cosi::arrival;

	BOOST_AUTO( e1, boost::make_shared< ArrivalProcess< OneEvent<> > >() );
	e1->name = "event 1";
	e1->eventTime = 5.0;
	e1->nextEvent = doPrint;

	BOOST_AUTO( e2, boost::make_shared< ArrivalProcess< OneEvent<> > >() );
	e2->name = "event 2";
	e2->eventTime = 10.0;
	e2->nextEvent = doPrint2;

	ArrivalProcess< Compound< OneEvent<> > > proc;
	proc.name = "aggregate";
	proc._components.push_back( e1 );
	proc._components.push_back( e2 );

	double t = 0;
	std::cout << "next update: " << ( t = nextEventTime( t, proc ) ) << "\n";
	executeNextEvent( proc );
	std::cout << "next update: " << ( t = nextEventTime( t, proc ) ) << "\n";
	executeNextEvent( proc );
	std::cout << "next update: " << ( t = nextEventTime( t, proc ) ) << "\n";
	
}

