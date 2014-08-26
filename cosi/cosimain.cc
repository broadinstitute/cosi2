//
// Program: cosimain
//
// The top-level executable for the cosi simulator.
// Parses command-line args, and invokes the cosi library to run the simulation.
//

#include <cstdlib>
#include <iostream>
#include <boost/exception/all.hpp>
#include <cosi/cositop.h>

int main( int argc, char *argv[] ) {
	try {
		cosi::CoSiMain cosiMain;
		return cosiMain.cosi_main( argc, argv );
	} catch( const boost::exception& e ) {
		std::cerr << "cosi error: " << boost::diagnostic_information( e ) << std::endl;
		return EXIT_FAILURE;
	}
}

