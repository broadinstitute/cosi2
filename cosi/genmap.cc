/* $Id: recomb.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

/*
  File: genmap.cc

	Maintains the genetic map (the recombination rate at each genomic position), and provides a method to choose
	a recombination point according to this map.
*/

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <ios>
#include <fstream>
#include <vector>
#include <map>
#include <boost/next_prior.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/fstream.hpp>
#include <cosi/genmap.h>
#include <cosi/utils.h>

namespace cosi {

#define ForEach BOOST_FOREACH  

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::map;
using util::chkCond;

GenMap::GenMap( const boost::filesystem::path& fname, const len_bp_t length_, ploc_bp_diff_t genMapShift_ ):
	rec_recombrate( 0.0 ),
	rec_length( length_ ) {
	try {
	  boost::filesystem::ifstream f( fname );
		readFrom( f, genMapShift_ );
	} catch( boost::exception& e ) {
		BOOST_THROW_EXCEPTION( cosi_io_error()
													 << boost::errinfo_file_name( fname.string() ) <<
													 error_msg( "Error reading genetic map file" ) );
	}
	catch( const std::ios_base::failure& e ) {
		BOOST_THROW_EXCEPTION( cosi_io_error()
													 << boost::errinfo_file_name( fname.string() ) <<
													 error_msg( "Error opening genetic map file" ) );
	}	
}

GenMap::GenMap( istream& f, const len_bp_t length_, ploc_bp_diff_t genMapShift_ ):
	rec_recombrate( 0.0 ),
	rec_length( length_ ) {
	readFrom( f, genMapShift_ );
}

void GenMap::readFrom( istream& recombfp, ploc_bp_diff_t genMapShift_ ) {

	boost::io::ios_exception_saver save_exceptions( recombfp );
	recombfp.exceptions( ios::failbit | ios::badbit );

	try {
		long start;
		double rate;

		loc2cumRate.clear();
		cumRate2loc.clear();

		util::SumKeeper<glen_cM_t> cumRateSoFar;
	
		long lastStart = 0;
		cosi_double lastRate = 0.0;
		loc2cumRate.addPt( MIN_PLOC, MIN_GLOC );
		cumRate2loc.addPt( MIN_GLOC, MIN_PLOC );

		vector< pair< ploc_t, gloc_cM_t > > ploc2gloc;
		while ( True ) {

			std::string line;
			try {
				std::getline( recombfp, line );
			} catch( std::ios::failure f ) { break; }
			if ( sscanf(line.c_str(), "%ld %lf", &start, &rate) == 2  && start <= rec_length ) {
		
				if ( !ploc2gloc.empty() && start == lastStart ) {
					std::cerr << "warning: skipping duplicate genetic map point at " << start << "\n";
					continue;
				}
				start += genMapShift_;
				chkCond( ploc2gloc.empty() || start > lastStart, "genmap: genetic map locations must be increasing" );
				chkCond( rate > 0., "genmap: zero recomb rate not supported" );

				// ignore locations outside the simulated region; these would normally appear if genMapShift_ != 0
				if ( start < 0 ) continue;
				if ( start > rec_length ) break;
		
				ploc_t ploc( ((cosi_double)start) / rec_length );

				if ( ploc2gloc.empty() && ploc > MIN_PLOC )
					 // to ensure that every segment of positive physical length has positive genetic length,
					 // let the first entry in the genmap file specify the recomb rate to its left as well as
					 // to its right.  This makes it symmetric with the last entry uhich sets the distance from
					 // that point to MAX_PLOC.
					 lastRate = rate;

				cumRateSoFar.add( glen_cM_t( ( start - lastStart ) * lastRate ) );
				ploc2gloc.push_back( make_pair( ploc, ZERO_GLOC_CM + cumRateSoFar.getSum() ) );
													 
				lastStart = start;
				lastRate = rate;
			} else break;
		}

		{
			cumRateSoFar.add( glen_cM_t( ( rec_length - lastStart ) * lastRate ) );
			rec_recombrate = cumRateSoFar.getSum();
			ploc2gloc.push_back( make_pair( MAX_PLOC, ZERO_GLOC_CM + rec_recombrate ) );

			typedef pair< ploc_t, gloc_cM_t > pair_t;
			ForEach( pair_t it, ploc2gloc ) {
				gloc_t gloc( ToDouble( it.second ) / ToDouble( rec_recombrate ) );
				loc2cumRate.addPt( it.first, gloc );
				cumRate2loc.addPt( gloc, it.first );
			}
		}
	} catch( const std::ios_base::failure& e ) {
		throw cosi_io_error() << boost::errinfo_nested_exception( boost::current_exception() );
	}
}  // GenMap::readFrom()
  

}  // namespace cosi

