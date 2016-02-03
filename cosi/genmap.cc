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
#include <string>
#include <typeinfo>
//#include <ctime>
#include <boost/core/demangle.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/exception.hpp>
#include <boost/next_prior.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/exception/all.hpp>
#include <cosi/general/utils.h>
#include <cosi/general/math/cosirand.h>
#include <cosi/genmap.h>

namespace cosi {

// using std::cout;
// using std::endl;
// using std::ofstream;
// using std::ios;
// using std::map;
// using util::chkCond;

GenMap::GenMap( const boost::filesystem::path& fname, const len_bp_t length_ ):
	pd_range_len( length_ ) {
	try {
		try {
			boost::filesystem::ifstream f( fname );
			readFrom( f );
		} catch( const std::ifstream::failure& e ) {
			BOOST_THROW_EXCEPTION( cosi_io_error() << boost::errinfo_nested_exception( boost::copy_exception( e ) ) );
		} catch( const std::exception& e ) {
			BOOST_THROW_EXCEPTION( cosi_io_error() << boost::errinfo_nested_exception( boost::copy_exception( e ) ) );
		}	
	} catch( boost::exception& e ) {
		e << boost::errinfo_file_name( fname.string() )
			<< error_msg3( "Error reading genetic map file" );
		throw;
	}
}

GenMap::GenMap( std::istream& f, const len_bp_t length_ ):
	pd_range_len( length_ ) {
	readFrom( f );
}

void GenMap::readFrom( std::istream& recombfp ) {

	boost::io::ios_exception_saver save_exceptions( recombfp );
	recombfp.exceptions( std::istream::failbit | std::istream::badbit );
	
	pd_locs.clear();
	gd_locs.clear();
	
	//std::cerr.precision(15);

	int lineNum = 0;
	std::string line;
	try {
		try {
			util::SumKeeper<gd_orig_len_t> gd_orig_so_far;
			loc_bp_int_t lastStart = 0;
			double lastRate = -1;

			while ( recombfp ) {
				if ( recombfp.eof() ) break;
				line.clear();
				try {
					std::getline( recombfp, line );
				} catch( ::std::ios_base::failure const& f ) {
					//std::cerr << "genMap::readFrom - caught exception\n";
					break;
				}
				catch( std::exception const& e ) {
					// std::cerr << "genMap::readFrom - caught UNKNOWN exception of type " <<
					// 	 typeid( e ).name() << " demangled " << 
					// 	 ( boost::core::demangle( typeid( e ).name() ) ) << " ; exception is " << e.what() << "\n";
					break;
				}
				catch( ... ) {
					//					std::cerr << "genMap::readFrom - caught UNKNOWN exception\n";
					break;
				}
				//std::cerr << "line: " << line << "\n";
				++lineNum;
				loc_bp_int_t start;
				double rate;

				//std::cerr << "reding line " << line << "\n";
				int nread = sscanf(line.c_str(), "%ld %lf", &start, &rate);
				if ( nread != 2 )
					 BOOST_THROW_EXCEPTION( cosi_io_error() );
				
				if ( !pd_locs.empty() && start == lastStart ) {
					std::cerr << "warning: skipping duplicate genetic map point at " << start << "\n";
					continue;
				}
				if ( !( start > lastStart ) )
					 BOOST_THROW_EXCEPTION( cosi_io_error() <<
																	error_msg( "genetic map locations must be increasing" ) );
				if ( !( rate > 0 ) )
					 BOOST_THROW_EXCEPTION( cosi_io_error() <<
																	error_msg( "zero recomb rate not supported" ) );
			
				if ( pd_locs.empty() ) {
					pd_locs.push_back( static_cast<pd_orig_loc_t>( 0.0 ) );
					gd_locs.push_back( static_cast<gd_orig_loc_t>( 0.0 ) );
					lastStart = 0;
					lastRate = rate;
				}
				gd_orig_so_far.add( static_cast< gd_orig_len_t>( ( start - lastStart ) * lastRate ) );
				pd_locs.push_back( static_cast< pd_orig_loc_t>( start ) );
				gd_locs.push_back( static_cast< gd_orig_loc_t>( 0 ) + gd_orig_so_far.getSum() );
				//std::cerr << pd_locs.back() << "\t" << gd_locs.back() << "\n";
				
				lastStart = start;
				lastRate = rate;
			}  // while ( recombfp )

			if ( pd_range_len > lastStart ) {
				gd_orig_so_far.add( static_cast< gd_orig_len_t>( ( pd_range_len - lastStart ) * lastRate ) );
				pd_locs.push_back( static_cast< pd_orig_loc_t>( pd_range_len ) );
				gd_locs.push_back( static_cast< gd_orig_loc_t>( 0 ) + gd_orig_so_far.getSum() );
			}
			//std::cerr << pd_locs.back() << "\t" << gd_locs.back() << "\n";
			
		} catch( const std::istream::failure& e ) {
			BOOST_THROW_EXCEPTION(
				cosi_io_error()
				<< boost::errinfo_nested_exception( boost::copy_exception( e ) )
				);
		}
	} catch( boost::exception& e ) {
		e << error_bad_line( line ) << boost::errinfo_at_line( lineNum );
		throw;
	}
	setStart( static_cast<pd_orig_loc_t>( 0 ) );
}  // GenMap::readFrom()

// time_t tot_sec;
// long tot_nsec;
// long callCount;
// std::clock_t totClock;

// struct ShowTime {
// 	 ShowTime() { };
// 	 ~ShowTime() {
// //		 if ( btm ) delete btm;
// 		 std::cerr << "clockSec=" << ( ((double)totClock)/((double)CLOCKS_PER_SEC) ) << " tot_sec=" << tot_sec << " tot_nsec=" << tot_nsec << " callCount=" << callCount << "\n";
// 	 }
// } showTime;


void GenMap::setStart( pd_orig_loc_t start ) {
	// std::clock_t startTime = std::clock();
	// ++callCount;
	// struct timespec tm;
	// clock_gettime( CLOCK_REALTIME, &tm );
	BOOST_AUTO( pd_beg_it, boost::lower_bound( pd_locs, start ) );
	util::chkCond( pd_beg_it != pd_locs.end() && *pd_beg_it >= start );
	BOOST_AUTO( pd_end_it, std::upper_bound( boost::next( pd_beg_it ), pd_locs.end(), start + pd_range_len ) );
	util::chkCond( pd_beg_it != pd_end_it );

	BOOST_AUTO( gd_beg_it, gd_locs.begin() + ( pd_beg_it - pd_locs.begin() ) );
	BOOST_AUTO( gd_end_it, gd_beg_it + ( pd_end_it - pd_beg_it ) );

	pd_locs_range.clear();
	gd_locs_range.clear();

	if ( *pd_beg_it > start ) {
		pd_locs_range.push_back( start );
		gd_locs_range.push_back( util::interp( pd_locs, gd_locs, start ) );
	}

	pd_locs_range.insert( pd_locs_range.end(), pd_beg_it, pd_end_it );
	gd_locs_range.insert( gd_locs_range.end(), gd_beg_it, gd_end_it );

	if ( pd_locs_range.back() != start + pd_range_len ) {
		pd_locs_range.push_back( start + pd_range_len );
		gd_locs_range.push_back( util::interp( pd_locs, gd_locs, start + pd_range_len ) );
	}
	
	gd_range_len = gd_locs_range.back() - gd_locs_range.front();
	
	if ( gd_range_len <= static_cast< gd_orig_len_t >( 1e-16) )
		 BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg( "genetic map length is zero!" ) );

	// totClock += ( std::clock() - startTime );
	// struct timespec tm2;
	// clock_gettime( CLOCK_REALTIME, &tm2 );
	// tot_sec += ( tm2.tv_sec - tm.tv_sec );
	// tot_nsec += ( tm2.tv_nsec - tm.tv_nsec );
}


}  // namespace cosi

