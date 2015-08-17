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
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/exception.hpp>
#include <boost/next_prior.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/exception/all.hpp>
#include <cosi/utils.h>
#include <cosi/cosirand.h>
#include <cosi/genmap.h>

namespace cosi {

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::map;
using util::chkCond;

GenMap::GenMap( const boost::filesystem::path& fname, const len_bp_t length_, ploc_bp_diff_t genMapShift_,
								bool_t genmapRandomRegions_, RandGenP randGen ):
	rec_recombrate( 0.0 ),
	rec_length( length_ ) {
	try {
		try {
			boost::filesystem::ifstream f( fname );

			if ( genmapRandomRegions_ ) {
				// get length of file:
				f.seekg (0, f.end);
				istream::streampos stream_length = f.tellg();
				// find start of last line
				using std::ios_base;
				f.seekg( -1, ios_base::cur );
				if ( f.peek() == '\n' ) f.seekg( -1, ios_base::cur );
				while( f.tellg() > 0  &&  f.peek() != '\n' )
					 f.seekg( -1, ios_base::cur );
				loc_bp_int_t endPos = -1;
				f >> endPos;

				while( true ) {
					istream::streampos seekTo = static_cast< istream::streampos >( randGen->random_double() * stream_length );
					//std::cerr << "seekTo=" << seekTo << "\n";
					f.seekg( seekTo );
					std::string line;
					std::getline( f, line );
					std::getline( f, line );
					loc_bp_int_t begPos;
					double rate;
					if ( sscanf(line.c_str(), "%ld %lf", &begPos, &rate) != 2 )
						 BOOST_THROW_EXCEPTION( cosi_io_error()
																		<< boost::errinfo_file_name( fname.string() ) <<
																		error_msg( "Error reading genetic map file - could not parse line: " + line ) );

					if ( endPos - begPos >= length_ ) {
						genMapShift_ = -begPos;
						break;
					}
				}  // while( true )
			}  // if( genmapRandomRegions_ )
		
			readFrom( f, genMapShift_ );


		} catch( const std::ios_base::failure& e ) {
			BOOST_THROW_EXCEPTION( cosi_io_error() << boost::errinfo_nested_exception( boost::copy_exception( e ) ) );
		}			
	} catch( boost::exception& e ) {
		e << boost::errinfo_file_name( fname.string() )
			<< error_msg( "Error reading genetic map file" );
		throw;
	}
}

GenMap::GenMap( istream& f, const len_bp_t length_, ploc_bp_diff_t genMapShift_ ):
	rec_recombrate( 0.0 ),
	rec_length( length_ ) {
	readFrom( f, genMapShift_ );
}

void GenMap::readFrom( istream& recombfp, ploc_bp_diff_t /* genMapShift_ */ ) {

	boost::io::ios_exception_saver save_exceptions( recombfp );
	recombfp.exceptions( ios::failbit | ios::badbit );
	
	pd_locs.clear();
	gd_locs.clear();
	
	//std::cerr.precision(15);

	int lineNum = 0;
	try {
		try {
			util::SumKeeper<gd_orig_len_t> gd_orig_so_far;
			loc_bp_int_t lastStart = 0;
			double lastRate = -1;

			while ( recombfp ) {
				std::string line;
				try {
					std::getline( recombfp, line );
				} catch( std::ios::failure f ) { break; }
				//std::cerr << "line: " << line << "\n";
				++lineNum;
				loc_bp_int_t start;
				double rate;
			
				int nread = sscanf(line.c_str(), "%ld %lf", &start, &rate);
				if ( !pd_locs.empty() && start == lastStart ) {
					std::cerr << "warning: skipping duplicate genetic map point at " << start << "\n";
					continue;
				}
				if ( nread != 2 )
					 BOOST_THROW_EXCEPTION( cosi_io_error() <<
																	error_msg( "malformed genetic map file line " + line ) );
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

			gd_orig_so_far.add( static_cast< gd_orig_len_t>( ( rec_length - lastStart ) * lastRate ) );
			pd_locs.push_back( static_cast< pd_orig_loc_t>( rec_length ) );
			gd_locs.push_back( static_cast< gd_orig_loc_t>( 0 ) + gd_orig_so_far.getSum() );
			//std::cerr << pd_locs.back() << "\t" << gd_locs.back() << "\n";
			
		} catch( const std::ios_base::failure& e ) {
			BOOST_THROW_EXCEPTION(
				cosi_io_error()
				<< boost::errinfo_nested_exception( boost::copy_exception( e ) )
				);
		}
	} catch( boost::exception& e ) {
		e << boost::errinfo_at_line( lineNum );
		throw;
	}
	pd_locs_range = pd_locs;
	gd_locs_range = gd_locs;

	pd_range_len = rec_length;
//gd_range_len = pd2gd( rec_length ) - static_cast< gd_orig_loc_t>( 0 );
	rec_recombrate = gd_range_len = gd_locs.back() - static_cast< gd_orig_loc_t>( 0 );
}  // GenMap::readFrom()

#if 0 
void GenMap::readFrom( istream& recombfp, ploc_bp_diff_t genMapShift_ ) {

	boost::io::ios_exception_saver save_exceptions( recombfp );
	recombfp.exceptions( ios::failbit | ios::badbit );

	FILE *recombOut = NULL;
	if ( getenv( "COSI_SAVE_RECOMB" ) ) {
		extern int curSimNum;
		char fnbuf[1024];
		sprintf( fnbuf, "%s_%d.recom", getenv( "COSI_SAVE_RECOMB" ), curSimNum );
		recombOut = fopen( fnbuf, "wt" );
	}
	
	try {
		loc_bp_int_t start;
		double rate;

		loc2cumRate.clear();
		cumRate2loc.clear();

		util::SumKeeper<glen_cM_t> cumRateSoFar;
	
		loc_bp_int_t lastStart = 0;
		cosi_double lastRate = 0.0;
		loc2cumRate.addPt( MIN_PLOC, MIN_GLOC );
		cumRate2loc.addPt( MIN_GLOC, MIN_PLOC );


		vector< pair< ploc_t, gloc_cM_t > > ploc2gloc;
		while ( True ) {

			std::string line;
			try {
				std::getline( recombfp, line );
			} catch( std::ios::failure f ) { break; }
			//std::cerr << "line: " << line << "\n";
			if ( sscanf(line.c_str(), "%ld %lf", &start, &rate) == 2  && ( start + genMapShift_ ) <= rec_length ) {
		
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

				if ( recombOut ) {
					if ( lastStart == 0 && start > 1 ) 
						 fprintf( recombOut, "1\t%e\n", rate );
					fprintf( recombOut, "%ld\t%e\n", start, rate );
				}
		
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
			BOOST_FOREACH( pair_t it, ploc2gloc ) {
				gloc_t gloc( ToDouble( it.second ) / ToDouble( rec_recombrate ) );
				loc2cumRate.addPt( it.first, gloc );
				cumRate2loc.addPt( gloc, it.first );
			}
		}
	} catch( const std::ios_base::failure& e ) {
		throw cosi_io_error() << boost::errinfo_nested_exception( boost::current_exception() );
	}
	if ( recombOut ) fclose( recombOut );
	
}  // GenMap::readFrom()
#endif // #if 0 old readFrom  

}  // namespace cosi

