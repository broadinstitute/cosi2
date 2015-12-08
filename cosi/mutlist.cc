/* $Id: mutlist.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ios>
#include <boost/make_shared.hpp>
#include <boost/range.hpp>
#include <boost/next_prior.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>
#include <cosi/utils.h>
#include <cosi/mutlist.h>
#include <cosi/stats.h>
#include <cosi/genmap.h>

namespace cosi {

#define ForEach BOOST_FOREACH
  
using std::cout;
using std::ios;
using std::endl;
using std::ofstream;
using std::make_pair;
using std::string;
namespace lam = boost::lambda;
namespace ran = boost::range;

Mut::Mut() { }
Mut::Mut( loc_t loc_, leafset_p leaves_, genid gen_, popid popName_  ): loc( loc_ ), leaves( (leafset_p )leaves_ ),
																																				gen( gen_ ), popName( popName_ ) { }

void Mutlist::addMut( loc_t loc, leafset_p leaves, genid gen, popid popName ) {
	assert( !frozen );
  muts.push_back( Mut( loc, leaves, gen, popName ) );
}

// const Mut& Mutlist::find( loc_t loc ) const {
// 	auto it = muts.find( loc );
// 	assert( it != muts.end() );
// 	return it->second;
// }

struct CmpIntLoc {
	 len_bp_t length;

	 CmpIntLoc( len_bp_t length_ ): length( length_ ) {}
	 bool_t operator()( const Mut& m1, const Mut& m2 ) { return m1.get_int_loc( length ) == m2.get_int_loc( length ); }
};

void Mutlist::freeze( bool_t inf_sites, len_bp_t length ) {
	if ( frozen ) return;

	std::sort( muts.begin(), muts.end() );

  int i = 0;
  for ( iterator mi = muts.begin(); mi != muts.end(); mi++, i++ )
		 mi->mutIdOrig = i;

  if ( !inf_sites )
		 // Remove mutations that fall on the same integer (basepair) coordinate.
		 muts.erase(std::unique(muts.begin(), muts.end(), CmpIntLoc( length )),  muts.end() );

	////////
	i = 0;
  for ( iterator mi = muts.begin(); mi != muts.end(); mi++, i++ )
		 mi->mutId = i;

	////////
#ifndef COSI_LEAFSET_SIZEONLY
#ifndef COSI_FREQONLY	
	leaf2muts.resize( ToInt( leafset_get_max_leaf_id() ) );
	for ( const_iterator it = muts.begin(); it != muts.end(); it++ ) {
		vector< leaf_id_t > leavesVec;
		leafset_get_leaves( it->leaves, std::back_inserter( leavesVec ) );
		ForEach( leaf_id_t leaf_it, leavesVec )
			 leaf2muts[ ToInt( leaf_it ) ].push_back( it );
	}
#endif	
#endif	
  ///////

#if 0
	//
	// Record the genetic map position of each mutation
	//
	
  cM_t genDistShift = 0.0;

  iterator pmi;
  for ( iterator mi = muts.begin(); mi != muts.end(); pmi = mi, mi++ ) {
		Mut& m = mi->second;
		if ( mi != muts.begin() ) {
			const Mut& pm = pmi->second;
			assert( m.loc > pm.loc );
			assert( m.locGd > pm.locGd );
			len_bp_t gapLen = get_phys_len( m.loc - pm.loc ) * genMap->recomb_get_length();
			double factor = ( gapLen > 200000 ) ? 0.0 : ( ( gapLen > 20000 ) ? (20000.0 / ((double)gapLen)) : 1.0 );
			cM_t genDist = m.locGd - pm.locGd;
			genDistShift += (1.0f - factor) * genDist;
			m.locGdAdj -= genDistShift;
		}
  }
#endif	

	frozen = True;
}  // Mutlist::freeze()

#if 0
vector< Mut >
Mutlist::getMutsInDir( loc_t coreLoc, dir_t dir ) const {
	chkCond( False, "unimpl" );
  namespace ad = boost::adaptors;

  const_iterator it = muts.lower_bound( coreLoc );
	assert( it != muts.end() );

  vector< Mut > result;
#if 0	
  if ( dir == DIR_L )
		 boost::push_back( result,  make_pair<const_iterator,const_iterator>( muts.begin(), boost::next(it) ) | ad::map_values | ad::reversed );
  else {
		while ( it != muts.end() ) {
			result.push_back( it->second );
			it++;
		}
		//	 push_back( result, make_pair( it, muts.end() ) | map_values );
  }
#endif	
  return result;
}

const Mut& Mutlist::getMut( loc_t loc ) const {
  flat_map< loc_t, Mut >::const_iterator it = muts.find( loc );
	if ( it == muts.end() )
		 PRINT2( "loc not found", loc );
  assert(  it != muts.end() );
  assert( it->second.loc == loc && it->first == loc );
  return it->second;
}
#endif



// Method: loadFromMs
// Load a Mutlist from the output of the ms simulator
MutlistP Mutlist::loadFromMs( istream& is ) {

	// ms output for one simulation looks like:

	//
	// //
	// segsites: 4
	// positions: 0.0110 0.0765 0.6557 0.7571
	// 0010
	// 0100
	// 0000
	// 1001
	//

	// (ms output for a set of sims also includes a header, followed by a block for each sim;
	// here we only parse one block for one sim.)

	MutlistP mutlist = boost::make_shared<Mutlist>();

	boost::io::ios_exception_saver save_is( is );
	is.exceptions( std::ios_base::eofbit | std::ios_base::failbit | std::ios_base::badbit );

	using std::getline;
	string line;
	do { getline( is, line ); } while( !util::startsWith( line, "//" ) );

	string word;
	is >> word;
	chkCond( word == "segsites:", "missing segsites line in ms format" );
	unsigned int nmuts;
	is >> nmuts;

	is >> word;
	chkCond( word == "positions:", "missing positions line in ms format" );
	vector< ploc_t > mutLocs;
	for ( unsigned int i = 0; i < nmuts; i++ ) {
		ploc_t ploc;
		is >> ploc;
		mutLocs.push_back( ploc );
	}
	chkCond( mutLocs.size() == nmuts, "ms format error: number of mutations does not equal segsites" );

	vector< string > chroms;
	while ( True ) {
		line.clear();
		if ( is ) {
			try { getline( is, line ); }
			catch( const istream::failure& f ) { }
			bool_t okLine = !util::isSpace( line );
			if ( okLine ) {
				chkCond( line.size() == nmuts, "Mutlist::loadFromMs: mismatch between nsites and chromosome line length" );
				chroms.push_back( line );
			} else if ( !chroms.empty() ) break;
		} else break;
	}

	chkCond( !chroms.empty(), "ms format error: no chroms" );

	leafset_set_max_leaf_id( static_cast< leaf_id_t >( chroms.size() ) );

	for ( unsigned int mutId = 0; mutId < nmuts; mutId++ ) {
		leafset_p leafset( LEAFSET_NULL );
		leaf_id_t leafId( 0 );
		ForEach( string chrom, chroms ) {
			if ( chrom[ mutId ] == '1' ) {
				leafset = leafset_union( leafset, make_singleton_leafset( leafId ) );
			}
			leafId++;
		}
		chkCond( !leafset_is_empty( leafset ), "ms format error: a mutation has zero derived allele freq" );
		mutlist->addMut( loc_t( mutLocs[ mutId ] ), leafset, genid( 0.0 ), popid( 1 ) );
	}  // for each mut

	mutlist->freeze( /* inf_sites= */ True, /* regionLen= */ 0 /* unused for inf_sites=True */ );

	return mutlist;
}

// Method: print_haps_ms
//
// Write the output of one simulation in ms format to the specified stream.
// Note that the standard format also includes one header describing all the simulation replicas:
// the command line of the simulator, the random seeds, the total number of replicas,
// and the number of samples in each replica.  Those are written in <CoSiMain::cosi_main()>;
// print_haps_ms only writes a single simulation (replica).
//
// Params:
//
//   strm - the stream to which the ms-format simulation is written
//   sampleSizes - for each sampled population, the number of sampled haplotypes in it
//   treeStatsHook - the hook object used to gather requested tree statistics; may be NULL
//   outputMutGens - whether to write out the time of each mutation
//   recombLocs - if non-NULL, vector (unsorted) of recombination locations
//   outputPrecision - number of decimal places in the output
//
void Mutlist::print_haps_ms( ostream& strm, const vector< nchroms_t >& sampleSizes,
														 TreeStatsHookP treeStatsHook, bool_t outputMutGens,
														 const vector< loc_t > *recombLocs,
														 bool_t outputMutGlocs,
														 GenMapP genMap,
														 int outputPrecision,
#ifndef COSI_NO_CPU_TIMER
														 boost::timer::cpu_timer
#else														 
														 void
#endif														 
														 *cpuTimer,
														 const genid *endGen ) const {

	boost::io::ios_precision_saver strm_precision_saver( strm );
	strm.precision( outputPrecision );

	const Mutlist *mutlist = this;

  size_t nmuts = mutlist->size();

	if ( treeStatsHook.get() ) {
		treeStatsHook->printStats();
	}
#ifndef COSI_NO_CPU_TIMER	
	if ( cpuTimer ) {
		boost::timer::cpu_times elapsed = cpuTimer->elapsed();
		strm << "stat time " << (static_cast<double>( elapsed.user + elapsed.system ) / 1e8 ) << "\n";
	}
#endif	
	if ( endGen ) strm << "stat endGen " << *endGen << "\n";
	if ( outputMutGlocs )
		 strm << "region_len_cM: " << ( genMap->getRegionRecombRateAbs() * 100. ) << "\n";
			 
	strm << "segsites: " << nmuts << "\n";
	
	boost::scoped_array<char> line( new char[ nmuts+10 ] );
	memset( line.get(), 0, nmuts+3 );
	line.get()[ nmuts ] = '\n';

	if ( !mutlist->getMuts().empty() ) {

		strm << "positions:";
		ForEach( const Mut& m, mutlist->getMuts() ) strm << " " << get_loc( m.loc );
		strm << "\n";

		if ( outputMutGlocs ) {
			strm << "positions_genMap:";
			ForEach( const Mut& m, mutlist->getMuts() ) strm << " " << genMap->getGdPos( m.loc );
			strm << "\n";
		}

		if ( outputMutGens ) {
			strm << "muttimes:";
			ForEach( const Mut& m, mutlist->getMuts() ) strm << " " << m.gen;
			strm << "\n";
		}

		if ( recombLocs ) {
			vector< loc_t > recombLocsSorted( *recombLocs );
			std::sort( recombLocsSorted.begin(), recombLocsSorted.end() );
			strm << "recomblocs:";
			ForEach( const loc_t& recombLoc, recombLocsSorted )
				 strm << " " << get_ploc( recombLoc );
			strm << "\n";
		}
		//PRINT( "wrote header" );

		leaf_id_t leaf( 0 );
		for (size_t ipop = 0; ipop < sampleSizes.size(); ipop++) {
			if (sampleSizes[ipop] > nchroms_t(0) ) {
				leaf_id_t popEndLeaf = leaf + sampleSizes[ipop];

				for (; leaf < popEndLeaf; leaf++) {
					memset( line.get(), '0', nmuts );

					const vector< Mutlist::const_iterator >& leafMuts = mutlist->getLeafMuts( leaf );
					ForEach( Mutlist::const_iterator m, leafMuts ) {
					  chkCond( 0 <= m->mutId );
					  chkCond( m->mutId < int(nmuts) );
					
					  line.get()[ m->mutId ] = '1';
					}
					strm << line.get();
				}
			}
		}  // for each pop

		//PRINT( "wrote haps" );

		
	}
	strm << "\n";
}  // print_haps_ms



// End class impl: Mutlist 


//
// Class impl: MutProcessor_AddToMutlist
//

MutProcessor_AddToMutlist::~MutProcessor_AddToMutlist() {}
void MutProcessor_AddToMutlist::processMut(loc_t loc, leafset_p leaves, genid gen, popid popName) {
	mutlist->addMut( loc, leaves, gen, popName );
}
void MutProcessor_AddToMutlist::postprocess() {
	
}

// Class impl: MutProcessor_AddToMutlist_WithAscertainment
MutProcessor_AddToMutlist_WithAscertainment::~MutProcessor_AddToMutlist_WithAscertainment() { }

void MutProcessor_AddToMutlist_WithAscertainment::processMut(loc_t loc, leafset_p leaves, genid gen, popid popName) {
	if ( ( leafset_size( leaves ) != nchroms_t(1) ) || ( random_double() > dropSingletonsFrac ) )
		 PARENT::processMut( loc, leaves, gen, popName );
}


}  // namespace cosi
