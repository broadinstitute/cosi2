//
// Header: ehh.h
//
// Code to compute EHH and EHH-based statistics (iHS, deltaIHH, XP-EHH).
// See http://www.sciencemag.org/content/327/5967/883.abstract
//

#ifndef __INCLUDE_COSI_EHH_H
#define __INCLUDE_COSI_EHH_H

#include <cstdlib>
#include <vector>
#include <boost/range/value_type.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <cosi/defs.h>
#include <cosi/leafset.h>
#include <cosi/mutlist.h>

namespace cosi {

using namespace std;

namespace ehh {
// Type: ehh_t
// Represents an EHH value: of pairs of chromosomes
// with given allele of a given SNP, what fraction are identical
// out from the given SNP to a given genetic distance?
typedef frac_t ehh_t;
  
typedef double ihh_t;

struct iHSTestConfig {
	 ehh_t cutoffEHH;
	 ehh_t startAllEHH;
	 freq_t minDerivedFreq;

	 iHSTestConfig():
		 cutoffEHH( 0.05 ), startAllEHH( 1.00 ), minDerivedFreq( 0.0 ) { }
};

// FuncP: computeEHH
// Compute EHH for the given set of leaves, starting from the given core and going
// in a given direction.
vector<ehh_t> computeEHH( const vector<Mut>& mutlist,
													leafset_p coreLeafset,
													ehh_t ehhCutoff = 0.05 );


template <typename T> inline
T npairs_twice( T n ) { return n * ( n-1 ); }

// Func: computeEHH
// Compute EHH for the given set of leaves, starting from the given core and going
// in a given direction.
template <typename LeafsetRangeType>
vector<ehh_t> computeEHH( const LeafsetRangeType& leafsets,
													typename boost::range_value<LeafsetRangeType>::type coreLeafset,
													ehh_t ehhCutoff ) {
  using boost::adaptors::transformed;
  using boost::adaptors::filtered;
	
	typedef typename boost::range_value<LeafsetRangeType>::type leafset_type;

	vector<ehh_t> result;
	//result.reserve( boost::size( leafsets ) );
	result.push_back( 1.0 );
	vector< leafset_type  > classes;
	classes.push_back( coreLeafset );
	nchroms_t coreLeafsetCardinality = leafset_size( coreLeafset );
	//PRINT( leafset_size( coreLeafset ) );
	if ( coreLeafsetCardinality >= 2 ) {

		//size_t mutNum = 0;
		BOOST_FOREACH( leafset_type leafset, leafsets ) {
			vector< leafset_type  > newClasses;
			
			BOOST_FOREACH( leafset_type  cls, classes ) {
				leafset_p inters[2] = { leafset_intersection( cls, leafset ),
																leafset_difference( cls, leafset ) };
				if ( leafset_is_empty( inters[0] ) || leafset_is_empty( inters[1] ) ) {
					newClasses.push_back( cls );
				} else {
					for ( int i = 0; i < 2; i++ )
						 if ( leafset_size( inters[i] ) > 1 ) newClasses.push_back( inters[i] );
				}
			}
			if ( newClasses.empty() ) break;
			classes = newClasses;
			
			unsigned ehhSum = boost::accumulate( classes | transformed( leafset_size ) | transformed( npairs_twice<nchroms_t> ),
																					 0 );
			
      ehh_t ehhHere = util::getFrac( ehhSum, npairs_twice( coreLeafsetCardinality ) );
			if ( ehhHere < ehhCutoff ) break;
			result.push_back( ehhHere );
		}
		
		if ( result.size() > 1 ) result.pop_back();
	}
	//PRINT( result.size() );
	return result;
}  // computeEHH

#if 0
// Func: computeIHH
// Compute IHH for the given set of leaves, starting from the given core and going
// in a given direction.
template <typename LeafsetRangeType>
ihh_t computeIHH( const LeafsetRangeType& leafsets,
									//loc_t coreSNP, dir_t dir,
								  typename boost::range_value<LeafsetRangeType>::type coreLeafset,
									typename boost::range_value<LeafsetRangeType>::type popLeafset,
									ehh_t ehhCutoff
	) {
	namespace lam = ::boost::lambda;
	namespace ad = ::boost::adaptors;
	namespace ran = ::boost::range;
	using lam::bind;
	using lam::_1;
	
	//vector<Mut> muts = mutlist->getMutsInDir( coreSNP, dir );
	// vector<Mut> muts( muts_ );
	// ran::remove_erase_if( muts,  );
	// vector<Mut> sideMuts( muts.begin()+1, muts.end() );
	
	double oneOverN = 1.0 / static_cast<double>( leafset_size( popLeafset ) );
	vector<ehh_t> ehh = computeEHH( leafsets | filtered( ! bind( &Mut::isMonomorphicIn, _1, popLeafset ) ), coreLeafset, ehhCutoff );

	vector<glen_t> gdists;
	BOOST_FOREACH( const Mut& m, muts )
		 gdists.push_back( get_gloc( m.loc ) - gloc_t( 0.0 ) );
	// push_back( gdists, muts | sliced( 0, ehh.size() ) | transformed( lam::bind( &Mut::loc, lam::_1 ) )
	// 					 | transformed( get_gloc ) );

	PRINT( ehh.size() );

	std::for_each( ehh.begin()+1, ehh.end(), lam::_1 += oneOverN );

	double cutoffPt = math::findWhere( ehh, ehhCutoff );

	if ( cutoffPt < 0 ) cutoffPt = ehh.size()-1;

	ran::for_each( ehh, lam::_1 -= ehhCutoff );
	return static_cast<ihh_t>( ToDouble(cosi_fabs( math::integrate( ehh, gdists, 0.0, cutoffPt ) ) ) );
//#endif	
}  // computeIHH()


ihh_t computeIHH( const vector<Mut>& muts,
									//loc_t coreSNP, dir_t dir,
									leafset_p coreLeafset,
									leafset_p popLeafset,
									ehh_t ehhCutoff = 0.05 );

#endif
}  // namespace ehh

}  // namespace cosi

#endif
// #ifndef __INCLUDE_EHH_H
