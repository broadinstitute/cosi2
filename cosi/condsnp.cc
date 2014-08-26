
#include <algorithm>
#include <iterator>
#include <iostream>
#include <boost/phoenix/core/value.hpp>
#include <boost/phoenix/core/argument.hpp>
#include <boost/phoenix/bind.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <cosi/utils.h>
#include <cosi/generalmath.h>
#include <cosi/genmap.h>
#include <cosi/mutlist.h>
#include <cosi/demography.h>
#include <cosi/ehh.h>
#include <cosi/condsnp.h>

#if 0
namespace std {

template <>
class numeric_limits< cosi::gloc_t > {
	 typedef cosi::gloc_t T;
	 typedef cosi::cosi_double ValT;
	 typedef numeric_limits<ValT> NL;
public:
	 static const bool is_specialized = true;
	 static T min() throw() { return T( NL::min() ); }
	 static T max() throw() { return T( NL::max() ); }
	 static const int  digits = NL::digits;
	 static const int  digits10 = NL::digits10;
	 static const bool is_signed = NL::is_signed;
	 static const bool is_integer = NL::is_integer;
	 static const bool is_exact = NL::is_exact;
	 static const int radix = NL::radix;
	 static T epsilon() throw() { return T( NL::epsilon() ); }
	 static T round_error() throw() { return T( NL::round_error() ); }

	 static const int  min_exponent = NL::min_exponent;
	 static const int  min_exponent10 = NL::min_exponent10;
	 static const int  max_exponent = NL::max_exponent;
	 static const int  max_exponent10 = NL::max_exponent10;

	 static const bool has_infinity = NL::has_infinity;
	 static const bool has_quiet_NaN = NL::has_quiet_NaN;
	 static const bool has_signaling_NaN = NL::has_signaling_NaN;
	 static const float_denorm_style has_denorm = NL::has_denorm;
	 static const bool has_denorm_loss = NL::has_denorm_loss;
	 
	 static T infinity() throw() { return T( NL::infinity() ); }
	 static T quiet_NaN() throw() { return T( NL::quiet_NaN() ); }
	 static T signaling_NaN() throw() { return T( NL::signaling_NaN() ); }
	 static T denorm_min() throw() { return T( NL::denorm_min() ); }

	 static const bool is_iec559 = NL::is_iec559; 
	 static const bool is_bounded = NL::is_bounded;
	 static const bool is_modulo = NL::is_modulo;

	 static const bool traps = NL::traps;
	 static const bool tinyness_before = NL::tinyness_before;
	 static const float_round_style round_style = NL::round_style;
};

}  // namespace std
#endif

#undef NDEBUG

namespace cosi {

namespace math {

//
// Class: StepwiseFn
//
// An stepwise function: stores a vector of (x,f(x)) pairs and
// can compute the values of f at new points by simple linear interpolation.
//
template <typename TDomain = double, typename TRange = double>
class StepwiseFn /*: public AbstractFunction<TDomain, TRange>*/ {
public:
	 // Method: clear
	 // Remove all (x,f(x)) pairs from this function.
	 void clear() { x2y.clear(); }

	 // Method: addPt
	 // Store an (x,f(x)) pair specifying the value of the function at a given point 'x'.
	 // Points must be added in order of increasing 'x'.
	 void addPt( TDomain x, TRange y ) {
		 assert( x2y.empty() || x > x2y.rbegin()->first );
		 if ( !x2y.empty() && x2y.rbegin()->second == y ) return;
		 x2y.insert( x2y.end(), make_pair( x, y ) );
	 }

	 // Method: eval
	 // Evaluate the function at the specified point 'x', interpolating as necessary
	 // between the neighboring values.
	 TRange operator()( TDomain x ) const {
		 assert( !x2y.empty() );
		 assert( x >= x2y.begin()->first );
		 assert( x <= x2y.rbegin()->first );
		 BOOST_AUTO_TPL(it, x2y.lower_bound( x ) );
		 assert( it != x2y.end() );
		 if ( it->first == x ) return it->second;
		 assert( it != x2y.begin() );
		 BOOST_AUTO_TPL( it_p, boost::prior( it ) );
		 assert( it_p->first < x );
		 assert( it->first > x );
		 return it_p->second; 
	 }

	 virtual TRange eval( TDomain x ) const { return (*this)( x ); }

	 size_t getNumPts() const { return x2y.size(); }
	 
	 void save( filename_t filename ) const {
		 std::ofstream out( filename.c_str() );
		 out.exceptions( std::ios_base::failbit | std::ios_base::badbit );
		 out.precision( 20 );
		 out << "x\ty\n";
		 for ( x2y_const_iterator it = x2y.begin(); it != x2y.end(); it++ )
				out << it->first << "\t" << it->second << "\n";
	 }
	 
	 void load( filename_t filename ) {
		 std::ifstream inFile( filename.c_str() );
		 inFile.exceptions( std::ios_base::failbit | std::ios_base::badbit );
		 clear();
		 size_t lineNum = 0;
		 while ( inFile ) {
			 std::string line;
			 std::getline( inFile, line );
			 if ( lineNum == 0 ) continue;
			 std::istringstream lineStr( line );
			 lineStr.exceptions( std::ios_base::failbit | std::ios_base::badbit );
			 TDomain x;
			 TRange y;
			 lineStr >> x >> y;
			 addPt( x, y );
			 lineNum++;
		 }
	 }
	 
	 typedef std::map<TDomain,TRange> x2y_type;

	 x2y_type getMap() const { return x2y; }

private:

	 typedef typename x2y_type::const_iterator x2y_const_iterator;
	 
	 // Field: x2y
	 // Map from point x to the value f(x), for the points at which the function value
	 // is explicitly specified.
	 x2y_type x2y;

	friend class boost::serialization::access;
	 template <class Archive> void serialize( Archive& ar, const unsigned int /* version */ ) {
		 ar & BOOST_SERIALIZATION_NVP( x2y );
	 }
	 
};  // class StepwiseFn



}  // namespace math


namespace {
static gens_t treeTimeTot( 0.0 ) ;
static gens_t treeTimeMatching( 0.0 );
static nsims_t simsTot( 0 );
static nsims_t simsMatching( 0 );
static frac_t fracMatching( 0.0 );

static nsims_t simsAbove( 0 );

struct CondStatsPrinter {
	 CondStatsPrinter() { }
	 ~CondStatsPrinter() {
		 // PRINT5( fracMatching, simsTot, simsMatching, ((double)simsMatching) / ((double)simsTot),
		 // 				 treeTimeMatching / treeTimeTot );
		 using util::getFrac;
		 PRINT3( simsAbove, simsTot, getFrac( simsAbove, simsTot ) );
	 }
} instance;

}

bool CondSnpDef::sampleFreqsMatch( const std::map< popid, nchroms_t >& pop2Counts ) const {
// #ifdef COSI_DEV_PRINT	
// 	PRINT( pop2Counts.size() );
// 	for( std::map<popid,nchroms_t>::const_iterator it2 = pop2Counts.begin(); it2 != pop2Counts.end(); it2++ ){
// 		PRINT2( it2->first, it2->second );
// 	}
// #endif	
			 
	for ( pop2cond_map_t::const_iterator it = pop2cond.begin(); it != pop2cond.end(); it++ ) {
//		using namespace util;
//		PRINT( it->first );
		if ( !( it->second( util::map_get( pop2Counts, it->first, 0 ) ) ) ) return false;
	}
	return true;
}


// Virtual method: handle_add_edge
//
// Called after adding an edge to the ARG.  The new edge may be the result of a coalescence,
// a recombination or a gene conversion.
//
// Params:
//
//   nodeId_moreRecent - <nodeid> of the more recent (lower gen) node of the edge
//   nodeId_lessRecent - <nodeid> of the less recent (higher gen) node of the edge.
//   genId_moreRecent - <genid> of the more recent (lower gen) node of the edge
//   genId_lessRecent - <genid> of the less recent (higher gen) node of the edge.
//   seglist - the segments inherited along the edge.  NOTE: this seglist may be destroyed
//       after the call, so make a copy if you need to save it.
//       
void CondSnpMgr::handle_add_edge( nodeid /* nodeId_moreRecent */,
																	nodeid /*nodeId_lessRecent */,
																	genid genId_moreRecent,
																	genid genId_lessRecent, const Seglist *seglist,
																	edge_kind_t /*edgeKind*/ ) {
	const seglist::Seg *seg = seglist_find_seg( seglist, condSnpDef.getLoc() );
	if ( seg ) {
		assert( seg && seg->contains( condSnpDef.getLoc( ) ) );
		gens_t edgeLen = genId_lessRecent - genId_moreRecent; 
		condSnpTreeTimeTot += edgeLen;
		std::map< popid, nchroms_t > popCounts = demography->getLeafsetPopCounts( seg->getLeafset() );
		if ( condSnpDef.sampleFreqsMatch( popCounts ) ) {
			 condSnpTreeTimeMatching += edgeLen;
			 savedLeafsets.push_back( seg->leafset );
			 savedEdgesCumLens.push_back( edgeLen + ( savedEdgesCumLens.empty() ? gens_t(0.0) :
																								savedEdgesCumLens.back() ) );
		} else nwrongEdges++;
	}
}

CondSnpMgr::CondSnpMgr( DemographyP demography_, CondSnpDef condSnpDef_):
	demography( demography_ ), condSnpDef( condSnpDef_ ),
	condSnpTreeTimeTot( 0.0 ), condSnpTreeTimeMatching( 0.0 ), nwrongEdges( 0 ) {
}


CondSnpMgr::~CondSnpMgr() {
	
}

void CondSnpMgr::printResults( MutlistP mutlist, GenMapP genMap ) const {
//	std::cerr << "cond time tot: " << condSnpTreeTimeTot << " matching: " << condSnpTreeTimeMatching <<
//		 "\n";
	using std::vector;
	using std::map;
	using std::back_inserter;
	using std::cerr;
	using std::ifstream;
	using std::ofstream;
	using std::ios_base;
	using boost::adaptors::transformed;
	using boost::adaptors::filtered;
	using boost::phoenix::bind;
	using boost::phoenix::arg_names::arg1;
	using boost::phoenix::val;
	using ehh::ehh_t;
	using ehh::computeEHH;
	using math::StepwiseFn;

	static StepwiseFn<glen_t, ehh_t> loaded;

	static map<glen_t, ehh_t> loadedMap;
	static glen_t maxDist;
	if ( loadedMap.empty() ) {
		std::string loadFrom( getenv( "LOADEHH" ) ? getenv( "LOADEHH" ) : "ehhfn.xml" );
		ifstream ifs( loadFrom.c_str() );
		ifs.exceptions( ios_base::failbit | ios_base::badbit );
		boost::archive::xml_iarchive ia( ifs );
		ia >> BOOST_SERIALIZATION_NVP( loaded );
		loadedMap = loaded.getMap();
		maxDist = loadedMap.rbegin()->first;
		PRINT3( loadFrom, loadedMap.size(), maxDist );
	}
	
	const vector<Mut>& muts = mutlist->getMuts();

	vector<Mut> causalMut;
	boost::range::copy( muts | filtered( bind( &Mut::loc, arg1 ) == val( condSnpDef.getLoc() ) ), back_inserter( causalMut ) );
	loc_t causalMutLoc( condSnpDef.getLoc() );
	leafset_p causalMutLeafset;
	if ( !causalMut.empty() ) causalMutLeafset = this->sel_leaves;
	else {

		if ( !savedLeafsets.empty() ) {
			gens_t breakpt( demography->getRandGen()->random_double() * ToDouble( savedEdgesCumLens.back() ) );
			PRINT4( savedLeafsets.empty() ? gens_t(0) : savedEdgesCumLens.back(), breakpt, savedLeafsets.size(),
							nwrongEdges );
			for ( size_t i = 0; i < savedLeafsets.size() && !causalMutLeafset; i++ )
				 if ( breakpt < savedEdgesCumLens[ i ] )
						causalMutLeafset = savedLeafsets[ i ];
		}

		
	}



	if ( !causalMutLeafset ) { PRINT4( "----NO-LEAFSET----", savedLeafsets.empty(), nwrongEdges, simsTot ); }
	
	if ( causalMutLeafset ) {
		std::map< popid, nchroms_t > popCnts( demography->getLeafsetPopCounts( causalMutLeafset ) );
		using util::map_get;
		PRINT3( map_get( popCnts, popid(1), 0 ), map_get( popCnts, popid(4), 0 ), map_get( popCnts, popid(5), 0 ) );

		vector<Mut> mutsR;
		boost::range::copy( muts | filtered( bind( &Mut::loc, arg1 ) > val( condSnpDef.getLoc() ) ), back_inserter( mutsR ) );

		vector<ehh_t> ehhVals = computeEHH( mutsR | transformed( bind( &Mut::leaves, arg1 ) ),
																				causalMutLeafset,
																				/* ehhCutoff= */ 0.20 );

		gloc_t causalGloc = genMap->getGdPos( causalMutLoc );
		StepwiseFn<glen_t, ehh_t> ehhFn;
		ehhFn.addPt( glen_t( 0.0 ), 1.0 );
		//PRINT3( mutsR.size(), ehhVals.size(), causalGloc );
		bool isLarger = true;
		for ( size_t i = 0; i < mutsR.size() && i+1 < ehhVals.size(); i++ ) {
			gloc_t mut_gloc = genMap->getGdPos( mutsR[ i ].loc );
			glen_t mut_gdist = mut_gloc - causalGloc;

			if ( mut_gdist >= maxDist ) break;
			ehh_t loaded_ehh_here = loaded( mut_gdist );
			if ( loaded_ehh_here < .1 ) break;
				 
			if ( isLarger && ehhVals[ i+1 ] < loaded( mut_gdist ) ) {
				PRINT9( simsAbove, simsTot, i, mutsR[i].loc, mut_gloc, causalGloc, mut_gdist, ehhVals[ i+1 ], loaded_ehh_here );
				isLarger = false;
				break;
			}
			
			ehhFn.addPt(  mut_gdist, ehhVals[ i+1 ] );
		}
		if ( isLarger ) simsAbove++;
		if ( isLarger ) { PRINT5( isLarger, loaded.getNumPts(), ehhFn.getNumPts(), simsAbove, simsTot ); }
		if ( isLarger ) { PRINT( "*****************" ); }

		if ( false ) {
			ostringstream fname;
			fname << "saveehh" << simsTot << ".xml";
			std::string fnameStr( fname.str() );
			ofstream ofs( fnameStr.c_str() );
			ofs.exceptions( ios_base::failbit | ios_base::badbit );
			boost::archive::xml_oarchive oa( ofs );
			oa << BOOST_SERIALIZATION_NVP( ehhFn );
		}
		
		// std::cerr << "ehhVals: mutsR.size=" << mutsR.size() << " ehhVals.size=" << ehhVals.size() << " vals=";
		// std::copy( ehhVals.begin(), ehhVals.end(), std::ostream_iterator<ehh_t>( std::cerr, ", " ) );
		// std::cerr << "\n";
	}

	//PRINT( ehhVals );
	
	simsTot++;
	treeTimeTot += this->condSnpTreeTimeTot;
	treeTimeMatching += this->condSnpTreeTimeMatching;
	if ( condSnpTreeTimeMatching > gens_t(0.0) ) {
		simsMatching++;
		fracMatching += ( condSnpTreeTimeMatching / condSnpTreeTimeTot );
	}
}


}  // namespace cosi
