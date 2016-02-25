//
// * Computation of additional stats
//
// A framework for computing population genetics statistics for each simulation,
// directly within the simulator, rather than piping simulation output to
// a separate program.

// ** details

#ifndef COSI_INCLUDE_MORESTATS_H
#define COSI_INCLUDE_MORESTATS_H

#include <vector>
#include <string>
#include <utility>
#include <boost/mpl/if.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/mutlist.h>
#include <cosi/genmap.h>

namespace cosi {

namespace morestats {

// * Class SimData - data from one simulation
class SimData {
	 
};  // class SimData

namespace acc = boost::accumulators;

class CustomStats {
	 typedef acc::accumulator_set<double, acc::stats< acc::tag::sum_kahan, acc::tag::mean,
																										acc::tag::variance > > acc_t;

	 DemographyP dem;
	 GenMapP genMap;

	 std::vector< acc_t > fst;

	 void processSNP( leafset_p leafset, loc_t loc ) {

		 if ( !leafset_is_empty( leafset ) ) {
			 std::stack<const leafset_t *> s;
			 s.push( leafset.get() );
			 while( !s.empty() ) {
				 const leafset_t *p = s.top();
				 s.pop();
				 if ( p->leafId != NULL_LEAF_ID ) {
					 p->leafId;
				 } else {
					 s.push( p->childA.get() );
					 s.push( p->childB.get() );
				 }
			 }
		 }
	 }


		 
	 }
	 
	 
	 
};  // class CustomStats


// * Class Stats - abstract base class for stats computers; each instance of a stats computer computes
//      one ore more stats for each simulation.
class Stats {

protected:
	 
// ** Method do_getStatNames - returns the names of the relevant stats
	 virtual void do_getStatNames( std::vector< std::string >& ) = 0;

	 
// ** Method do_computeStats - returns the values of the relevant stats.
// Input params:
//    - muts :: the list of mutations, sorted by location.
	 virtual void do_computeStats( std::vector< Mut > const& muts, std::vector< double >& stats ) = 0;
};  // class Stat

// * Class LD - computes LD for pairs of SNPs
//
// Template params:
//     - use_plen :: if true, measure distance between SNPs in physical distance; if false, in genetic distance.
template <bool use_plen>
class LD: public Stats {
// ** private
	 GenMapP genMap;

	 
// *** Type: sep_len_t - the type for the separation between SNPs
	 typedef typename boost::mpl::if_c< use_plen, plen_t, glen_t >::type sep_len_t;
// *** Type: sep_loc_t - the type for the location of a SNP
	 typedef typename boost::mpl::if_c< use_plen, ploc_t, gloc_t >::type sep_loc_t;

	 
// *** Method: get_mut_loc - given a mutation's physical location, convert it to sep_loc_t
	 sep_loc_t
	 get_mut_loc( typename boost::enable_if_c< use_plen, loc_t >::type m ) { return get_ploc( m ); }

	 sep_loc_t
	 get_mut_loc( typename boost::disable_if_c< use_plen, loc_t >::type m ) { return genMap->getGdPos( m ); }
	 
protected:
// ** protected
	 
	 virtual void do_getStatNames( std::vector< std::string >& statNames ) {
		 statNames.push_back( std::string( "Dprime" ) );
		 statNames.push_back( std::string( "r2" ) );
	 }
	 virtual void do_computeStats( std::vector< Mut > const& muts, std::vector< double >& stats ) {
		 typedef typename std::vector< Mut >::iterator it_t;
		 it_t b = muts.begin();
		 it_t e = b;
		 while ( true ) {
			 bool foundSome = false;
			 while ( e != muts.end() && get_mut_loc( e ) - get_mut_loc( b ) < sepRange.first ) ++e;

			 for ( it_t t = e; t != muts.end() && get_mut_loc( t ) - get_mut_loc( b ) < sepRange.second; ++t ) {
				 foundSome = true;
			 }

			 if ( !foundSome ) break;
		 }
	 }

private:
	 std::pair< sep_len_t, sep_len_t > sepRange;
};


}  // namespace morestats


}  // namespace cosi


#endif  // #ifndef COSI_INCLUDE_MORESTATS_H
