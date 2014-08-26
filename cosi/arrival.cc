
#include <map>
#include <utility>
#include <iostream>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>

#include <cosi/utils.h>
#include <cosi/arrival.h>

namespace cosi {

void testArrival() {
	using std::map;
	using std::make_pair;
	
	RandGenP randGen = boost::make_shared<RandGen>();

	map<genid, rate_t> m;
	m.insert( make_pair( genid(0), rate_t(.001) ) );
	m.insert( make_pair( genid(4999), rate_t(.001) ) );
	m.insert( make_pair( genid(5000), rate_t(.002) ) );
	m.insert( make_pair( genid(100000), rate_t(.002) ) );

	InhomogeneousPoissonProcess ipp( m );

	namespace acc = boost::accumulators;
	typedef acc::accumulator_set<double, acc::stats< acc::tag::sum_kahan, acc::tag::mean, acc::tag::variance > > acc_t;

	acc_t counts;
	acc_t wtimes;

	for ( int k = 0; k < 1000000; k++ ) {
		genid g(0);
		int cnt = 0;
		while ( g < genid(10000) ) {
			genid last_g = g;
			g = ipp.getTimeOfNextEvent( randGen, g, factor_t(1.0) );
			wtimes( ToDouble( g - last_g ) );
			if ( genid(5000) < g  &&  g < genid( 10000 ) ) {
				//PRINT2( i, g );
				cnt ++;
			}
		}
		counts( double(cnt ) );
	}
	PRINT2( acc::mean( counts ), acc::variance( counts ) );
	PRINT2( acc::mean( wtimes ), sqrt( acc::variance( wtimes ) ) );
	
}

}  // namespace cosi
