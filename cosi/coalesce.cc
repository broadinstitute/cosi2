/* $Id: coalesce.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#include <algorithm>
#include <map>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <cosi/defs.h>
#include <cosi/coalesce.h>
#include <cosi/demography.h>
#include <cosi/node.h>

namespace cosi {

extern double poisPrec;
extern unsigned int poisMaxSteps;

namespace coal {

extern size_t lastIntersCount;

Coalesce::Coalesce( DemographyP demography_ ):
	demography( demography_ ),
	lastrate(0.0), nonhomogPopIdx( -1 ) {
}

	
/* 
 * Calculates the coalesce rate according to the number of
 * nodes in each population. 
 */
double 
Coalesce::coalesce_get_rate (void) const
{
	int numpops = demography->dg_get_num_pops();
	int i;
	double rate = 0;
	nchroms_t numnodes;
	nchroms_t popsize;

	popRates.clear();
	
	for (i = 0; i < numpops; i++) {
		popRates.push_back( 0 );
		Pop *popptr = demography->dg_get_pop_by_index (i);
		if ( popptr->getCoalArrivalProcess() ) continue;

		numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
		popsize = demography->dg_get_pop_size_by_index (i);
		if (numnodes > 1  /*&& popsize > 0*/ ) {
			prob_t coalRate = util::getFrac( numnodes * (numnodes - 1),
																			 4 * std::max( popsize, 1 ) );
#ifdef COSI_SUPPORT_COALAPX

			if ( popptr->restrictingCoalescence() )
				 coalRate = util::getFrac( popptr->getNumCoalesceableChromPairs(), 2.0 * std::max( popsize, 1 ) );
			
#endif  // #ifdef COSI_SUPPORT_COALAPX			
			
#ifdef COSI_DEV			
			demography->dg_get_pop_by_index( i )->setCoalRate( coalRate );
#endif
			popRates.back() = coalRate;
			rate += coalRate;
		}
	}
	
	lastrate = rate;
	return rate;
}

pop_idx_t
Coalesce::coalesce_pick_popindex () const
{
	return dist( *this->getRandGen(), boost::random::discrete_distribution<>::param_type( popRates ) );
}  // Coalesce::coalesce_pick_popindex()

gens_t Coalesce::coalesce_get_wait_time_nonhomog( genid gen, gens_t maxWaitTime ) const {
	this->nonhomogPopIdx = -1;
	int numpops = demography->dg_get_num_pops();
	
	for (pop_idx_t i = 0; i < numpops; i++) {
		Pop *popptr = demography->dg_get_pop_by_index (i);
		Pop::coal_arrival_process_type_ptr proc = popptr->getCoalArrivalProcess();
		// PRINT5( "getTimeNonhomog", i, numpops, demography->dg_get_num_nodes_in_pop_by_index (i),
		// 				popptr->restrictingCoalescence() );
		if ( proc ) {
			nchroms_t numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);

			
			nchromPairs_t npairs = numnodes * (numnodes-1) / 2;

#ifdef COSI_SUPPORT_COALAPX

			if ( popptr->restrictingCoalescence() )
				 npairs = popptr->getNumCoalesceableChromPairs();
			
#endif  // #ifdef COSI_SUPPORT_COALAPX
			if ( npairs > 0 ) {

			
				genid nextCoalTime = proc.nextArrivalTime( gen, gen + maxWaitTime,
																									 static_cast<double>( npairs ),
																									 *getRandGen(),
																									 /* eps= */ poisPrec,
																									 /* maxSteps= */ poisMaxSteps );
				//PRINT7( "waitNonhomog", i, gen, gen+maxWaitTime, maxWaitTime, npairs, nextCoalTime );
			
				if ( nextCoalTime < gen + maxWaitTime ) {
					maxWaitTime = nextCoalTime - gen;
					this->nonhomogPopIdx = i;
				}
			}
		}
	}
	//PRINT4( "returning", maxWaitTime, nonhomogPopIdx, gen );
	return this->nonhomogPopIdx == -1 ? gens_t( std::numeric_limits<double>::infinity() ) : maxWaitTime;
}
pop_idx_t Coalesce::coalesce_pick_popindex_nonhomog() const {
	return this->nonhomogPopIdx;
}

}  // namespace coal

}  // namespace cosi
