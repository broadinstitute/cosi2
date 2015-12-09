/* $Id: simulator.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cosi/historical.h>
#include <cosi/migrate.h>
#include <cosi/recomb.h>
#include <cosi/geneconversion.h>
#include <cosi/demography.h>
#include <cosi/coalesce.h>
#include <cosi/simulator.h>
#include <cosi/leafset.h>
#include <cosi/genmap.h>
#include <cosi/mutlist.h>
#include <cosi/mutate.h>
#include <cosi/sweep3.h>

namespace cosi {

Simulator::Simulator( DemographyP demography_, GenMapP genMap_ ):
	demography( demography_ ), genMap( genMap_ ),
	coalesce_rate(0.0),
	migrate_rate(0.0),
	recombination_rate(0.0),
	geneconv_rate(0.0),
	poisson_rate(0.0),
	use_what( use_none )
{
}

/* To execute the coalescent simulator:
 * 1. Calculate the times to the next historical and poisson events.
 * 2. Choose the closest event, and execute it.
 * 3. Repeat until sim_complete().
 * 
 * Important assumption:
 * Historical events are appropriately spaced such that
 * historical_event_time will never be negative unless
 * there are no historical events left.
 *
 */
genid
Simulator::sim_execute (void) 
{
	if ( getenv( "COSI_NEWSIM" ) ) {
		genid INF(1e30);
		genid gen( 0. );
		while( !demography->dg_done_coalescent() ) {
			genid nextEvtTime = nextEventTime( arrProcs, gen, INF, *getRandGen() );
			//std::cerr << "gen=" << gen << " nextEvt=" << nextEvtTime << "\n";
			if ( nextEvtTime >= INF ) break;

			executeNextEvent( arrProcs, nextEvtTime, *getRandGen() );
			gen = nextEvtTime;
		}
		return gen;
	}
	
	genid gen = ZERO_GEN;
	gens_t historical_event_time;
	gens_t poisson_event_time;	
	int coal = 0;
	bool_t complete_flag = False;	

	while (complete_flag == 0) {

	  coal = 0;
	  historical_event_time = sim_get_time_till_next_hist_event (gen);
	  poisson_event_time = sim_get_time_till_next_pois_event ( gen, historical_event_time );
	  
#ifdef COSI_DEV_PRINT
			PRINT7( gen, recombination_rate, coalesce_rate, migrate_rate, geneconv_rate, poisson_event_time, historical_event_time );
		 if ( !cosi::util::noDbgPrint ) {

			 for ( int i = 0; i < demography->dg_get_num_pops(); i++ )
					std::cerr << "|" << i << " " << demography->dg_get_pop_name_by_index( i ) << " " << demography->dg_get_pop_size_by_index( i ) << " " << demography->dg_get_pop_by_index( i )->pop_get_num_nodes() ;
			 std::cerr << "\n";
		 }
#endif

	  if ( is_null( historical_event_time )
				 || poisson_event_time < historical_event_time) {
			//PRINT2( "poisson: advance by", poisson_event_time );
	    gen += poisson_event_time;
	    coal = sim_do_poisson (gen);
	    if (coal) {complete_flag = demography->dg_done_coalescent();}
	  }
	  else {
			//PRINT2( "historical: advance by", historical_event_time );
	    gen += historical_event_time;
			//PRINT2( "gen bef", gen );
	    gen = histEvents->historical_event_execute(gen);
			//PRINT2( "gen aft", gen );
	    complete_flag = demography->dg_done_coalescent(); 
	  }

	}
	histEvents->processSimEnd( gen );
	return gen;
}

void Simulator::sim_setCoalesce( CoalesceP coalesce_ ) {
	coalesce = coalesce_;
}


/*****************************************************************/
/* INTERNAL FUNCTIONS */

gens_t
Simulator::sim_get_time_till_next_hist_event (genid gen) 
{
	return histEvents->historical_get_time_till_next_event (gen);
}

gens_t
Simulator::sim_get_time_till_next_pois_event (genid gen, gens_t maxWaitTime) 
{
	gens_t time_till_homog_poisson = gens_t( poisson_get_next (sim_get_poisson_rate( gen )) );
	gens_t waitTimeBound = time_till_homog_poisson;
	if ( !is_null( maxWaitTime) && maxWaitTime < time_till_homog_poisson )
		 waitTimeBound = maxWaitTime;
	gens_t time_till_nonhomog_poisson = coalesce->coalesce_get_wait_time_nonhomog( gen, waitTimeBound );
	//PRINT2( time_till_homog_poisson, time_till_nonhomog_poisson );
	gens_t result;
	if ( time_till_homog_poisson < time_till_nonhomog_poisson ) {
		use_what = use_homog;
		waitTimeBound = time_till_homog_poisson;
	} else {
		use_what = use_coal;
		waitTimeBound = time_till_nonhomog_poisson;
	}

	gens_t time_till_side_recomb = sweep3::getSideRecombWaitTime( gen, waitTimeBound );
	if ( time_till_side_recomb < waitTimeBound ) {
		use_what = use_side_recomb;
		waitTimeBound = time_till_side_recomb;
	}
	
	return waitTimeBound;
}

prob_t
Simulator::sim_get_poisson_rate( genid gen )
{
	coalesce_rate = coalesce->coalesce_get_rate();
	migrate_rate = ToDouble( migrate->migrate_get_all_nodes_rate( gen ) );
	recombination_rate = ToDouble( recomb->getAllNodesRecombRate() * genMap->getRegionRecombRateAbs() );
	geneconv_rate = ToDouble( geneConversion->getAllNodesGeneConvRate() * genMap->getRegionRecombRateAbs() );

	poisson_rate = (double) (coalesce_rate + migrate_rate + recombination_rate + geneconv_rate);
	return poisson_rate;
}

int 
Simulator::sim_do_poisson (genid gen) 
{
  int did_coal = 0;
	if ( use_what == use_coal ) {
			int popindex = coalesce->coalesce_pick_popindex_nonhomog();
			//if ( gen > genid(400) ) PRINT3( gen, coalesce_rate, popindex );
			assert( popindex >= 0 );
			//PRINT( "NONHOMOG POISSON" );
			demography->dg_coalesce_by_index (popindex, gen);
			did_coal = 1;
	} else if ( use_what == use_side_recomb ) {
		sweep3::sideRecombExecute();
	} else {

		prob_t randdouble = random_double();
		int popindex;

		if (randdouble < recombination_rate / poisson_rate) {
			recomb->recomb_execute( gen, randdouble / ( recombination_rate / poisson_rate ) );
		}
		else if (randdouble < (recombination_rate + migrate_rate) / poisson_rate) {
			migrate->migrate_execute(gen);
		}
		else if (randdouble < (recombination_rate + migrate_rate + coalesce_rate) / poisson_rate) {
			popindex = coalesce->coalesce_pick_popindex();
			//if ( gen > genid(400) ) PRINT3( gen, coalesce_rate, popindex );
			demography->dg_coalesce_by_index (popindex, gen);
			did_coal = 1;
		}
		else {
			geneConversion->gc_execute(gen,
																 ( randdouble - (recombination_rate + migrate_rate + coalesce_rate) / poisson_rate ) /
																 ( 1.0 - (recombination_rate + migrate_rate + coalesce_rate) / poisson_rate ) );
		}

		
	}
  return did_coal;
}

bool_t
Simulator::sim_complete (void) const {
	return demography->dg_done_coalescent();

}
 
}  // namespace cosi
