/* $Id: simulator.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_SIMULATOR_H
#define __INCLUDE_COSI_SIMULATOR_H
#include <cstdio>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/cosirand.h>

namespace cosi {

//
// Class: Simulator
//
// The "inner loop" of the simulations.  Runs the actual backwards simulation.
//
class Simulator: public HasRandGen {
public:
	 Simulator( DemographyP demography_, GenMapP genMap_ );

	 // MethodP: sim_execute
	 // Run the simulation.
	 genid sim_execute (void);

	 // Setup methods
	 
	 void sim_setGeneConversion( GeneConversionP geneConversion_ ) { geneConversion = geneConversion_; }
	 void sim_setMigrate( MigrateP migrate_ ) { migrate = migrate_; }
	 void sim_setHistEvents( HistEventsP histEvents_ ) { histEvents = histEvents_; }
	 void sim_setRecomb( RecombP recomb_ ) { recomb = recomb_; }
	 void sim_setCoalesce( CoalesceP coalesce_ );

private:
	 DemographyP demography;
	 GenMapP genMap;
	 prob_t coalesce_rate;
	 prob_t migrate_rate;
	 prob_t recombination_rate;
	 prob_t geneconv_rate;
	 prob_t poisson_rate;
	 
	 RecombP recomb;
	 GeneConversionP geneConversion;
	 MigrateP migrate;
	 HistEventsP histEvents;
	 MutateP mutate;
	 CoalesceP coalesce;

	 prob_t sim_get_poisson_rate( genid gen );
	 int sim_do_poisson (genid gen);
	 bool_t sim_complete (void) const;
	 gens_t sim_get_time_till_next_hist_event (genid gen);
	 gens_t sim_get_time_till_next_pois_event (genid gen, gens_t maxWaitTime);

	 enum use_what_t { use_homog, use_coal, use_side_recomb, use_none };
	 use_what_t use_what;
};

}  // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_SIMULATOR_H
