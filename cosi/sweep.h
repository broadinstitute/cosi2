/* $Id: sweep.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_SWEEP_H
#define __INCLUDE_COSI_SWEEP_H

#include <string>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/nodefwd.h>
#include <cosi/cosirand.h>
#include <cosi/generalmath.h>

namespace cosi {

//
// Class: Sweep
//
// Implements a selective sweep (currently, on a single allele in a single population).
// The key method is <sweep_execute()>.
//
class Sweep: public HasRandGen {
public:

	 //
	 // Method group: Initialization/setup
	 //
	 
	 Sweep( DemographyP demography_ );
	 
	 // MethodP: set_trajFN
	 // Set the name of the file from which to load the frequency trajectory of the selected allele
	 // in the selected pop, for use during sweep simulation.
	 void setTrajFN( filename_t trajFN_ );

	 void sweep_set_fix_freq( bool_t );
	 void sweep_setGenMap( GenMapP genMap_ );
	 void sweep_setRecomb( RecombP recomb_ );
	 void sweep_setGeneConversion( GeneConversionP geneConversion_ );
	 void sweep_setMutate( MutateP mutate_ );
	 void sweep_setNodePool( node::NodePoolP nodePool_ );
	 void set_deltaTfactor( factor_t deltaTfactor_ ) { deltaTfactor = deltaTfactor_; }

	 freq_t getFinalFreq() const { return finalFreq; }
	 loc_t getSelLoc() const { return selLoc; }

	 void setVerbose( bool_t verbose_ ) { verbose = verbose_; }
	 bool_t getVerbose() const { return verbose; }

	 void setTrajOutFN( filename_t trajOutFN_ ) { trajOutFN = trajOutFN_; }

	 //
	 // Method group: executing the sweep
	 //

	 //
	 // MethodP: sweep_execute
	 //
	 // Execute the selected sweep on one selected allele in one population,
	 // with a deterministically chosen selected allele trajectory.
	 //
	 // Params:
	 //
	 //   popname - the population in which the sweep happens
	 //   selCoeff - the selection coefficient.  We assume the heterozygous effect is .5.
	 //   gen - the generation when the sweep ends (moving forward) or starts (moving pastward).
	 //   sel_pos - the location of the selected allele within the simulated region
	 //   final_sel_freq - the frequency to which the selected allele rises by the end of the sweep
	 //      (moving forward), i.e. at time 'gen'.
	 //
	 genid sweep_execute(popid popname, gensInv_t selCoeff, genid gen, loc_t sel_pos, 
											 freq_t final_sel_freq);

	 //
	 // MethodP: sweep_stochastic_execute
	 //
	 // Execute selected sweep with a stochastically chosen selected allele trajectory
	 // (work-in-progress).
	 //
	 //
	 // Params:
	 //
	 //   popname - the population in which the sweep happens
	 //   selCoeff - the selection coefficient.  We assume the heterozygous effect is .5.
	 //   gen - the generation when the sweep ends (moving forward) or starts (moving pastward).
	 //   sel_pos - the location of the selected allele within the simulated region
	 //   final_sel_freq - the frequency to which the selected allele rises by the end of the sweep
	 //      (moving forward), i.e. at time 'gen'.
	 //
	 genid sweep_stochastic_execute(popid popname, gensInv_t selCoeff, double heterozygousEffect,
																	genid gen, loc_t sel_pos, 
																	freq_t final_sel_freq);
	 
	 
private:

	 //
	 // Field group: Infrastructure parts
	 //
	 // Other major singleton objects implementing parts of the simulator
	 //
	 
	 // Field: demography
	 // The <Demography> object containing the state of the simulation.
	 // When the sweep is executed, the state of the simulation in this object is updated.
	 DemographyP demography;

	 // Field: genMap
	 // The genetic map of the simulated region.
	 GenMapP genMap;

	 // Field: recomb
	 // The object that executes recombination events.
	 RecombP recomb;

	 // Field: geneConversion
	 // The object that executes gene conversion events.
	 GeneConversionP geneConversion;

	 // Field: mutate
	 // The object that handles the recording of mutations generated
	 // during the simulation.
	 MutateP mutate;

	 // Field: nodePool
	 // The object that keeps track of all the <Nodes> in the simulation.
	 node::NodePoolP nodePool;

	 //
	 // Field group: Sweep params
	 //
	 // Values affecting how the sweep is simulated.
	 // 

	 // Field: finalFreq
	 // The frequency to which the causal allele rises in the selected population by the end
	 // of the sweep (going forward).
	 freq_t finalFreq;

	 // Field: selLoc
	 // Location of the causal allele within the simulated region.
	 loc_t selLoc;

	 // Field: fix_freq
	 // If true, the frequency of the selected allele in the sample is set precisely to
	 // its frequency in the population undergoing sweep.  If false (the default), frequency of the selected
	 // allele in the sample is obtained by doing a binomial sampling based on the population frequency
	 // of the allele.
	 bool_t fix_freq;

	 // Field: traj_file
	 // If non-empty, name of file specifying the frequency trajectory of the causal allele in the
	 // population under sweep.
	 filename_t trajFN;

	 // Field: deltaTfactor
	 // Factor by which to multiply the step during sweep simulation.
	 factor_t deltaTfactor;

	 //
	 // Field group: Output specifications
	 //
	 // Control details of how sweep output it produced.
	 //

	 // Field: trajOutFN
	 // File to which we write out the trajectory of the causal allele
	 filename_t trajOutFN;

	 // Field: verbose
	 // Whether to print verbose status messages during the sweep.
	 bool_t verbose;

	 // Field group: Fields used during simulation of the sweep

	 //
	 // Fields: Non-selected event rates
   // 
	 // Rates of various events in all populations EXCEPT the one under sweep.
	 //
	 //   coalesce_rate_nonsel - total coalescence rate in non-sweep pops
	 //   migrate_rate_nonsel  - total migration rate in non-sweep pops
	 //   recombination_rate_nonsel - total recombination rate in non-sweep pops
	 //   geneconv_rate_nonsel - total gene conversion rate in non-sweep pops
	 //   poisson_rate_rate - total rate of <Poisson events> in non-sweep pops
	 //
	 
	 double coalesce_rate_nonsel;
	 double migrate_rate_nonsel;
	 double recombination_rate_nonsel;
	 double geneconv_rate_nonsel;
	 double poisson_rate_nonsel;

	 // Field: selAlleleFreqTraj
	 // Frequency trajectory of the selected allele, as an interpolated
	 // function from generation to the frequency of the selected allele in
	 // the selected pop at that generation.
	 math::InterpFn<genid,freq_t> selAlleleFreqTraj;


	 //////////////////////
	 // Private methods //
	 ////////////////////

	 // MethodP: sweep_load_traj
	 // Load the causal allele trajectory from <trajFN> into <selAlleleFreqTraj>, if
	 // trajFN was given.
	 void sweep_load_traj();
	 
	 int sw_coalesce(int nnode, Node **nodes, Pop *popptr, genid t);
	 double sw_get_poisson_rate_nonsel(popid sel_popname);
	 double sw_coalesce_get_rate_nonsel (popid sel_popname) const;
	 double sw_recomb_get_rate_nonsel (popid sel_popname) const;
	 double sw_gc_get_rate_nonsel (popid sel_popname) const;
	 void sw_do_poisson_nonsel(popid sel_popname, genid gen);
	 int sw_trajectory(popid popname, double s);
	 void sweep_sim_freq_traj( genid start_gen );

	 int sw_coalesce_pick_popindex (popid sel_popname);
	 int sw_recomb_pick_popindex(popid sel_popname);

	 
};  // class Sweep

// End class: Sweep

/////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace cosi


#endif
// #ifndef __INCLUDE_COSI_SWEEP_H
