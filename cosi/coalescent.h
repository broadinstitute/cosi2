#ifndef __INCLUDE_COSI_COALESCENT_H
#define __INCLUDE_COSI_COALESCENT_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/parameter.hpp>
#include <cosi/defs.h>
#include <cosi/cosirand.h>
#include <cosi/stats.h>
#include <cosi/decls.h>
#include <cosi/nodefwd.h>

namespace cosi {


//
// Class: CoSi
//
// The top-level simulator object.  Orchestrates all simulation tasks:
// parsing of the demographic model from a parameter file,
// setting up and connecting together parts of the simulator, running the simulation,
// reporting the results.  End users who use cosi as a library should need to interact mainly
// with this object.  See <cosimain.cc> for a usage example.
//

class CoSi: public HasRandGen {
public:

	 CoSi();
	 ~CoSi();

	 // Group: Before the simulation

	 // MethodP: setUpSim
	 // Set up the simulator.  This allocates the various sub-objects representing parts of the simulator,
	 // and connects them together properly.
	 void setUpSim( filename_t paramfile, RandGenP randGenToUse_ = RandGenP() );

	 // Method: setMutProcessor
	 // Set the object that will process mutations as they are generated.  Call after <setUpSim()>.
	 void setMutProcessor( MutProcessorP mutProcessor_ );

	 void set_segfp( FILE *segfp_ ) { segfp = segfp_; }
	 void set_logfp( FILE *logfp_ ) { logfp = logfp_; }
	 void set_verbose( bool_t verbose_ ) { verbose = verbose_; }

	 // Method: set_trajFN
	 // Set the name of the file from which to load the frequency trajectory of the selected allele
	 // in the selected pop, for use during sweep simulation.
	 void set_trajFN( string trajFN_ ) { trajFN = trajFN_; }

	 // Method: set_deltaTfactor
	 // Set the factor by which to multiply the step during sweep simulation.
	 void set_deltaTfactor( factor_t deltaTfactor_ ) { deltaTfactor = deltaTfactor_; }

#ifdef COSI_SUPPORT_COALAPX	 
	 void set_maxCoalDist( plen_t maxCoalDist_ ) { maxCoalDist = maxCoalDist_; }
	 void set_maxCoalDistCvxHull( bool maxCoalDistCvxHull_ ) { maxCoalDistCvxHull = maxCoalDistCvxHull_; }
#endif	 
	 
	 void set_trajOutFN( filename_t trajOutFN_ ) { trajOutFN = trajOutFN_; }
	 void set_outputTreeStats( bool_t outputTreeStats_ ) { outputTreeStats = outputTreeStats_; }
	 void set_outputMutGens( bool_t outputMutGens_ ) { outputMutGens = outputMutGens_; }
	 bool_t get_outputMutGens() const { return outputMutGens; }
	 void set_outputRecombLocs( bool_t outputRecombLocs_ ) { outputRecombLocs = outputRecombLocs_; }
	 bool_t get_outputRecombLocs() const { return outputRecombLocs; }

	 void set_genMapShift( ploc_bp_diff_t genMapShift_ ) { this->genMapShift = genMapShift_; }
	 void set_sweepFracSample( bool sweepFracSample_ ) { this->sweepFracSample = sweepFracSample_; }
	 void set_condSnpDef( CondSnpDefP condSnpDef_ ) { this->condSnpDef = condSnpDef_; }
	 void set_recombfileFN( filename_t recombfileFN_ ) { this->recombfileFN = recombfileFN_; }
	 void set_outputARGedges( bool outputARGedges_ ) { outputARGedges = outputARGedges_; }
	 void set_genmapRandomRegions( bool genmapRandomRegions_ ) { genmapRandomRegions = genmapRandomRegions_; }
	 
	 //
	 // Group: Running the simulation
	 //

	 // MethodP: runSim
	 // Run the simulation.
	 // Returns: the generation at which simulation finished.
	 genid runSim();

	 //
	 // Group: After the simulation
	 //

	 ParamFileReaderP getParams() const { return params; }
	 SweepP getSweep() const { return sweep; }
	 DemographyP getDemography() const { return demography; }
	 MutateP getMutate() const { return mutate; }
	 RecombP getRecomb() const { return recomb; }
	 TreeStatsHookP getTreeStatsHook() const { return treeStatsHook; }
	 RecombRecorderP getRecombRecorder() const { return recombRecorder; }
	 CondSnpMgrP getCondSnpMgr() const { return condSnpMgr; }
	 GenMapP getGenMap() const { return genMap; }
	 HistEventsP getHistEvents() { return histEvents; }

	 void write_ms_flags( std::ostream& ) const;

protected:

	 // Field: params
	 // Parameters read from the parameter file.
	 ParamFileReaderP params;

	 // Field: genMap
	 // The genetic map of the simulated region.
	 GenMapP genMap;

	 // Field: recomb
	 // Keeps track of recombination probabilities and executes recombination.
	 RecombP recomb;

	 // Field: geneConversion
	 // Keeps track of gene conversion probabilities and executes gene conversions.
	 GeneConversionP geneConversion;

	 // Field: migrate
	 // Keeps track of migration probabilities and executes migrations.
	 MigrateP migrate;

	 // Field: coalesce
	 // Keeps track of coalescence probabilities and executes coalescences.
	 CoalesceP coalesce;

	 // Field: mutate
	 // Handles the placement of mutations on branches of the ARG.
	 MutateP mutate;

	 // Field: histEvents
	 // Manages historical events.
	 HistEventsP histEvents;

	 // Field: sweep
	 // Manages selective sweeps (this field may be moved out of here).
	 SweepP sweep;

	 // Field: nodePool
	 // Managers the <nodes>, and implements <poisson events> (recombinations,
	 // coalescences, gene conversions) on the nodes.
	 node::NodePoolP nodePool;

	 // Field: demography
	 // The current state of a simulation.
	 DemographyP demography;

	 // Field: simulator
	 // Runs the central loop of a simulation.
	 SimulatorP simulator;

	 // Field: mutProcessor
	 // Handles mutations generated during the simulation.
	 MutProcessorP mutProcessor;

	 // Field: hooks
	 // The list of hooks (callbacks) currently active.  Callbacks
	 // can be registered to be called when one of a fixed set of things
	 // happen in the simulator.  
	 HooksP hooks;

	 // Field: randGenToUse
	 // If non-NULL, use this random generator and ignore any random seed
	 // from the config file.
	 RandGenP randGenToUse;

	 // Field: segfp
	 // A logging file.
	 FILE *segfp;

	 // Field: logfp
	 // A logging file.
	 FILE *logfp;

	 // Field: verbose
	 // Whether to output extra debug information.
	 bool_t verbose;

	 // Field: trajOutFN
	 // Name of file to which sweep frequency trajectory should be written.
	 filename_t trajOutFN;

	 // Field: trajFN
	 // The name of the file from which to load the frequency trajectory of the selected allele
	 // in the selected pop, for use during sweep simulation.
	 filename_t trajFN;

	 // Field: deltaTfactor
	 // Factor by which to multiply the step during sweep simulation.
	 factor_t deltaTfactor;

	 // Field: outputTreeStats
	 // Whether to output tree statistics for each simulation.
	 bool_t outputTreeStats;

	 // Field: outputMutGens
	 // Whether to output mutation times.
	 bool_t outputMutGens;

	 // Field: outputRecombLocs
	 // Whether to output recombination locations
	 bool_t outputRecombLocs;

	 // Field: treeStatsHook
	 // If outputting tree stats, the hook object that gathered the stats.
	 TreeStatsHookP treeStatsHook;

	 // Field: recombRecorder
	 // If recording recomb locs, the recorder hook.
	 RecombRecorderP recombRecorder;

#ifdef COSI_SUPPORT_COALAPX	 
	 // Field: maxCoalDist
	 // Max distance between coalescing segments
	 plen_t maxCoalDist;

	 // Field: maxCoalDistCvxHull
	 // Whether max distance between coalescing segments is computed using convex hull.
	 bool_t maxCoalDistCvxHull;
#endif

	 // Field: genMapShift
	 // A shift applied to the genetic map: added to all physical positions in the genetic map file.
	 ploc_bp_diff_t genMapShift;

	 bool sweepFracSample;

	 // Field: condSnpDef
	 // If conditioning the simulation on the presence of a SNP at a given loc with given freq(s),
	 // the condition.
	 CondSnpDefP condSnpDef;

	 CondSnpMgrP condSnpMgr;

	 // Field: recombfileFN
	 // If non-empty, overrides the recomb file specified in the paramfile
	 filename_t recombfileFN;

	 // Field: outputARGedges
	 // Whether to print ARG edges to stdout.
	 bool_t outputARGedges;

	 // ** Field: genmapRandomRegions
	 // Whether to take for each simulation a different random region from the genetic map.
	 bool_t genmapRandomRegions;

};  // class CoSi


  

}  // namespace cosi  

#endif
// #ifndef __INCLUDE_COSI_COALESCENT_H
