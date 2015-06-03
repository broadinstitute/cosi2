
// * Header cositop.h - top-level 

#ifndef __COSI_INCLUDE_COSITOP_H
#define __COSI_INCLUDE_COSITOP_H

#include <string>
#include <boost/optional.hpp>
#include <cosi/coalescent.h>
#include <cosi/output.h>
#include <cosi/cosirand.h>
#include <cosi/condsnp.h>

namespace cosi {

using std::string;

//
// Class: CoSiMain
//
// The top-level CoSi program.  Handles parsing of command-line, and invoking
// the simulator.
//
class CoSiMain {
public:
	 CoSiMain();
	 int cosi_main( int argc, char *argv[] );
	 
private:
	 filename_t paramfile, recombfileFN, segfile, logfile, outfilebase, trajOutFN;

	 // Field: trajFN
	 // If non-empty, the name of the file from which to load the trajectory of the causal allele
	 // for use during selected <Sweep>.
	 string trajFN;

	 // Field: deltaTfactor
	 // Factor by which to multiply deltaT (step size) during sweep simulation.
	 factor_t deltaTfactor;
	 
	 // Field: msOutput
	 // Whether to produce output in ms simulator format, written to stdout,
	 // or in cosi format, written to separate .pos and .hap files for each population.
	 bool_t msOutput;

	 // Field: outputMutGens
	 // Whether to output mutation times.
	 bool_t outputMutGens;

	 // Field: outputRecombLocs
	 // Whether to output recombination locations
	 bool_t outputRecombLocs;
	 
	 FILE *segfp, *logfp;
	 len_bp_int_t len_length;
	 bool_t sim_only;
	 unsigned int showProgress;
	 bool_t verbose;
	 bool_t treeSizeOnly;
	 int nsims;
	 bool_t showNumRecombs;

	 // Field: outputTreeStats
	 // Output tree statistics, when writing out ms format.
	 bool_t outputTreeStats;

	 // Field: outputPrecision
	 // Number of decimal places in the output
	 int outputPrecision;

	 // Field: outputMutContextsFor
	 // Mutations for which to output the context
	 vector< loc_bp_int_t > outputMutContextsFor;

	 // Field: randSeed
	 // Random seed to be used.
	 RandGen::rseed_t randSeed;

#ifdef COSI_SUPPORT_COALAPX	 
	 double maxCoalDist;
	 bool maxCoalDistCvxHull;
#endif

	 // Field: genMapShift
	 // A shift applied to the genetic map: added to all physical positions in the genetic map file.
	 ploc_bp_diff_t genMapShift;

	 // Field: sweepFracSample
	 // When simulating sweeps, set the derived allele sample freq at end of sweep to exactly
	 // the pop freq, instead of doing a binomial sampling.
	 bool sweepFracSample;

	 // Field: condSnpDef
	 // If conditioning the simulation on the presence of a SNP at a given loc with given freq(s),
	 // the condition.
	 boost::optional<CondSnpDef> condSnpDef;

	 // Field: outputSimTimes
	 // If true, for each simulation output the time it ook.
	 bool outputSimTimes;

	 // Field: outputEndGens
	 // If true, for each simulation output the generation at which it ended.
	 bool outputEndGens;

	 // Field: stopAfterMinutes
	 // If positive, stop simulation after this many minutes (rather than after a given number has been simulated)
	 double stopAfterMinutes;

	 // Field: outputARGedges
	 // Whether to print ARG edges to stdout.
	 bool_t outputARGedges;

	 // Field: freqsOnly
	 // Output allele freqs only.
	 bool_t freqsOnly;

	 // Field: dropSingletonsFrac
	 // Drop this fraction of SNPs which appear on only one chrom across all pops
	 frac_t dropSingletonsFrac;

	 // ** Field: genmapRandomRegions
	 // Whether to take for each simulation a different random region from the genetic map.
	 bool_t genmapRandomRegions;

	 // ** Field: outputPopInfo
	 // Whether to output population info (in ms output mode).
	 bool_t outputPopInfo;
	 
	 // ** Field: outputGenMap
	 // Whether to output the genetic map (in ms output mode).
	 bool_t outputGenMap;
	 
	 int parse_args( int argc, char *argv[] );
	 static void printCompileOptions();
};  // class CoSiMain


}  // namespace cosi

#endif  // __COSI_INCLUDE_COSITOP_H
