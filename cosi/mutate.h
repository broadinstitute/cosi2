/* $Id: mutate.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */
#ifndef __INCLUDE_COSI_MUTATE_H
#define __INCLUDE_COSI_MUTATE_H
#include <cstdio>
#include <iostream>
#include <cosi/general/math/cosirand.h>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>
#include <cosi/seglistfwd.h>

namespace cosi {


// Abstract class: MutProcessor
//
// Receives mutations generated during the simulation, and decides what
// to do with them.
//
// See <MutProcessor_AddToMutlist>.
class MutProcessor {
public:
	 MutProcessor() {}
	 virtual ~MutProcessor();

	 // Virtual method: processMut
	 // Called for each mutation encountered during the simulation.  Default implementation
	 // does nothing.
	 virtual void processMut(loc_t loc, leafset_p leaves, genid gen, popid popName);

	 // Virtual method: postprocess
	 // Called after the simulation has finished, to let the object to postprocessing
	 // tasks (e.g. sort the mutations by location).  The default implementation does nothing.
	 virtual void postprocess();
};

//
// Class: Mutate
//
// Code for putting mutations on branches of the <ARG>.
// The mutations that are generated are sent to a <MutProcessor>.
//
class Mutate: public HasRandGen {
public:
	 typedef seglist::Seglist Seglist;
	 
	 Mutate( RandGenP randGen_, prob_per_bp_per_gen_t mu_, len_bp_t chrom_len_bp_, bool_t treeSizeOnly_ = False );

	 // MethodP: mutate_put_muts_on_seglist
	 // Put mutations on an edge of the ARG.
	 // Params:
	 //   seglist - the list of segs inherited along the ARG edge
	 //   edge_gens - the length of the ARG edge.
	 void mutate_put_muts_on_seglist( const Seglist *seglist, gens_t edge_gens,
																		genid edgeMinGen, popid popName )
#ifdef COSI_DISABLE_MUTS
			{}
#else	 
		 ;
#endif	 

	 // Method: setMutProcessor
	 // Sets the object to receive mutations as they are generated during the simulation.
	 void setMutProcessor( MutProcessorP mutProcessor_ );
	 
	 // MethodP: mutate_print_leafset
	 // Add a mutation at a specified location, inherited by the specified leaves.
	 void mutate_print_leafset (loc_t loc, leafset_p leaves, genid gen, popid popName);
	 
	 prob_per_len_per_gen_t getMu() const { return mutate_mu; }
	 len_bp_t getLength() const { return mutate_chrom_len_bp; }

	 void setFreqsOnly( bool_t freqsOnly_ ) { freqsOnly = freqsOnly_; }
	 void setDemography( const Demography *demography_ ) { demography = demography_; }
	 void writeTreeSize() const;

private:
	 // Field: mutProcessor
	 // The object that receives the mutations generated during the simulation.
	 MutProcessorP mutProcessor;

	 // Field: mutate_mu
	 // The mutation rate, as the probability of mutation per generation per region length.
	 // Currently we assume that the mutation rate does not depend on location on the chromosome.
	 prob_per_len_per_gen_t mutate_mu;

   /* Field: mutate_chrom_len_bp */
	 /* The length of the chromosome, in basepairs, stored locally in this module. */
	 len_bp_t mutate_chrom_len_bp;

	 /*
		* Static var: wait_left
		*
		* When placing mutations on edges of the ARG, keeps the poisson "wait time" left
		* until the next mutation.  (Of course the waiting times are memoryless, but
		* keeping this waiting time saves us from having to simulate additional Poisson-distributed
		* waiting times.)
		*/
	 mutime_t wait_left;

	 // Field: waitInitialized
	 // Whether we have chosen the first wait time, until the first mutation.
	 bool_t waitInitialized;

	 bool_t freqsOnly;
	 gens_t totTreeSize;
	 const Demography *demography;
	 std::vector<nchroms_t> edgeDerCounts;
	 std::vector<gens_t> edgeLens;
	 
};  // class Mutate

}  // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_MUTATE_H


