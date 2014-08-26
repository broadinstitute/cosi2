/* $Id: output.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

//
// Header: output.h
//
// Routines for outputting the result of a simulation.
//


#ifndef __INCLUDE_COSI_OUTPUT_H
#define __INCLUDE_COSI_OUTPUT_H

#include <utility>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <cosi/defs.h>
#include <cosi/mutcontext.h>

namespace cosi {

class Demography;
typedef boost::shared_ptr<Demography> DemographyP;
class Mutlist;
typedef boost::shared_ptr<Mutlist> MutlistP;

typedef boost::multi_array< std::pair<loc_t,loc_t>, 2 > recomblessNeighborhoods_t;

// Method: print_haps
//
// Write the output of one simulation in cosi format.
//
// Params:
//
//    demography - the <Demography> object describing what populations are there, size of sample from each pop, etc
//    filebase - the common filename prefix of the output files.  For each population with <popid> of 'i', two files will be
//      written: a 'filebase'.pos-i file giving SNP positions and a 'filebase'.hap-i file giving the haplotypes.
//    length - the length of the simulated region, in bases
//    mutlist - the list of mutations
//    inf_sites - whether to use an infinite-sites model (if not, of mutations falling on the same integer coordinate only
//      the first is kept)
void print_haps( DemographyP demography, const string& filebase, len_bp_int_t length, MutlistP mutlist, bool_t inf_sites );

void print_mut_contexts( DemographyP demography, const string& filebase, len_bp_int_t length,
												 const mutcontext::mutContexts_t& mutContexts );

}  // namespace cosi

#endif // __INCLUDE_COSI_OUTPUT_H
