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
#include <string>
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
void print_haps( DemographyP demography, const string& filebase, len_bp_int_t length, MutlistP mutlist,
								 bool_t inf_sites, int outputPrecision );

void print_mut_contexts( DemographyP demography, const string& filebase, len_bp_int_t length,
												 const mutcontext::mutContexts_t& mutContexts );

class ARGOutputHook: public Hook {
public:
	 // Virtual method: handle_add_edge
	 //
	 // Called after adding an edge to the ARG.  The new edge may be the result of a coalescence,
	 // a recombination or a gene conversion.
	 //
	 // Params:
	 //
	 //   nodeId_moreRecent - <nodeid> of the more recent (lower gen) node of the edge
	 //   nodeId_lessRecent - <nodeid> of the less recent (higher gen) node of the edge.
	 //   genId_moreRecent - <genid> of the more recent (lower gen) node of the edge
	 //   genId_lessRecent - <genid> of the less recent (higher gen) node of the edge.
	 //   seglist - the segments inherited along the edge.  NOTE: this seglist may be destroyed
	 //       after the call, so make a copy if you need to save it.
	 //       
	 virtual void handle_add_edge( nodeid nodeId_moreRecent,
																 nodeid nodeId_lessRecent,
																 genid genId_moreRecent,
																 genid genId_lessRecent, const Seglist *seglist,
																 edge_kind_t edgeKind );
	 virtual ~ARGOutputHook();
	 
private:
	 static char edgeKind2code( edge_kind_t edgeKind );
	 
}; // class ARGOutputHook

}  // namespace cosi

#endif // __INCLUDE_COSI_OUTPUT_H
