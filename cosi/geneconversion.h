/* $Id: geneconversion.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_GENECONVERSION_H
#define __INCLUDE_COSI_GENECONVERSION_H

#include <cstdlib>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/cosirand.h>
#include <cosi/arrproc2.h>

namespace cosi {

using std::string;
namespace node { class Node; }
using node::Node;

//
// Class: Gene Conversion
//
// Code for coordinating the execution of gene conversions.
// There are two main tasks: determining the probability of
// a gene conversion occurring given the current state of the simulation,
// and executing a gene conversion once it has been chosen to occur.
//
// The probability of a gene conversion occurring in each <Node> is actually
// kept by <NodePool>.  The transformation of node <Seglists> resulting from
// gene conversion is implemented in <Seglist>.  The updating of the
// list of nodes in each population resulting from gene conversion is
// implemented in <Demography>.  The GeneConversion class coordinates all that.
//
// The gene conversion model is adapted from the msHOT simulator
// (http://bioinformatics.oxfordjournals.org/content/23/4/520.abstract).
// Briefly, the probability of a gene conversion initiation at a given <loc>
// is proportional to the recombination rate at that loc times a specified constant factor.
// Once a gene conversion is initiated, it spreads independently left and right with length
// sampled from a geometric distribution with a given mean.  Thus, the total length
// of a gene conversion tract is the sum of two geometric distributions.  Note that only
// the initiation point of the gene conversion is chosen according to the genetic map;
// the length of the left and right parts of the gene conversion tract, and therefore the actual boundaries
// of the tract, depend only on the gene conversion length parameter.
//
// See also:
//
//   <Node::gcRate>, <NodePool::node_gc>, <Demography::dg_gc>.
//
class GeneConversion: boost::noncopyable, public HasRandGen {
public:
	 // Enum: GCModel
	 // Choice of gene conversion model.  In all cases, we choose an origination point of the
	 // gene conversion based on the genetic map, and then determine the boundaries of the gene
	 // conversion tract.  The models differ in how we determine the boundaries:
	 //
	 //    RIGHT_FROM_ORIGIN - go right from the origin.  This is the model used by the ms simulator.
	 //    LEFT_AND_RIGHT_FROM_ORIGIN - sample two lengths (left and right) independently.  This is
	 //       the model used by the msHOT simulator.
	 //    LEFT_OR_RIGHT_FROM_ORIGIN - choose a direction (left or right with equal probability),
	 //       and go in that direction from origin.
	 enum GCModel { GCM_RIGHT_FROM_ORIGIN, GCM_LEFT_AND_RIGHT_FROM_ORIGIN, GCM_LEFT_OR_RIGHT_FROM_ORIGIN };

	 static GCModel parseGCModel( const string& gcModel );
	 
	 // Constructor:: GeneConversion
	 //
	 // Params:
	 //
	 //   demography - the <Demography> object
	 //   genMap - the genetic map of the simulated region. 
	 //   gc2recombRateRatio - ratio of gene conversion rate to crossover recombination rate.
	 //   gcMeanTractLenBp - mean length of gene conversion tracts, in basepairs
	 //   regionLen - length of simulated region, in basepairs.
	 GeneConversion( RandGenP randGen_,
									 DemographyP demography_,
									 GenMapP genMap_,
									 len_bp_int_t regionLen_,
									 factor_t gc2recombRateRatio_,
									 len_bp_int_t gcMeanTractLenBp_,
									 len_bp_int_t gcMinTractLenBp_,
									 GCModel gcModel_ );
	 ~GeneConversion();
	 
	 // Method: getAllNodesGeneConvRate
	 // Returns the probability, per generation, of a gene conversion that splits some segment in some <Node>.
	 glen_t getAllNodesGeneConvRate() const;
	 
	 glen_cM_t getRegionGeneConvRate() const;
	 
	 void gc_execute (genid gen, int popindex, loc_t *loc1, loc_t *loc2, Node**nodes_out);

	 //
	 // Method: gc_execute
	 //
	 // Execute a gene conversion.  See description of gene conversion model in the class
	 // documentation for <GeneConversion>.
	 //
	 // Params:
	 //
	 //    gen - the current generation of the simulation.  The two new nodes resulting from
	 //      the gene conversion will have this as their <Node::gen>.
	 //    frac - value indicating in what node, and at what location in that node, to
	 //      initiate the gene conversion.  More specifically: for each Node we keep the
	 //      probability of initiating a gene conversion within that node; frac is a
	 //      fraction of the total sum of these per-node probabilities.  
	 //
	 void gc_execute (genid gen, frac_t frac);

	 typedef
	 arrival2::ArrivalProcess< genid, arrival2::Stoch< RandGen, arrival2::Poisson< math::Const<>, double > > >
	 geneconv_processes_type;

	 boost::shared_ptr< geneconv_processes_type > createGeneConvProcesses();

private:
	 // Field: demography
	 // The demography object; also gives access to the <NodePool>.
	 DemographyP demography;
	 
	 // Field: genMap
	 // The genetic map of the simulated region.
	 GenMapP genMap;
	 
	 // Field: regionLen
	 // Length of the simulated region, in basepairs.
	 len_bp_int_t regionLen;
	 
	 // Field: gc2recombRateRatio
	 // Ratio of gene conversion rate to crossover recombination rate.
	 factor_t gc2recombRateRatio;
	 
	 // Field: gcMeanTractLenHalf
	 // Half of mean length of gene conversion tract.
	 len_bp_int_t gcMeanTractLenHalf;

	 // Field: gcMinTractLenHalf
	 // Half of min length of gene conversion tract.
	 len_bp_int_t gcMinTractLenHalf;

	 // Field: gcModel
	 // The gene conversion model.  See <GCModel>.
	 GCModel gcModel;
	 

	 // Field: gcTractMargin
	 // The 99th percentile of the distance from gc origination point to tract boundary.
	 plen_t gcTractMargin;

	 // Field: gcTractLenDistr
	 // The distribution from which we draw the lengths of the left and right parts
	 // of a gene conversion tract.
	 boost::random::geometric_distribution<len_bp_int_t,prob_t> gcTractLenDistr;

	 //
	 // Private methods
	 //

	 // Method: sampleTractLen
	 // Sample the length of one side of a gene conversion tract.
	 plen_t sampleTractLen() const;
	 
};  // class GeneConversion

//
// Topic: gene conversion origin
//
// In our gene conversion model (adapted from the simulator msHOT), the place where a
// gene conversion originates.  We model a gene conversion by picking a gene conversion
// origin (according to the gene conversion rate at each point which is in turn proportional
// to the recombination rate according to a specified proportionality factor), then choose
// the left and right endpoints of the gene conversion tract by choosing the length from
// a geometric distribution.  In the physical reality, of course, gene conversion originates
// at one endpoint of the tract and ends at another, so "gene conversion origin" should refer
// to an endpoint of the tract.  But in our model the term refers to the point from which
// we go left and right to find the endpoints of the gene conversiont tract.

}  // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_RECOMB_H
