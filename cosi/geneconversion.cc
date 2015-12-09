/* $Id: geneconversion.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

/*
	File: geneconversion.cc

	Code for coordinating the execution of gene conversions: determining the probability
	of gene conversion, and invoking the actual execution.
*/

#include <limits>
#include <stdexcept>
#include <boost/math/distributions/geometric.hpp>
#include <cosi/geneconversion.h>
#include <cosi/utils.h>
#include <cosi/genmap.h>
#include <cosi/demography.h>
#include <cosi/node.h>
#include <cosi/pop.h>

namespace cosi {

//
// Class implementation: GeneConversion
//

GeneConversion::GCModel GeneConversion::parseGCModel( const string& gcModel ) {
	if ( gcModel == "right_from_origin" ) return GCM_RIGHT_FROM_ORIGIN;
	if ( gcModel == "left_and_right_from_origin" ) return GCM_LEFT_AND_RIGHT_FROM_ORIGIN;
	if ( gcModel == "left_or_right_from_origin" ) return GCM_LEFT_OR_RIGHT_FROM_ORIGIN;
	throw std::domain_error( "invalid gcModel: " + gcModel );
}


// Constructor:: GeneConversion
//
// Params:
//
//   demography - the <Demography> object
//   genMap - the genetic map of the simulated region. 
//   gc2recombRateRatio - ratio of gene conversion rate to crossover recombination rate.
//   gcMeanTractLenBp - mean length of gene conversion tracts, in basepairs
//   regionLen - length of simulated region, in basepairs.
GeneConversion::GeneConversion( RandGenP randGen_, DemographyP demography_, GenMapP genMap_, len_bp_int_t regionLen_,
																factor_t gc2recombRateRatio_,
																len_bp_int_t gcMeanTractLenBp_, len_bp_int_t gcMinTractLenBp_, GCModel gcModel_ ):
	HasRandGen( randGen_ ),
	demography( demography_ ), genMap( genMap_ ), regionLen( regionLen_ ),
	gc2recombRateRatio( gc2recombRateRatio_ ),
	gcMeanTractLenHalf( gcMeanTractLenBp_ / ( gcModel_ == GCM_LEFT_AND_RIGHT_FROM_ORIGIN ? 2 : 1 ) ), 
	gcMinTractLenHalf( gcMinTractLenBp_ / ( gcModel_ == GCM_LEFT_AND_RIGHT_FROM_ORIGIN ? 2 : 1 ) ),
	gcModel( gcModel_ ),
	gcTractLenDistr( 1. / ToDouble( gcMeanTractLenHalf+1 ) )
{
	assert( equal_eps( ToDouble( boost::math::mean<cosi_double>( boost::math::geometric_distribution<cosi_double>( gcTractLenDistr.p() ) ) ),
										 ToDouble( gcMeanTractLenHalf ) ) );
}

GeneConversion::~GeneConversion() {
}


// Method: getAllNodesGeneConvRate
// Returns the probability, per generation, of a gene conversion that splits some segment in some <Node>.
glen_t GeneConversion::getAllNodesGeneConvRate() const {
	return demography->getNodePool()->getAllNodesGeneConvRate() * gc2recombRateRatio;
}

glen_cM_t GeneConversion::getRegionGeneConvRate() const {
	return genMap->getRegionRecombRateAbs() * gc2recombRateRatio;
}

// Method: sampleTractLen
// Sample the length of one side of a gene conversion tract.
plen_t GeneConversion::sampleTractLen() const {
	len_bp_int_t tractLenBp;
	do {
		tractLenBp = gcTractLenDistr( *getRandGen() );
	} while( tractLenBp < gcMinTractLenHalf );
	return plen_t( ToDouble( tractLenBp ) / ToDouble( regionLen ) );
}

void GeneConversion::gc_execute (genid gen, int popindex, loc_t *loc1, loc_t *loc2, Node**nodes_out) {
	Pop *popptr = demography->dg_get_pop_by_index( popindex );
  int nodeindex = (int) (random_double() * popptr->pop_get_num_nodes());
  Node *n = popptr->pop_get_node (nodeindex);

	loc_t gcOrigin = genMap->getLoc( gloc_t( random_double() ) );
	*loc1 = genMap->sub( gcOrigin, sampleTractLen() );
	*loc2 = genMap->add( gcOrigin, sampleTractLen() );
  demography->dg_gc( n, gen, *loc1, *loc2, nodes_out );
}  // GeneConversion::gc_execute()

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
void GeneConversion::gc_execute (genid gen, frac_t frac) {
	loc_t loc1, loc2;
	gloc_t gcOrigin;
	Node *n = demography->getNodePool()->findGcOrigin( frac, &gcOrigin );
	loc_t gcOriginLoc = genMap->getLoc( gcOrigin );
	if ( gcModel == GCM_LEFT_AND_RIGHT_FROM_ORIGIN ) {
		plen_t leftLen = sampleTractLen();
		plen_t rightLen = sampleTractLen();
		loc1 = genMap->sub( gcOriginLoc, leftLen );
		loc2 = genMap->add( gcOriginLoc, rightLen );
	} else if ( gcModel == GCM_RIGHT_FROM_ORIGIN ) {
		loc1 = gcOriginLoc;
		loc2 = genMap->add( gcOriginLoc, sampleTractLen() );
	} else if ( gcModel == GCM_LEFT_OR_RIGHT_FROM_ORIGIN ) {
		if ( random_bit() ) {
			loc1 = gcOriginLoc;
			loc2 = genMap->add( gcOriginLoc, sampleTractLen() );
		} else {
			loc2 = gcOriginLoc;
			loc1 = genMap->sub( gcOriginLoc, sampleTractLen() );
		}
	} else throw std::logic_error( "invalid gcModel" );
	assert( MIN_LOC <= loc1 && loc1 < loc2 && loc2 <= MAX_LOC ); 
  demography->dg_gc( n, gen, loc1, loc2, (Node **)NULL );
}

}  // namespace cosi

#include <cosi/arrproc2.h>

namespace cosi {

class GeneConvProcessDef: public arrival2::ArrivalProcessDef< genid, RandGen, double > {

	 GeneConversion *geneConv;

public:

	 GeneConvProcessDef( GeneConversion *geneConv_ ):
		 geneConv( geneConv_ ) {
	 }

	 virtual double getRateFactor() const { return ToDouble( geneConv->getAllNodesGeneConvRate() ); }
	 virtual void executeEvent( genid gen, RandGen& ) {
		 geneConv->gc_execute( gen, geneConv->random_double() );
	 }
};  // class GeneConvProcessDef


boost::shared_ptr< GeneConversion::geneconv_processes_type >
GeneConversion::createGeneConvProcesses() {
	return boost::make_shared<geneconv_processes_type>(
			 math::fn_const<genid>( gensInv_t( ToDouble( genMap->getRegionRecombRateAbs() ) ) ),
			 boost::make_shared<GeneConvProcessDef>( this ) );
}  // createGeneConvProcesses

}  // namespace cosi
