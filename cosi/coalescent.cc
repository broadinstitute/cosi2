/* $Id: coalescent.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

//
// File: coalescent.cc
//
// Defines implementation of <CoSi>, the top-level simulator object.
// <CoSi::setUpSim()> allocates and wires together the various objects implementing different
// parts of simulator behavior; <CoSi::runSim()> then runs the simulation.
//

//#define COSI_DEV_PRINT

#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cosi/hooks.h>
#include <cosi/node.h>
#include <cosi/coalescent.h>
#include <cosi/mutlist.h>
#include <cosi/file.h>
#include <cosi/sweep.h>
#include <cosi/demography.h>
#include <cosi/simulator.h>
#include <cosi/mutate.h>
#include <cosi/seglist.h>
#include <cosi/historical.h>
#include <cosi/recomb.h>
#include <cosi/geneconversion.h>
#include <cosi/coalesce.h>
#include <cosi/recomb.h>
#include <cosi/genmap.h>
#include <cosi/migrate.h>
#include <cosi/utils.h>
#include <cosi/stats.h>
#include <cosi/sweep1.h>
#include <cosi/sweep2.h>
#include <cosi/sweep3.h>
#include <cosi/condsnp.h>
#include <cosi/output.h>

namespace cosi {

//
// Class impl: CoSi
//
// Implementation of class CoSi, the top-level object representing the simulator.
//

CoSi::CoSi(): segfp( NULL ), logfp( NULL ), verbose( False ),
							deltaTfactor( 1.0 ),
							outputTreeStats( False ), outputMutGens( False ),
							outputRecombLocs( False )
#ifdef COSI_SUPPORT_COALAPX							
						, maxCoalDist( plen_t( 1.0 ) ), maxCoalDistCvxHull( False )
#endif
						, genMapShift( 0 ), sweepFracSample( False ), outputARGedges( False ),
							genmapRandomRegions( False )
{
	seglist::seglist_init_module();
}

CoSi::~CoSi() {
}

// Method: setUpSim
// Set up the simulator.  This allocates the various sub-objects representing parts of the simulator,
// and connects them together properly.
//
// There are several singleton objects handling different aspects of the simulation.
// (Singleton as in, one per CoSi object).  They need to share some state and refer to each other.
// This method allocates them all and connects them to each other as needed.
void CoSi::setUpSim( filename_t paramfile, RandGenP randGenToUse_, GenMapP genMapToUse_ ) {
	using boost::make_shared;

	hooks = make_shared<Hooks>();

	demography = make_shared<Demography>();
	demography->dg_set_logfile(logfp);
	demography->setVerbose( verbose );
	demography->setHooks( hooks );

	nodePool = make_shared<node::NodePool>();
	nodePool->setHooks( hooks );

	if ( outputTreeStats )
		hooks->addHook( treeStatsHook = make_shared<TreeStatsHook>() );

	if ( condSnpDef ) {
		hooks->addHook( condSnpMgr = make_shared<CondSnpMgr>( demography, *condSnpDef ) );
	}

	if ( outputARGedges ) {
		hooks->addHook( make_shared<ARGOutputHook>() );
	}
	
	//if ( treeSizeOnly ) nodePool->setOutputMuts( False );
	demography->dg_setNodePool( nodePool );
			 
	params = make_shared<ParamFileReader>( demography, randGenToUse_ );
	params->set_recombfileFN( this->recombfileFN );

	params->file_read(paramfile, segfp);
	setRandGen( params->getRandGen() );

	// RandGenP rgen;
	// if ( randGenToUse_ ) rgen = randGenToUse_;
	// else if ( params->isSeeded() ) {
	// 	unsigned long rseed = params->getRandSeed();
	// 	rgen = make_shared<RandGen>( rseed );
	// } else rgen = make_shared<RandGen>();
 //	setRandGen( rgen );
	demography->setRandGen( getRandGen() );
#ifdef COSI_SUPPORT_COALAPX
	demography->setMaxCoalDist( this->maxCoalDist );
	demography->setMaxCoalDistCvxHull( this->maxCoalDistCvxHull );
#endif	
	nodePool->setRandGen( getRandGen() );
	params->getHistEvents()->setRandGen( getRandGen() );

	if ( genMapToUse_ ) {
		genMap = genMapToUse_;
		if ( genmapRandomRegions )
			 genMap->setStart( static_cast<len_bp_t>( getRandGen()->random_idx( genMapToUse_->recomb_get_map_length() - params->getLength() ) ) );
	}
	else {
		 genMap = boost::make_shared<GenMap>( params->get_recombfileFN(), params->getLength() );
		 if ( genMapShift > 0 ) genMap->setStart( genMapShift );
	}


	nodePool->setGenMap( genMap );
	if ( params->getGeneConv2RecombRateRatio() == ZERO_FACTOR ) nodePool->setEnableGeneConv( False );

	recomb = make_shared<Recomb>( demography, genMap );
	recomb->setRandGen( getRandGen() );
	recomb->setIgnoreRecombsInPop( params->getIgnoreRecombsInPop() );

	sweep = make_shared<Sweep>( demography );
	sweep->setRandGen( getRandGen() );
	sweep->setVerbose( verbose );
	params->getHistEvents()->historical_setSweep( sweep );
	sweep->sweep_setNodePool( nodePool );
	sweep->setTrajFN( trajFN );
	sweep1::sweep1_setTrajFN( trajFN );
	sweep3::sweep3_setTrajFN( trajFN );
	sweep->setTrajOutFN( trajOutFN );
	sweep->set_deltaTfactor( deltaTfactor );
	sweep1::sweep1_set_deltaTfactor( deltaTfactor );

	sweep2::sweep2_set_sweepFracSample( this->sweepFracSample );

	simulator = make_shared<Simulator>( demography, genMap );
	simulator->setRandGen( getRandGen() );
			 
	simulator->sim_setRecomb( recomb );
	simulator->sim_setHistEvents( params->getHistEvents() );

	coalesce = make_shared<coal::Coalesce>( demography );
	simulator->sim_setCoalesce( coalesce );
	coalesce->setRandGen( getRandGen() );

	{
		len_bp_int_t len = params->getLength();
		factor_t geneConv2RecombRateRatio = params->getGeneConv2RecombRateRatio();
		len_bp_int_t geneConversionMeanTractLength = params->getGeneConversionMeanTractLength();
		len_bp_int_t geneConversionMinTractLength = params->getGeneConversionMinTractLength();
		GeneConversion::GCModel geneConversionModel = params->getGeneConversionModel();
		geneConversion = make_shared<GeneConversion>( getRandGen(), demography, genMap, len,
																									geneConv2RecombRateRatio,
																									geneConversionMeanTractLength,
																									geneConversionMinTractLength,
																									geneConversionModel );
	}

	simulator->sim_setGeneConversion( geneConversion );

	sweep->sweep_setGeneConversion( geneConversion );
	
	migrate = make_shared<Migrate>( demography );
	migrate->setRandGen( getRandGen() );
	migrate->setHooks( hooks );
	params->getHistEvents()->historical_setMigrate( migrate );
	simulator->sim_setMigrate( migrate );
	
	demography->dg_complete_initialization();

	BOOST_AUTO( traj, params->get_pop2sizeTraj() );
	for( BOOST_AUTO( it, traj.begin() ); it != traj.end(); it++ ) {
		Pop *pop = demography->dg_get_pop_by_name( it->first );
		if ( !pop ) throw std::runtime_error( "trajectory specified for unknown population" );
		pop->setSizeTraj( it->second, genid( 0.0 ) );
	}

	if ( getenv( "COSI_NEWSIM" ) ) {
		BaseModelP baseModel = params->getBaseModel();
		for( BOOST_AUTO( it, baseModel->popInfos.begin() ); it != baseModel->popInfos.end(); it++ ) {
			Pop *pop = demography->dg_get_pop_by_name( it->first );
			if ( !pop ) throw std::runtime_error( "trajectory specified for unknown population" );
			BaseModel::PopInfo const& popInfo = it->second;
			pop->setCoalRateFn( popInfo.coalRateFn, genid( 0.0 ) );
			BOOST_AUTO( base, ( math::Function< genid, double, math::Const<> >( 2.0 ) * popInfo.popSizeFn ) );
			PRINT( base );
			// BOOST_AUTO( coalRateFn,
			// 						(math::Function< genid, double, math::Const<> >( .5 ) *
			// 						 math::pow( base,
			// 												math::fn_const<double>(-1.) ) ) );
			// PRINT( coalRateFn );
			// PRINT2( it->first, coalRateFn );
		}
	}
	
	mutate.reset( new Mutate( getRandGen(), params->getMu(), params->getLength() ) );
	demography->dg_setMutate( mutate );
	
	sweep->sweep_setGenMap( genMap );
	sweep->sweep_setRecomb( recomb );
	sweep->sweep_setMutate( mutate );
	//			 sweep->sweep_set_fix_freq( True );
	
	if ( outputRecombLocs ) {
		recombRecorder = make_shared<RecombRecorder>();
		hooks->addHook( recombRecorder );
	}

}  // CoSi::setUpSim()

void CoSi::setMutProcessor( MutProcessorP mutProcessor_ )  { mutate->setMutProcessor( mutProcessor_ ); }

genid CoSi::runSim() {
	return simulator->sim_execute();
}  // runSim()

// End class impl: CoSi

}  // namespace cosi


