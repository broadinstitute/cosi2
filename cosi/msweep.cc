
// * File: msweep.cc - implementation of selective sweeps

//#define COSI_DEV_PRINT

#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/numeric.hpp>
#include <boost/array.hpp>
#include <boost/assign/std/map.hpp>
//#include <boost/units/detail/utility.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/exception/all.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/exception.hpp>
#include <boost/tokenizer.hpp>
#include <cosi/general/utils.h>
#include <cosi/general/math/cosirand.h>
#include <cosi/general/math/generalmath.h>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/genmap.h>
#include <cosi/basemodel.h>
#include <cosi/node.h>
#include <cosi/hooks.h>
#include <cosi/mutlist.h>
#include <cosi/demography.h>
#include <cosi/leavesinfo.h>
#include <cosi/msweep.h>
#include <cosi/module.h>

namespace cosi {

// ** struct: SweepInfo - definition of a selected sweep
struct SweepInfo {
	 
// *** Field: selGen - the generation at which the selected allele is born
	 genid selGen;
// *** Field: selCoeff - the selection coefficient.
	 double selCoeff;
// *** Field: the location of the selected allele in the simulated region (should permit location outside?)
	 loc_t selPos;
// *** Field: selPop - the population in which the selected allele is born.
	 popid selPop;
	 util::ValRange<freq_t> final_sel_freq;
// *** Field: selBegPop - the population in which selection begins
	popid selBegPop;
// *** Field: selBegGen - the time at which selection begins
	genid selBegGen;

	 SweepInfo(): selGen( NULL_GEN ), selCoeff( 0.0 ), selPos( 0.0 ),
								selPop( NULL_POPID ), selBegPop( NULL_POPID ), selBegGen( NULL_GEN ) { }
	 SweepInfo( genid selGen_, double selCoeff_, loc_t selPos_, popid selPop_,
							util::ValRange<freq_t> final_sel_freq_, popid selBegPop_, genid selBegGen_ ):
		 selGen( selGen_ ), selCoeff( selCoeff_ ), selPos( selPos_ ), selPop( selPop_ ),
		 final_sel_freq( final_sel_freq_ ), selBegPop( selBegPop_ ), selBegGen( selBegGen_ ) { }
};  // struct SweepInfo

void setSweepInfo( BaseModel& baseModel,
									 genid selGen, double selCoeff, loc_t selPos, popid selPop,
									 util::ValRange<freq_t> final_sel_freq, popid selBegPop, genid selBegGen ) {
	baseModel.sweepInfo = boost::make_shared<SweepInfo>( selGen, selCoeff, selPos, selPop, final_sel_freq, selBegPop, selBegGen );
}

// ** class MSweep - implementation of selective sweep in multiple populations
//
// Note: some code here is adapted from simuPOP simulator by Bo Peng et al ; see
// http://simupop.sourceforge.net/manual_svn/build/userGuide_ch7_sec2.html	 
class MSweep: public Module {
public:

// *** Type: pop_traj_t - frequency trajectory of an allele in one pop	 
	 typedef math::Function< genid, freq_t, math::Piecewise< math::Const<> > > pop_traj_t;

// *** Type: mpop_traj_t - frequency trajectory of an allele in multiple pops
	 typedef std::map< popid, pop_traj_t > mpop_traj_t;

	 DemographyP demography;

// *** Field: mtraj - the frequency trajectory of the selected allele in each pop
	boost::shared_ptr<mpop_traj_t> mtraj;
	 BaseModelP baseModel;
	 BaseModelP sweepModel;

	 RandGenP randGen;
	 
	 leafset_p selLeaves;

	 std::map< popid, popid > pop2sib;
	 std::set< popid > selPops;

	 COSI_DECL(SweepHook);
	 SweepHookP sweepHook;

	 MSweep( DemographyP demography_, BaseModelP baseModel_, RandGenP randGen_ ):
		 demography( demography_ ), baseModel( baseModel_ ), randGen( randGen_ ) {

		 selLeaves = make_empty_leafset();
		
		 if ( baseModel->sweepInfo ) {
			 cosi_using2( std::map, boost::assign::insert );

			 SweepInfo const& sweepInfo = *baseModel->sweepInfo;
			 
			 map<popid, map<genotype_t,double> > fits;
			 map<popid, freq_t> begFreqs;
			 map<popid, util::ValRange<freq_t> > endFreqs;

			 cosi_for_map( pop, popInfo, baseModel->popInfos ) {
				 double s;
				 if ( pop != sweepInfo.selPop ) {
					 begFreqs[ pop ] = 0.;
					 endFreqs[ pop ] = util::make_val_range( 0., 1. );
					 s = 0;
				 } else {
					 begFreqs[ pop ] = 1. / ( 2. * ToDouble( popInfo.popSizeFn( sweepInfo.selGen ) ) );
					 endFreqs[ pop ] = sweepInfo.final_sel_freq;
					 s = sweepInfo.selCoeff;
				 }
				 insert( fits[ pop ] )( GT_AA, 1. + s )( GT_Aa, 1. + 0.5*s )( GT_aa, 1. );
			 } cosi_end_for;

			 size_t maxAttempts = 1000000;
			 if ( getenv( "COSI_MAXATTEMPTS" ) ) {
				 try { maxAttempts = boost::lexical_cast<size_t>( getenv( "COSI_MAXATTEMPTS" ) ); }
				 catch( const boost::bad_lexical_cast& ) { BOOST_THROW_EXCEPTION( cosi_error() <<
																																					error_msg( "invalid COSI_MAXATTEMPTS" ) ); }
			 }

			 if ( getenv( "COSI_LOAD_TRAJ" ) ) {
				 this->mtraj = this->loadTraj( getenv( "COSI_LOAD_TRAJ" ) );
			 } else
					this->mtraj = this->simulateTrajFwd( baseModel, fits, sweepInfo.selGen,
																							 begFreqs, endFreqs,
																							 sweepInfo.selBegPop, sweepInfo.selBegGen, sweepInfo.selCoeff,
																							 *randGen, maxAttempts );

			 if(0){
				 std::cerr << "traj  is:\n";
				 gens_t STEP(1);
				 for ( genid gen = sweepInfo.selGen; gen >= genid(0); gen -= STEP ) {
					 std::cerr << gen;
					 cosi_for_map_values( traj, *mtraj ) {
						 std::cerr << "\t" << traj( gen );
					 } cosi_end_for;
					 std::cerr << "\n";
				 }
			 }

			 using util::operator<<;
			 if ( getenv( "COSI_SAVE_TRAJ") ) {

				 static int simNum = 0;
				 ++simNum;

				 std::ofstream f;

				 f.exceptions( std::ios::failbit | std::ios::badbit );
				 std::ofstream::openmode mode = std::ofstream::out;
				 if ( simNum > 1 ) mode = std::ofstream::out | std::ofstream::app;
				 f.open( getenv( "COSI_SAVE_TRAJ"), mode );
				 f.precision(12);

				 if ( simNum == 1 ) {
					 f << "sim\tgen";
					 cosi_for_map_keys( pop, *mtraj ) {
						 f << "\tselfreq_" << pop;
					 } cosi_end_for;
					 f << "\n";
				 }

				 if ( getenv("COSI_SAVE_TRAJ_SELBEG_ONLY") ) {
						 f << simNum << "\t" << sweepInfo.selBegGen;
						 cosi_for_map_values( traj, *mtraj ) {
							 f << "\t" << traj( sweepInfo.selBegGen );
						 } cosi_end_for;
						 f << "\n";
				 } else {

					 gens_t STEP(1);
					 for ( genid gen = sweepInfo.selGen; gen >= genid(0); gen -= STEP ) {
						 f << simNum << "\t" << gen;
						 cosi_for_map_values( traj, *mtraj ) {
							 f << "\t" << traj( gen );
						 } cosi_end_for;
						 f << "\n";
					 }
				 }
			 }


			 if ( getenv( "COSI_SAVE_SWEEP_INFO") ) {

				 static int simNum = 0;
				 ++simNum;

				 std::ofstream f;

				 f.exceptions( std::ios::failbit | std::ios::badbit );
				 std::ofstream::openmode mode = std::ofstream::out;
				 if ( simNum > 1 ) mode = std::ofstream::out | std::ofstream::app;
				 f.open( getenv( "COSI_SAVE_SWEEP_INFO"), mode );
				 f.precision(12);

				 // if ( simNum == 1 ) {
				 // 	 f << "sim\tgen";
				 // 	 cosi_for_map_keys( pop, *mtraj ) {
				 // 		 f << "\tselfreq_" << pop;
				 // 	 } cosi_end_for;
				 // 	 f << "\n";
				 // }

				 {
					 f << simNum << "\t" 
						 << sweepInfo.selPop << "\t" 
						 << sweepInfo.selGen << "\t"
						 << sweepInfo.selBegPop << "\t"
						 << sweepInfo.selBegGen << "\t"
						 << sweepInfo.selCoeff << "\t"
					 cosi_for_map( pop, traj, *mtraj ) {
						 if ( pop == sweepInfo.selBegPop ) 
							 f << traj( genid(0) );
					 } cosi_end_for;
					 f << "\n";
				 }
			 }

			 // std::cerr.precision(8);
			 // std::cerr << "got  traj: " << *this->mtraj << "\n";

			 // if ( this->trajOnly ) {
			 // 	 cosi_for_map( pop, traj, mtraj ) {
					
			 // 	 } cosi_end_for;
			 // }
			
			 this->sweepModel = this->makeSweepModel( this->baseModel, *this->mtraj,
																								sweepInfo.selPos, this->pop2sib, this->randGen );
			 demography->getNodePool()->setSelPos( sweepInfo.selPos );
		 }  // if selCoeff nonzero
		 else {
				this->sweepModel = baseModel;
				//std::cerr << "GOT NO SWEEP\n";
		 }

		 //std::cerr << "baseModel=" << *baseModel << "\n";
		 //std::cerr << "sweepModel=" << *sweepModel << "\n";
	 }  // MSweep()

	 class SweepHook: public Hook {

			MSweep *msweep;
			loc_t selPos;

	 public:			
			SweepHook( MSweep *msweep_, loc_t selPos_ ): msweep( msweep_ ), selPos( selPos_ ) { }
			
			virtual void handle_recomb( Node *node1, Node *node2, loc_t loc, genid curGen ) {
				msweep->determineAlleleAtSelPos( selPos < loc ? node2 : node1, curGen );
			}
			virtual void handle_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2,
															genid curGen ) {
				msweep->determineAlleleAtSelPos( loc1 <= selPos && selPos <= loc2 ? node2 : node1,
																				 curGen );
			}
			
	 };  // class SweepHook

	 void determineAlleleAtSelPos( Node *node, genid gen ) {
		 cosi_using2( util::at, util::STLContains );
		 
		 if ( node ) {
			 popid pop_this = node->getPop()->pop_get_name();
			 BaseModel::PopInfo const& popInfo_this = at( sweepModel->popInfos, pop_this );
			 popid pop_othr = at( pop2sib, pop_this );
			 BaseModel::PopInfo const& popInfo_othr = at( sweepModel->popInfos, pop_othr );
			 nchroms_float_t popsize_this = popInfo_this.popSizeFn( gen );
			 nchroms_float_t popsize_othr = popInfo_othr.popSizeFn( gen );

			 nchroms_float_t popsize_sel, popsize_uns;
			 Pop *curPop = demography->dg_get_pop_by_name( pop_this );
			 Pop *selPop, *unsPop;
			 if ( STLContains( selPops, pop_this ) ) {
				 popsize_sel = popsize_this;
				 popsize_uns = popsize_othr;
				 selPop = demography->dg_get_pop_by_name( pop_this );
				 unsPop = demography->dg_get_pop_by_name( pop_othr );
			 } else {
				 popsize_sel = popsize_othr;
				 popsize_uns = popsize_this;
				 selPop = demography->dg_get_pop_by_name( pop_othr );
				 unsPop = demography->dg_get_pop_by_name( pop_this );
			 }
			 
			 freq_t testFreq = popsize_sel  / ( popsize_sel + popsize_uns );
			 prob_t pval = randGen->random_double();
			 Pop * shouldBeInPop = ( pval < testFreq ) ? selPop : unsPop;
			 if ( shouldBeInPop != curPop ) {
				 curPop->pop_remove_node( node );
				 shouldBeInPop->pop_add_node( node );
			 }
		 }  // if ( node ) 
	 } // void determineAlleleAtSelPos_( Node *node, genid gen )

	 
	 
// *** Method: makeSweepModel - given an original BaseModel, construct a BaseModel for simulating sweeps.
// 	  For each original pop, the sweep model has two pops: one for chroms carrying the selected allele,
//    and one for those carrying the unselected allele.  
	 BaseModelP makeSweepModel( boost::shared_ptr<const BaseModel> baseModel,
															mpop_traj_t const& pop2freqSelFn, loc_t selPos,
															std::map< popid, popid >& pop2sib, 
															RandGenP randGen_
															
		 ) {
		 BaseModelP sweepModel = boost::make_shared<BaseModel>();
		 
		 int nextPopId = ToInt( baseModel->popInfos.rbegin()->first ) + 1;
		 
		 using namespace math;
		 using util::at;
//		 using boost::units::simplify_typename;
		 
		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 popid selPop( nextPopId++ );
			 popid unsPop( pi->first );
			 pop2sib[ selPop ] = unsPop;
			 pop2sib[ unsPop ] = selPop;

			 demography->dg_create_pop( selPop, demography->dg_get_pop_by_name( unsPop )->get_label() +
																	std::string( "_sel" ), genid(0.) );
			 nchroms_t sampSz = demography->find_pop_request( unsPop )->members;
			 using boost::algorithm::clamp;
			 //std::cerr << "unsPop=" << unsPop << " fn=" << util::at( pop2freqSelFn, unsPop ) << "\n";

			 freq_t selFreq = at( pop2freqSelFn, unsPop )( genid( 0. ) );
			 //std::cerr << "unsPop=" << unsPop << " selPop=" << selPop << " selFreq=" << selFreq << "\n";
			 boost::random::binomial_distribution<nchroms_t> bdist( sampSz, selFreq );

			 nchroms_t selSampleSize = bdist(*randGen_);

			 // std::cerr << "makeSweepModel: selPop=" << selPop << " unsPop=" << unsPop <<
			 // 		" sampSz=" << sampSz << " frac=" << util::at( pop2freqSelFn, unsPop )( genid( 0. ) ) << "\n";
			 
			 demography->dg_populate_by_name( selPop,
																				selSampleSize );
			 nchroms_t unsSampleSize( sampSz - selSampleSize );
			 demography->find_pop_request( unsPop )->members = unsSampleSize;
			 // std::cerr << "selPop=" << selPop << " unsPop=" << unsPop << " sampSz=" << sampSz << " selFreq=" << selFreq << 
			 // 	 " selSampleSize=" << selSampleSize << " unsSampleSize=" << unsSampleSize << "\n";
		 }
		 sweepHook = boost::make_shared<SweepHook>( this, selPos );
		 demography->addHook( sweepHook );
		 
		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 //Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
			 BOOST_AUTO( const& popInfo, pi->second );

			 popid unsPop( pi->first );
			 popid selPop( pop2sib[ unsPop ] );

			 BOOST_AUTO( & popInfoSel, sweepModel->popInfos[ selPop ] );
			 BOOST_AUTO( & popInfoUns, sweepModel->popInfos[ unsPop ] );

			 pop2sib.insert( std::make_pair( selPop, unsPop ) );
			 pop2sib.insert( std::make_pair( unsPop, selPop ) );

			 selPops.insert( selPop );

			 BOOST_AUTO( const& freqSelFn, util::at( pop2freqSelFn, unsPop ) );
			 for( BOOST_AUTO( it, freqSelFn.getPieces().rbegin() ); it != freqSelFn.getPieces().rend(); ++it )  {
				 genid gen = it->first;
				 freq_t freqSel = it->second( gen );
				 freq_t freqUns = freq_t(1.0) - freqSel;
				 popsize_float_t sz = popInfo.popSizeFn( gen );
				 popsize_float_t szSel = freqSel * sz;
				 popsize_float_t szUns = freqUns * sz;
				 popInfoSel.setSizeFrom( gen, szSel );
				 popInfoUns.setSizeFrom( gen, szUns );
			 }

			 genid selBegGen = freqSelFn.getPieces().begin()->first;
			 //std::cerr << "selBegGen=" << selBegGen << "\n";

			 {
				 BOOST_AUTO( lb, popInfo.popSizeFn.getPieces().lower_bound( selBegGen ) );
				 popInfoUns.popSizeFn.getPieces().insert( popInfo.popSizeFn.getPieces().begin(), lb );
				 popInfoUns.popSizeFn.getPieces().insert( std::make_pair( selBegGen, lb->second ) );
			 }

			 {
				 BOOST_AUTO( lb, popInfo.coalRateFn.getPieces().lower_bound( selBegGen ) );
				 popInfoUns.coalRateFn.getPieces().insert( popInfo.coalRateFn.getPieces().begin(), lb );
				 popInfoUns.coalRateFn.getPieces().insert( std::make_pair( selBegGen, lb->second ) );
			 }

			 popInfoUns.migrRateTo = popInfo.migrRateTo;
			 for( BOOST_AUTO( migr_it, popInfo.migrRateTo.begin() ); migr_it != popInfo.migrRateTo.end(); ++migr_it )
					popInfoSel.migrRateTo[ pop2sib[ migr_it->first ] ] = migr_it->second;

			 set( popInfoSel.migrRateTo[ unsPop ], selBegGen, prob_per_chrom_per_gen_t( 1 ) );
		 }

		 return sweepModel;
	 } // makeSweepModel


	 typedef double fitness_t;

// *** Func: simulateTrajFwd - simulate the trajectory of a selected allele in a set of populations.
//
// Params:
//   baseModel - the BaseModel on which we're simulating	 
//   fit - the absolute fitness of AA, Aa, aa where A is the selected allele
//   begGen - the generation at which the simulation begins.
//   begFreqs - for each pop, the freq of A at begGen
//   endFreqs - for each pop, the acceptable final freqs at generation 0
//   randGen - the uniform random number generator
//
//  Returns:
//   for each pop, the trajectory of the selected allele in that pop.
//
//  Note: some code adapted from simuPOP by Bo Peng et al
	 template <typename URNG>
	 static boost::shared_ptr<mpop_traj_t>
	 simulateTrajFwd( boost::shared_ptr<const BaseModel> baseModel,
										std::map<popid, std::map<genotype_t,double> > pop2fits,
										genid begGen, std::map<popid,freq_t> begFreqs,
										std::map<popid, util::ValRange<freq_t> > endFreqs,
										popid selBegPop, genid selBegGen, double selCoeff,
										URNG& urng,
										size_t maxAttempts = 1000000 ) {
		 cosi_using5( util::at, std::map, boost::assign::insert, boost::adaptors::map_keys,
									boost::range::push_back );
		 
		 boost::shared_ptr<mpop_traj_t> pop2freqSelFn( new mpop_traj_t );

		 std::vector<popid> pops;
		 push_back( pops, baseModel->popInfos | map_keys );
		 
		 // map<popid, map<genotype_t,double> >	pop2s;
		 // BOOST_FOREACH( popid pop, pops ) {
		 // 	 map<genotype_t,double> const& fit = at( fits, pop );
		 // 	 insert( pop2s[ pop ] )
		 // 			(GT_AA, 0.)
		 // 			(GT_Aa, at(fit,GT_Aa) / at(fit,GT_AA) - 1.)
		 // 			(GT_aa, at(fit,GT_aa) / at(fit,GT_AA) - 1.);
		 // }
		 
		 //std::cerr << "s=(" << s[0] << "," << s[1] << "," << s[2] << "\n";

		 std::ostringstream msgs;

		 bool foundTrajectory = false;
		 while( !foundTrajectory && maxAttempts-- >= 1 ) {
			 msgs.str("");
			 msgs.clear();
			 //std::cerr << "-------------\n";
			 pop2freqSelFn->clear();

			 //if ( !(maxAttempts % 100000) ) std::cerr << "attempts left=" << maxAttempts << "\n";
			 
			 BOOST_AUTO( freqs, begFreqs );
			 gens_t STEP(1.);
			 bool trajFailed = false;
			 bool tr = false; //( maxAttempts == ( getenv( "COSI_TRAJ_ID" ) ? ((size_t)atol( getenv( "COSI_TRAJ_ID" ) ) ) : 0 ) ) ;
			 for( genid gen = begGen; gen-STEP >= genid(0) && !trajFailed; gen -= STEP ) {
				 if ( tr ) PRINT2( "mygen", gen );
				 genid gen_next = gen - STEP;  // generations are numbered into the past, with present time being generation 0;
				 //                               we use this convention for both fwd and bwd sims.
				 bool haveNonZero = false;
				 BOOST_FOREACH( popid pop, pops ) {
					 if ( tr ) PRINT4( gen, gen_next, pop, freqs[pop] );
					 chk_freq( freqs[pop] );
					 
					 // Record the current freq of sel allele in this pop
					 if ( tr ) std::cerr << "befset: " << (*pop2freqSelFn)[ pop ] << "\n";
					 set( (*pop2freqSelFn)[ pop ], gen, freqs[ pop ] );
					 if ( tr ) std::cerr << "aftset: " << (*pop2freqSelFn)[ pop ] << "\n";
					 if ( tr ) std::cerr << "justset: gen=" << gen << " pop=" << pop << " freq=" << freqs[pop] << " f=" << (*pop2freqSelFn)[pop](gen) << "\n";
					 cosi_chk(  (*pop2freqSelFn)[pop](gen) == freqs[ pop ], "bug in piecewise fns" );
					 
					 if ( freqs[ pop ] > 0 ) haveNonZero = true;
				 }  // for each pop
				 if ( !haveNonZero ) { 
					 trajFailed = true; 
				 }
				 else 
					 { // if !trajFailed
						 // for each pop, determine n_A and n_a after random mating and selection, but before migration
						 std::map<popid, freq_t> b4mig_p_A;
						 cosi_for_map_keys( pop, pop2fits ) {
							 double s = ( (pop == selBegPop) && (gen <= selBegGen) ) ? selCoeff : 0.;
							 double fit_AA = 1.0 + s, fit_Aa = 1.0 + 0.5*s, fit_aa = 1.0;
							 //	 double fit_AA = at( fits, GT_AA ), fit_Aa = at( fits, GT_Aa ), fit_aa = at( fits, GT_aa );
							 freq_t p_A = freqs[ pop ];
							 freq_t p_a = 1. - p_A;
							 double amt_A = p_A * p_A * fit_AA + p_A * p_a * fit_Aa;
							 double amt_a = p_a * p_a * fit_aa + p_A * p_a * fit_Aa;
							 b4mig_p_A[ pop ] = amt_A / ( amt_A + amt_a );
							 if ( tr ) PRINT9( pop, fit_AA, fit_Aa, fit_aa, p_A, p_a, amt_A, amt_a, b4mig_p_A[pop] );
						 } cosi_end_for;

						 // now model migration.
						 // code below influenced by  msms simulator by Ewing and Hermisson, http://bioinformatics.oxfordjournals.org/content/suppl/2010/06/20/btq322.DC1/InternalManual.pdf
						 cosi_for_map( pop, popInfo, baseModel->popInfos ) {
							 frac_t nonMigFrac(1.0);
							 freq_t p_A(0);
							 cosi_for_map( srcPop, migrRateFn, popInfo.migrRateTo ) {
								 frac_t migFrac_from_srcPop = migrRateFn( gen_next ) * STEP * nchroms_float_t(1);
								 p_A += migFrac_from_srcPop * at( b4mig_p_A, srcPop );
								 nonMigFrac -= migFrac_from_srcPop;
							 } cosi_end_for;
							 p_A += nonMigFrac * at( b4mig_p_A, pop );
							 
							 // genetic drift
							 nchroms_float_t N = popInfo.popSizeFn( gen_next );
							 boost::random::binomial_distribution<nchroms_t> bdist( 2 * nchroms_t( ToDouble( N ) ), p_A );
							 nchroms_t nsel_next_gen = bdist( urng );
							 if ( tr ) PRINT5( 2*N, p_A, nsel_next_gen, bdist.param(), bdist );
							 freqs[ pop ] = nchroms_float_t( ToDouble( nsel_next_gen ) ) / ( 2 * N);
							 if ( tr ) std::cerr << "gen=" << gen << " pop=" << pop << " 2*N=" << (2*N) << " p_A=" << p_A  << 
													 " nsel_next_gen=" << nsel_next_gen << " nonMigFrac=" << nonMigFrac << " p_A'=" << freqs[pop] << "\n";

						 } cosi_end_for;  // cosi_for_map( pop, popInfo, baseModel->popInfos )
					 } // if !trajFailed

					 // // Find the sel freq in pop at time gen-1
					 // freqs[ pop ] =
					 // 		getNextXt( freqs[ pop ],
					 // 							 at( baseModel->popInfos, pop ).popSizeFn( gen ), at( pop2s, pop ),
					 // 							 urng );
				 
				 // migrations; note that the direction is reversed for the fwd vs the bwd sim.
			 } // for each gen

			 if ( !trajFailed ) {
				 bool freqWrong = false;
				 BOOST_FOREACH( popid pop, pops ) {
					 //msgs << "EF pop=" << pop << " EF=" << freqs[pop] << "\n";
					 if ( endFreqs[ pop ]( freqs[ pop ] ) )
						 set( (*pop2freqSelFn)[ pop ], genid(0.), freqs[ pop ] );
					 else
						 freqWrong = true;
				 }
				 if ( !freqWrong ) {
					 foundTrajectory = true;
					 using util::operator<<;
					 using math::operator<<;
					 //std::cerr << "got traj! maxAttempts=" << maxAttempts << " freqs=" << freqs << "\n";

					 // BOOST_FOREACH( popid pop, pops ) {
					 // 	 std::cerr << "EF pop=" << pop << " EF=" << freqs[pop] << " fn=" << at(*pop2freqSelFn, pop) <<  "\n";
					 // }

					 //std::cerr << "GOT!\n" << msgs.str() << "\n";
				 }
			 }  // if ( !trajFailed )
		 }  // while traj not found
		 if ( foundTrajectory )
			 return pop2freqSelFn;
		 else
			 BOOST_THROW_EXCEPTION( cosi::cosi_error() << error_msg( "no trajectory found within given number of attempts" ) );
	 }  // simulateTrajFwd

	 static bool readTsvLine( std::istream& s, std::vector< std::string >& vec, size_t& lineNo,
														std::streampos *savePos = NULL ) {
		 cosi_using2(boost::tokenizer,boost::char_separator);
		 typedef tokenizer< char_separator<char> > Tokenizer;
		 std::string line;
		 vec.clear();
		 if ( savePos ) *savePos = s.tellg();
		 if ( !s || !std::getline( s, line ) ) return false;
		 else {
			 //std::cerr << "got line: " << line << " lineNo was " << lineNo << "\n";
			 ++lineNo;
			 char_separator<char> sep("\t");
			 Tokenizer tok(line,sep);
			 vec.assign(tok.begin(),tok.end());
			 return true;
		 }
	 }

	 static boost::shared_ptr<mpop_traj_t>
	 loadTraj( filename_t fname ) {
		 cosi_using5(std::vector,std::string,util::chkCond,boost::algorithm::starts_with,boost::filesystem::ifstream);
		 cosi_using2(boost::lexical_cast,boost::bad_lexical_cast);
		 size_t lineNo = 0;
		 try {
			 try {
				 static ifstream is;
				 static unsigned simId = 0;
				 static vector<popid> col2pop(2,NULL_POPID);
				 static boost::shared_ptr<mpop_traj_t> pop2freqSelFn;

				 bool repeatTraj = (bool)getenv( "COSI_REPEAT_TRAJ" );


				 ++simId;
				 lineNo=0;
				 if ( simId == 1 ) {
					 is.open( fname );
					 vector<string> titles;
					 readTsvLine( is, titles, lineNo );
					 //std::cerr << "got " << titles.size() << "titles.\n";
					 if ( !( titles.size() > 2 ) ) BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg("too few columns" ) );
					 if ( !( titles[0]=="sim" ) ) BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg("first col not sim" ) );
					 if ( !( titles[1]=="gen" ) ) BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg("2nd col not gen" ) );

					 //chkCond( (titles.size() > 2) && (titles[0]=="sim") && titles[1]=="gen", "bad traj file header" );
					 for( size_t i=2; i<titles.size(); ++i ) {
						 chkCond( starts_with( titles[i], "selfreq_" ), "bad traj file header" );
						 col2pop.push_back( lexical_cast<popid>( titles[i].substr( strlen("selfreq_") ) ) );
						 //std::cerr << "col2pop.size=" << col2pop.size() << " back=" << col2pop.back() << "\n";
					 }
				 } else {
					 if ( repeatTraj ) {
						 //std::cerr << "returning repeatTraj\n";
						 return pop2freqSelFn;
					 }
				 }
				 pop2freqSelFn.reset( new mpop_traj_t );
				 vector<string> vals;
				 genid lastGen(NULL_GEN);
				 bool isFirst = true;
				 std::streampos is_pos(0);
				 std::string simIdStr = lexical_cast<string>( simId );
				 unsigned linesRead = 0;
				 while ( readTsvLine( is, vals, lineNo, &is_pos ) ) {
					 chkCond( vals.size() == col2pop.size(), "bad traj file line" );
					 if ( vals[0] != simIdStr ) BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg("wrong sim id") );
					 genid gen = lexical_cast<genid>( vals[1] );
					 //std::cerr << "gen=" << gen << " vals[2]=" << vals[2] << "\n";
					 for( size_t i=2; i<vals.size(); ++i )
							set( (*pop2freqSelFn)[ col2pop[i] ], gen, lexical_cast<freq_t>(vals[i]) );
					 chkCond( isFirst || (gen < lastGen), "generations do not decrease" );
					 
					 isFirst = false;
					 lastGen = gen;
					 ++linesRead;
					 if ( gen < genid(1) ) {
						 if ( gen != genid(0) ) {
							 for( size_t i=2; i<vals.size(); ++i )
									set( (*pop2freqSelFn)[ col2pop[i] ], genid(0), lexical_cast<freq_t>(vals[i]) );
						 }
						 break;
					 }
				 }
				 if ( linesRead < 1 ) BOOST_THROW_EXCEPTION( cosi_io_error() << error_msg("traj too short" ) );
				 
				 return pop2freqSelFn;
			 } cosi_pkg_exception3( cosi_io_error(), ifstream::failure, bad_lexical_cast, std::exception );
		 } catch( boost::exception& e ) {
			 e << boost::errinfo_file_name( fname.string() )
				 << boost::errinfo_at_line( lineNo )
				 << error_stage( "reading trajectory file" );
			 throw;
		 }
		 
	 }


	 // Fn: defineParams - defines command-line params relevant to this module
	 virtual void defineParams( boost::program_options::options_description& opts ) {
		 namespace po = boost::program_options;
		 po::options_description sweep_opts( "Sweep options" );

		 sweep_opts.add_options()
				( "sweep.max-fwd-traj-attempts", po::value(&maxFwdTrajAttempts)->default_value(1000000),
					"max number of attempts to generate a forward trajectory during sweep simulations" );
		 
		 opts.add( sweep_opts );
	 }
	 
private:
	 boost::uint32_t maxFwdTrajAttempts;
	 
// *** Private methods
	 
// **** Fn: getNextXt - given the freq of selected allele in a pop at gen g, determine the freq at gen g-1.
// Params:
//    x - freq of causal allele at time g
//    Nt - pop size (number of diploids) at time g
//    s - the relative fitness of AA, Aa, aa (where A is the selected allele).
//    urng - a uniform random number generator	 
// Returns:
//    frequency of selected allele at ge-1
// Note: code adopted from Bo Peng's simuPOP simulator	 
	 template <typename URNG>
	 static freq_t getNextXt( freq_t x, popsize_float_t Nt, std::map<genotype_t,double> const& s, URNG& urng ) {
		 using util::at;
     // if the selected allele has already been either lost or fixed in the population, it stays that way
		 if ( x == 0 || x == 1 ) return x;
		 
		 double s_Aa = at( s, GT_Aa );
		 double s_aa = at( s, GT_aa );
		 // with s1 and s2 on hand, calculate freq at the next generation
		 double num = x * (1. + s_aa * x + s_Aa * (1. - x));
		 double denom = (1. + s_aa * x * x + 2 * s_Aa * x * (1. - x));
		 freq_t y =  num/ denom;
		 // y is obtained, is the expected allele frequency for the next generation t+1
		 boost::random::binomial_distribution<nchroms_t> bdist( 2 * nchroms_t( ToDouble( Nt ) ), y );
		 nchroms_t nsel_next_gen = bdist( urng );
		 // std::cerr << "x=" << x << " Nt=" << Nt << " num=" << num << " denom=" << denom << " y=" << y
		 // 					 << " nsel_next_gen=" << nsel_next_gen << "\n";

		 return ToDouble( nsel_next_gen ) / ToDouble( 2*Nt );
	 }
};  // class MSweep

MSweepP make_MSweep( DemographyP demography, BaseModelP baseModel, RandGenP randGen ) {
	return boost::make_shared<MSweep>( demography, baseModel, randGen );
}

ModuleP as_Module( MSweepP msweep ) { return msweep; }

BaseModelP getSweepModel( MSweepP msweep ) { return msweep->sweepModel; }

LeavesInfoP
computeLeavesInfo( MSweepP msweep ) {
	using util::STLContains;
	LeavesInfoP leavesInfo;
	BaseModelP sweepModel = msweep->sweepModel;
	DemographyP demography = msweep->demography;
	if ( msweep->baseModel->sweepInfo ) {
		leavesInfo = boost::make_shared< LeavesInfo >();
		std::vector< leafset_p > const& pop2leaves = demography->get_pop2leaves();
		std::vector< bool > leavesUsed( leafset_get_max_leaf_id() );
		cosi_for_map_keys( pop, sweepModel->popInfos ) {
			//std::cerr << "computeLeafOrder: pop=" << pop << "\n";
			if ( STLContains( msweep->selPops, pop ) ) {
				popid pop_sel = pop;
				popid pop_uns = util::at( msweep->pop2sib, pop_sel );
				leafset_p leaves_sel = pop2leaves[ demography->dg_get_pop_index_by_name( pop_sel ) ];
				leafset_p leaves_uns = pop2leaves[ demography->dg_get_pop_index_by_name( pop_uns ) ];
				//std::cerr << "leaves_sel=" << leaves_sel << " leaves_uns=" << leaves_uns << "\n";
				msweep->selLeaves = leafset_union( msweep->selLeaves, leaves_sel );
			
				nchroms_t nleavesBef = leavesInfo->leafOrder.size();
				COSI_FOR_LEAFSET( leaves_uns, leaf, {
						chkCond( (0<=leaf) && (((size_t)leaf)<leavesUsed.size()) && !leavesUsed[ leaf ], "leaf collision" );
						leavesInfo->leafOrder.push_back( leaf );
						leavesUsed[ leaf ] = true;
					});
				COSI_FOR_LEAFSET( leaves_sel, leaf, {
						chkCond(  (0<=leaf) && (((size_t)leaf)<leavesUsed.size()) && !leavesUsed[ leaf ], "leaf collision" );
						leavesInfo->leafOrder.push_back( leaf );
						leavesUsed[ leaf ] = true;
					});
				leavesInfo->sampleSizes.push_back( leavesInfo->leafOrder.size() - nleavesBef );
				leavesInfo->popNames.push_back( pop_uns );
			}  // if ( STLContains( msweep->selPops, pop ) )
			else {
				//std::cerr << "pop not in selPops\n";
			}
		} cosi_end_for;  // cosi_for_map( popInfo, sweepModel->popInfos )
		chkCond( boost::accumulate( leavesUsed, true, std::logical_and<bool>() ), "some leaf unused" );
	}
	return leavesInfo;
}  // computeLeafOrder()

void addSelMut( MSweepP msweep, MutlistP muts ) {
	// std::cerr << "addSelMut: adding\n";
	// std::cerr << "sweepInfo null? " << ( !msweep->baseModel || !msweep->baseModel->sweepInfo ) << "\n";
	// std::cerr << "selLeaves size: " << leafset_size( msweep->selLeaves ) << "\n";
	if ( msweep->baseModel->sweepInfo && !leafset_is_empty( msweep->selLeaves ) ) {
		SweepInfo const& si = *msweep->baseModel->sweepInfo;
		muts->addMut( si.selPos, msweep->selLeaves, si.selGen, si.selPop );
	} else {
		//std::cerr << "addSelMut: NOT adding\n";
	}
}



}  // namespace cosi

