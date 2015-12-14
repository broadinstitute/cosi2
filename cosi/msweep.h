#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

//#define COSI_DEV_PRINT

#include <map>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <boost/range/adaptor/map.hpp>  
//#include <boost/units/detail/utility.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/exception/all.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/clamp.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/utils.h>
#include <cosi/cosirand.h>
#include <cosi/generalmath.h>
#include <cosi/basemodel.h>
#include <cosi/node.h>
#include <cosi/hooks.h>

namespace cosi {

// * class MSweep - implementation of selective sweep in multiple populations
//
// Note: some code here is adapted from simuPOP simulator by Bo Peng et al ; see
// http://simupop.sourceforge.net/manual_svn/build/userGuide_ch7_sec2.html	 
class MSweep {
public:

// ** Type: pop_traj_t - frequency trajectory of an allele in one pop	 
	 typedef math::Function< genid, freq_t, math::Piecewise< math::Const<> > > pop_traj_t;

// ** Type: mpop_traj_t - frequency trajectory of an allele in multiple pops
	 typedef std::map< popid, pop_traj_t > mpop_traj_t;
	 

	 class SweepHook: public Hook {
			DemographyP demography;
			BaseModelP sweepModel;
			RandGenP randGen;
			loc_t selPos;

	 public:
			SweepHook( DemographyP demography_, BaseModelP sweepModel_, RandGenP randGen_, loc_t selPos_ ):
				demography( demography_ ), sweepModel( sweepModel_ ), randGen( randGen_ ), selPos( selPos_ )  { }
			
			virtual void handle_recomb( Node *node1, Node *node2, loc_t loc, genid curGen ) {
				determineAlleleAtSelPos_( selPos < loc ? node2 : node1, curGen );
			}
			virtual void handle_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2,
															genid curGen ) {
				determineAlleleAtSelPos_( loc1 <= selPos && selPos <= loc2 ? node2 : node1,
																	curGen );
			}

	 private:

			void determineAlleleAtSelPos_( Node *node, genid gen ) {
				using util::at;
				if ( node ) {
					popid pop_this = node->getPop()->pop_get_name();
					BaseModel::PopInfo const& popInfo_this = at( sweepModel->popInfos, pop_this );
					popid pop_othr = at( sweepModel->pop2sib, pop_this );
					BaseModel::PopInfo const& popInfo_othr = at( sweepModel->popInfos, pop_othr );
					nchroms_float_t popsize_this = popInfo_this.popSizeFn( gen );
					nchroms_float_t popsize_othr = popInfo_othr.popSizeFn( gen );

					nchroms_float_t popsize_sel, popsize_uns;
					Pop *curPop = demography->dg_get_pop_by_name( pop_this );
					Pop *selPop, *unsPop;
					if ( popInfo_this.isSelPop ) {
						popsize_sel = popsize_this;
						popsize_uns = popsize_othr;
						selPop = demography->dg_get_pop_by_name( pop_this );
						unsPop = demography->dg_get_pop_by_name( pop_othr );
					} else {
						popsize_sel = popsize_this;
						popsize_uns = popsize_othr;
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
	 };  // class SweepHook

	 
	 
// ** Method: getSweepModel - given an original BaseModel, construct a BaseModel for simulating sweeps.
// 	  For each original pop, the sweep model has two pops: one for chroms carrying the selected allele,
//    and one for those carrying the unselected allele.  
BaseModelP getSweepModel( boost::shared_ptr<const BaseModel> baseModel,
													mpop_traj_t const& pop2freqSelFn, DemographyP demography, loc_t selPos ) {
		 BaseModelP sweepModel = boost::make_shared<BaseModel>();

		 int nextPopId = ToInt( baseModel->popInfos.rbegin()->first ) + 1;
		 
		 using namespace math;
//		 using boost::units::simplify_typename;
		 
		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 popid selPop( nextPopId++ );
			 popid unsPop( pi->first );
			 sweepModel->pop2sib[ selPop ] = unsPop;
			 sweepModel->pop2sib[ unsPop ] = selPop;

			 demography->dg_create_pop( selPop, demography->dg_get_pop_by_name( unsPop )->get_label() +
																	std::string( "_sel" ), genid(0.) );
			 nchroms_t sampSz = demography->find_pop_request( unsPop )->members;
			 using boost::algorithm::clamp;
			 nchroms_t selSampleSize( clamp( int( util::at( pop2freqSelFn, unsPop )( genid( 0. ) ) *
																						ToDouble( sampSz ) ), 0, sampSz ) );
			 
			 demography->dg_populate_by_name( selPop,
																				selSampleSize );
			 nchroms_t unsSampleSize( sampSz - selSampleSize );
			 demography->find_pop_request( unsPop )->members = unsSampleSize;
			 //std::cerr << "selSampleSize=" << selSampleSize << " unsSampleSize=" << unsSampleSize << "\n";
		 }
		 boost::shared_ptr<Hook> hookPtr( new SweepHook( demography, sweepModel, demography->getRandGen(),
																										 selPos ) );
		 demography->addHook( hookPtr );
		 
		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 //Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
			 BOOST_AUTO( const& popInfo, pi->second );

			 popid unsPop( pi->first );
			 popid selPop( sweepModel->pop2sib[ unsPop ] );

			 BOOST_AUTO( & popInfoSel, sweepModel->popInfos[ selPop ] );
			 BOOST_AUTO( & popInfoUns, sweepModel->popInfos[ unsPop ] );

			 sweepModel->pop2sib.insert( std::make_pair( selPop, unsPop ) );
			 sweepModel->pop2sib.insert( std::make_pair( unsPop, selPop ) );

			 popInfoSel.isSelPop = true;
			 popInfoUns.isSelPop = false;

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
					popInfoSel.migrRateTo[ sweepModel->pop2sib[ migr_it->first ] ] = migr_it->second;

			 set( popInfoSel.migrRateTo[ unsPop ], selBegGen, prob_per_chrom_per_gen_t( .99999 ) );
		 }

		 return sweepModel;
	 } // getSweepModel


	 // adapted from simuPOP

	 typedef double fitness_t;

// ** Func: simulateTrajFwd - simulate the trajectory of a selected allele in a set of populations.
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
	 template <typename URNG>
	 mpop_traj_t
	 simulateTrajFwd( boost::shared_ptr<const BaseModel> baseModel, const double fit[3],
										genid begGen, std::map<popid,freq_t> begFreqs,
										std::map<popid, std::pair<freq_t,freq_t> > endFreqs,
										URNG& urng,
										size_t maxAttempts = 1000000 ) {
		 using util::at;
		 
		 mpop_traj_t pop2freqSelFn;
		 
		 double s[3] = { 0., fit[1] / fit[0] - 1., fit[2] / fit[0] - 1. };
		 
		 //std::cerr << "s=(" << s[0] << "," << s[1] << "," << s[2] << "\n";
		 using boost::adaptors::map_keys;
		 using boost::range::push_back;
		 
		 std::vector<popid> pops;
		 push_back( pops, baseModel->popInfos | map_keys );
		 
		 bool found = false;
		 while( !found && maxAttempts-- >= 1 ) {
			 
			 BOOST_AUTO( freqs, begFreqs );
			 namespace rng = boost::range;
			 std::cerr.precision(8);
			 for( genid gen = begGen; gen > genid(0); gen -= gens_t(1) ) {
				 // record the current freqs
				 bool haveNonZero = false;
				 BOOST_FOREACH( popid pop, pops ) {
					 //PRINT3( gen, pop, freqs[pop] );
					 set( pop2freqSelFn[ pop ], gen, freqs[ pop ] );
					 freqs[ pop ] = getNextXt( freqs[ pop ], at( baseModel->popInfos, pop ).popSizeFn( gen ), s, urng );
					 if ( freqs[ pop ] > 0 ) haveNonZero = true;
				 }
				 if ( !haveNonZero ) break;
				 
				 // migrations; note that the direction is reversed for the fwd vs the bwd sim.
				 cosi_for_map( dstPop, popInfoDst, baseModel->popInfos ) {
					 cosi_for_map( srcPop, migrRateFn, popInfoDst.migrRateTo ) {
						 if ( freqs[ srcPop ] > 0 ) {
							 // find the number of chroms bearing the selected allele in the source population
							 // what if fixed or lost?
							 // check also for small pop size.
							 BaseModel::PopInfo const& popInfoSrc = at( baseModel->popInfos, srcPop );
							 nchroms_t NtSrc( 2 * ToDouble( popInfoSrc.popSizeFn( gen ) ) );
							 nchroms_t NtDst( 2 * ToDouble( popInfoDst.popSizeFn( gen ) ) );
							 nchroms_t nSrcSel( 2 * NtSrc * freqs[ srcPop ] );
							 nchroms_t nDstSel( 2 * NtDst * freqs[ dstPop ] );
							 double migrRate = ToDouble(  migrRateFn( gen ) );
							 if ( nSrcSel > nchroms_t(1) && migrRate > 0. ) {
								 
								 nchroms_t nSrcUns = 2 * NtSrc - nSrcSel;
								 
								 // so, for each src chrom, there's a chance of migration; the number of migrators
								 // is then binomial.  note that we're treating this as a collection of haploids;
								 // more correctly might be to assume it's randomly mixing and see how many diploids
								 // of each type migrate?  on the other hand these are probs of haploid chroms migrating?
								 // so, if there's a fixed chance of a diploid _individual_ migrating, then
								 // a chrom's chance of migrating is : there are p*p chance of being in an aa indiv.
								 // would we get the same migration rates?
								 
								 using boost::random::binomial_distribution;
								 binomial_distribution<nchroms_t> bdistSel( nSrcSel, migrRate );
								 binomial_distribution<nchroms_t> bdistUns( nSrcUns, migrRate );
								 nchroms_t nMigSel = bdistSel( urng );
                 nchroms_t nMigUns = bdistUns( urng );
								 
								 //PRINT9( gen, srcPop, dstPop, NtSrc, NtDst, nSrcSel, nDstSel, nMigSel, nMigUns );
								 
                 nchroms_t nMig = nMigSel + nMigUns;
								 freqs[ srcPop ] = double( nSrcSel - nMigSel ) / double( 2*NtSrc - nMig );
								 freqs[ dstPop ] = double( nDstSel + nMigSel ) / double( 2*NtDst + nMig );
							 }  // if ( nSrcSel > nchroms_t(1) && migrRate > 0. )
						 } // if ( freqs[ srcPop ] > 0 )
					 } cosi_end_for;  // cosi_for_map( srcPop, migrRateFn, popInfo.migrRateTo )
				 } cosi_end_for;  // cosi_for_map( trgPop, popInfo, baseModel->popInfos )
			 } // for each gen

			 bool freqWrong = false;
			 BOOST_FOREACH( popid pop, pops ) {
				 if ( ( endFreqs[ pop ].first <= freqs[ pop ] ) && ( freqs[ pop ] <= endFreqs[ pop ].second ) )
						set( pop2freqSelFn[ pop ], genid(0.), freqs[ pop ] );
				 else
						freqWrong = true;
			 }
			 if ( !freqWrong )
					return pop2freqSelFn;
		 }
		 BOOST_THROW_EXCEPTION( cosi::cosi_error() << error_msg( "no traj found within given number of attempts" ) );
	 }  // simulateTrajFwd
	 
private:
// ** Private methods
	 
// *** PrivMethod: getNextXt - given the current freq of selected allele in a pop, compute next-generation freq.
// Params:
//    x - current freq on causal allele 
//    Nt - current pop size (number of diploids).
//    s - the relative fitness of AA, Aa, aa (where A is the selected allele).
// Returns:
//    frequency of selected allele in the next generation.
	 template <typename URNG>
	 static freq_t getNextXt( freq_t x, popsize_float_t Nt, const double s[3], URNG& urng ) {
     // if current allele freq in subpop sp at locus loc has already been 0 or 1,
     // set it to be 0 or 1 for next gens
		 if ( x == 0 || x == 1 ) return x;
		 freq_t s1 = s[1];
		 freq_t s2 = s[2];
		 // with s1 and s2 on hand, calculate freq at the next generation
		 double num = x * (1. + s2 * x + s1 * (1. - x));
		 double denom = (1. + s2 * x * x + 2 * s1 * x * (1. - x));
		 freq_t y =  num/ denom;
		 // y is obtained, is the expected allele frequency for the next generation t+1
		 boost::random::binomial_distribution<nchroms_t> bdist( 2 * nchroms_t( ToDouble( Nt ) ), y );
		 nchroms_t nsel_next_gen = bdist( urng );
		 // std::cerr << "x=" << x << " Nt=" << Nt << " num=" << num << " denom=" << denom << " y=" << y
		 // 					 << " nsel_next_gen=" << nsel_next_gen << "\n";

		 return ToDouble( nsel_next_gen ) / ToDouble( 2*Nt );
	 }

	 
};  // class MSweep

}  // namespace cosi



#endif // #ifndef COSI_INCLUDE_MSWEEP_H
