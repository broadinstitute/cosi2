#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

#define COSI_DEV_PRINT

#include <map>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <boost/units/detail/utility.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/exception/all.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/utils.h>
#include <cosi/cosirand.h>
#include <cosi/generalmath.h>
#include <cosi/basemodel.h>

namespace cosi {

// * class MSweep - implementation of selective sweep in multiple populations
//
// Note: some code here is adapted from simuPOP simulator by Bo Peng et al ; see
// http://simupop.sourceforge.net/manual_svn/build/userGuide_ch7_sec2.html	 
class MSweep {
public:
	 
// ** Method: getSweepModel - given an original BaseModel, construct a BaseModel for simulating sweeps.
// 	  For each original pop, the sweep model has two pops: one for chroms carrying the selected allele,
//    and one for those carrying the unselected allele.  
	 BaseModelP getSweepModel( boost::shared_ptr<const BaseModel> baseModel,
														 std::map< popid, math::Function< genid, freq_t, math::Piecewise< math::Const<> > > >
														 pop2freqSelFn ) {
		 BaseModelP sweepModel = boost::make_shared<BaseModel>();
		 
		 int nextPopId = ToInt( baseModel->popInfos.rbegin()->first ) + 1;

		 using namespace math;
		 using boost::units::simplify_typename;

		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 popid selPop( nextPopId++ );
			 popid unsPop( pi->first );
			 sweepModel->pop2sib[ selPop ] = unsPop;
			 sweepModel->pop2sib[ unsPop ] = selPop;
		 }

		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 //Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
			 BOOST_AUTO( const& popInfo, pi->second );

			 popid unsPop( pi->first );
			 popid selPop( sweepModel->pop2sib[ unsPop ] );

			 BOOST_AUTO( & popInfoSel, sweepModel->popInfos[ selPop ] );
			 BOOST_AUTO( & popInfoUns, sweepModel->popInfos[ unsPop ] );

			 popInfoSel.sweepSibPop = unsPop;
			 popInfoUns.sweepSibPop = selPop;

			 BOOST_AUTO( const& freqSelFn, pop2freqSelFn[ unsPop ] );
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
		 }

		 return sweepModel;
	 } // getSweepModel


	 // adapted from simuPOP

	 typedef double fitness_t;

// ** Type: pop_traj_t - frequency trajectory of an allele in one pop	 
	 typedef math::Function< genid, freq_t, math::Piecewise< math::Const<> > > pop_traj_t;

// ** Type: mpop_traj_t - frequency trajectory of an allele in multiple pops
	 typedef std::map< popid, pop_traj_t > mpop_traj_t;

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
	 simulateTrajFwd( BaseModelP baseModel, const double fit[3],
										genid begGen, std::map<popid,freq_t> begFreqs,
										std::map<popid, std::pair<freq_t,freq_t> > endFreqs,
										URNG& urng,
										size_t maxAttempts = 10000 ) {
		 mpop_traj_t pop2freqSelFn;

		 double s[3] = { 0., fit[1] / fit[0] - 1., fit[2] / fit[0] - 1. };

		 bool found = false;
		 while( !found && maxAttempts-- >= 1 ) {

			 BOOST_AUTO( freqs, begFreqs );
			 namespace rng = boost::range;
			 for( genid gen = begGen; gen > genid(0); gen -= gens_t(1) ) {
				 // record the current freqs
				 cosi_for_map( pop, popInfo, baseModel->popInfos ) {
					 set( pop2freqSelFn[ pop ], gen, freqs[ pop ] );
					 freqs[ pop ] = getNextXt( freqs[ pop ], baseModel->popInfos[ pop ].popSizeFn( gen ), s, urng );
				 } cosi_end_for;
			 }

			 bool freqWrong = false;
			 cosi_for_map( pop, popInfo, baseModel->popInfos ) {
				 if ( ( endFreqs[ pop ].first <= freqs[ pop ] ) && ( freqs[ pop ] <= endFreqs[ pop ].second ) )
						set( pop2freqSelFn[ pop ], genid(0.), freqs[ pop ] );
				 else
						freqWrong = true;
			 } cosi_end_for;
			 if ( !freqWrong )
					return pop2freqSelFn;
		 }
		 BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "no traj found within given number of attempts" ) );
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
		 std::cerr << "x=" << x << " Nt=" << Nt << " num=" << num << " denom=" << denom << " y=" << y << "\n";
		 // y is obtained, is the expected allele frequency for the next generation t+1
		 boost::random::binomial_distribution<nchroms_t> bdist( 2 * nchroms_t( ToDouble( Nt ) ), y );
		 nchroms_t nsel_next_gen = bdist( urng );

		 return ToDouble( nsel_next_gen ) / ToDouble( 2*Nt );
	 }
};  // class MSweep

}  // namespace cosi



#endif // #ifndef COSI_INCLUDE_MSWEEP_H
