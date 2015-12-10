#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

#define COSI_DEV_PRINT

#include <map>
#include <utility>
#include <iostream>
#include <boost/units/detail/utility.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/generalmath.h>
#include <cosi/basemodel.h>

namespace cosi {

class MSweep {
public:
	 BaseModelP getSweepModel( boost::shared_ptr<const BaseModel> baseModel,
														 std::map< popid, math::Function< genid, freq_t, math::Piecewise< math::Const<> > > >
														 pop2freqSelFn ) {
		 BaseModelP sweepModel = boost::make_shared<BaseModel>();
		 
		 int nextPopId = ToInt( baseModel->popInfos.rbegin()->first ) + 1;

		 using namespace math;
		 using boost::units::simplify_typename;

		 for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					pi != baseModel->popInfos.end(); ++pi ) {
			 //Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
			 BOOST_AUTO( const& popInfo, pi->second );

			 popid selPop( nextPopId++ );
			 popid unsPop( pi->first );

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

			 // std::cerr << "freqUns=" << freqUns << "\n";
			 // PRINT2( freqUns, popInfo.popSizeFn );
			 // PRINT( freqUns * popInfo.popSizeFn );
			 // PRINT( simplify_typename( freqUns * popInfo.popSizeFn ) );
			 
			 //popInfoUns.popSizeFn = freqUns * popInfo.popSizeFn;
			 // popInfoUns.coalRateFn = (cval( 1. ) / freqUns) * popInfo.coalRateFn;
			 // popInfoSel.popSizeFn = freqSel * popInfo.popSizeFn;
			 // popInfoSel.coalRateFn = (cval( 1. ) / freqSel) * popInfo.coalRateFn;
			 
			 // for( BOOST_AUTO( mi, popInfo.migrRateTo.begin() );
			 // 			mi != popInfo.migrRateTo.end(); ++mi ) {
			 // 	 Pop *dstPop = demography->dg_get_pop_by_name( mi->first );
			 
			 
			 // }
		 }

		 return sweepModel;
	 } // getSweepModel
	 
};  // MSweep

}  // namespace cosi



#endif // #ifndef COSI_INCLUDE_MSWEEP_H
