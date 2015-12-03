#ifndef INCLUDE_COSI_BASEMODEL_H
#define INCLUDE_COSI_BASEMODEL_H

#include <utility>
#include <vector>
#include <map>
#include <iostream>
#include <boost/typeof/typeof.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/generalmath.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {

// * Class BaseModel - demographic model expressed in the most basic terms.
// 
// A representation of the demographic model in terms of just a few primitive elements:
//   - a set of pops
//   - for each pop, its size at each generation
//   - for each pop pair, the migration rate between the pop pair at each time point.
struct BaseModel {


// *** Class PopInfo - information about one population
	 struct PopInfo {
			
			// Field: popSizeFn - the size of this pop at each gen
			math::Function< genid, popsize_float_t, math::Piecewise< math::Any<> > > popSizeFn;

			// Field: coalRateIntegralFn - the integral of the coalescence rate, equal to
			// indefiniteIntegral( 1/2*popSizeFn )
			math::Function< genid, popsizeInv_float_integral_t, math::Piecewise< math::Any<> > > coalRateIntegralFn;

			// Field: migrRateTo - for each target pop, migration rate from this pop to the target pop.
			std::map< popid, math::Function< genid, prob_per_chrom_per_gen_t,
																			 math::Piecewise< math::Any<> > > > migrRateTo;

			template <typename TSpec>
			void setSizeFrom( genid fromGen, math::Function< genid, popsize_float_t, TSpec> const& f ) {
				using namespace math;
				Function< genid, double, math::Const<> > c_1( 1.0 ) ;
				Function< genid, double, math::Const<> > c_2( 2.0 ) ;
				std::cerr << "setSizeFrom: f=" << f << "\n";
				std::cerr << "setSizeFrom: coalrate=" << ( c_1 / ( c_2 * f ) ) << "\n";
				popSizeFn.getPieces()[ fromGen ] = f;
				// coalRateIntegralFn.getPieces()[ fromGen ] =
				// 	 indefiniteIntegral( Function< genid, double, Const<> >( 0.5 ) *
				// 											 pow( f, math::Function< genid, double, math::Const<> >( -1. ) ) );
			}

			void setSizeFrom( genid fromGen, popsize_float_t sz ) {
				setSizeFrom( fromGen, math::fn_const< genid >( sz ) );
			}

			bool isSizeSetFrom( genid gen ) const { return popSizeFn.getPieces().count( gen ) > 0; }

			void setMigrRate( popid dstPop, genid fromGen, prob_per_chrom_per_gen_t rate ) {
				BOOST_AUTO( &pieces, migrRateTo[ dstPop ].getPieces() );
				if ( pieces.empty() )
					 pieces[ ZERO_GEN ] =
							math::fn_const< genid >( prob_per_chrom_per_gen_t( 0. ) );
				pieces[ fromGen ] = math::fn_const< genid >( rate );
			}
			
	 };  // struct PopInfo
	 
	 // Field: popInfos - map from pop name to <PopInfo> for  that pop.
	 std::map< popid, PopInfo > popInfos;
};  // struct BaseModel

inline
std::ostream& operator<<( std::ostream& s, BaseModel const& baseModel ) {
	s << "[BaseModel:\n";
	for( BOOST_AUTO( pi, baseModel.popInfos.begin() );
			 pi != baseModel.popInfos.end(); ++pi ) {
		s << "pop " << pi->first << "\n";
		s << "  size fn         : " << pi->second.popSizeFn << "\n";
		s << "  size fn integral: " << indefiniteIntegral( pi->second.popSizeFn ) << "\n";
	}
	s << "]\n";
	return s;
}

//std::map< genid, ModelState > getBaseModel( DemographyP demography, HistEventsP histEvents );

}  // namespace cosi

BOOST_TYPEOF_REGISTER_TYPE(cosi::BaseModel)
BOOST_TYPEOF_REGISTER_TYPE(cosi::BaseModel::PopInfo)

#endif  // #ifndef INCLUDE_COSI_BASEMODEL_H
