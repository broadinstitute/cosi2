//
// Header: coalrate.h
//
// Computation of coalescence rate, especially in the face of varying population size.
//

#ifndef COSI_INCLUDE_COALRATE_H
#define COSI_INCLUDE_COALRATE_H

#include <boost/call_traits.hpp>
#include <cosi/defs.h>
#include <cosi/generalmath.h>

namespace cosi {

namespace math {

struct SweepPopSizeTraj;
template <> struct IsFunctionSpec<SweepPopSizeTraj>: public boost::true_type {};

struct SweepPopSizeComplementTraj;
template <> struct IsFunctionSpec<SweepPopSizeComplementTraj>: public boost::true_type {};

struct SweepCoalRateIntegral;
template <> struct IsFunctionSpec<SweepCoalRateIntegral>: public boost::true_type {};

struct SweepCoalRateComplementIntegral;
template <> struct IsFunctionSpec<SweepCoalRateComplementIntegral>: public boost::true_type {};

 namespace detail {
 typedef AreaType<genid,popsizeInv_float_t>::type genid_times_popsizeInv_float_t;
 }
 
template <>
class Function< /* TDomain= */ genid, /* TRange= */ popsize_float_t, /* TSpec= */ SweepPopSizeTraj> {

public:
	 typedef genid argument_type;
	 typedef popsize_float_t result_type;
	 typedef SweepPopSizeTraj spec_type;

	 Function( popsize_float_t basePopSize_,
						 gensInv_t sel_coeff_,
						 freq_t epsilon_,
						 gens_t end_shift_,
						 genid tend_ ):
		 basePopSize( basePopSize_ ), sel_coeff( sel_coeff_ ), epsilon( epsilon_ ),
		 end_shift( end_shift_ ), tend( tend_ ) { }
	 
	 result_type operator()( genid t ) const {
		 return basePopSize * ( epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*(t + end_shift - tend))) );
	 }

   friend ostream& operator<<( ostream& s, const Function& f ) {
		 s << "SweepPopSizeTraj[basePopSize=" << f.basePopSize << ", sel_coeff=" << f.sel_coeff <<
				", epsilon=" << f.epsilon << ", end_shift=" << f.end_shift << ", tend=" << f.tend << "]";
		 return s;
	 }
	 
	 
private:
	 popsize_float_t basePopSize;
	 gensInv_t sel_coeff;
	 freq_t epsilon;
	 gens_t end_shift;
	 genid tend;

	 friend class Function<genid, detail::genid_times_popsizeInv_float_t, SweepCoalRateIntegral >;
};																											 

template <>
class Function< /* TDomain= */ genid, /* TRange= */ popsize_float_t, /* TSpec= */ SweepPopSizeComplementTraj> {

public:
	 typedef genid argument_type;
	 typedef popsize_float_t result_type;
	 typedef SweepPopSizeComplementTraj spec_type;

	 Function( popsize_float_t basePopSize_,
						 gensInv_t sel_coeff_,
						 freq_t epsilon_,
						 gens_t end_shift_,
						 genid tend_ ):
		 basePopSize( basePopSize_ ), sel_coeff( sel_coeff_ ), epsilon( epsilon_ ),
		 end_shift( end_shift_ ), tend( tend_ ) { }
	 
	 result_type operator()( genid t ) const {
		 return basePopSize *
				( 1.0 -
					( epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*(t + end_shift - tend))) ) );
	 }
	 
   friend ostream& operator<<( ostream& s, const Function& f ) {
		 s << "SweepPopSizeComplementTraj[basePopSize=" << f.basePopSize << ", sel_coeff=" << f.sel_coeff <<
				", epsilon=" << f.epsilon << ", end_shift=" << f.end_shift << ", tend=" << f.tend << "]";
		 return s;
	 }
	 
private:
	 popsize_float_t basePopSize;
	 gensInv_t sel_coeff;
	 freq_t epsilon;
	 gens_t end_shift;
	 genid tend;

	 friend class Function<genid, detail::genid_times_popsizeInv_float_t,
												 SweepCoalRateComplementIntegral >;
};

template <>
class Function< /* TDomain= */ genid, /* TRange= */ detail::genid_times_popsizeInv_float_t,
															 /* TSpec= */ SweepCoalRateIntegral> {

public:
	 typedef genid argument_type;
	 typedef detail::genid_times_popsizeInv_float_t result_type;
	 typedef SweepCoalRateIntegral spec_type;

	 
	 Function( const Function<genid, popsize_float_t, SweepPopSizeTraj>& f ):
		 basePopSize( f.basePopSize ), selCoeff( f.sel_coeff ), epsilon( f.epsilon ),
		 end_shift( f.end_shift ), tend( f.tend ) { }

	 Function( popsize_float_t basePopSize_,
						 gensInv_t selCoeff_,
						 freq_t epsilon_,
						 gens_t end_shift_,
						 genid tend_ ):
		 basePopSize( basePopSize_ ), selCoeff( selCoeff_ ), epsilon( epsilon_ ),
		 end_shift( end_shift_ ), tend( tend_ ) { }

	 
	 result_type operator()( genid x ) const {
		 //return basePopSize * ( epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*ToDouble(t + end_shift - tend))) );
		 // return ( 1.0 / ( 2 * basePopSize * epsilon ) ) * ( epsilon * ToDouble( x ) +
		 // 																										(epsilon - 1.0 ) *
		 // 																										exp( selCoeff *
		 // 																												 ( selCoeff + ToDouble( x ) ) ) / selCoeff );
		 double num = ( (epsilon - 1) * exp( selCoeff * ( x + end_shift - tend  )  )  );
		 gens_t g = num / selCoeff;
		 gens_t num2 = ( epsilon * x - g - ZERO_GEN);
		 popsize_float_t denom = ( 2 * basePopSize * epsilon );

		 return num2 / denom;
		 
	 }
	 
   friend ostream& operator<<( ostream& s, const Function& f ) {
		 s << "SweepPopRateIntegral[basePopSize=" << f.basePopSize << ", selCoeff=" << f.selCoeff <<
				", epsilon=" << f.epsilon << ", end_shift=" << f.end_shift << ", tend=" << f.tend << "]";
		 return s;
	 }
	 
private:
	 popsize_float_t basePopSize;
	 gensInv_t selCoeff;
	 freq_t epsilon;
	 gens_t end_shift;
	 genid tend;

	 
};

template <>
class Function< /* TDomain= */ genid, /* TRange= */ detail::genid_times_popsizeInv_float_t,
															 /* TSpec= */ SweepCoalRateComplementIntegral> {

public:
	 typedef genid argument_type;
	 typedef detail::genid_times_popsizeInv_float_t result_type;
	 typedef SweepCoalRateComplementIntegral spec_type;
	 
	 Function( const Function<genid, popsize_float_t, SweepPopSizeComplementTraj>& f ):
		 basePopSize( f.basePopSize ), selCoeff( f.sel_coeff ), epsilon( f.epsilon ),
		 end_shift( f.end_shift ), tend( f.tend ) { }

	 Function( popsize_float_t basePopSize_,
						 gensInv_t selCoeff_,
						 freq_t epsilon_,
						 gens_t end_shift_,
						 genid tend_ ):
		 basePopSize( basePopSize_ ), selCoeff( selCoeff_ ), epsilon( epsilon_ ),
		 end_shift( end_shift_ ), tend( tend_ ) { }

	 
	 result_type operator()( genid x ) const {
		 //return basePopSize * ( epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*ToDouble(t + end_shift - tend))) );
		 // return ( 1.0 / ( 2 * basePopSize * epsilon ) ) * ( epsilon * ToDouble( x ) +
		 // 																										(epsilon - 1.0 ) *
		 // 																										exp( selCoeff *
		 // 																												 ( selCoeff + ToDouble( x ) ) ) / selCoeff );
		 // double num = ( (epsilon - 1) * exp( selCoeff * ( x + end_shift - tend  )  )  );
		 // gens_t g = num / selCoeff;
		 // gens_t num2 = ( epsilon * x - g - ZERO_GEN);
		 // popsize_float_t denom = ( 2 * basePopSize * epsilon );

		 // return num2 / denom;


		 BOOST_AUTO_TPL( numer1, (  epsilon * exp( -selCoeff * ( x + end_shift - tend )   )  ) );

		 BOOST_AUTO_TPL( denom1, (  (epsilon-1) * selCoeff ) );

		 BOOST_AUTO_TPL( ratio1, ( numer1
															 /   denom1    ) );

		 BOOST_AUTO_TPL( sum1, ( x + ratio1 ) );
		 
		 BOOST_AUTO_TPL( denom2, ( 2 * basePopSize ) );

		 return  ( sum1 - ZERO_GEN ) / denom2;
		 
	 }
	 
   friend ostream& operator<<( ostream& s, const Function& f ) {
		 s << "SweepPopRateComplementIntegral[basePopSize=" << f.basePopSize << ", selCoeff=" << f.selCoeff <<
				", epsilon=" << f.epsilon << ", end_shift=" << f.end_shift << ", tend=" << f.tend << "]";
		 return s;
	 }
	 
private:
	 popsize_float_t basePopSize;
	 gensInv_t selCoeff;
	 freq_t epsilon;
	 gens_t end_shift;
	 genid tend;
};



template <typename TPopSizeTrajSpec>
struct result_of_coalRateFunction {
	 
	 typedef
	 BOOST_TYPEOF_TPL((
											boost::declval< Function<genid,double,Const<> > >() /
											( boost::declval< Function<genid,double, Const<> > >( )
												* boost::declval< Function<genid, popsize_float_t, TPopSizeTrajSpec> >() )))
	 type;
};

template <typename TPopSizeTrajSpec>
typename result_of_coalRateFunction<TPopSizeTrajSpec>::type
coalRateFunction( const Function<genid, popsize_float_t, TPopSizeTrajSpec>& popSizeTraj ) {
	return
		 Function<genid,double,Const<> >( 1.0 ) /
							( Function<genid,double, Const<> >( 2.0 )
								* popSizeTraj );
}

namespace detail {
	typedef result_of_coalRateFunction<SweepPopSizeTraj>::type::result_type SweepPopSizeTraj_result_type;
	typedef result_of_coalRateFunction<SweepPopSizeTraj>::type::spec_type SweepPopSizeTraj_spec_type;
}
 
template <> struct result_of_indefiniteIntegral<genid,
	detail::SweepPopSizeTraj_result_type,
	detail::SweepPopSizeTraj_spec_type > {

	typedef detail::genid_times_popsizeInv_float_t atype;
	typedef const Function< genid, atype, SweepCoalRateIntegral>
		type;
};
																								
result_of_indefiniteIntegral<genid, popsizeInv_float_t,
	detail::SweepPopSizeTraj_spec_type >::type
inline
indefiniteIntegral( result_of_coalRateFunction<SweepPopSizeTraj>::type const&
										 f ) {
 	return Function< genid, detail::genid_times_popsizeInv_float_t, SweepCoalRateIntegral>
		 (f.second().second());
}

 namespace detail {
	 typedef result_of_coalRateFunction<SweepPopSizeComplementTraj>::type::result_type SweepPopSizeComplementTraj_result_type;
	 typedef result_of_coalRateFunction<SweepPopSizeComplementTraj>::type::spec_type SweepPopSizeComplementTraj_spec_type;
 }

template <> struct result_of_indefiniteIntegral<genid,
	detail::SweepPopSizeComplementTraj_result_type,
	detail::SweepPopSizeComplementTraj_spec_type > {
	 
	 typedef const Function< genid, detail::genid_times_popsizeInv_float_t, SweepCoalRateComplementIntegral>
	 type;
};
																								
result_of_indefiniteIntegral<genid, popsizeInv_float_t,
	detail::SweepPopSizeComplementTraj_spec_type >::type
inline
indefiniteIntegral( result_of_coalRateFunction<SweepPopSizeComplementTraj>::type const&
										 f ) {
 	return Function< genid, detail::genid_times_popsizeInv_float_t, SweepCoalRateComplementIntegral>
		 (f.second().second());
}

																												 

}  // namespace math



}  // namespace cosi


#endif // #ifndef COSI_INCLUDE_COALRATE_H
