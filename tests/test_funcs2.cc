//
// Test driver: test_funcs
//
// Tests the generic-functions functionality in generaicmath.h
//

#define BOOST_TEST_MODULE general_funcs test

#include <utility>
#include <boost/test/included/unit_test.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cosi/general/math/generalmath2.h>
#include <cosi/coalrate.h>
#include <cosi/defs.h>
#if 0
using cosi::math::Function;
using cosi::math::Const;
using cosi::math::CompileTime;
using cosi::math::Piecewise;
using cosi::math::RunTime;
using cosi::math::definiteIntegral;
using cosi::math::indefiniteIntegral;
using cosi::genid;
using cosi::popsize_float_t;
using std::make_pair;
using boost::numeric::bounds;
using boost::math::isnan;

BOOST_AUTO_TEST_SUITE( test_math_Function )

BOOST_AUTO_TEST_CASE( const_funcs_return_their_const_value ) {
	Function< double, double, Const< CompileTime< 3 > > > f_const_3;
	BOOST_CHECK_EQUAL( eval( f_const_3, 100), 3 );
	
	Function< double, double, Const<> > f_const_rt( 239 );
	BOOST_CHECK_EQUAL( eval( f_const_rt, 0), 239 );
	
}

struct piecewise_func_fixture {
	Function< double, double, Piecewise< Const<> > > f_piecewise;
};

BOOST_FIXTURE_TEST_SUITE( piecewise_functions, piecewise_func_fixture )

BOOST_AUTO_TEST_CASE( piecewise_funcs_test ) {

	BOOST_CHECK( (isnan)( f_piecewise( 3 ) ) );
	f_piecewise.getPieces().insert( std::make_pair( 0.0, Function< double, double, Const<> >( 1.0 ) ) );
	BOOST_CHECK( (isnan)( f_piecewise( -1e-100 ) ) );
	BOOST_CHECK_EQUAL( f_piecewise( 0 ), 1.0 );
	BOOST_CHECK_EQUAL( (f_piecewise( 1e-100 )), 1.0 );

	typedef std::numeric_limits<double> dbl_lim;

	const double dbl_lowest = bounds<double>::lowest();
	const double dbl_highest = bounds<double>::highest();
	const double dbl_smallest = bounds<double>::smallest();
	const double dbl_nan = dbl_lim::quiet_NaN();
	const double dbl_inf = dbl_lim::infinity();

	std::cerr.precision( 20 );
	PRINT4( dbl_lowest, dbl_highest, dbl_smallest, dbl_nan );
	
	BOOST_CHECK_EQUAL( f_piecewise( dbl_highest ), 1.0 );
	double ten(10.0);
	f_piecewise.getPieces().insert( make_pair( ten, Function< double, double, Const<> >( 15.0 ) ) );
	BOOST_CHECK_EQUAL( f_piecewise( ten ), 15.0 );
	BOOST_CHECK_EQUAL( f_piecewise( ten + dbl_smallest ), 15.0 );
	BOOST_CHECK_EQUAL( f_piecewise( 9.999 ), 1.0 );
	double bef10 = static_cast<double>( ten ) - static_cast<double>( 1e-15 );
	PRINT( bef10 );
	BOOST_CHECK( bef10 < ten );
	BOOST_CHECK_EQUAL( f_piecewise( bef10 ), 1.0 );
	BOOST_CHECK_EQUAL( f_piecewise( dbl_highest ), 15.0 );
	BOOST_CHECK( (boost::math::isnan)( f_piecewise( dbl_lowest ) ) );

	BOOST_CHECK_EQUAL( definiteIntegral( f_piecewise, 0.0, 2.0 ), 2.0 );
	// PRINT( definiteIntegral( f_piecewise, 0.0, 9.0 ) );
	// PRINT( definiteIntegral( f_piecewise, 0.0, 9.99999 ) );
	BOOST_CHECK_EQUAL( definiteIntegral( f_piecewise, 0.0, 10.0 ), 10 );
	//PRINT( definiteIntegral( f_piecewise, 10.0, 11.0 ) );
	BOOST_CHECK_EQUAL( definiteIntegral( f_piecewise, 0.0, 11.0 ), 25.0 );

	f_piecewise.getPieces().insert( make_pair( dbl_lowest, Function<double,double,Const<> >( 0.0 ) ) );
	BOOST_CHECK_EQUAL( definiteIntegral( f_piecewise, dbl_lowest * .99, 3.0 ), 3.0 );

	BOOST_CHECK_EQUAL( eval( Function< double, double, Const<> >( 10.0 ) * f_piecewise, 10.0 ), 150.0 );

	PRINT2( f_piecewise, (Function< double, double, Const<> >( 10.0 ) * f_piecewise) );
	PRINT( indefiniteIntegral( Function< double, double, Const<> >( 1.0 ) /
														 ( Function< double, double, Const<> >( 10.0 ) * f_piecewise ) ) );


	Function< genid, popsize_float_t, Piecewise< Const<> > > f_ps;
	f_ps.getPieces().insert( std::make_pair( genid(0.0), popsize_float_t( 10000.0 ) ) );
	f_ps.getPieces().insert( std::make_pair( genid(10.0), popsize_float_t( 5000.0 ) ) );
	f_ps.getPieces().insert( std::make_pair( genid(20.0), popsize_float_t( 50000.0 ) ) );

	BOOST_AUTO( crf, coalRateFunction( f_ps ) );
	PRINT( crf );
	PRINT( indefiniteIntegral( crf ) );
}


BOOST_AUTO_TEST_SUITE_END() // piecewise_functions

BOOST_AUTO_TEST_SUITE_END()  // test_math_Function
#endif
