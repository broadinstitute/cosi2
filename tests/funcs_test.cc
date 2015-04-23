
#ifdef NDEBUG
#undef NDEBUG
#endif
#ifndef COSI_DEV_PRINT
#define COSI_DEV_PRINT
#endif

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <boost/typeof/std/utility.hpp>
#include <boost/typeof/std/vector.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <cosi/typedval.h>
#include <cosi/generalmath.h>
#include <cosi/coalrate.h>

namespace cosi {


template <typename TFunc>
void testFunc( std::string msg, const TFunc& f ) {
	using namespace std;
	using namespace math;

	BOOST_CONCEPT_ASSERT((FunctionConcept<TFunc>));

	PRINT9( msg, eval( f, 1 ), eval( f, 100 ), eval( f, 1000 ), eval( f, 5 ), eval( f, 10 ), eval( f, 15 ),
					eval( f, 20 ), eval ( f, 25 ) );
}

template <typename T>
std::vector<T> makeVec( const T& a1 ) { std::vector<T> v; v.push_back( a1 ); return v; } 
template <typename T>
std::vector<T> makeVec( const T& a1, const T& a2 ) { std::vector<T> v; v.push_back( a1 );
	v.push_back( a2 ); return v; } 
template <typename T>
std::vector<T> makeVec( const T& a1, const T& a2, const T& a3 ) { std::vector<T> v; v.push_back( a1 );
	v.push_back( a2 ); v.push_back( a3 ); return v; } 

void testArrival() {
	using std::map;
	using std::make_pair;
	using std::cerr;
	using std::endl;
	using namespace math;

	Function<genid, popsize_float_t, SweepPopSizeTraj> sweepPopSizeTraj( popsize_float_t( 10000.0 ),
																																			 gensInv_t( .02 ),
																																			 freq_t( 1.0 / (2*10000.0) ),
																																			 gens_t( 0.0 ),
																																			 genid( 100.0 ) );

	BOOST_AUTO(coalRateFunc, coalRateFunction( sweepPopSizeTraj ) );
	PRINT( typeid( DomainType< BOOST_TYPEOF( coalRateFunc ) >::type ).name() );
	PRINT( typeid( RangeType< BOOST_TYPEOF( coalRateFunc ) >::type ).name() );
	//PRINT( typeid( SpecType< BOOST_TYPEOF( coalRateFunc ) >::type ).name() );

	//(void)coalRateFunc;
	BOOST_AUTO( ha, indefiniteIntegral( coalRateFunc ) );
	PRINT( typeid( DomainType< BOOST_TYPEOF( ha ) >::type ).name() );
	PRINT( typeid( RangeType< BOOST_TYPEOF( ha ) >::type ).name() );
	//PRINT( typeid( SpecType< BOOST_TYPEOF( ha ) >::type ).name() );
	BOOST_AUTO( nhpp, makeNonHomogeneousPoissonProcess( coalRateFunc, genid( 0.0 ), "test process" ) );
	PRINT( typeid( nhpp ).name() );

	{
		BOOST_AUTO( ha2, indefiniteIntegral( Function<genid, popsizeInv_float_t, Const<> >( popsizeInv_float_t(1.0) ) - coalRateFunc ) );
		PRINT( typeid( DomainType< BOOST_TYPEOF( ha2 ) >::type ).name() );
		PRINT( typeid( RangeType< BOOST_TYPEOF( ha2 ) >::type ).name() );
		//PRINT( typeid( SpecType< BOOST_TYPEOF( ha ) >::type ).name() );
		typedef BOOST_TYPEOF( coalRateFunc ) coalRateFunc_type;
		typedef coalRateFunc_type::result_type coalRateFunc_range_type;
		//BOOST_AUTO( nhpp2, ( makeNonHomogeneousPoissonProcess( Function<genid, popsizeInv_float_t, Const<> >( popsizeInv_float_t( 1.0 ) ) - coalRateFunc ) ) );
		//PRINT( typeid( nhpp2 ).name() );


		
	}


	Function<double,double,Const<> > f239( .1 );
	BOOST_AUTO( piece_pair, std::make_pair( 1.0, f239 ) );
	std::vector< BOOST_TYPEOF( piece_pair ) > piece_pairs;
	piece_pairs.push_back( piece_pair );
	Function<double,double,Piecewise< Const<> > >
		 rateFunc( piece_pairs );

	typedef double time_type;
	typedef double rate_type;

	Function<double,double,Const<> > f239a( .2 );
	rateFunc.getPieces().insert( std::make_pair( 0.0, f239a ) );
	
  BOOST_AUTO(proc, ( makeNonHomogeneousPoissonProcess( rateFunc, 0.0, "test process 2" ) ) );
	PRINT3( rateFunc, proc, proc.getRateFunctionIntegral() );
//ArrivalProcess<Time, Poisson< double, NonHomogeneous< BOOST_TYPEOF( 

	boost::random::mt19937 rng;
	
	
	// RandGenP randGen = boost::make_shared<RandGen>();

	// map<genid, rate_t> m;
	// m.insert( make_pair( genid(0), rate_t(.001) ) );
	// m.insert( make_pair( genid(4999), rate_t(.001) ) );
	// m.insert( make_pair( genid(5000), rate_t(.002) ) );
	// m.insert( make_pair( genid(100000), rate_t(.002) ) );

	// InhomogeneousPoissonProcess ipp( m );

	namespace acc = boost::accumulators;
	typedef acc::accumulator_set<double, acc::stats< acc::tag::sum_kahan, acc::tag::mean, acc::tag::variance > > acc_t;

	acc_t counts;
	acc_t wtimes;

	typedef double genid;

	time_type maxTime = 10000;

	for ( int k = 0; k < 1000; k++ ) {
		if ( !( k % 100 ) ) { PRINT( k ); }
		genid g(0);
		int cnt = 0;
		while ( g < maxTime ) {
			genid last_g = g;
			g = proc.nextArrivalTime( /* fromTime= */ g,
																/* maxTime= */ maxTime,
																/* rateFactor= */ 1.0,
																/* randGen= */ rng,
																/* eps= */ 1e-5 );
			wtimes( g - last_g );
			if ( g < maxTime ) {
				if ( genid(5000) < g ) {
					//PRINT2( i, g );
					cnt ++;
				}
			}
		}
		counts( double(cnt ) );
	}
	PRINT2( acc::mean( counts ), acc::variance( counts ) );
	PRINT2( acc::mean( wtimes ), sqrt( acc::variance( wtimes ) ) );
	
}



void test_funcs() {
	using namespace std;
	using namespace math;

	Function< double, int, Const< RunTime > > fr( 239 );

	
	Function< double, int, Const< CompileTime<10> > > fc;

	testFunc( "rtime", fr );
	testFunc( "ctime", fc );

	testFunc( "rintegral", indefiniteIntegral( fr ) );
	testFunc( "rintegrallin", indefiniteIntegral( fr * fr ) );
	//testFunc( "rintegra2", indefiniteIntegral( fc ) );


	// Function< double, int, X_To<0> > fc0;
	// testFunc( "rintegral3", indefiniteIntegral( fc0 ) );


	testFunc( "rintegral_fc1", indefiniteIntegral( fr * fr ) );

	Function< double, double, X_To<1> > fx;
	Function< double, double, Const<> > fxc( 2.0 );
	testFunc( "rintegral_fcx", indefiniteIntegral( fxc * fx ) );
	testFunc( "rintegral_fcx2", indefiniteIntegral( fx * fxc ) );

	PRINT( definiteIntegral( fxc, 0.0, 1.0 ) );
	PRINT( eval( integralFromPoint( fxc * fx, .5 ), .5 ) );
	PRINT( sizeof( integralFromPoint( fxc * fx, 0.0 ), .5 ) );	
	
	//testFunc( "rintegral_fcxa", indefiniteIntegral( fr * fx ) );
	//testFunc( "rintegral_fcxb", indefiniteIntegral( fx * fa ) );
	

	Function< double, int, Piecewise< Const<> > > f( makeVec( std::make_pair( 10.0, Function< double, int, Const<> >( 111 ) ),
																														std::make_pair( 20.0, Function< double, int, Const<> >( 222 ) ) ) ) ;
	// f.addPiece( 10.0, 111 );
	// f.addPiece( 20.0, 222 );
	testFunc( "piecewise", f );

	BOOST_AUTO( ff, fc * f );
	testFunc( "multiplied", ff );

	BOOST_AUTO( fline, make_line_through( 0.0, 0.0, 10.0, 5.0 ) );
	PRINT( eval( fline, 1.0 ) );
	PRINT( sizeof( fline ) );
	PRINT( sizeof(( Function<double,int,Const<> >( 3 ) )) );
	PRINT( sizeof(( Function<double,int,X_To<1> >( ) ) ));
	PRINT( sizeof(( Function<double,int,Const<> >( 3 ) ) * Function<double,int,X_To<1> >( ) ));
	PRINT( sizeof(( boost::compressed_pair< Function<double,int,Const<> >,
									Function<double,int,X_To<1> > >() )));
	PRINT( sizeof(( boost::compressed_pair< Function<double,int,X_To<1> >,
									Function<double,int,Const<> > >() )));
	
	PRINT( sizeof(( boost::compressed_pair< Function<double,int,X_To<1> >, Function<double,int,X_To<1> > >())));
	
	PRINT( sizeof(( Function<double,int,Const<> >( 3 ) ) * Function<double,int,X_To<1> >( )
								+ Function<double,int,Const<> >( 4 ) ) );
	
	PRINT(( ::boost::is_empty< Function<double,int,X_To<1> > >::value ));

	BOOST_AUTO( flinemin, -fline );
	PRINT( eval( flinemin, 2 ) );

	Function< double, double, Any<> > f_any( Function<double,double,Const<> >( 3.0 ) );
	PRINT( eval( f_any, 100.0 ) );
	f_any.reset( integralFromPoint( fxc * fx, 0.0 ) );
	PRINT( eval( f_any, .5 ) );

	PRINT( evalInverse( Function<double,double,X_To<2> >( 1.0 ), 0.0, 100000.0, 250.0, 1e-6 ) );

	testArrival();

	cerr << "tests done\n";
}

void test_funcs2() {
	using namespace std;
	using namespace math;
	
	Function< double, double, Const<> > f_const_7( 7.0 );
	Function< double, double, Const< CompileTime< 3 > > > f_const_3;
	Function< double, double, X_To<1> > f_x( 2.0 );
	
	BOOST_AUTO( f, ( f_const_7 * ( f_x - f_const_3 ) ) );
	
	BOOST_AUTO( e, exp_( f ) );
	
	BOOST_AUTO( i, indefiniteIntegral( e ) );
	
	
	typedef BOOST_TYPEOF( e ) e_type;
	
	typedef e_type::unop_spec_type e_spec;
	typedef result_of_differentiate< double, double, e_spec >::type r_d;
	
	PRINT( e );
	PRINT(( boost::is_same< SpecType< r_d >::type, Const<> >::type::value ));
	PRINT( i );

	//PRINT( typeid( e_spec ).name() );

//	typedef result_of_differentiate<double,double, SpecType<   >::spec_type s_t;
}

}  // namespace cosi

int main( int /*argc*/, char ** /*argv */ ) {
	cosi::test_funcs2();
	return EXIT_SUCCESS;
}
