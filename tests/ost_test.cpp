
#undef NDEBUG

#include <ctime>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <list>
#include <algorithm>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <boost/type_traits/has_nothrow_constructor.hpp>
#include <boost/concept_check.hpp>
#include <boost/phoenix/bind/bind_member_variable.hpp>
#define private public
#include <cosi/order_statistics.hpp>

namespace cosi {
namespace util {

struct PairWeight {
	 size_t a;
	 double b;

	 PairWeight(): a(0), b(0) {}
	 PairWeight( size_t a_, double b_ ): a( a_ ), b( b_ ) {}
	 PairWeight& operator+=( const PairWeight& other ) { a += other.a; b += other.b; return *this; }
	 PairWeight& operator-=( const PairWeight& other ) { a -= other.a; b -= other.b; return *this; }
	 friend PairWeight operator+( const PairWeight& pw1, const PairWeight& pw2 ) { return PairWeight( pw1.a + pw2.a, pw1.b + pw2.b ); }
	 friend PairWeight operator-( const PairWeight& pw1, const PairWeight& pw2 ) { return PairWeight( pw1.a - pw2.a, pw1.b - pw2.b ); }



	 friend PairWeight operator*( size_t factor, const PairWeight& pw ) {
		 return PairWeight( factor * pw.a, factor * pw.b );
	 }
	 friend PairWeight operator*( const PairWeight& pw, size_t factor ) {
		 return PairWeight( factor * pw.a, factor * pw.b );
	 }
	 
	 friend bool operator==( const PairWeight& pw1, const PairWeight& pw2 ) { return pw1.a == pw2.a && pw1.b == pw2.b; }
	 friend bool operator!=( const PairWeight& pw1, const PairWeight& pw2 ) { return !( pw1 == pw2 ); }
};

std::ostream& operator<<( std::ostream& s, const PairWeight& p ) { s << "(" << p.a << "," << p.b << ")"; return s; }


//////////

typedef size_t loc_t;

template <
	class PropObject,
	class PropType, 
	PropType PropObject::* Prop
	>
struct PropReader
{
	 typedef PropObject argument_type;
	 typedef PropType result_type;

	 typedef PropType value_type;

	 PropType operator()( const PropObject& p ) const { return p.*Prop; }
	 PropType& operator()( PropObject& p ) const { return p.*Prop; }
};

typedef std::pair< loc_t, size_t> p_t;

void showType( const char *msg, boost::mpl::true_ ) { std::cerr << "TRUE TYPE! " << msg << std::endl; }
void showType( const char *msg, boost::mpl::false_ ) { std::cerr << "FALSE TYPE! " << msg << std::endl; }

int ost_test() {

	

	using std::cerr;
	using std::endl;
	using std::pair;
	using std::make_pair;
	using std::vector;
	using std::random_shuffle;

		


//	typedef boost::order_statistics_tree< pair< loc_t, int >, CmpFirst, WeightExtractor > ost_t;
//	typedef ost_t::iterator ost_iter_t;


	typedef pair< loc_t, size_t > pair_t;

	namespace ph = boost::phoenix;
	
	typedef order_statistics_tree< pair_t , std::less<loc_t>, /* KeyExtractor= */ PropReader< pair_t, loc_t, &pair_t::first >,
																 /* WeightExtractor= */ PropReader< pair_t, loc_t, &pair_t::second > > ost_t;
	typedef ost_t::iterator ost_iter_t;

	showType( "ost_iter_t", boost::has_nothrow_constructor< ost_iter_t >() );
	showType( "ost_t", boost::has_nothrow_constructor< ost_t >() );

	BOOST_CONCEPT_ASSERT((boost::AssociativeContainer<ost_t>));


	{
		typedef pair< double, int > pair_t;

		typedef order_statistics_tree< pair_t , std::less<double>, /* KeyExtractor= */ PropReader< pair_t, double, &pair_t::first >,
																 /* WeightExtractor= */ PropReader< pair_t, int, &pair_t::second > > ost_t;


		unsigned int seed = time( NULL );
		cerr << "random seed: " << seed << endl;
		srand( seed );
		
		ost_t t;

		int wsum = 0;
		const int N = 1007;
		for ( int i = 0; i < N; i++ ) {
			int w = i * 100;
			t.insert( make_pair( ((double)( i + ( rand() % 355 ) ) ), w ) );
			 wsum += w;
			 assert( t.totalWeight() == wsum );
		}

		for ( int z = 0; z < 10000; z++ ) {
			assert( ((int)t.size()) == N );
			ost_t::iterator it = t[ rand() % N ];
			//int wsum_bef = wsum;
			wsum -= it.weightActual();
			// PRINT8( z, t.totalWeight(), it.position(), it->first, it->second, it.weightActual(), wsum,
			// 				wsum_bef );
			
			//t.printSVG( "/home/unix/ilya/public_html/be" );
			t.erase( it );
			//t.printSVG( "/home/unix/ilya/public_html/ae" );
			t.check();
			//PRINT2( t.totalWeight(), wsum );
			assert( t.totalWeight() == wsum );
			double k = rand() % 30;
			int v = rand() % 333;
			t.insert( make_pair( k, v ) );
			t.check();
			wsum += v;
			assert( t.totalWeight() == wsum );

			t.addWeightFast( t.begin(), t.end(), 3 );
			t.check();
			wsum += 3 * N;
			//PRINT2( t.totalWeight(), wsum );
			assert( t.totalWeight() == wsum );

			int m = rand() % N;
			int ww = rand() % 3737;
			t.addWeightFast( t.begin(), t[ m ], ww );
			t.check();
			t.addWeightFast( t[ m ], t.end(), ww );
			t.check();
			t.addWeightFast( t.begin(), t.end(), -ww );
			t.check();
			
			assert( t.totalWeight() == wsum );
		}

		
		//		t.printSVG( "/home/unix/ilya/public_html/cb" );
		t.check();
		t.addWeightFast( t[3], t[6], 1 );
		t.check();
		//		t.printSVG( "/home/unix/ilya/public_html/ca" );

		exit( 0);
	}

	typedef order_statistics_tree< double > ost_dbl_t;
	ost_dbl_t dblTree;
	dblTree.insert( 1.234 );
	PRINT2( sizeof( detail::DummyWeight ), sizeof( ost_dbl_t::tree_node ) );

	unsigned int seed = time( NULL );
	cerr << "random seed: " << seed << endl;
	srand( seed );

	{
		typedef pair< loc_t, PairWeight > pair_t;
		typedef order_statistics_tree< pair_t , std::less<loc_t>, /* KeyExtractor= */ PropReader< pair_t, loc_t, &pair_t::first >,
																	 /* WeightExtractor= */ PropReader< pair_t, PairWeight, &pair_t::second > > ost_t;
		typedef ost_t::iterator ost_iter_t;
		
		ost_t t;
		for ( size_t i = 0; i < 30; i++ )
			 t.insert( make_pair( i, PairWeight( 2*i, static_cast<double>( 3*i ) ) ) );
		ost_iter_t i = t.find( 4 );
		t.check();
		t.printSVG( "/home/unix/ilya/public_html/cb" );
		i.addWeight( 500.0, PropReader< PairWeight, double, &PairWeight::b >() );
		t.check();
		t.printSVG( "/home/unix/ilya/public_html/ca" );
	}
	
	for ( int k = 0; k <3000000; k++ ) {
		ost_t t;

		const size_t N = rand() % 500;

		vector<size_t> v;
		vector<size_t> inds;
		vector<size_t> weights;
		for ( size_t i = 0; i < N; i++ ) {
			v.push_back( i );
			inds.push_back( i );
			weights.push_back( rand() % 373737 );
		}
		random_shuffle( inds.begin(), inds.end() );
		random_shuffle( weights.begin(), weights.end() );
		for ( size_t i = 0; i < N; i++ ) {
			t.insert( make_pair( v[ inds[i] ], weights[ inds[i] ] ) );
			t.check();
			//t.printSVG( "/home/unix/ilya/public_html/t" + ToString( k ) + "_" + ToString( i ) );
		}

		size_t n_ok = 0;
		loc_t targetWeight = 0;
		for ( size_t i = 0; i < N; i++ ) {
			ost_iter_t it = t.find( i );

			//PRINT6( i, v[i], *it, it.position(), it.cumWeight(), targetWeight );
			if ( i != v[i] || it.position() != i || it.cumulWeight() != targetWeight )
				 cerr << "CHEGO??" << endl;
			else
				 n_ok++;
			targetWeight += weights[ i ];
		}
		if ( !( k % 1000 ) ) cerr << "k=" << k << " n_ok=" << n_ok << endl;

		random_shuffle( inds.begin(), inds.end() );
		size_t totWeight = t.totalWeight();
		for ( size_t i = 0; i < N; i++ ) {
			//cerr << "deleting " << i << endl;
			assert( t.size() == N - i );
			//t.printSVG( "/home/unix/ilya/public_html/bef" );
			
			t.erase( v[ inds[i] ] );
			//cerr << "deleted " << i << endl;
			//t.printSVG( "/home/unix/ilya/public_html/aft" );
			t.check();
			totWeight -= weights[ inds[i] ];
			assert( t.size() == N - i - 1 );
			assert( totWeight == t.totalWeight() );
			//cerr << "checked deletion of " << i << endl;
		}
		
	}  // for each experiment

#if 0
	

	t.printDot( "rb0.dot" );
	ost_iter_t v1 = t.insert( make_pair( loc_t( .5 ), 27 ) ).first;
	t.printDot( "rb1.dot" );
	ost_iter_t v2 = t.insert( make_pair( loc_t( .6 ), 17 ) ).first;
	t.printDot( "rb2.dot" );
	ost_iter_t v3 = t.insert( make_pair( loc_t( .7 ), 7 ) ).first;
	t.printDot( "rb3.dot" );
	PRINT3( *v1, *v2, *v3 );
	PRINT4( v1.position(), v2.position(), v3.position(), t.end().position() );
	PRINT4( v1.cumulWeight(), v2.cumulWeight(), v3.cumulWeight(), t.end().cumulWeight() );
	PRINT( t.lower_bound( make_pair( loc_t( .5 ), 0 ) ).position() );
	PRINT( t.upper_bound( make_pair( loc_t( .5 ), 0 ) ).position() );

	t.printDot( "rb.dot" );
	
	t.erase( v2 );
	PRINT2( v1.position(), v3.position() );

#endif	
	return 0;
}

} // namespace util
} // namespace cosi

int main( int /* argc */, char ** /* argv */ ) { return cosi::util::ost_test(); }
