//
// * Program: sample_stats_extra
//
// Compute various summary statistics for a set of simulations in ms output format.
// Some extensions of ms format are accepted.  Which statistics are gathered, can be configured
// on the command line.
//
// For each simulation, one line of statistics is output.
//
// Adapted from the ms package by Richard Hudson.
// See http://home.uchicago.edu/rhudson1/source/mksamples.html 
//
// Usage:
//
// sample_stats_extra MARGIN afsMax (dmin dmax)*
//
// Args:
//
//    MARGIN - margin of simulated region to ignore; this fraction from the start and from the end
//       of the region will be ignored, stats will be gathered only from the middle portion of the region.
//    afsMax - when gathering the allele frequency spectrum, gather the fraction of SNPs having 1, 2, ..., afsMax
//       chroms with derived allele
//    dmin dmax - gather parameters of the distribution of LD for pairs of SNPs at physical separation of between dmin and dmax
//       (specified as fraction of total region length)
//

// ** includes
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <stdint.h>
#include <limits>
#include <stdexcept>
#include <exception>
#include <sstream>
#include <ios>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>
#include <set>
#include <map>
#include <string>
#include <locale>
#include <limits>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/assign.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/swap.hpp>
#include <boost/multi_array.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/numeric.hpp>
#include <boost/concept/assert.hpp>
#include <boost/range/numeric.hpp>
#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/array.hpp>
#include <boost/cstdint.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/next_prior.hpp>
#include <boost/program_options.hpp>
#include <boost/integer_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/declval.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>
#include <boost/phoenix/bind/bind_member_variable.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/core/value.hpp>
#include <boost/phoenix/core/argument.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/regex.hpp>
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
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/array.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/throw_exception.hpp>
#include <boost/exception/all.hpp>
//#include <Sequence/SimParams.hpp>

#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wmissing-declarations"
#endif

namespace cosi {

using std::vector;
using std::set;
using std::ostream;
using std::istream;
using std::ostream_iterator;
using std::ostringstream;
using std::istringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios_base;
using std::pair;
using std::make_pair;
using std::locale;
using std::ctype;
using std::string;

typedef std::string filename_t;

// ** utils


//
// cosi exceptions
//

struct cosi_error: virtual std::exception, virtual boost::exception { };
struct cosi_io_error: virtual cosi_error { };

typedef boost::error_info<struct tag_errno_code,int> errno_code;
typedef boost::error_info<struct tag_error_msg,std::string> error_msg;

double tajd(int, int, double) ;

string Date( );

typedef double cosi_double;
inline cosi_double ToDouble( const cosi_double& x ) { return x; }
vector<string> Split( string s, const string& separators = " \t" );
void Split( const string& s, const string& sep, string& part1, string& part2 );
void Split( const string& s, const string& sep, string& part1, string& part2, string& part3 );
string Join( const string& sep, const vector<string>& strings );
string Join( const string& sep, const string& s1, const string& s2 );
string Join( const string& sep, const string& s1, const string& s2, const string& s3 );
string Join( const string& sep, const string& s1, const string& s2, const string& s3, const string& s4 );
string Join( const string& sep, const string& s1, const string& s2, const string& s3, const string& s4, const string& s5 );
string MakeValidIdentifier( const string& s );
double Frac( int a, int b );

class SNPCond;
ostream& operator<<( ostream& os, const SNPCond& snpCond );


#define ForEach BOOST_FOREACH

// The following PRINT macros are intended primarily for debugging.

#define PRCORE(X) #X " = " << (X)
#define PRINT(X) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << endl
#define PRINT2(X, Y) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << endl
#define PRINT3(X, Y, Z) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
	<< PRCORE(Z) << endl
#define PRINT4(X, Y, Z, W) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
	<< PRCORE(Z) << ", " << PRCORE(W) << endl
#define PRINT5(X, Y, Z, W, T) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
	<< PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << endl
#define PRINT6(X, Y, Z, W, T, U) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
	<< PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "				\
	<< PRCORE(U) << endl
#define PRINT7(X, Y, Z, W, T, U, V) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) \
  << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", " \
  << PRCORE(U) << ", " << PRCORE(V) << endl
#define PRINT10(X, Y, Z, W, T, U, V, A, B, C) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) \
  << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", " \
  << PRCORE(U) << ", " << PRCORE(V) << PRCORE(A) << ", " << PRCORE(B) << ", " << PRCORE(C) << endl


string Date( )
{    char nowstr[80];
	time_t nowbin;
	const struct tm *nowstruct;
	static bool locale_has_been_set = false;
	if( ! locale_has_been_set ) {
		(void)setlocale(LC_ALL, "");
		locale_has_been_set = true;
	}
	if (time(&nowbin) == (time_t) - 1) return "(date unavailable - time failed)";
	nowstruct = localtime(&nowbin);
	if (strftime(nowstr, 80, "%a %b %d %H:%M:%S %Y", nowstruct) == (size_t) 0)
		 return "(date unavailable - strftime failed)";
	return string(nowstr);
}

// ** Class SumKeeper

template <typename ValT = cosi_double, typename CountT = size_t>
class SumKeeper {
public:
	 SumKeeper() { clear(); }

	 typedef ValT value_t;
	 typedef CountT count_t;

	 // Method: add
	 // Add a value to the sum keeper.
	 void add( ValT x ) {
		 if ( (boost::math::isnan)( ToDouble( x ) ) ) numNaNs++;
		 else if ( (boost::math::isinf)( ToDouble( x ) ) ) numInfs++;
		 else {
			 numVals++;
		
			 int i = 0;
			 for ( typename vector<ValT>::const_iterator yi = partials.begin(); yi != partials.end(); yi++ ) {
				 ValT y = *yi;

				 if ( ::fabs( ToDouble( x ) ) < ::fabs( ToDouble( y ) ) ) boost::swap( x, y );
				 ValT hi = x + y;
				 ValT lo = y - ( hi - x );
				 if ( ToDouble( lo ) != 0.0 ) partials[ i++ ] = lo;
				 x = hi;
			 }
			 partials.erase( partials.begin()+i, partials.end() );
			 partials.push_back( x );
		 }
	 }

	 // Method: add
	 // Add a range of values to this SumKeeper.
	 template <class ValRange>
	 void add( const ValRange& valRange ) { ForEach( ValT val, valRange ) add( valRange ); }

	 SumKeeper<ValT,CountT>& operator+=( ValT x ) { add( x ); return *this; }
	 SumKeeper<ValT,CountT>& operator+=( const SumKeeper<ValT,CountT>& sk  ) {
		 add( sk.getSum() );
		 numVals += sk.numVals;
		 numNaNs += sk.numNaNs;
		 numInfs += sk.numInfs;
		 return *this;
	 }

	 ValT getSum() const { return boost::accumulate( partials, ValT(0.0) ); }
	 CountT getNumVals() const { return numVals; }
	 CountT getNumNaNs() const { return numNaNs; }
	 CountT getNumInfs() const { return numInfs; }

	 ValT getMean() const { return numVals > 0 ? ( getSum() / ValT(numVals) ) : std::numeric_limits<ValT>::quiet_NaN(); }

	 void clear() {
		 partials.clear();
		 numVals = 0;
		 numNaNs = 0;
		 numInfs = 0;
	 }
  
private:
	 vector<ValT> partials;

	 CountT numVals;
	 CountT numNaNs;
	 CountT numInfs;

#if 0	 
	 friend std::ostream& operator<<( std::ostream& s, const SumKeeper<ValT,CountT>& k ) {
		 s << "SumKeeper[ " << k.getSum() << " " << k.getNumVals() << " " << k.getNumNaNs() << " "
			 << k.getNumInfs() << " ]";
		 return s;
	 }

	 friend std::istream& operator>>( std::istream& s, SumKeeper<ValT,CountT>& k ) {
		 std::string txt;
		 s >> txt;
		 if ( txt != "SumKeeper{" ) throw std::ios::failure( "SumKeeper extractor: missing header" );
		 k.clear();
		 ValT theSum;
		 s >> theSum >> k.numVals >> k.numNaNs >> k.numInfs >> txt;
		 if ( txt != "}" ) throw std::ios::failure( "SumKeeper extractor: missing closing paren" );
		 k.partials.push_back( theSum );
		 return s;
	 }
#endif	 
  
};  // class SumKeeper

// ** Class StatKeeper
//
// Keeps a running sum and sum-of-squares without round-off errors,
// allowing accurate computation of the mean and stddev.
//
template <typename ValT = cosi_double, typename CountT = size_t>
class StatKeeper {
public:
	 typedef ValT value_type;
	 typedef CountT count_type;
	 
	 StatKeeper() { clear(); }
	 StatKeeper( const std::string& name_ ): name( name_ ) { clear(); }

	 void add( ValT x ) { sum.add( x ); sumSq.add( x * x ); }

	 void clear() { sum.clear(); sumSq.clear(); }

	 ValT getSum() const { return sum.getSum(); }
	 ValT getSumSq() const { return sumSq.getSum(); }
	 CountT getNumVals() const { return sum.getNumVals(); }
	 CountT getNumNaNs() const { return sum.getNumNaNs(); }
	 CountT getNumInfs() const { return sum.getNumInfs(); }

	 // Method: add
	 // Add a range of values to this SumKeeper.
	 template <class ValRange>
	 void add( const ValRange& valRange ) { ForEach( ValT val, valRange ) add( val ); }

	 void operator()( value_type val ) { this->add( val ); }

	 ValT getMean() const { return sum.getMean(); }
	 ValT getStd() const {
		 ValT meanSoFar = getMean();
		 return std::sqrt( sumSq.getMean() - ( meanSoFar * meanSoFar ) );
	 }
	 const std::string& getName() const { return name; }
	 
private:
	 // Fields:
	 //
	 //   sum - sum of values passed to <add()>
	 //   sumSq - sum of squares of values passed to <add()>
	 SumKeeper<ValT,CountT> sum, sumSq;
	 std::string name;
};  // class StatKeeper
	

//////////////////////////////////////////////////////

// ** more utils

double nucdiv(int, int, char **);
double hfay(int, int, char **);
double thetah(int, int, char **);
void
biggerlist(int nsam, unsigned nmax, char **list );

void chk_helper( bool cond, const char *file, int line, const char *msg );
void chk_helper( bool cond, const char *file, int line, const char *msg ) {
	if ( !cond ) {
		ostringstream s;
		s << "ERROR: file=" << file << " line=" << line << " msg=" << msg << endl;
		throw std::logic_error( s.str() );
	}
}

#define chk( cond ) chk_helper( cond, __FILE__, __LINE__, #cond )

char **
cmatrix(int nsam,int len);
void free_cmatrix(char **cm, int nsam );

int maxsites = 1000 ;

template <typename T>
inline T getVal( string str ) {
	istringstream s( str.c_str() );
	s.exceptions( std::ios::failbit | std::ios::badbit );
	T val;
	s >> val;
	chk( s.eof() );
	return val;
}

template <typename TA, typename TB>
inline ostream& operator<<( ostream& s, const pair<TA,TB>& p ) {
	s << "(" << p.first << "," << p.second << ")";
	return s;
}

template <class T>
inline ostream& operator<<( ostream& s, const vector<T>& v) {
	s << "[";

	for(typename std::vector<T>::const_iterator it = v.begin();
			it != v.end();
			it++) {
		s<<(*it)<<", ";
	}
	s << "]";
	return s;
}

// *** string utils

//
// Function: Split
//
// Split a string based on separators.
//
vector<string> Split( string s, const string& separators ) {
	typedef boost::tokenizer<boost::char_separator<char> > 
		 tokenizer;
  boost::char_separator<char> sep(separators.c_str());
  tokenizer tokens(s, sep);
	vector<string> result( tokens.begin(), tokens.end() );
	return result;
}

void Split( const string& s, const string& sep, string& part1, string& part2 ) {
	vector<string> parts( Split( s, sep ) );
	if ( parts.size() > 2 ) throw std::logic_error( "too many parts" );
	part1.clear();
	part2.clear();
	if ( parts.size() >= 1 ) part1 = parts[0];
	if ( parts.size() >= 2 ) part2 = parts[1];
}

void Split( const string& s, const string& sep, string& part1, string& part2, string& part3 ) {
	vector<string> parts( Split( s, sep ) );
	if ( parts.size() > 3 ) throw std::logic_error( "too many parts" );
	part1.clear();
	part2.clear();
	part3.clear();
	if ( parts.size() >= 1 ) part1 = parts[0];
	if ( parts.size() >= 2 ) part2 = parts[1];
	if ( parts.size() >= 3 ) part3 = parts[2];
}

//
// Function: Join
//
// Joint strings with a separator
//
string Join( const string& sep, const vector<string>& strings ) {
	vector<string> strings_f;
	ForEach( string s, strings ) if ( !s.empty() ) strings_f.push_back( s );
	ostringstream out;
	copy( strings_f.begin(), strings_f.end(), ostream_iterator<string>( out, sep.c_str() ) );
	return out.str();
}

template <typename T> vector<T> MakeVec() { return vector<T>(); }
template <typename T> vector<T> MakeVec( const T& v1 ) { return vector<T>( 1, v1 ); }
template <typename T> vector<T> MakeVec( const T& v1, const T& v2 ) {
	vector<T> v( MakeVec( v1 ) );
	v.push_back( v2 );
	return v;
}
template <typename T> vector<T> MakeVec( const T& v1, const T& v2, const T& v3 ) {
	vector<T> v( MakeVec( v1, v2 ) );
	v.push_back( v3 );
	return v;
}
template <typename T> vector<T> MakeVec( const T& v1, const T& v2, const T& v3, const T& v4 ) {
	vector<T> v( MakeVec( v1, v2, v3 ) );
	v.push_back( v4 );
	return v;
}
template <typename T> vector<T> MakeVec( const T& v1, const T& v2, const T& v3, const T& v4, const T& v5 ) {
	vector<T> v( MakeVec( v1, v2, v3, v4 ) );
	v.push_back( v5 );
	return v;
}

string Join( const string& sep, const string& s1, const string& s2 ) { return Join( sep, MakeVec( s1, s2 ) ); }
string Join( const string& sep, const string& s1, const string& s2, const string& s3 ) { return Join( sep, MakeVec( s1, s2, s3 ) ); }
string Join( const string& sep, const string& s1, const string& s2, const string& s3, const string& s4 ) {
	return Join( sep, MakeVec( s1, s2, s3, s4 ) );
}
string Join( const string& sep, const string& s1, const string& s2, const string& s3, const string& s4, const string& s5 ) {
	return Join( sep, MakeVec( s1, s2, s3, s4, s5 ) );
}


//
// Function template: ToString
//
// A general method for converting arbitrary types to a human-readable
// string representation.
//
template <typename T>
inline string ToString( const T& v ) {
	std::ostringstream s;
	s << v;
	return s.str();
}

// Func: MakeValidIdentifier
// Replace any characters in the string that are not valid in a variable name, with underscores.
// Sequences of invalid chars are replaced by a single underscore.
// If the first character is not a letter, prepend an underscore.
string MakeValidIdentifier( const string& s ) {
	static const boost::regex invalidChar( "\\W+" );
	locale loc;
	string result( boost::regex_replace( s, invalidChar, "_" ) );
	if ( result.empty() || !std::use_facet< ctype<char> >( loc ).is( ctype<char>::alpha, result.at( 0 ) ) )
		 result = string("_") + result;
	return result;
}

double Frac( int a, int b ) { return double(a) / double(b); }

template <typename T> inline int isize( const vector<T>& v ) { return int( v.size() ); }

//
// End section: General utils
//

// ** logical types

// Logical type: nchroms_t
// A count of chromosomes.
typedef int nchroms_t;

// Logical type: nsnps_t
// A count of SNPs.
typedef int nsnps_t;

// Logical type: loc_t
// A location of a SNP (represented as a fraction of the simulated region).
typedef double loc_t;

// Logical type: gloc_t
// A genetic position of a SNP
typedef double gloc_t;

// Logical type: genid
// A generation (a time point).
typedef double genid;

// Logical type: prob_t
// A probability.
typedef double prob_t;
typedef double freq_t;

// Logical type: popid
// Name of a population
typedef int popid;

typedef size_t popNum_t;

template <typename T>
T get_null() { throw std::logic_error( "unimpl" ); return T(); }

template <>
inline loc_t get_null<loc_t>() { return std::numeric_limits<loc_t>::quiet_NaN(); }
template <>
inline int get_null<int>() { return -1; }

template <typename T>
inline bool is_null( const T& val ) { return val == get_null<T>(); }

template <>
inline bool is_null<double>( const double& loc ) { return (boost::math::isnan)( loc ); }

// Logical type: len_t
// A length of a segment of the simulated region; a difference of two loc_t's.
typedef double len_t;

// Logical type: gens_t
// Length of a time interval, in generations.
typedef double gens_t;

//const loc_t NULL_LEN = std::numeric_limits<loc_t>::quiet_NaN();

// Logical type: snp_id_t;
// Identifies a specific SNP.
typedef int snp_id_t;

// Var: derCounts
// Number of derived alleles of each SNP.  Indexed by <snp_id>.
vector< nchroms_t > derCounts;


vector< nchroms_t > derCounts_DIND_der, derCounts_DIND_anc;


// Var: posit
// The location of each SNP
loc_t *posit;

// Var: posit_gloc - genetic positions of SNPs
gloc_t *posit_gloc;

bool compute_LD( int site_A, int site_B, int nsam, char **list,
								 double *r2, double *Dprime );
bool compute_LD( snp_id_t site_A, snp_id_t site_B, nchroms_t bsam, nchroms_t nsam, char **list,
								 double *r2, double *Dprime );

// ** Class ValRange
//
// A range of values.
//
template <typename T>
class ValRange {
public:
	 ValRange( T lwr = get_null<T>(), T upr = get_null<T>() ): bounds( make_pair( lwr, upr ) ) {}
	 ValRange( string rangeDef ) { init( rangeDef ); }

	 void init( string rangeDef ) {
		 if ( !rangeDef.empty() ) {

			 typedef boost::tokenizer<boost::char_separator<char> > 
					tokenizer;
			 boost::char_separator<char> sep("-");
			 tokenizer tokens(rangeDef, sep);
			 vector<string> result( tokens.begin(), tokens.end() );

			 chk( result.size() == 1  ||  result.size() == 2 );
			 if ( result.size() == 1 ) {
				 bounds.first = bounds.second = boost::lexical_cast<T>( result.at(0) );
			 } else {
				 bounds.first = boost::lexical_cast<T>( result.at(0) );
				 bounds.second = boost::lexical_cast<T>( result.at(1) );
			 }
		 }
	 }

	 bool operator() ( T val ) const { return ( is_null( bounds.first ) || val >= bounds.first ) &&
				( is_null( bounds.second ) || val <= bounds.second ); }

	 T getMin() const { return bounds.first; }
	 T getMax() const { return bounds.second; }
	 void setMin( T min_ ) { bounds.first = min_; }
	 void setMax( T max_ ) { bounds.second = max_; }
	 
private:
	 // Field: bounds
	 // The pair of bounds of the range values
	 pair<T,T> bounds;

};  // class ValRange

template <typename T>
istream& operator>>( istream& is, ValRange<T>& valRange );

template <typename T>
istream& operator>>( istream& is, ValRange<T>& valRange ) {
	string valRangeDef;
	is >> valRangeDef;
	valRange.init( valRangeDef );
	
	return is;
}


template <typename T>
ostream& operator<<( ostream& os, const ValRange<T>& valRange ) {
	if ( !is_null( valRange.getMin() ) ) os << valRange.getMin();
	if ( !is_null( valRange.getMax() ) && !( valRange.getMin() == valRange.getMax() ) ) os << "-" << valRange.getMax();
	
	return os;
}

//
// ** Class SNPCond
//
// A condition on SNPs.
//
class SNPCond {
public:

	 void init( string locRangeDef ) { locRange.init( locRangeDef ); }

	 bool operator() ( snp_id_t snpId ) const {
		 return locRange( posit[ snpId ] );
	 }

	 string getColSfx() const { return ToString( locRange ); }

	 const ValRange<loc_t>& getLocRange() const { return locRange; }

private:
	 // Field: locRange
	 // Range of locations of SNPs satisfying this condition.
	 ValRange<loc_t> locRange;
	 
};  // class SNPFilter

ostream& operator<<( ostream& os, const SNPCond& snpCond ) {
	os << snpCond.getColSfx();
	return os;
}


//
// ** Class AFS
//
// Allele frequency spectrum for a specified subset of SNPs.
//
class AFS {
public:

	 // class: Bin
	 // A bin of the allele frequency spectrum.
	 class Bin {
	 public:
			Bin(): nsnps( 0 ) {}

			Bin( string binDef ): nsnps( 0 ) { derCountRange.init( binDef ); }
			// BIN( nchroms_t binMin, nchroms_t binMax ):
			// 	derCountRange( binMin, binMax )

			nchroms_t getMinDerCount() const { return derCountRange.getMin(); }
			nchroms_t getMaxDerCount() const { return derCountRange.getMax(); }

			const ValRange<nchroms_t>& getDerCountRange() const { return derCountRange; }

			void clear() { nsnps = 0; }
			void increment() { nsnps++; }
			nsnps_t getSnpCount() const { return nsnps; }

			string getColSfx() const { return ToString( derCountRange ); }

	 private:
			// Field: derCountRange
			// Range of derived-allele counts for SNPs going into this bin.
			ValRange<nchroms_t> derCountRange;
			
			// Field: nsnps
			// Number of SNPs in this bin
			nsnps_t nsnps;
			
	 };  // class Bin

	 // Method: init
	 // Initialize an AFS from a string definition.
	 // Definition looks like:
	 // 1,2,3-4
	 // 1,2,3-4;30000-40000
	 // 20-360:20 means bins of size 20
	 // The required part before semicolon specifies the bins.
	 // The optional part after the semicolon specifies the range of locations from which we take SNPs.
	 void init( const string& afsDef ) {

		 string binsDef, snpCondDef;
		 Split( afsDef, "/",  binsDef, snpCondDef );
		 snpCond.init( snpCondDef );
		 vector<string> binDefs = Split( binsDef, "," );
		 ForEach( string binDef, binDefs ) {
			 string binRange, binStep;
			 Split( binDef, ":", binRange, binStep );

			 if ( binStep.empty() ) {
				 int binId = bins.size();
			 
				 bins.push_back( Bin( binRange ) );
				 for ( nchroms_t nchroms = bins[ binId ].getMinDerCount(); nchroms <= bins[ binId ].getMaxDerCount(); nchroms++ ) {
					 if ( nchroms+1 >= isize( derCount2bin ) ) derCount2bin.resize( nchroms+2, -1 );
					 derCount2bin[ nchroms ] = binId;
				 }
			 } else {
				 // !binStep.empty()
				 // nchroms_t binStep = getVal<nchroms_t>( binStep );
				 // ValRange<nchroms_t> rng( binRange );
				 // nchroms_t rmin = rng.getMin();
				 // nchroms_t rmax = rmin + binStep;
				 // nchroms_t rend = rng.getMax();

				 // while ( rmax <= rend ) {
					 
				 // }
				 
			 }
		 }  // for binDef in binDefs
	 }  // init()

	 void setSampleSize( nchroms_t nsam ) {
		 if ( nsam+1 >= isize( derCount2bin ) ) derCount2bin.resize( nsam+2, -1 );
	 }

	 void clear() { ForEach( Bin& bin, bins ) bin.clear(); }

	 void processSNP( snp_id_t snpId ) {
		 if ( snpCond( snpId ) ) {
			 nchroms_t derCount = derCounts[ snpId ];
			 chk( 0 <= derCount && derCount < int(derCount2bin.size()) );
			 int bin = derCount2bin[ derCount ];
			 if ( bin != -1 ) {
				 chk( 0 <= bin && bin < int(bins.size()) );
				 bins[ bin ].increment();
			 }
		 }
	 }

	 void writeHeadings( ostream& s ) const {
		 ForEach( const Bin& bin, bins )
				s << "\t" << MakeValidIdentifier( Join( "_", "afs", snpCond.getColSfx(), bin.getColSfx() ) );
	 }

	 void writeData( ostream& s, nsnps_t totSnps ) {
		 ForEach( const Bin& bin, bins ) {
			 s << "\t" << Frac( bin.getSnpCount(), totSnps); 
		 }
	 }

	 const vector<Bin>& getBins() const { return bins; }
	 const SNPCond& getSNPCond() const { return snpCond; }

	 const vector< int >& getDerCount2bin() const { return derCount2bin; }
   
	 
private:
	 // Field: bins
	 // All the bins of this AFS
	 vector< Bin > bins;

	 // Field: derCount2bin
	 // For each possible derived allele count, the corresponding bin in <bins>,
	 // or NULL if there is no bin for this count.
	 vector< int > derCount2bin;

	 // Field: snpCond
	 // Condition on SNPs going into this AFS
	 SNPCond snpCond;
	 
};  // class AFS

istream& operator>>( istream& is, AFS& afs );

istream& operator>>( istream& is, AFS& afs ) {
	string afsDef;
	is >> afsDef;
	afs.init( afsDef );
	
	return is;
}

ostream& operator<<( ostream& os, const AFS::Bin& bin );
ostream& operator<<( ostream& os, const AFS::Bin& bin ) {
	os << bin.getDerCountRange();
	return os;
}

ostream& operator<<( ostream& os, const AFS& afs );
ostream& operator<<( ostream& os, const AFS& afs ) {
	os << "[AFS: " << afs.getSNPCond() << " / " << afs.getBins() << "]";
	return os;
}

//
// ** Class SNPPairCond
//
// A condition on pairs of SNPs.
//
class SNPPairCond {
public:
	 
	 // Method: init
	 // Initialize this SNP pair condition from a string definition.
	 // Syntax:
	 // .2-.3/snpcond1/snpcond2
	 void init( const string& def ) {
		 string lenRangeDef, snpCond1def, snpCond2def;
		 Split( def, "/", lenRangeDef, snpCond1def, snpCond2def );
		 lenRange.init( lenRangeDef );
		 snpConds[0].init( snpCond1def );
		 snpConds[1].init( snpCond2def );
	 }
	 
	 bool operator() ( snp_id_t snp1, snp_id_t snp2 ) const {
		 return lenRange( fabs( posit[ snp2 ] - posit[ snp1 ] ) ) &&
				( ( snpConds[0]( snp1 ) && snpConds[1]( snp2 ) ) ||
					( snpConds[0]( snp2 ) && snpConds[1]( snp1 ) ) )
					;
	 }
	 
	 const ValRange<len_t>& getLenRange() const { return lenRange; }
	 const SNPCond& getSnpCond1() const { return snpConds[0]; }
	 const SNPCond& getSnpCond2() const { return snpConds[1]; }
	 
private:
	 // Field: lenRange
	 // Range physical distances for SNP pairs in the set.
	 ValRange<len_t> lenRange; 
	 
	 // Field: snpConds
	 // Conditions that must be satisfied by the two SNPs.
	 SNPCond snpConds[2];
	 
};  // class SNPPairCond

ostream& operator<<( ostream& os, const SNPPairCond& pairCond );
ostream& operator<<( ostream& os, const SNPPairCond& pairCond ) {
	os << Join( "_", ToString( pairCond.getLenRange() ), ToString( pairCond.getSnpCond1() ), ToString( pairCond.getSnpCond2() ) );
	return os;
}

namespace acc = boost::accumulators;
// using acc::acumulator_set;
// using acc::tag;
// using acc::stats;
// using acc::mean;
// using acc::variance;
// using acc::sum_kahan;
//using acc::extract;

typedef acc::accumulator_set<double, acc::stats< acc::tag::sum_kahan, acc::tag::mean,
																								 acc::tag::variance > > acc_t;

// typedef acc::accumulator_set<double, acc::stats< acc::tag::sum_kahan, acc::tag::mean, acc::tag::variance,
// 																								 acc::tag::extended_p_square > > acc2_t;

// typedef boost::shared_ptr< acc2_t > acc2_p;

std::vector< popid > popNames;
std::vector< nchroms_t > sampleSizes;
std::vector< nchroms_t > sampleStarts;
bool perPopStats = false;


//
// ** Class LD
//
// A subset of the possible SNP pairs, and LD info for it.
//
class LD {
	 
public:
	 
	 void init( string def ) { pairCond.init( def ); }
	 
	 void clear() {
		 r2_stats.clear();
		 r2_stats.resize( popNames.size() );
		 Dprime_stats.clear();
		 Dprime_stats.resize( popNames.size() );
	 }
	 
	 const SNPPairCond& getPairCond() const { return pairCond; }
	 
	 void writeHeadings( ostream& s ) const {
		 vector<string> statNames = MakeVec<string>( "mean", "var" );

		 for ( size_t popNum = 0; popNum < popNames.size(); ++popNum ) 
				for ( int which = 0; which < 2; which++ ) {
					ForEach( string statName, statNames ) {
						string accName = ( which == 0 ? "r2" : "Dprime" );
						s << "\t" << MakeValidIdentifier( Join( "_", "ld", ToString( pairCond ), accName, statName ) );
						if ( perPopStats ) s << popNames[ popNum ];
					}
				}
	 }
	 
	 void writeData( ostream& s ) const {
		 for ( size_t popNum = 0; popNum < popNames.size(); ++popNum ) 
				for ( int which = 0; which < 2; which++ ) {
					const acc_t *accum = ( which == 0 ? &r2_stats[popNum] : &Dprime_stats[popNum] );
					s << "\t" << acc::mean( *accum ); 
					s << "\t" << acc::variance( *accum ); 
					// s << "\t" << acc::extract::count( *accum ); 
					// s << "\t" << acc::sum_kahan( *accum ); 
				}
	 }
	 
	 const SNPPairCond& getSNPPairCond() const { return pairCond; }

	 void processSNPPair( snp_id_t snp1, snp_id_t snp2, char **snps ) {
		 // if ( ( fabs( posit[snp1] - .5 ) < 1e-5  ||
		 // 				fabs( posit[snp2] - .5 ) < 1e-5 ) )
		 if ( pairCond( snp1, snp2 ) ) {
			 //PRINT3( popNames.size(), sampleStarts.size(), sampleSizes.size() );
			 for ( size_t popNum = 0; popNum < sampleStarts.size(); ++popNum ) {
				 double r2, Dprime;
				 if ( compute_LD( snp1, snp2, sampleStarts[ popNum ], sampleSizes[ popNum ], snps, &r2, &Dprime ) ) {
					 //PRINT7( popNum, snp1, snp2, sampleStarts[popNum], sampleSizes[popNum], r2, Dprime );
					 r2_stats[ popNum ]( r2 );
					 Dprime_stats[ popNum ]( Dprime );
				 }
			 }
		 }
	 }
	 
private:
	 // Field: pairCond
	 // Condition that must be satisfied by SNP pairs in the set
	 SNPPairCond pairCond;

	 // Field: r2_stats
	 // Statistics on r2 for SNP pairs in the set.
	 std::vector< acc_t > r2_stats;

	 // Field: Dprime_stats
	 // Statistics on D' for SNP pairs in the set.
	 std::vector< acc_t > Dprime_stats;


//	 static const vector<prob_t> quantile_probs;
	 
};  // class LD


//const vector<prob_t> LD::quantile_probs = boost::assign::list_of(.1)(.5)(.9);

istream& operator>>( istream& is, LD& ld );
istream& operator>>( istream& is, LD& ld ) {
	string ldDef;
	is >> ldDef;
	ld.init( ldDef );
	return is;
}

ostream& operator<<( ostream& os, const LD& ld );
ostream& operator<<( ostream& os, const LD& ld ) {
	os << ld.getPairCond();
	return os;
}

// ** Class TreeEdge
class TreeEdge {
public:
	 TreeEdge( const char *line ) {
				chk( sscanf( line, "ARG %lf %lf %lf %lf %d", &beg, &end, &genMoreRecent, &genLessRecent, &leafsetSize ) == 5 );
	 }

	 bool containsLoc( loc_t loc ) const { return beg <= loc && loc <= end; }
	 
	 loc_t beg;
	 loc_t end;
	 genid genMoreRecent;
	 genid genLessRecent;
	 nchroms_t leafsetSize;
};

//
// ** Class TreeEdgeSet
//
// Gather stats for a collection of tree edges
//
class TreeEdgeSet {
public:
	 typedef boost::function<bool (const TreeEdge& treeEdge)> treeEdgePred_t;
	 
	 TreeEdgeSet( treeEdgePred_t pred_, string predName_ ): pred( pred_ ), predName( predName_ ) { clear(); }

	 void processTreeEdge( const TreeEdge& treeEdge ) {
		 if ( pred( treeEdge ) ) {

			 len_t segLen = treeEdge.end - treeEdge.beg;
			 gens_t edgeGens = treeEdge.genLessRecent - treeEdge.genMoreRecent;
			 
			 mutimes( segLen * edgeGens );
			 lens( segLen );
			 gens( edgeGens );
			 
			 maxTreeHeight = std::max( maxTreeHeight, treeEdge.genLessRecent );
		 }
	 }

	 void writeHeadings( ostream& s ) const {
		 s << "\t" << "treeEdge_" << predName << "_n";
		 vector<string> statNames = MakeVec<string>( "mutime", "len", "gens" );
		 vector<string> summNames = MakeVec<string>( "sum", "mean", "var" );
		 ForEach( string statName, statNames ) {
				ForEach( string summName, summNames ) {
					s << "\t" << MakeValidIdentifier( Join( "_", "treeEdge", predName, statName, summName ) );
				}
		 }
		 s << "\t" << "treeEdge_" << predName << "_maxTreeHeight";
	 }

	 void writeData( ostream& s ) const {
		 s << "\t" << acc::count( mutimes );

		s << "\t" << acc::sum_kahan( mutimes ) << "\t" << acc::mean( mutimes ) << "\t" << acc::variance( mutimes );
		s << "\t" << acc::sum_kahan( lens ) << "\t" << acc::mean( lens ) << "\t" << acc::variance( lens );
		s << "\t" << acc::sum_kahan( gens ) << "\t" << acc::mean( gens ) << "\t" << acc::variance( gens );

		s << "\t" << maxTreeHeight;
	 }

	 void clear() {
		 mutimes = acc_t();
		 lens = acc_t();
		 gens = acc_t();
		 maxTreeHeight = 0;
	 }

private:
	 treeEdgePred_t pred;
	 string predName;
	 
	 acc_t mutimes, lens, gens;

	 double maxTreeHeight;
	 
};  // class TreeEdgeSet

/////////////////////////////
// * Section: Histograms
/////////////////////////////

//
// ** Concept: Binner
//
// A mapping from values to bins.
//
template <class X>
struct BinnerConcept: public boost::Assignable<X>, boost::UnaryFunction<X,
																																				// X maps values of type X::value_type
																																				// to bin ids of type X::bin_id_type
																																				typename X::bin_id_type, typename X::value_type > {
	 typedef typename X::value_type value_type;
	 typedef typename X::bin_id_type bin_id_type;

	 BOOST_CONCEPT_USAGE(BinnerConcept) {
		 const X& binnerConst = this->binner;
		 binId = binnerConst( val );
		 
		 nbins = binnerConst.getNumBins();
		 val = binnerConst.getBinBeg( binId );
		 val = binnerConst.getBinEnd( binId );
	 }

private:
	 X binner;
	 value_type val;
	 bin_id_type binId;
	 size_t nbins;
};  // struct BinnerConcept

// Metafunction: DiffType
// Returns the type of the difference of values of the given type.
template <typename TVal> struct DiffType {
	 typedef BOOST_TYPEOF_TPL( boost::declval<TVal>() - boost::declval<TVal>()) type;
};

// ** Class: UniformBinner
// Maps values to bins by subdividing a specified range into equal-sized bins, and providing
// underflow and overflow bins for values below and above the range, respectively.
template <typename ValT = double, typename BinIdType = size_t>
class UniformBinner {
public:
	 typedef ValT value_type;
	 typedef BinIdType bin_id_type;
	 typedef typename DiffType<ValT>::type value_diff_type;

	 BOOST_CONCEPT_ASSERT((boost::LessThanComparable<value_type>));
	 BOOST_CONCEPT_ASSERT((boost::Convertible<value_diff_type,double>));
	 BOOST_MPL_ASSERT((boost::is_integral<bin_id_type>));
	 BOOST_MPL_ASSERT((boost::is_arithmetic<value_type>));

	 UniformBinner( ValT rangeMin_, ValT rangeMax_, size_t nbins_ ):
		 rangeMin( rangeMin_ ), rangeMax( rangeMax_ ), nbins( nbins_ + 2 ),
		 binSize( static_cast<double>( rangeMax_ - rangeMin_ ) / nbins_ ) { }

	 bin_id_type operator()( ValT val_ ) const {
		 return ( val_ < this->rangeMin ) ? 0 :
				( ( val_ > this->rangeMax ) ? nbins-1 :
					( 1 + static_cast<bin_id_type>( static_cast<double>( val_ - this->rangeMin ) / this->binSize ) ) );
	 }

	 size_t getNumBins() const { return nbins; }
	 ValT getBinBeg( size_t binId ) const { return binId == 0 ? boost::numeric::bounds<ValT>::lowest() : this->rangeMin + (binId-1) * binSize; }
	 ValT getBinEnd( size_t binId ) const { return binId == nbins-1 ? boost::numeric::bounds<ValT>::highest() : this->rangeMin + binId * binSize; }

private:
	 ValT rangeMin, rangeMax;
	 size_t nbins;
	 double binSize;
};  // class UniformBinner

BOOST_CONCEPT_ASSERT((BinnerConcept< UniformBinner<> >));

//
// ** Class: Histogram1D
//
// A one-dimensional histogram.
//
// Template params:
//
//   BinnerT - maps values to bins.
//   StatKeeperT - keeps stats for values within each bin
template <typename BinnerT = UniformBinner<>, typename StatKeeperT = StatKeeper< typename BinnerT::value_type > >
class Histogram1D {
public:
	 BOOST_CONCEPT_ASSERT((BinnerConcept<BinnerT>));
	 
	 typedef BinnerT binner_type;
	 typedef StatKeeperT stat_keeper_type;
	 typedef typename stat_keeper_type::value_type value_type;

	 Histogram1D( const binner_type& binner_ ):
		 binner( binner_ ), statKeepers( binner_.getNumBins() ) { }

	 void operator()( value_type val_ ) {
				statKeepers[ binner( val_ ) ]( val_ );
	 }
	 void operator()( typename binner_type::value_type binnerVal_, value_type val_ ) {
				statKeepers[ binner( binnerVal_ ) ]( val_ );
	 }

	 static void writeHeaders( ostream& s ) {
		 s << "binBeg\tbinEnd\tcount\tsum\tsumSq\tnumNaNs\tnumInfs\n";
	 }
	 
	 void save( ostream& s, bool dropEmptyBins = false ) const {
		 for ( size_t binId = 0; binId < statKeepers.size(); binId++ ) {
			 const stat_keeper_type& sk = statKeepers[ binId ];
			 if ( !dropEmptyBins || sk.getNumVals() > 0 ) 
					s << binner.getBinBeg( binId ) << "\t" << binner.getBinEnd( binId ) << "\t"
						<< sk.getNumVals() << "\t" << sk.getSum() << "\t" <<
						 sk.getSumSq() << "\t" << sk.getNumNaNs() << "\t" <<
						 sk.getNumInfs() << 
						 "\n";
		 }
	 }

	 void save( filename_t f, bool dropEmptyBins = false ) const {
		 std::ofstream s( f );
		 s.exceptions( std::ios::failbit | std::ios::badbit );
		 writeHeaders( s );
		 save( s, dropEmptyBins );
	 }

private:
	 binner_type binner;
	 std::vector< stat_keeper_type > statKeepers;
};  // class Histogram1D


//
// ** Class: Histogram2D
// A two-dimensional histogram
template <typename Binner1T, typename Binner2T, typename count_type = size_t>
class Histogram2D {
public:
	 typedef Binner1T binner1_type;
	 typedef Binner2T binner2_type;
	 typedef typename Binner1T::value_type value1_type;
	 typedef typename Binner2T::value_type value2_type;

	 Histogram2D( const Binner1T& binner1_, const Binner2T& binner2_ ):
		 binner1( binner1_ ), binner2( binner2_ ), counts( binner1_.getNumBins() * binner2_.getNumBins(), 0 ) { }

	 void operator()( value1_type val1_, value2_type val2_ ) {
		 if ( !( (boost::math::isnan)( ToDouble( val1_ ) ) ) &&
					!( (boost::math::isnan)( ToDouble( val2_ ) ) ) )
				counts[ flatBinId( binner1( val1_ ), binner2( val2_ ) ) ]++;
	 }

	 void save( filename_t f ) const {
		 std::ofstream s( f );
		 s.exceptions( std::ios::failbit | std::ios::badbit );
		 s << "binId1\tbinBeg1\tbinEnd1\tbinId2\tbinBeg2\tbinEnd2\tcount\tflatBinId\n";
		 for ( size_t binId1 = 0; binId1 < binner1.getNumBins(); binId1++ )
				for ( size_t binId2 = 0; binId2 < binner2.getNumBins(); binId2++ )
					 s <<
							binId1 << "\t" << binner1.getBinBeg( binId1 ) << "\t" << binner1.getBinEnd( binId1 ) << "\t" <<
							binId2 << "\t" << binner2.getBinBeg( binId2 ) << "\t" << binner2.getBinEnd( binId2 ) << "\t" <<
							
							counts[ flatBinId( binId1, binId2 ) ] << "\t" <<
							flatBinId( binId1, binId2 ) << "\n";
	 }

private:

	 size_t flatBinId( size_t id1_, size_t id2_ ) const { return id1_ * binner2.getNumBins() + id2_; }

	 binner1_type binner1;
	 binner2_type binner2;
	 std::vector< count_type > counts;
};  // class Histogram2D

void histogramTest() {
	Histogram1D< UniformBinner<double> > h( UniformBinner<double>( 0.0, 10.0, 5 ) );
	h( 2 );
	h( 3 );
	h( 2 );
	h( 0 );
	h( -3 );
	h( 100 );
	h( 8.3 ); 
	h.save( "h1.tsv" );

	Histogram2D< UniformBinner<double>,
							 UniformBinner<double> > h2( UniformBinner<double>( 0.0, 10.0, 5 ),
																					 UniformBinner<double>( 1000, 10000, 15 ) );

	h2( 3, 4 );
	h2( 100, 1000000 );
	h2.save( "h2.tsv" );
	
}

// * fst computation

// ** Function ~compute_fst_for_one_snp~ - compute a Weir-Hill estimator of Fst for one SNP
//
// Params:
//
//   - nai,naj :: counts of the two alleles in pop i,j
//
// Returns: estimator of Fst, or NAN if can't be computed
//
double compute_fst_for_one_snp( const nchroms_t nai[2], const nchroms_t naj[2] ) {
	// adapted from jvitti's code
	if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) return NAN;
	nchroms_t ni = nai[0] + nai[1], nj = naj[0] + naj[1];
	if (ni == 0 || nj == 0) return NAN;
	
	// Weir-Hill estimator
	freq_t p[2] = { (double) nai[0] / ni, (double) naj[0] / nj };
	double pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
	double nic = ni - (double) ni * ni / (ni + nj);
	double njc = nj - (double) nj * nj / (ni + nj);
	nchroms_t nc = nic + njc;
	double msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
	double msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
	double num = msp - msg;
	double denom = msp + (nc - 1) * msg;

	return denom == 0 ? static_cast<double>( NAN ) : ( num / denom );
}  // double compute_fst_for_one_snp( const nchroms_t nai[2], const nchroms_t naj[2] )


//////////////////////////////////////////

bool StartsWith( const char *s, const char *pfx ) { return strlen( s ) >= strlen( pfx ) && !strncmp( s, pfx, strlen( pfx ) ); }

int  segsub( int nsam, int segsites, char **list ) ;

// * Main

int sample_stats_main(int argc, char *argv[])
{
// ** var decls
	
	using std::map;
	using std::string;
	
	int nsam, i,  howmany  ;
	const size_t MAX_LINE_LEN = 65000;
	char **list, line[MAX_LINE_LEN+1];
	filename_t simFile, outFile;
	FILE *pfin ;
	nsnps_t   segsites;
	int count;
	double pi , h, th  ,prob ;
	char dum[20], astr[100] ;
	int precision(8);

	vector< pair< double, double > > snpPairDists;

	len_t MARGIN(0.0);
	
	cerr.precision(8);

 	StatKeeper<> piStats, nsitesStats, dindStats;
	vector< StatKeeper<> > statKeepers;

	namespace po = boost::program_options;

	vector<AFS> AFSs;
	vector<LD> LDs;
	vector<TreeEdgeSet> treeEdgeSets;

	po::positional_options_description positional_opts;

	bool includeTreeStats = false;
	bool summaryOnly = false;
	bool ld_use_cM = false;
	int progressEvery = 0;

	ValRange<nchroms_t> chromRange( 0, 100000 );

	//double regionLen_bp(-1);
	double regionLen_cM(-1);

	loc_t dindLoc = 0.0;
	snp_id_t dindLocIdx = -1;

	string globalAFS_fname, globalLD_r2_fname, globalLD_Dprime_fname;
	vector< nsnps_t > globalLD_snpCountDists;

	string ldSepsDef;
	vector< nsnps_t >ldSeps;

	typedef int nsims_t;

	nsims_t maxSims = 0;

	boost::random::mt19937::result_type rng_seed;

	bool addQuantiles = false;

// ** program options
	
	po::options_description desc("Allowed options");
	desc.add_options()
		 ("help,h", "produce help message")

		 ("simfile,i", po::value(&simFile), "name of file containing simulator output (default is to read from stdin)")
		 ("outfile,o", po::value(&outFile), "name of output file (default is to write to stdout)")
			 
		 ("precision,p", po::value(&precision)->default_value(8), "specify precision of output")
		 
		 ("margin,m", po::value(&MARGIN)->default_value(.1), "ignore part of simulated region within this margin of ends of region")
		 ("chromRange,c", po::value( &chromRange ), "restrict chroms to this range" )
		 
		 ("afs,a", po::value(&AFSs)->composing(), "compute allele frequency spectrum")
		 ("ld,l", po::value(&LDs)->composing(), "compute LD stats for these definitions")
		 ("ld-use-cM", po::bool_switch(&ld_use_cM), "use cM for LD snp distances")

		 ("ld-seps", po::value(&ldSepsDef), "snp count separation for LD stats" )
		 
		 ("DIND,D", po::value(&dindLoc), "compute DIND for this loc" )
		 
		 ("tree,t", po::bool_switch(&includeTreeStats), "include tree stats")
		 
		("per-pop-stats", po::bool_switch(&perPopStats), "break down stats by pop")
		("summary-only,s", po::bool_switch(&summaryOnly), "print summary only")
		 ("global-afs", po::value(&globalAFS_fname), "save global AFS to this file")
		 ("global-ld-r2", po::value(&globalLD_r2_fname), "save global LD histogram for r^2 to this file")
		 ("global-ld-Dprime", po::value(&globalLD_Dprime_fname), "save global LD histogram for Dprime to this file")
		 ("global-ld-dists-snp-count", po::value( &globalLD_snpCountDists )->composing(), "for LD computation, use SNP pairs at this distance (in # of snps between them)" )
		 //("seed", po::value( &rng_seed )->default_value( 373737 ), "seed for random choices" )
		 ("max-sims", po::value( &maxSims )->default_value( 0 ), "max sims to consider" )
		 ("quantiles,q", po::bool_switch( &addQuantiles), "add quantiles" )
		 
		("progress-every,g", po::value(&progressEvery), "print progress every N sims")
		 ;

	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_opts).run(), vm);
		po::notify(vm);    
	} 
	catch(std::exception& e) {
		cerr << e.what() << "\n";
		return EXIT_FAILURE;
	}

	if ( vm.count( "help" ) ) { cerr << desc; return EXIT_FAILURE; }

	vector<string> ldSepParts = Split( ldSepsDef, "," );
	BOOST_FOREACH( string ldSepPart, ldSepParts )
		 ldSeps.push_back( getVal<int>( ldSepPart ) );
	boost::sort( ldSeps );

	//boost::random::mt19937 rng( rng_seed == 0 ? time( NULL ) : rng_seed );

	std::sort( globalLD_snpCountDists.begin(), globalLD_snpCountDists.end() );
	
	vector< Histogram1D<> > globalLD_r2, globalLD_Dprime;
	BOOST_FOREACH( nsnps_t dist, globalLD_snpCountDists ) {
		double dist_dbl = dist;
		globalLD_r2.push_back( Histogram1D<>( UniformBinner<>( dist_dbl - .5, dist_dbl + .5, /* nbins= */ 1 ) ) );
		globalLD_Dprime.push_back( Histogram1D<>( UniformBinner<>( dist_dbl - .5, dist_dbl + .5, /* nbins= */ 1 ) ) );
	}
	
	if ( includeTreeStats ) {
		namespace phx = boost::phoenix;
		using phx::placeholders::arg1;
		using phx::bind;
		using phx::val;
		
		treeEdgeSets.push_back( TreeEdgeSet( val(true), "all" ) );
		treeEdgeSets.push_back( TreeEdgeSet( bind( &TreeEdge::containsLoc, arg1, .5 ), "mid" ) );
		treeEdgeSets.push_back( TreeEdgeSet( bind( &TreeEdge::genMoreRecent, arg1 ) == 0, "bot" ) );
		treeEdgeSets.push_back( TreeEdgeSet( bind( &TreeEdge::genMoreRecent, arg1 ) == 0 && bind( &TreeEdge::containsLoc, arg1, .5 ), "botmid" ) );
	}

	// chk( ( argc & 1 ) );

	// double MARGIN = getVal<double>( argv[1] );
	// int afsMax = getVal<int>( argv[2] );
	
	// for ( int ni = 3; ni < argc; ni+=2 )
	// 	snpPairDists.push_back( make_pair( getVal<double>( argv[ni] ), getVal<double>( argv[ni+1] ) ) );


/* read in first two lines of output  (parameters and seed) */
  pfin = simFile.empty() ? stdin : fopen( simFile.c_str(), "rt" );
	if ( !pfin ) throw std::logic_error( "sample_stats: could not open file " + simFile );

	std::ofstream outStrm;
	outStrm.exceptions( std::ios::failbit | std::ios::badbit );
	if ( !outFile.empty() ) outStrm.open( outFile.c_str() );
	ostream& fout = outFile.empty() ? cout : outStrm;

	fout.precision( precision );
	
  chk( fgets( line, MAX_LINE_LEN, pfin) );
	cerr << line;
  chk( sscanf(line," %s  %d %d", dum,  &nsam, &howmany) == 3 );

	if ( maxSims > 0 && howmany > maxSims ) howmany = maxSims;

	nchroms_t nsam_orig = nsam;
	chromRange.setMax( std::min( chromRange.getMax(), nsam-1 ) );

	nsam = chromRange.getMax() - chromRange.getMin() + 1;
	PRINT3( chromRange, nsam_orig, nsam );

	popNames.push_back( popid(1) );
	sampleStarts.push_back( 0 );
	sampleSizes.push_back( nsam );
	
  chk( fgets( line, MAX_LINE_LEN, pfin) );
	
	if ( boost::algorithm::starts_with( line, "pops" ) ) {
		popNames.clear();
		sampleStarts.clear();
		sampleSizes.clear();
		std::istringstream is( line );
		is.exceptions( std::ios::failbit | std::ios::badbit );

		std::string popsLineHeader;

		try {
		
			size_t npops;
			is >> popsLineHeader >> npops;

			nchroms_t sampleStart = 0;
			for ( size_t popNum = 0; popNum < npops; ++popNum ) {
				popid popName;
				nchroms_t sampleSize;
				is >> popName >> sampleSize;
				popNames.push_back( popName );
				sampleSizes.push_back( sampleSize );
				sampleStarts.push_back( sampleStart );
				sampleStart += sampleSize;
				cerr << "got pop name " << popName << " size " << sampleSize << "\n";
			}
		} catch( std::ios_base::failure const& e ) {
			 BOOST_THROW_EXCEPTION( cosi_io_error()
															<< error_msg( "error parsing population info line" ) );
		}
		chk( fgets( line, MAX_LINE_LEN, pfin) );
	}
	//cerr << line;
	chk( boost::accumulate( sampleSizes, 0 ) == nsam );

  list = cmatrix(nsam,maxsites+1);
	char **trimmed_list = ( char **)malloc( nsam * sizeof( char * ) );
  posit = (loc_t *)malloc( maxsites*sizeof( loc_t ) ) ;
  posit_gloc = (gloc_t *)malloc( maxsites*sizeof( gloc_t ) ) ;

	chk( trimmed_list && posit && posit_gloc );

	//boost::shared_ptr< vector< boost::shared_ptr< string > > > statNames = boost::make_shared< vector< boost::shared_ptr< string > > >();

	//vector< nsnps_t > globalAFS( nsam + 1, 0 );

	Histogram1D<> globalAFS( UniformBinner<>( 0.0, 1.0, nsam ) );

// ** main loop: for each sim

// *** misc

	PRINT( LDs.size() );

  count=0;
	while( howmany-count++ ) {

		if ( progressEvery > 0 && !( (count-1) % progressEvery ) )
			cerr << "sim " << (count-1) << " of " << howmany << " (" << static_cast<int>( static_cast<double>(count-1) / static_cast<double>( howmany ) * 100 ) << "%)" << endl;

		vector< acc_t > ld_r2( ldSeps.size() ), ld_Dprime( ldSeps.size() );
		vector< vector< double > > ld_r2_vals( ldSeps.size() ), ld_Dprime_vals( ldSeps.size() );

		boost::array<prob_t,5> probs = { .10, .25, .50, .75, .90 };

		// for ( size_t i = 0; i < ldSeps.size(); i++ ) {
		// 	ld_r2.push_back( boost::make_shared< acc2_t >( acc::extended_p_square::probabilities = probs ) );
		// 	ld_Dprime.push_back( boost::make_shared< acc2_t >( acc::extended_p_square::probabilities = probs ) );
		// }

		//
		// Find the beginning of one simulation.  Simulations begin with
		// a double forward slash at the beginning of a line.
		//
		do { chk( fgets( line, MAX_LINE_LEN, pfin) ); }
		while ( !StartsWith( line, "//" ) );
		if ( StartsWith( line, "// cosi-early-exit" ) ) {
			std::cerr << "sample_stats_extra: simulation stopped early after " << (count-1) << " of " << howmany << " sims\n";
			break;
		}

		//
		// Read statistics gathered for the simulation, and/or ARG edges
		// output for this simulation.
		//


		int statIdx = 0;
		vector<double> statVals;
		set<loc_t> recombLocs, recombBegs, recombEnds;

		ForEach( TreeEdgeSet& treeEdgeSet, treeEdgeSets ) treeEdgeSet.clear();

// *** process simulation headers
		
		while( true ) {
			chk( fgets( line, MAX_LINE_LEN, pfin) );

			if ( StartsWith( line, "segsites:" ) )
				 break;
			else if( StartsWith( line, "prob:" ) ) {
				chk( sscanf( line, "  prob: %lf", &prob ) == 1 );
			}
			else if( boost::algorithm::starts_with( line, "region_len_cM" ) ) {
				chk( sscanf( line, "region_len_cM: %lf", &regionLen_cM ) == 1 );
			}
// **** process ARG edges
			else if ( StartsWith( line, "ARG" ) ) {

				TreeEdge treeEdge( line );

				ForEach( TreeEdgeSet& treeEdgeSet, treeEdgeSets ) treeEdgeSet.processTreeEdge( treeEdge );
				
				recombLocs.insert( treeEdge.beg );
				recombLocs.insert( treeEdge.end );
				recombBegs.insert( treeEdge.beg );
				recombEnds.insert( treeEdge.end );
			}
// **** process generic stats
			else if ( StartsWith( line, "stat " ) ) {
				string dummy, statName;
				double statVal;
				istringstream istrm( line );
				istrm.exceptions( std::ios::failbit | std::ios::badbit );
				istrm >> dummy >> statName >> statVal;
				//cerr << "statName=" << statName << " statVal=" << statVal << endl;
				// if ( count == 1 )
				// 	statNames->push_back( boost::make_shared<string>( statName.c_str() ) );
				// else
				// 	chk( *( statNames->at( statIdx ) ) == statName );
				statVals.push_back( statVal );
				if ( count == 1 ) statKeepers.push_back( StatKeeper<>( statName ) );
				else chk( statKeepers[ statIdx ].getName() == statName );
				statKeepers[ statIdx++ ].add( statVal );
			} else
				 throw std::logic_error( string("unknown input line: ") + string( line ) );
		}
// *** read SNP locs
		chk( sscanf( line, "  segsites: %d", &segsites ) == 1 );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (loc_t *)realloc( posit, maxsites*sizeof( loc_t ) ) ;
			chk( posit );
			posit_gloc = (gloc_t *)realloc( posit_gloc, maxsites*sizeof( gloc_t ) ) ;
			chk( posit_gloc );
			biggerlist(nsam,maxsites, list) ;
		}

		boost::scoped_array<char> list_dummy( new char[ segsites + 1 ] );

		int first_snp = 0, last_snp = segsites-1, trimmed_segsites = segsites;
		const loc_t *trimmed_posit = posit;
		const gloc_t *trimmed_posit_gloc = posit_gloc;

// *** read SNP data
		{
			if ( segsites > 0 ) {
				chk( fscanf(pfin," %s", astr) == 1 );
				chk( !strcmp( astr, "positions:" ) );
				for( i=0; i<segsites ; i++) {
					chk( fscanf(pfin," %lf",posit+i) == 1 );
					chk(0.0 <= *(posit+i) && *(posit+i) <= 1.0 );
				}

				if ( regionLen_cM > 0 ) {
					chk( fscanf(pfin," %s", astr) == 1 );
					chk( !strcmp( astr, "positions_genMap:" ) );
					for( i=0; i<segsites ; i++) {
						chk( fscanf(pfin," %lf",posit_gloc+i) == 1 );
						chk(0.0 <= *(posit_gloc+i) && *(posit_gloc+i) <= 1.0 );
					}
					
				}
			}
			first_snp = 0;
			while( first_snp < segsites && posit[ first_snp ] < MARGIN ) first_snp++;
			last_snp = segsites-1;
			while( last_snp > first_snp && posit[ last_snp ] > ( 1.0 - MARGIN ) ) last_snp--;

			//chk( 0 <= first_snp && first_snp < last_snp && last_snp < segsites );
			
			trimmed_segsites = last_snp - first_snp + 1;
			
			trimmed_posit = posit + first_snp;
			trimmed_posit_gloc = posit_gloc + first_snp;
			derCounts.clear();
			derCounts.resize( trimmed_segsites, 0 );

			if ( dindLoc > 0 ) {
				derCounts_DIND_der.clear();
				derCounts_DIND_der.resize( trimmed_segsites, 0 );
				derCounts_DIND_anc.clear();
				derCounts_DIND_anc.resize( trimmed_segsites, 0 );
				const loc_t *dindLocPtr = std::find( trimmed_posit, trimmed_posit + trimmed_segsites, dindLoc );
				chk( dindLocPtr < trimmed_posit + trimmed_segsites );
				dindLocIdx = dindLocPtr - trimmed_posit;
				PRINT3( dindLocIdx, *dindLocPtr, trimmed_posit[ dindLocIdx ] );
			}

			if ( segsites > 0 ) {
				ostringstream fmtStream;
				fmtStream << " " << "%" << ( segsites+1 ) << "s" << std::ends;
				string fmt = fmtStream.str();
				
				for( i=0; i<nsam_orig;i++) {
					//PRINT4( "reading sample ", i, "of", nsam_orig );
					memset( list_dummy.get(), 0, segsites+1 );
					int nread = fscanf(pfin, fmt.c_str(), list_dummy.get() );
					//PRINT4( i, nread, strlen(list_dummy), segsites );
					chk( nread == 1 );
					chk( int( strlen( list_dummy.get() ) ) == segsites );


					//PRINT2( i, chromRange( i ) );
					if ( !chromRange( i ) ) continue;
					strcpy( list[i], list_dummy.get() );
					
					trimmed_list[ i ] = list[ i ] + first_snp;
					
					size_t num = trimmed_segsites;
					const char *p = trimmed_list[ i ];
					while( const char *p_new = (const char *)memchr( p, '1', num ) ) {
						chk( *p_new == '1' );
						chk( size_t( p_new - p ) <= num );
						num -= ( p_new - p );
						int snpId = p_new - trimmed_list[i];
						chk( 0 <= snpId && snpId < trimmed_segsites );
						derCounts[ snpId ]++;
						if ( dindLoc > 0 ) {
							( trimmed_list[i][ dindLocIdx ] == '1' ? derCounts_DIND_der : derCounts_DIND_anc )[ snpId ]++;
						}
						p = p_new+1;
						if ( num >= 1 ) num--;
					}
				}  // read each sample
				
				//	 list[i][ last_snp+1 ] = '\0';
			}  // if ( segsites > 0 )
		}  // if ( segsites > 0 )

		

// *** analyse sample ( do stuff with segsites and list)

		if ( !globalAFS_fname.empty() ) {
			ForEach( nchroms_t thisSnpCount, derCounts ) globalAFS( ((double)thisSnpCount) / ((double)nsam), 1 );
		}

		if ( !ldSeps.empty() ) {
			size_t highestDistIdx = ldSeps.size();
			for ( snp_id_t snp1 = 0; snp1 < trimmed_segsites; snp1++ ) {
				for ( size_t distIdx = 0; distIdx < highestDistIdx; distIdx++ ) {
					nsnps_t dist = ldSeps[ distIdx ];
					snp_id_t snp2 = snp1 + dist;
					if ( snp2 >= trimmed_segsites-2 )
						 highestDistIdx = distIdx;
					else {
						double r2, Dprime;
						if ( compute_LD( snp1, snp2, nsam, trimmed_list, &r2, &Dprime ) ) {
							ld_r2[ distIdx ]( r2 );
							ld_Dprime[ distIdx ]( Dprime );
							if ( addQuantiles ) {
								ld_r2_vals[ distIdx ].push_back( r2 );
								ld_Dprime_vals[ distIdx ].push_back( Dprime );
							}
						}
					}
				}
			}
		}
		

		if ( !globalLD_snpCountDists.empty() ) {
			size_t highestDistIdx = globalLD_snpCountDists.size();
			for ( snp_id_t snp1 = 0; snp1 < trimmed_segsites; snp1++ ) {
				for ( size_t distIdx = 0; distIdx < highestDistIdx; distIdx++ ) {
					nsnps_t dist = globalLD_snpCountDists[ distIdx ];
					snp_id_t snp2 = snp1 + dist;
					if ( snp2 >= trimmed_segsites-2 )
						 highestDistIdx = distIdx;
					else {
						double r2, Dprime;
						if ( compute_LD( snp1, snp2, nsam, trimmed_list, &r2, &Dprime ) ) {
							globalLD_r2[ distIdx ]( dist, r2 );
							globalLD_Dprime[ distIdx ]( dist, Dprime );
						}
					}
				}
			}
		}

// *** if first sim, initialize stats keepers 
		
		if ( count == 1 && !summaryOnly ) {
//			fout << "pi\tss\tD\ttheta\tH\tnrecombLocs\tnrecombBegs\tnrecombEnds";
			fout << "pi\tss\tD\ttheta\tH\tnucdiv";
			BOOST_FOREACH( const StatKeeper<>& sk, statKeepers )
				fout << "\t" << sk.getName();

			ForEach( const TreeEdgeSet& treeEdgeSet, treeEdgeSets ) {
				treeEdgeSet.writeHeadings( fout );
			}

			ForEach( AFS& afs, AFSs ) {
				afs.writeHeadings( fout );
				afs.setSampleSize( nsam );
			}

			for ( size_t popNum1 = 0; popNum1 < popNames.size(); ++popNum1 )
				for ( size_t popNum2 = popNum1+1; popNum2 < popNames.size(); ++popNum2 )
					fout << "\tfst_" << popNames[ popNum1 ] << "_" << popNames[ popNum2 ];

			BOOST_FOREACH( nsnps_t ldSep, ldSeps ) {
				 fout
						<< "\tld_sep" << ldSep << "_r2_mean"
						<< "\tld_sep" << ldSep << "_r2_std"
						<< "\tld_sep" << ldSep << "_Dprime_mean"
						<< "\tld_sep" << ldSep << "_Dprime_std";
				 if ( addQuantiles ) {
					 BOOST_FOREACH( prob_t p, probs ) {
						 fout
								<< "\tld_sep" << ldSep << "_r2_q" << p
								<< "\tld_sep" << ldSep << "_Dprime_q" << p;
					 }
				 }
			}
			
			ForEach( const LD& ld, LDs ) ld.writeHeadings( fout );
			fout << "\n";
		}  // if first sim

// *** compute basic stats
		
		double tajd_val( 0.0 );

		pi = nucdiv(nsam, trimmed_segsites, trimmed_list) ;
		h = hfay(nsam, trimmed_segsites, trimmed_list) ;
		th = thetah(nsam, trimmed_segsites, trimmed_list) ;
		tajd_val = tajd(nsam,trimmed_segsites,pi);

		double dind_val = 0.0;
		if ( dindLoc > 0 ) {
			double sumDer = 0.0, sumAnc = 0.0;
			nchroms_t nder = derCounts[ dindLocIdx ];
			nchroms_t nanc = nsam - nder;
			PRINT3( nsam, nder, nanc );
			for ( snp_id_t  i = 0; i < trimmed_segsites; i++ ) {
				nchroms_t derCountsHere = derCounts_DIND_der[ i ];
				sumDer += derCountsHere * ( nder - derCountsHere );
				nchroms_t ancCountsHere = derCounts_DIND_anc[ i ];
				sumAnc += ancCountsHere * ( nanc - ancCountsHere );
			}
			PRINT2( sumDer, sumAnc );
			sumDer /= ( nder * ( nder-1.0 ) / 2.0 );
			sumAnc /= ( nanc * ( nanc-1.0 ) / 2.0 );
			dind_val = sumDer / sumAnc;
		}

		if ( !summaryOnly )
			fout << pi << "\t" << trimmed_segsites << "\t" << tajd_val << "\t"
					 << th << "\t"  << h << "\t" << Frac( pi, trimmed_segsites ); // << "\t" << recombLocs.size() << "\t" << recombBegs.size() << "\t" << recombEnds.size();

		piStats.add( pi );
		nsitesStats.add( static_cast< cosi_double >( trimmed_segsites ) );
		dindStats.add( dind_val );

		if ( !summaryOnly ) {
			ForEach( double statVal, statVals ) fout << "\t" << statVal;

			ForEach( const TreeEdgeSet& treeEdgeSet, treeEdgeSets ) treeEdgeSet.writeData( fout );
		}

// *** compute AFS

		ForEach( AFS& afs, AFSs ) {
			afs.clear();
			for ( snp_id_t ii = 0; ii < trimmed_segsites; ii++ )
				 afs.processSNP( ii );
			if ( !summaryOnly ) afs.writeData( fout, trimmed_segsites );
		}

// *** compute FSTs
		{
			size_t numPopPairs = popNames.size() * ( popNames.size() - 1 ) / 2; 
			std::vector<acc_t> FSTs( numPopPairs );
			for ( snp_id_t snp = 0; snp < trimmed_segsites; snp++ ) {

				static vector< nchroms_t > popDerCounts;
				popDerCounts.clear();
				popDerCounts.resize( popNames.size() );
				for ( size_t popNum = 0; popNum < popNames.size(); ++popNum ) {
					nchroms_t popStart = sampleStarts[ popNum ];
					for ( nchroms_t chrom = 0; chrom < sampleSizes[ popNum ]; ++chrom ) {
						if ( trimmed_list[ snp ][ popStart + chrom ] == '1' )
							 ++popDerCounts[ popNum ];
					}
				}

				std::vector<acc_t>::iterator popPairIter = FSTs.begin();
				for ( size_t popNum1 = 0; popNum1 < popNames.size(); ++popNum1 ) {
					for ( size_t popNum2 = popNum1+1; popNum2 < popNames.size(); ++popNum2 ) {
						nchroms_t nai[2] = { popDerCounts[ popNum1 ], sampleSizes[ popNum1 ] - popDerCounts[ popNum1 ] };
						nchroms_t naj[2] = { popDerCounts[ popNum2 ], sampleSizes[ popNum2 ] - popDerCounts[ popNum2 ] };
						//PRINT4( nai[0], nai[1], naj[0], naj[1] );
						double fstHere = compute_fst_for_one_snp( nai, naj );
						if ( !(boost::math::isnan)( fstHere ) ) {
							(*popPairIter)( fstHere );
						}
						++popPairIter;
					}
				}
			}

			std::vector<acc_t>::const_iterator popPairIter = FSTs.begin();
			for ( size_t popNum1 = 0; popNum1 < popNames.size(); ++popNum1 )
				for ( size_t popNum2 = popNum1+1; popNum2 < popNames.size(); ++popNum2 )
					fout << "\t" << acc::mean( *popPairIter++ );
		}
		

		{
			for ( size_t ii = 0; ii < ldSeps.size(); ++ii ) {
				fout
					 << "\t" << acc::mean( ld_r2[ ii ] )
					 << "\t" << sqrt( acc::variance( ld_r2[ ii ] ) )
					 << "\t" << acc::mean( ld_Dprime[ ii ] )
					 << "\t" << sqrt( acc::variance( ld_Dprime[ ii ] ) );


				 if ( addQuantiles ) {
					 boost::sort( ld_r2_vals[ ii ] );
					 boost::sort( ld_Dprime_vals[ ii ] );
					 // {
					 // 	 std::ostringstream dprimeFNstrm;
					 // 	 dprimeFNstrm << "dprime_sep" << ii << "_sim" << count;
					 // 	 std::string dprimeFNstring( dprimeFNstrm.str() );
					 // 	 std::ofstream dprime( dprimeFNstring.c_str() );
					 // 	 dprime << "dprime\n";
					 // 	 BOOST_FOREACH( double d, ld_Dprime_vals[ ii ] )
					 // 			dprime << d << "\n";
					 // }
					 BOOST_FOREACH( prob_t p, probs ) {
						 
						 fout
								<< "\t" << ld_r2_vals[ ii ][ boost::algorithm::clamp<size_t>( static_cast< size_t >( p * static_cast< double >( ld_r2_vals[ ii ].size() ) ),
																																							0, ld_r2_vals[ ii ].size()-1 ) ]
								<< "\t" << ld_Dprime_vals[ ii ][ boost::algorithm::clamp<size_t>( static_cast< size_t >( p * static_cast< double >( ld_Dprime_vals[ ii ].size() ) ),
																																									0, ld_Dprime_vals[ ii ].size()-1 ) ];

					 }
				 }
				
			}
		}

// *** compute LD stats
		//PRINT( LDs.size() );
		ForEach ( LD& ld, LDs ) {
			//PRINT2( "computing LD", ld );
			ld.clear();
			ValRange<len_t> lenRange = ld.getSNPPairCond().getLenRange();
			const double *pos = ( !ld_use_cM ? trimmed_posit : trimmed_posit_gloc );
			int snp2 = 1;
			for ( int snp1 = 0; snp1 < trimmed_segsites; snp1++ ) {
				while( snp2 < trimmed_segsites && ( snp2 <= snp1 || pos[snp2]-pos[snp1] < lenRange.getMin() ) )
					 snp2++;
				int snp2_tmp = snp2;
				while( snp2_tmp < trimmed_segsites && pos[snp2_tmp]-pos[snp1] < lenRange.getMax() ) {
					chk( snp1 != snp2_tmp );
					ld.processSNPPair( snp1, snp2_tmp, trimmed_list );
					snp2_tmp++;
				}
			}
			if ( !summaryOnly ) ld.writeData( fout );
		}

		if ( !summaryOnly ) fout << "\n";
		
  }  // while( howmany-count++ )
	
// * post-loop processing

	if ( !globalAFS_fname.empty() ) globalAFS.save( globalAFS_fname, /* dropEmptyBins= */ false );
	if ( !globalLD_r2_fname.empty() ) {
		std::ofstream s( globalLD_r2_fname.c_str() );
		s.exceptions( std::ios::failbit | std::ios::badbit );

		
		Histogram1D<>::writeHeaders( s );
		BOOST_FOREACH( const Histogram1D<>& h, globalLD_r2 )
			 h.save( s, /* dropEmptyBins= */ true );
	}
	if ( !globalLD_Dprime_fname.empty() ) {
		std::ofstream s( globalLD_Dprime_fname.c_str() );
		s.exceptions( std::ios::failbit | std::ios::badbit );
		
		Histogram1D<>::writeHeaders( s );
		BOOST_FOREACH( const Histogram1D<>& h, globalLD_Dprime )
			 h.save( s, /* dropEmptyBins= */ true );
	}

	cerr << "\n";
	cerr << "pi.mean=" << piStats.getMean() << " pi.std=" << piStats.getStd() << " pi.count=" << piStats.getNumVals() << "\n";
	cerr << "nsites.mean=" << nsitesStats.getMean() << " nsites.std=" << nsitesStats.getStd() << " nsites.count=" << nsitesStats.getNumVals() << "\n";
	if ( dindLoc > 0 ) 
		 cerr << "dind.mean=" << dindStats.getMean() << " dind.std=" << dindStats.getStd() << " dind.count=" << dindStats.getNumVals() << "\n";	
	//PRINT2( statKeepers.size(), statNames->size() );
	//cerr << "got " << statKeepers.size() << " stats." << endl;
	BOOST_FOREACH( const StatKeeper<>& sk, statKeepers ) {
		string sname = sk.getName();
		cerr << sname << ".mean=" << sk.getMean() << " " << sname << ".std=" << sk.getStd() <<
			" " << sname << ".count=" << sk.getNumVals() << "\n";
	}
	cerr << "\n";

	free_cmatrix( list, nsam );
	free( trimmed_list );
	free( posit_gloc );
	free( posit );
	return EXIT_SUCCESS;
}

// * util impls

/* allocates space for gametes (character strings) */
char **
cmatrix(int nsam,int len)
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) ) {
	   perror("alloc error in cmatrix") ;
		 chk( m );
	}
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) ))) {
			 perror("alloc error in cmatric. 2");
			 chk( m[i] );
		}
	}
	return( m );
}

void free_cmatrix(char **cm, nchroms_t nsam ) {
	for ( int i = 0; i < nsam; i++ )
		 free( cm[ i ] );
	free( cm );
}

void
biggerlist(int nsam, unsigned nmax, char **list )
{
	int i;

	maxsites = nmax  ;
	for( i=0; i<nsam; i++){
		chk( ( list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ) );
	}
}                        


double
nucdiv( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(1.0 -p1)*nnm1 ;
	}
	return( pi ) ;
}

/*   thetah - pi   */
double
hfay( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
	}
	return( -pi ) ;
}

/* Fay's theta_H  */
double
thetah( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd /*, nnm1*/  ;

	pi = 0.0 ;

	nd = nsam;
	//nnm1 = nd/(nd-1.0) ;
	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list) ;
		pi += p1*p1 ; 
	}
	return( pi*2.0/( nd*(nd-1.0) )  ) ;
}


int
frequency( char allele,int site,int nsam,  char **list);

void get_freqs( int site_i, int site_j, int nsam, char **list,
								double *p_i, double *p_j, double *p_ij );
int
frequency( char allele,int site,int nsam,  char **list)
{
	int i, count=0;
	for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}

void get_freqs( nsnps_t site_A, nsnps_t site_B, nchroms_t nsam, char **list,
								freq_t *p_1, freq_t *q_1, freq_t *x_11 ) {

	nchroms_t n_11 = 0;
	for ( nchroms_t s = 0; s < nsam; s++ ) {
		bool A_der = ( list[s][site_A] == '1' );
		bool B_der = ( list[s][site_B] == '1' );
		if ( A_der && B_der ) n_11++;
	}
	nchroms_t n_A_1 = derCounts[ site_A ];
	nchroms_t n_B_1 = derCounts[ site_B ];

	*p_1 = double( n_A_1 ) / double( nsam );
  *q_1 = double( n_B_1 ) / double( nsam );
	*x_11 = double( n_11 ) / double( nsam ); 
}

bool compute_LD( int site_A, int site_B, nchroms_t nsam, char **list,
								 double *r2, double *Dprime ) {
	freq_t p_1, q_1, x_11;
	assert( r2 && Dprime && list );
	get_freqs( site_A, site_B, nsam, list, &p_1, &q_1, &x_11 );
	assert( 0.0 <= p_1 && p_1 <= 1.0 );
	assert( 0.0 <= q_1 && q_1 <= 1.0 );
	assert( x_11 <= p_1 && x_11 <= q_1 );

	freq_t p_2 = 1.0 - p_1;
	freq_t q_2 = 1.0 - q_1;

	
	if ( p_1 < 1e-10 || p_2 < 1e-10 || q_1 < 1e-10 || q_2 < 1e-10 ) {
		*Dprime = std::numeric_limits<cosi_double>::quiet_NaN();
		*r2 = std::numeric_limits<cosi_double>::quiet_NaN();
		return false;
	} else {

		double D = x_11 - p_1 * q_1;
		double D_max = ( D < 0 ) ?  std::min( p_1 * q_1, p_2 * q_2  )  : std::min( p_1 * q_2, p_2 * q_1 );
		
		*Dprime = std::abs( D / D_max );
		*r2 = ( D * D ) / ( p_1 * p_2 * q_1 * q_2 );
		//PRINT10( site_A, site_B, nsam, p_1, q_1, x_11, D, D_max, *r2, *Dprime );
#ifndef NDEBUG		
		double eps = 1e-6;
		assert( -eps <= *Dprime && *Dprime <= ( 1.0 + eps ) );
		assert( -eps <= *r2 && *r2 <= ( 1.0 + eps ) );
#endif		
		return true;
	}
}

void get_freqs( nsnps_t site_A, nsnps_t site_B, nchroms_t bsam, nchroms_t nsam, char **list,
								freq_t *p_1, freq_t *q_1, freq_t *x_11 ) {

	nchroms_t n_11 = 0, n_A_1 = 0, n_B_1 = 0;
	for ( nchroms_t s = bsam; s < bsam+nsam; s++ ) {
		bool A_der = ( list[s][site_A] == '1' );
		bool B_der = ( list[s][site_B] == '1' );
		if ( A_der ) ++n_A_1;
		if ( B_der ) ++n_B_1;
		if ( A_der && B_der ) n_11++;
	}

	*p_1 = double( n_A_1 ) / double( nsam );
  *q_1 = double( n_B_1 ) / double( nsam );
	*x_11 = double( n_11 ) / double( nsam ); 
}

bool compute_LD( int site_A, int site_B, nchroms_t bsam, nchroms_t nsam, char **list,
								 double *r2, double *Dprime ) {
	freq_t p_1, q_1, x_11;
	assert( r2 && Dprime && list );
	get_freqs( site_A, site_B, bsam, nsam, list, &p_1, &q_1, &x_11 );
	assert( 0.0 <= p_1 && p_1 <= 1.0 );
	assert( 0.0 <= q_1 && q_1 <= 1.0 );
	assert( x_11 <= p_1 && x_11 <= q_1 );

	freq_t p_2 = 1.0 - p_1;
	freq_t q_2 = 1.0 - q_1;

	
	if ( p_1 < 1e-10 || p_2 < 1e-10 || q_1 < 1e-10 || q_2 < 1e-10 ) {
		*Dprime = std::numeric_limits<cosi_double>::quiet_NaN();
		*r2 = std::numeric_limits<cosi_double>::quiet_NaN();
		return false;
	} else {

		double D = x_11 - p_1 * q_1;
		double D_max = ( D < 0 ) ?  std::min( p_1 * q_1, p_2 * q_2  )  : std::min( p_1 * q_2, p_2 * q_1 );
		
		*Dprime = std::abs( D / D_max );
		*r2 = ( D * D ) / ( p_1 * p_2 * q_1 * q_2 );
		//PRINT10( site_A, site_B, nsam, p_1, q_1, x_11, D, D_max, *r2, *Dprime );
#ifndef NDEBUG		
		double eps = 1e-6;
		assert( -eps <= *Dprime && *Dprime <= ( 1.0 + eps ) );
		assert( -eps <= *r2 && *r2 <= ( 1.0 + eps ) );
#endif		
		return true;
	}
}



int
segsub( int nsub, int segsites, char **list )
{
	int i, count = 0 , c1 ;
	int frequency( char, int, int, char**) ;

	for(i=0; i < segsites ; i++){
	  c1 = frequency('1',i,nsub, list);
	  if( ( c1 > 0 ) && ( c1 <nsub )  ) count++;
	}
	return( count ) ;
}

	
}  // namespace ms

// * top-level main
int main( int argc, char *argv[] ) {
	int exitCode = EXIT_FAILURE;
	try { 
		exitCode = cosi::sample_stats_main( argc, argv );
	} catch( const boost::exception& e ) {
		std::cerr << "cosi error: " << boost::diagnostic_information( e ) << std::endl;
	}
#ifdef COSI_DEV_PRINT	
	std::cerr << "exit code " << exitCode << std::endl;
#endif
	return exitCode;
}

#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic pop
#endif
#endif

