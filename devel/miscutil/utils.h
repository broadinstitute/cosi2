/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright (2005) by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/

//
// Header: utils.h
//
// Miscelaneous utilities, not directly related to simulation or genetics.
//

#ifndef __INCLUDE_COSI_UTILS_H
#define __INCLUDE_COSI_UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cctype>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <regex>
#include <arpa/inet.h>
#include <boost/swap.hpp>
#include <boost/concept/assert.hpp>
#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/numeric.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/next_prior.hpp>
#include <boost/foreach.hpp>
//#include <cosi/defs.h>

namespace miscutil {

using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;

using std::map;
using std::vector;
using std::make_pair;
using std::pair;
using std::string;

#ifdef COSI_LONG_DOUBLE
typedef long double cosi_double;
#else
typedef double cosi_double;
#endif


/* Type: bool_t */
/* A Boolean value.  */
typedef bool bool_t;

const bool_t True = true;
const bool_t False = false;

typedef int idx_t;
typedef double frac_t;


#ifdef COSI_VALGRIND
#include "valgrind.h"
#include "memcheck.h"
#else
#define VALGRIND_PRINTF_BACKTRACE(...) printf(__VA_ARGS__)
#define VALGRIND_DO_LEAK_CHECK  
#endif

#define COSI_LEAK_CHECK do { static int chk_count = 0; printf( "leak check %d at %s:%d\n", chk_count++, __FILE__, __LINE__ ); VALGRIND_DO_LEAK_CHECK; } while(0)  

const char *DateStr();
void dbg(const char *fmt, ...);

/// FuncProt: chk
/// Check that a pointer is non-NULL.  If it is non-NULL, return the pointer.
/// If it is NULL, cause a fatal error with the specified error message.
void *chk(const void *p, const char *fmt, ...);

template <class T>
boost::shared_ptr<T> chk( boost::shared_ptr<T> p, const char *fmt, ... ) {
  va_list ap;

  if ( p.get() ) return p;

  printf( "Error: null ptr - %s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
  throw std::logic_error( "error" );
	return p;
}


#define Chk(p) chk(p, #p)

/// FuncProt: chkCond
/// Check that a condition is true; abort with an error if not.
/// Same as assert(), except the check is always done regardless of -DNDEBUG status.
void chkCond(int cond, const char *fmt = "error", ...);

/// FuncProt: fopenChk
/// Open a file, with error-checking.
FILE *fopenChk( const char *fname, const char *mode );

#define BITS_PER_BYTE 8  

#define COSI_CALLOC(T,n) ( (T *)miscutil::chk(calloc( (n), sizeof( T ) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )
#define COSI_ALLOC(T) COSI_CALLOC(T,1)
#define COSI_CALLOC_VARSIZE(T,Tvar,n) ( (T *)miscutil::chk(calloc( 1, sizeof( T ) + (n)*sizeof(Tvar) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )

#define COSI_MALLOC(T,n) ( (T *)miscutil::chk(malloc( (n) * sizeof( T ) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )
#define COSI_MALLOC_VARSIZE(T,Tvar,n) ( (T *)miscutil::chk(malloc( sizeof( T ) + (n)*sizeof(Tvar) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )

#define COSI_REALLOC(p,T,n) do { p = ( (T *)miscutil::chk(realloc( p, (n) * sizeof( T ) ), \
																											"error allocating memory at %s:%d", __FILE__, __LINE__ )  ); } while( 0 )

#define BIN_WRITE_DOUBLE(X,F) do { fwrite(&(X), sizeof(X), 1, F); } while(0)
#define BIN_READ_DOUBLE(X,F) do { fread(&(X), sizeof(X), 1, F); } while(0)
  
#define BIN_WRITE(X,F) do { uint32_t v = htonl((uint32_t)X); fwrite( &v, sizeof(v), 1, F); } while(0)
#define BIN_READ(X,F) do { uint32_t v; fread( &v, sizeof(v), 1, F); X = (int)ntohl( v ); } while(0)

#define BIN_WRITE_2(X1,X2,F) do { BIN_WRITE(X1,F); BIN_WRITE(X2,F); } while(0) 
#define BIN_READ_2(X1,X2,F) do { BIN_READ(X1,F); BIN_READ(X2,F); } while(0)

#define BIN_WRITE_3(X1,X2,X3,F) do{ BIN_WRITE_2(X1,X2,F); BIN_WRITE(X3,F); } while(0)
#define BIN_READ_3(X1,X2,X3,F) do{ BIN_READ_2(X1,X2,F); BIN_READ(X3,F); } while(0)

#define BIN_WRITE_4(X1,X2,X3,X4,F) do{ BIN_WRITE_3(X1,X2,X3,F); BIN_WRITE(X4,F); } while(0)
#define BIN_READ_4(X1,X2,X3,X4,F) do{ BIN_READ_3(X1,X2,X3,F); BIN_READ(X4,F); } while(0)

#define BIN_WRITE_5(X1,X2,X3,X4,X5,F) do{ BIN_WRITE_4(X1,X2,X3,X4,F); BIN_WRITE(X5,F); } while(0)
#define BIN_READ_5(X1,X2,X3,X4,X5,F) do{ BIN_READ_4(X1,X2,X3,X4,F); BIN_READ(X5,F); } while(0)

long GetMemUsage( );

#define is_odd(x) (((x) & 1) > 0)   
#define is_even(x) (((x) & 1) == 0)
  
char *get_ith_token( const char *str, int i );
int get_ith_token_int( const char *str, int i );
double get_ith_token_double( const char *str, int i );
int index_of_token( const char *str, const char *token );

int comparePointers( const void *pp1, const void *pp2 );
void dontFreePtr( void *);

#ifndef HAVE_TDESTROY
void tdestroy (void *root, void (*free_node)(void *nodep));
#endif

int binarySearch( const cosi_double *a, int from, int to, cosi_double key );
  
bool_t random_bit(void);  

#define SWAP(T,x,y) do{ T tmp; tmp = x; x = y; y = tmp; } while(0)

char *cosi_strdup( const char *s );
char * cosi_strtok_r(char *s, const char *delim, char **last);

#define cosi_fwrite(ptr, size, nmemb, stream) cosi_fwrite_helper(ptr, size, nmemb, stream, #ptr, __FILE__, __LINE__ )
#define cosi_fread(ptr, size, nmemb, stream) cosi_fread_helper(ptr, size, nmemb, stream, #ptr, __FILE__, __LINE__ )
void cosi_fwrite_helper(const void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line);
void cosi_fread_helper(void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line);

#define cosi_free(x) do { free(x); x = NULL; } while(0)

#define cosi_fwrite_val(val, stream) cosi_fwrite_helper(&(val), sizeof(val), 1, stream, #val, __FILE__, __LINE__ )
#define cosi_fread_val(val, stream) cosi_fread_helper(&(val), sizeof(val), 1, stream, #val, __FILE__, __LINE__ )

int *make_random_permutation( size_t n );

void print_trace (void);

#ifdef NDEBUG
#define cosi_assert(cond)
#else  
#define cosi_assert(cond) do{ if( !(cond)) { print_trace(); assert(cond); } } while(0)
#endif

int ptr_id( const void *p );

typedef long obj_id_t;

obj_id_t get_obj_id(void);

struct timespec;  
  
unsigned long timespec_since( const struct timespec *start);

#define ForVec(T,x,v) for( vector< T >::const_iterator x = (v).begin(); x != (v).end(); x++ )
#define ForVecMut(T,x,v) for( vector< T >::iterator x = (v).begin(); x != (v).end(); x++ )
#define ForMap( KeyT, ValT, x, m ) for( map< KeyT, ValT >::const_iterator x = (m).begin(); x != (m).end(); x++ )

// Macro: For
// Loop over elements of a range.  Less generic but possibly more efficient than
// BOOST_FOREACH.
// Params:
//    x - name of the loop variable, which will be assigned an _iterator_ through the range.
//    v - the range.
#define For(x,v) for ( auto x = boost::begin(v); x != boost::end(v); x++ )

#define ForEach BOOST_FOREACH

template <typename T> struct IndexSorter {
	 const vector<T>& origVec;
	 IndexSorter( const vector<T>& origVec_ ): origVec( origVec_ ) { }
	 bool operator()( int i, int j ) { return origVec[i] < origVec[j]; }
};

cosi_double interpolate( const vector<cosi_double>& f, cosi_double i );
cosi_double findWhere( const vector<cosi_double>& f, cosi_double c );
cosi_double integrate( const vector<cosi_double>& f, const vector<cosi_double>& x, cosi_double from, cosi_double to );

template <class T>
struct refcounted {
public:
	 refcounted(): nrefs(0) { }

	 unsigned int nrefs;


	 //  friend template <> void intrusive_ptr_add_ref<T>( refcounted<T> * );
	 //  friend template <> void intrusive_ptr_release<T>( refcounted<T> * );
};

template <class T>
inline void intrusive_ptr_add_ref( refcounted<T> *r ) { r->nrefs++; }

template <class T>
inline void intrusive_ptr_release( refcounted<T> *r ) {
  if ( !( --r->nrefs ) )
		 delete ((T *)r);
}

std::string Date( );

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
#define PRINT8(V1, V2, V3, V4, V5, V6, V7, V8) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << endl
#define PRINT9(V1, V2, V3, V4, V5, V6, V7, V8, V9) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << endl
#define PRINT10(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << endl
#define PRINT11(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11) cerr << Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) << endl


template <typename T>
struct Identity {
	 T operator()( const T& x ) { return x; }
};

//
// Class: SumKeeper
//
// Keeps a running sum without round-off errors.  Also keeps count of
// NaN and infinite values.
//
// Adapted from: http://code.activestate.com/recipes/393090/
//
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
	 void add( const ValRange& v ) { For( i, v ) add( *i ); }

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
  
};


//
// Class: PartialSumTree
//
// For an array of objects each of which has an associated weight,
// efficiently support the following operations:
//
//    - changing an object's weight
//    - getting the total sum of all weights
//    - finding the longest prefix of the array such that the sum of
//      weights for objects in the prefix does not
//      exceed a specified fraction of the total sum of all weights
//
// Objects are addressed by their index in the array.  The weight
// of the object with index 0 is always 0. The array is
// expanded automatically when higher-indexed objects are accessed.
// The weights of all objects are initially zero.
//
// The implementation uses Fenwick trees:
// http://en.wikipedia.org/wiki/Fenwick_tree
// http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=binaryIndexedTrees
//
// Template params:
//
//    T - the type of value (weight) associated with each object.
//
template <typename T> 
class PartialSumTree {
public:
	 typedef T value_type;
	 
	 PartialSumTree(): partialSums( 9, T(0.0) )  { partialSums.reserve( 16385 ); }

	 // Method: add
	 // To the weight of the object at index 'itemIdx', add the value 'delta'
	 // (which can be positive or negative).  Note that the weight of the
	 // object with index 0 is always 0 and cannot be changed; thus, object
	 // indexing effectively starts at 1.
   void add( idx_t itemIdx, T delta ) {
		 assert( itemIdx > 0 );
		 ensureCapacity( itemIdx );
	   for (; itemIdx < isize( partialSums ); itemIdx += ( itemIdx & -itemIdx ) )
				partialSums[ itemIdx ] += delta;
   }

	 // Method: getSotalSum
	 // Returns the sum of weights of all objects.  Recall that the
	 // weights of all objects are initially zero until changed
	 // by the <add()> method, so objects whose weights were never changed
	 // count as weight zero.
	 T getTotalSum() const { return partialSums.back(); }

	 // Method: ensureCapacity
	 // Prepares this PartialSumTree for storing weights of objects with indices
	 // of at least 'size' by pre-allocating internal storage.  Like <std::vector::ensureCapacity>,
	 // affects only performance and not functionality.
	 void ensureCapacity( size_t size ) {
		 size_t oldSize = partialSums.size()-1;
		 if ( size+2 >= oldSize ) {
			 size_t newSize = oldSize * 2;
			 partialSums.resize( newSize+1, T( 0.0 ) );
			 partialSums[ newSize ] = partialSums[ oldSize ];
		 }
	 }

	 //
	 // Method: findCumulativeFraction
	 //
	 // Given a fraction of the sum of all weights, find the largest index i
	 // such that the sum of weights through i divided by <getTotalSum()>
	 // is less than that fraction.  Also returns the residue, such that
	 // adding the residue to the sum of weights through i and dividing by
	 // <getTotalSum()> equals 'cumulativeFraction'.
	 //
	 // Input params:
	 //
	 //     cumulativeFraction - specifies a fraction of <getTotalSum()>
	 //
	 // Returns:
	 //
	 //     the highest index i such that the sum of weights of elements 1..i
	 //     is less than cumulativeFraction*getTotalSum().   Also, in 'residue'
	 //     stores the value r such that the sum of weights of elements 1..i
	 //     plus r equals cumulativeFraction*getTotalSum().
	 //
	 idx_t findCumulativeFraction( frac_t cumulativeFraction, T *residue ) const {
	
		 idx_t searchRegionBase = 0;
		 int searchRegionHalfSize = partialSums.size() >> 1;
		 T amt = cumulativeFraction * getTotalSum();
		 T amtOrig = amt;
		 while( searchRegionHalfSize > 0 ) {
			 idx_t searchRegionMid = searchRegionBase + searchRegionHalfSize;
			 T amtInFirstHalf = partialSums[ searchRegionMid ];
			 if ( amt > amtInFirstHalf ) {
				 amt -= amtInFirstHalf;
				 searchRegionBase = searchRegionMid;
			 }
			 searchRegionHalfSize >>= 1;
		 }
		 idx_t itemIdx = searchRegionBase+1;
		 *residue = amt;
		 return itemIdx;
	 }  // findCumulativeFraction
	 
private:
	 // Field: partialSums
	 // Partial sums of sub-ranges of elements of the underlying vector of weights.
	 // (The vector of weights itself is not stored, though can be recovered from <partialSums>).
	 // See the references in the <PartialSumTree> class comment for details.
	 std::vector<T> partialSums;
};  // class PartialSumTree


//
// Class: InterpFn
//
// An interpolated function: stores a vector of (x,f(x)) pairs and
// can compute the values of f at new points by simple linear interpolation.
//
template <typename DomT = cosi_double, typename RanT = cosi_double>
class InterpFn {
public:
	 typedef DomT domain_type;
	 typedef RanT range_type;

	 // Method: clear
	 // Remove all (x,f(x)) pairs from this function.
	 void clear() { x2y.clear(); }

	 // Method: addPt
	 // Store an (x,f(x)) pair specifying the value of the function at a given point 'x'.
	 // Points must be added in order of increasing 'x'.
	 void addPt( DomT x, RanT y ) {
		 assert( x2y.empty() || x > boost::prior( x2y.end() )->first );
		 x2y.insert( x2y.end(), make_pair( x, y ) );
	 }

	 // Method: eval
	 // Evaluate the function at the specified point 'x', interpolating as necessary
	 // between the neighboring values.
	 RanT eval( DomT x ) const {
		 auto it = x2y.lower_bound( x );
		 assert( it != x2y.end() );
		 if ( it->first == x ) return it->second;
		 assert( it != x2y.begin() );
		 auto it_p = boost::prior( it );
		 cosi_double frac = ( ToDouble( x ) - ToDouble( it_p->first ) ) / ( ToDouble( it->first ) - ToDouble( it_p->first ) );
		 cosi_double result = ( ToDouble( it_p->second ) * (1.0-frac) + ToDouble( it->second ) * frac );
		 return RanT( result );
	 }

	 // Operator: function call
	 // Evaluate the function at the specified point 'x', interpolating as necessary
	 // between the neighboring values.
	 RanT operator()( DomT x ) const { return eval( x ); }

	 
#if 0
	 struct Pt {
			DomT x;
			RanT y;
			typename map<DomT,RanT>::const_iterator lwrBnd, uprBnd;
	 };

	 Pt evalPt( DomT x ) const {
		 auto it = x2y.lower_bound( x );
		 assert( it != x2y.end() );
		 Pt p;
		 p.x = x;
		 p.lwrBnd = it;
		 if ( it->first == x ) {
			 p.uprBnd = it;
			 p.y = it->second;
		 } else {
			 assert( it != x2y.begin() );
			 auto it_p = boost::prior( it );
			 p.uprBnd = it_p;
			 cosi_double frac = ( x - it_p->first ) / ( it->first - it_p->first );
			 cosi_double result = ( it_p->second * (1.0-frac) + it->second * frac );
			 p.y = result;
		 }
		 
		 return p;
	 }
	 
	 Pt evalPt( Pt from, Pt to, DomT x ) const {
		 auto it = lower_bound( from.uprBnd, to.lwrBnd, x );

		 assert( it != x2y.end() );
		 Pt p;
		 p.x = x;
		 p.lwrBnd = it;
		 if ( it->first == x ) {
			 p.uprBnd = it;
			 p.y = it->second;
		 } else {
			 assert( it != x2y.begin() );
			 auto it_p = boost::prior( it );
			 p.uprBnd = it_p;
			 cosi_double frac = ( x - it_p->first ) / ( it->first - it_p->first );
			 cosi_double result = ( it_p->second * (1.0-frac) + it->second * frac );
			 p.y = result;
		 }
		 
		 return p;
	 }
#endif

private:
	 // Field: x2y
	 // Map from point x to the value f(x), for the points at which the function value
	 // is explicitly specified.
	 map<DomT,RanT> x2y;
	 
};  // class InterpFn


// Function template: isize
// Return the size of a container as a signed integer.  Helps avoid "comparison between
// signed and unsigned integer" warnings.
template <class T> inline int isize( const T& container ) { return (int)container.size(); }

//
// Class: InterpBiFun
//
// An bijective interpolated function, which can map values both from
// domain to range and from range to domain.
//
template <typename DomT = cosi_double, typename RanT = cosi_double >
class InterpBiFun: public InterpFn<DomT, RanT> {
public:
	 typedef InterpFn<DomT,RanT> PARENT;
	 
	 void clear() { PARENT::clear(); inverse.clse(); }
	 void addPt( DomT x, RanT y ) { PARENT::addPt( x, y ); inverse.addPt( y, x ); }
	 DomT evalInv( RanT y ) const { return inverse( y ); }
	 
private:
	 InterpFn<RanT,DomT> inverse;
};

template <class T, typename FT, FT T::*fieldAddr>
struct CmpByField {
	 bool operator()( const FT& f, const T& x ) const { return f < x.*fieldAddr; }
	 bool operator()( const T& x, const FT& f ) const { return x.*fieldAddr < f; }
};

template <class Range1, class Range2, class Range3>
void save( const string& col1, const Range1& r1, const string& col2, const Range2& r2,
					 const string& col3, const Range3& r3,
					 const string& fname,
					 int precision ) {
	ofstream out( fname.c_str() );
	out.setf(ios::fixed,ios::floatfield);
	out.precision( precision );
	out << col1 << "\t" << col2 << "\t" << col3 << "\n";
	auto i1 = boost::begin( r1 );
	auto i2 = boost::begin( r2 );
	auto i3 = boost::begin( r3 );
	while( i1 != boost::end( r1 ) && i2 != boost::end( r2 ) && i3 != boost::end( r3 ) ) {
		out << *i1 << "\t" << *i2 << "\n";
		i1 = boost::next( i1 );
		i2 = boost::next( i2 );
		i3 = boost::next( i3 );
	}
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

inline string ToString( const string& s ) { return s; }

template <typename T1, typename T2>
string tabjoin( const T1& v1, const T2& v2 ) { return ToString( v1 ) + "\t" + ToString( v2 ); }
template <typename T1, typename T2, typename T3>
string tabjoin( const T1& v1, const T2& v2, const T3& v3 ) { return tabjoin( v1, v2 ) + "\t" + v3; }
template <typename T1, typename T2, typename T3, typename T4>
string tabjoin( const T1& v1, const T2& v2, const T3& v3, const T4& v4 ) { return tabjoin( v1, v2, v3 ) + "\t" + v4; }
template <typename T1, typename T2, typename T3, typename T4, typename T5>
string tabjoin( const T1& v1, const T2& v2, const T3& v3, const T4& v4, const T5& v5 ) { return tabjoin( v1, v2, v3, v4 ) + "\t" + v5; }


template <class Range1, class Range2, class Range3, class Range4>
void save( const string& col1, const Range1& r1, const string& col2, const Range2& r2,
					 const string& col3, const Range3& r3,
					 const string& col4, const Range4& r4,
					 const string& fname,
					 int precision ) {
	ofstream out( fname.c_str() );
	out.setf(ios::fixed,ios::floatfield);
	out.precision( precision );
	out << col1 << "\t" << col2 << "\t" << col3 << "\t" << col4 << "\n";
	auto i1 = boost::begin( r1 );
	auto i2 = boost::begin( r2 );
	auto i3 = boost::begin( r3 );
	auto i4 = boost::begin( r4 );
	while( i1 != boost::end( r1 ) && i2 != boost::end( r2 ) && i3 != boost::end( r3 ) ) {
		out << *i1 << "\t" << *i2 << "\t" << *i3 << "\t" << *i4 << "\n";
		i1 = boost::next( i1 );
		i2 = boost::next( i2 );
		i3 = boost::next( i3 );
		i4 = boost::next( i4 );
	}
}

// Function template: STLContains
// Test whether an STL container contains a given element.
template <class T, typename V>
bool_t STLContains( const T& container, const V& val ) { return container.find( val ) != container.end(); }

// Function template: map_get
// Return the value from a map given the key; *assumes the map has a mapping for this key*.
// If key is not in the map, result is undefined!
// Unlike accessing the map via array notation, this function treats the map as read-only.
template <typename K, typename V>
const V& map_get( const std::map<K,V>& m, const K& k ) { return m.find( k )->second; }

// Func: startsWith
// Test whether string 's' begins with string 't'.
inline bool_t startsWith( const string& s, const string& t ) {
	return s.size() >= t.size() && !s.compare( 0, t.size(), t );
}

template <class T>
ostream& operator<<( ostream& s, const vector<T>& v) {
	s << "[";
	std::copy( v.begin(), v.end(), std::ostream_iterator<T>( s, ", " ) );
	s << "]";
	return s;
}

inline bool_t isSpace( const string& s ) {
	For( it, s ) if ( !std::isspace( *it ) ) return False;
	return True;
}


}  // namespace cosi  

#endif
// #ifndef __INCLUDE_COSI_UTILS_H

