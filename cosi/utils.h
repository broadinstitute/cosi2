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
#include <limits>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/swap.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/and.hpp>
#include <boost/optional.hpp>
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
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/checked_delete.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/all.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/next_prior.hpp>
#include <cosi/utildefs.h>

#ifdef COSI_VALGRIND
#include "valgrind.h"
#include "memcheck.h"
#endif

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {

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

//#define COSI_LEAK_CHECK do { static int chk_count = 0; printf( "leak check %d at %s:%d\n", chk_count++, __FILE__, __LINE__ ); VALGRIND_DO_LEAK_CHECK; } while(0)

//
// cosi exceptions
//

struct cosi_error: virtual std::exception, virtual boost::exception { };
struct cosi_io_error: virtual cosi_error { };

typedef boost::error_info<struct tag_errno_code,int> errno_code;
typedef boost::error_info<struct tag_error_bad_line,std::string> error_bad_line;
typedef boost::error_info<struct tag_error_msg,std::string> error_msg;
typedef boost::error_info<struct tag_error_msg2,std::string> error_msg2;
typedef boost::error_info<struct tag_error_msg3,std::string> error_msg3;
typedef boost::error_info<struct tag_error_msg4,std::string> error_msg4;
typedef boost::error_info<struct tag_error_msg5,std::string> error_msg5;

namespace util {

const char *DateStr();
void dbg(const char *fmt, ...);

/// FuncProt: chk
/// Check that a pointer is non-NULL.  If it is non-NULL, return the pointer.
/// If it is NULL, cause a fatal error with the specified error message.
void *chk(const void *p, const char *fmt, ...);

template <typename T>
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


/// FuncProt: chk
/// Check that a pointeer is non-NULL.  If it is non-NULL, return the pointer.
/// If it is NULL, cause a fatal error with the specified error message.
template <typename T>
T *chk(T *p, const char *fmt, ...) {
  va_list ap;

  if ( p ) return p;

  printf( "Error: null ptr - %s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
  throw std::logic_error( "error" );
  return NULL;
}

template <typename T>
const T *chk(const T *p, const char *fmt, ...) {
  va_list ap;

  if ( p ) return p;

  printf( "Error: null ptr - %s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
  throw std::logic_error( "error" );
  return NULL;
}


#define Chk(p) chk(p, #p)

/// FuncProt: chkCond
/// Check that a condition is true; abort with an error if not.
/// Same as assert(), except the check is always done regardless of -DNDEBUG status.
void chkCond(int cond, const char *fmt = "error", ...);

/// FuncProt: fopenChk
/// Open a file, with error-checking.
FILE *fopenChk( const char *fname, const char *mode );

const int BITS_PER_BYTE = 8;  

#define COSI_CALLOC(T,n) ( (T *)cosi::chk(calloc( (n), sizeof( T ) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )
#define COSI_ALLOC(T) COSI_CALLOC(T,1)
#define COSI_CALLOC_VARSIZE(T,Tvar,n) ( (T *)cosi::chk(calloc( 1, sizeof( T ) + (n)*sizeof(Tvar) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )

#define COSI_MALLOC(T,n) ( (T *)cosi::chk(malloc( (n) * sizeof( T ) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )
#define COSI_MALLOC_VARSIZE(T,Tvar,n) ( (T *)cosi::chk(malloc( sizeof( T ) + (n)*sizeof(Tvar) ), "error allocating memory at %s:%d", __FILE__, __LINE__ )  )

#define COSI_REALLOC(p,T,n) do { p = ( (T *)cosi::chk(realloc( p, (n) * sizeof( T ) ), \
																											"error allocating memory at %s:%d", __FILE__, __LINE__ )  ); } while( 0 )


long GetMemUsage( );

int comparePointers( const void *pp1, const void *pp2 );
void dontFreePtr( void *);

#ifndef HAVE_TDESTROY
void tdestroy (void *root, void (*free_node)(void *nodep));
#endif

char *cosi_strdup( const char *s );
char * cosi_strtok_r(char *s, const char *delim, char **last);

template <typename T> inline void updateMin( T& curMin, T newVal ) { if ( newVal < curMin ) curMin = newVal; }
template <typename T> inline void updateMax( T& curMax, T newVal ) { if ( newVal > curMax ) curMax = newVal; }

#define cosi_fwrite(ptr, size, nmemb, stream) cosi_fwrite_helper(ptr, size, nmemb, stream, #ptr, __FILE__, __LINE__ )
#define cosi_fread(ptr, size, nmemb, stream) cosi_fread_helper(ptr, size, nmemb, stream, #ptr, __FILE__, __LINE__ )
void cosi_fwrite_helper(const void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line);
void cosi_fread_helper(void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line);

template <typename T>
inline void cosi_free(T& x) { free(x); x = NULL; }

template <typename T> inline
void cosi_fwrite_val_helper(const T& val, FILE *stream, const char *expr, const char *fname, int line) { cosi_fwrite_helper(&(val), sizeof(val), 1, stream, expr, fname, line ); }
template <typename T> inline
void cosi_fread_val_helper(T& val, FILE *stream, const char *fname, int line) { cosi_fread_helper(&(val), sizeof(val), 1, stream, fname, line ); }

#define cosi_fwrite_val(val, stream) cosi_fwrite_val_helper(val, stream, #val, __FILE__, __LINE__ )
#define cosi_fread_val(val, stream) cosi_fread_val_helper(val, stream, __FILE__, __LINE__ )

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
//#define For(x,v) for ( auto x = boost::begin(v); x != boost::end(v); x++ )

#define cosi_for_map( k, v, m ) for( BOOST_AUTO( it, m.begin() ); it != m.end(); ++ it ) { \
	BOOST_AUTO( const& k, it->first ); \
	BOOST_AUTO( const& v, it->second );

#define cosi_end_for }

#define ForEach BOOST_FOREACH

template <typename K, typename V>
V const& at( const std::map< K, V>& m, const K k ) {
	typename std::map< K, V>::const_iterator it = m.find( k );
	if ( it == m.end() ) BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "lookup in map failed" ) );
	else
		 return it->second;
}

template <typename T> struct IndexSorter {
	 const vector<T>& origVec;
	 IndexSorter( const vector<T>& origVec_ ): origVec( origVec_ ) { }
	 bool operator()( int i, int j ) { return origVec[i] < origVec[j]; }
};

template <class T>
class refcounted {
	 unsigned int nrefs;

protected:
	 refcounted(): nrefs(0) { }

public:
	 template <class U> friend void intrusive_ptr_add_ref( refcounted<U> * );
	 template <class U> friend void intrusive_ptr_release( refcounted<U> * );
};

template <class T>
inline void intrusive_ptr_add_ref( refcounted<T> *r ) { r->nrefs++; }

template <class T>
inline void intrusive_ptr_release( refcounted<T> *r ) {
  if ( !( --r->nrefs ) )
		 boost::checked_delete((T *)r);
}

std::string Date( );

// The following PRINT macros are intended primarily for debugging.

extern bool noDbgPrint;

#ifdef COSI_DEV_PRINT
#define PRCORE(X) #X " = " << (X)
#define PRINT(X) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << std::endl; } while(0)
#define PRINT2(X, Y) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << std::endl; } while(0)
#define PRINT3(X, Y, Z) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
																													<< PRCORE(Z) << std::endl; } while(0)
#define PRINT4(X, Y, Z, W) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
																														 << PRCORE(Z) << ", " << PRCORE(W) << std::endl; } while(0)
#define PRINT5(X, Y, Z, W, T) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
																																<< PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << std::endl; } while(0)
#define PRINT6(X, Y, Z, W, T, U) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) << ", " \
	<< PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "				\
																																	 << PRCORE(U) << std::endl; } while(0)
#define PRINT7(X, Y, Z, W, T, U, V) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " << PRCORE(X) << ", " << PRCORE(Y) \
  << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", " \
																																			<< PRCORE(U) << ", " << PRCORE(V) << std::endl; } while(0)
#define PRINT8(V1, V2, V3, V4, V5, V6, V7, V8) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
																																								 << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << std::endl; } while(0)
#define PRINT9(V1, V2, V3, V4, V5, V6, V7, V8, V9) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
																																										 << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << std::endl; } while(0)
#define PRINT10(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
																																													 << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << std::endl; } while(0)
#define PRINT11(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
																																																<< PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) << std::endl; } while(0)
#define PRINT12(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) \
																																																		 << ", " << PRCORE(V12) << std::endl; } while(0)
#define PRINT13(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13) do { if( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) \
																																																				 << ", " << PRCORE(V12) << ", " << PRCORE(V13) << std::endl; } while(0)
#define PRINT14(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14) do { if( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) \
																																																							<< ", " << PRCORE(V12) << ", " << PRCORE(V13) << ", " << PRCORE(V14) << std::endl; } while(0)
#define PRINT15(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) \
	<< ", " << PRCORE(V12) << ", " << PRCORE(V13) << ", " << PRCORE(V14) << ", " << PRCORE(V15) << std::endl
#define PRINT16(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16) do { if ( !cosi::util::noDbgPrint ) std::cerr << ::cosi::util::Date() << " at " << __FILE__ << ":" << __LINE__ << " " \
  << PRCORE(V1) << ", " << PRCORE(V2)																		\
  << ", " << PRCORE(V3) << ", " << PRCORE(V4) << ", " << PRCORE(V5) << ", "	\
  << PRCORE(V6) << ", " << PRCORE(V7) << ", " << PRCORE(V8) << ", " << PRCORE(V9) << ", " << PRCORE(V10) << ", " << PRCORE(V11) \
	<< ", " << PRCORE(V12) << ", " << PRCORE(V13) << ", " << PRCORE(V14) << ", " << PRCORE(V15) ", " << PRCORE(V16) << std::endl; } while(0)
#else  // #ifdef COSI_DEV_PRINT

#define PRINT(X)
#define PRINT2(X,Y)
#define PRINT3(X,Y,Z)
#define PRINT4(X,Y,Z,W)
#define PRINT5(X,Y,Z,W,T)
#define PRINT6(X,Y,Z,W,T,U)
#define PRINT7(X,Y,Z,W,T,U,V)
#define PRINT8(V1,V2,V3,V4,V5,V6,V7,V8)
#define PRINT9(V1,V2,V3,V4,V5,V6,V7,V8,V9)
#define PRINT10(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
#define PRINT11(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11)
#define PRINT12(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12)
#define PRINT13(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13)
#define PRINT14(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14)
#define PRINT15(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15)
#define PRINT16(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16)

#endif  // #ifdef COSI_DEV_PRINT


template <typename T>
struct Identity {
	 T operator()( const T& x ) { return x; }
};

// Function template: isize
// Return the size of a container as a signed integer.  Helps avoid "comparison between
// signed and unsigned integer" warnings.
template <class T> inline int isize( const T& container ) { return (int)container.size(); }
 
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

   ValT getMean() const { return static_cast< ValT >( numVals > 0 ? ( ToDouble( getSum() ) / cosi_double(numVals) ) : std::numeric_limits<cosi_double>::quiet_NaN() ); }

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
  
};  // class SumKeeper

//
// Class: StatKeeper
//
// Keeps a running sum and sum-of-squares without round-off errors,
// allowing accurate computation of the mean and stddev.
//
template <typename ValT = cosi_double, typename CountT = size_t>
class StatKeeper {
public:
	 StatKeeper() { clear(); }

   void add( ValT x ) { sum.add( x ); sumSq.add( static_cast< ValT >( ToDouble( x ) * ToDouble( x ) ) ); }

	 void clear() { sum.clear(); sumSq.clear(); }

	 ValT getSum() const { return sum.getSum(); }
	 CountT getNumVals() const { return sum.getNumVals(); }
	 CountT getNumNaNs() const { return sum.getNumNaNs(); }
	 CountT getNumInfs() const { return sum.getNumInfs(); }

	 // Method: add
	 // Add a range of values to this SumKeeper.
	 template <class ValRange>
	 void add( const ValRange& valRange ) { ForEach( ValT val, valRange ) add( val ); }

	 ValT getMean() const { return sum.getMean(); }
	 ValT getStd() const {
		 ValT meanSoFar = getMean();
		 return static_cast< ValT >( std::sqrt( ToDouble( sumSq.getMean() - ( static_cast< ValT >( ToDouble( meanSoFar ) * ToDouble( meanSoFar ) ) ) ) ) );
	 }
	 
private:
	 // Fields:
	 //
	 //   sum - sum of values passed to <add()>
	 //   sumSq - sum of squares of values passed to <add()>
	 SumKeeper<ValT,CountT> sum, sumSq;
};  // class StatKeeper


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
		 ensureCapacity( itemIdx+1 );
		 assert( itemIdx < static_cast<idx_t>( partialSums.size() ) );
	   for (; itemIdx < isize( partialSums ); itemIdx += ( itemIdx & -itemIdx ) ) {
				partialSums[ itemIdx ] += delta;
				assert( partialSums[ itemIdx ] >= T(-1e-10) );
		 }
   }

	 T readCumulative(idx_t idx) const {
		 T sum( 0.0 ) ;
		 while (idx > 0){
			 sum += partialSums[idx];
			 idx -= (idx & -idx);
		 }
		 return sum;
	 }

	 T read( idx_t idx ) const {
		 assert( idx > 0 );
		 const_cast< PartialSumTree * >(this)->ensureCapacity( idx+1 );

		 // copied from TopCoder http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=binaryIndexedTrees
		 T sum = partialSums[idx]; // sum will be decreased
		 if (idx > 0){ // special case
			 idx_t z = idx - (idx & -idx); // make z first
			 idx--; // idx is no important any more, so instead y, you can use idx
			 while (idx != z){ // at some iteration idx (y) will become z
				 sum -= partialSums[idx]; 
// substruct tree frequency which is between y and "the same path"
				 idx -= (idx & -idx);
			 }
		 }
		 assert( sum == readCumulative( idx ) - readCumulative( idx-1 ) );
		 return sum;
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
		 while ( partialSums.size() <= size+3 ) {
			 size_t oldSize = partialSums.size()-1;
			 if ( size+2 >= oldSize ) {
				 size_t newSize = oldSize * 2;
				 partialSums.resize( newSize+1, T( 0.0 ) );
				 partialSums[ newSize ] = partialSums[ oldSize ];
			 }
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
		 //T amtOrig = amt;
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

	 const std::vector<T>& getPartialSums() const { return partialSums; }
	 
private:
	 // Field: partialSums
	 // Partial sums of sub-ranges of elements of the underlying vector of weights.
	 // (The vector of weights itself is not stored, though can be recovered from <partialSums>).
	 // See the references in the <PartialSumTree> class comment for details.
	 std::vector<T> partialSums;
};  // class PartialSumTree

template <typename T1, typename T2>
inline
typename boost::enable_if< boost::mpl::and_< boost::is_convertible< T1, double >,
																						 boost::is_convertible< T2, double > >, double >::type
getFrac( T1 num, T2 denom ) { return static_cast<double>( num ) / static_cast<double>( denom ); }


struct unspecified;

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
	BOOST_AUTO_TPL( i1, boost::begin( r1 ) );
	BOOST_AUTO_TPL( i2, boost::begin( r2 ) );
	BOOST_AUTO_TPL( i3, boost::begin( r3 ) );
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
	BOOST_AUTO_TPL( i1, boost::begin( r1 ) );
	BOOST_AUTO_TPL( i2, boost::begin( r2 ) );
	BOOST_AUTO_TPL( i3, boost::begin( r3 ) );
	BOOST_AUTO_TPL( i4, boost::begin( r4 ) );
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
const V& map_get( const std::map<K,V>& m, const K& k ) {
	typename std::map<K,V>::const_iterator it = m.find( k );
	if ( it == m.end() ) throw std::logic_error( "map lookup error in map_get" );
	return it->second;
}

template <typename K, typename V>
const V& map_get( const std::map<K,V>& m, const K& k, const V& dflt ) {
	typename std::map<K,V>::const_iterator it = m.find( k );
	return ( it == m.end() ) ? dflt : it->second;
}


// Func: startsWith
// Test whether string 's' begins with string 't'.
inline bool_t startsWith( const string& s, const string& t ) {
	return s.size() >= t.size() && !s.compare( 0, t.size(), t );
}

template <class T>
ostream& operator<<( ostream& s, const vector<T>& v) {
	s << "[";
	bool first = true;
	for ( typename vector<T>::const_iterator it = v.begin(); it != v.end(); ++it ) {
		if ( !first ) s << ", ";
		using namespace std;
		s << *it;
		first = false;
	}
	//std::copy( v.begin(), v.end(), std::ostream_iterator<T>( s, ", " ) );
	s << "]";
	return s;
}

template <typename K, typename V>
ostream& operator<<( ostream& s, const std::map<K,V>& m) {
	s << "[";
	bool first = true;
	for ( typename std::map<K,V>::const_iterator it = m.begin(); it != m.end(); ++it ) {
		if ( !first ) s << ", ";
		using namespace std;
		s << "(" << it->first << "," << it->second << ")";
		first = false;
	}
	
	//std::copy( m.begin(), m.end(), std::ostream_iterator< typename std::map<K,V>::value_type >( s, ", " ) );
	s << "]";
	return s;
}

template <typename T1, typename T2>
ostream& operator<<( ostream& s, const std::pair<T1,T2>& p ) {
	s << "(" << p.first << "," << p.second << ")";
	return s;
}

template <typename T1, typename T2>
ostream& operator<<( ostream& s, const std::pair<const T1,T2>& p ) {
	s << "(" << p.first << "," << p.second << ")";
	return s;
}

inline bool_t isSpace( const string& s ) {
	ForEach( char c, s ) if ( !std::isspace( c ) ) return False;
	return True;
}

template <
	class PropObject,
	class PropType, 
	PropType PropObject::* Prop
	>
struct FieldRef
{
	 typedef PropObject argument_type;
	 typedef PropType result_type;

	 typedef PropType value_type;

	 PropType operator()( const PropObject& p ) const { return p.*Prop; }
	 PropType& operator()( PropObject& p ) const { return p.*Prop; }
};

namespace STLExt {

// Function: last_among_equal
// Returns an iterator pointing to the rightmost item in a range of equal items.
template <typename ForwardIterator>
ForwardIterator last_among_equal( ForwardIterator i, ForwardIterator end ) {
	while ( i != end ) {
		ForwardIterator n = boost::next( i );
		if ( n != end && *n == *i )
			 i = n;
		else
			 break;
	}
	return i;
}

// Function: first_among_equal
// Returns an iterator pointing to the leftmost item in a range of equal items.
template <typename BidiIterator>
BidiIterator first_among_equal( BidiIterator i, BidiIterator begin ) {
	while ( i != begin ) {
		BidiIterator p = boost::prior( i );
		if ( *p == *i )
			 i = p;
		else
			 break;
	}
	return i;
}


}

template <typename T> inline
T countPairs( T n ) { return n * ( n-1 ) / 2 ; }

#if 0
class EventCounter {
public:
	 void count() 
private:
	 unsigned long long nevents;
};
#endif

//
// Class: ValRange
//
// A range of values.
//
template <typename T>
class ValRange {
public:
	 ValRange() { }
	 explicit ValRange( string rangeDef ) { init( rangeDef ); }
	 ValRange( T min_, T max_ ): bounds( min_, max_ ) { }

	 void init( string rangeDef ) {
		 if ( !rangeDef.empty() ) {

			 typedef boost::tokenizer<boost::char_separator<char> > 
					tokenizer;
			 boost::char_separator<char> sep("-");
			 tokenizer tokens(rangeDef, sep);
			 vector<string> result( tokens.begin(), tokens.end() );

			 if ( result.size() == 1 ) {
				 bounds.first = bounds.second = boost::lexical_cast<T>( result.at(0) );
			 } else if ( result.size() == 2 && !result.at(0).empty() && !result.at(0).empty() ) {
				 try { 
					 bounds.first = boost::lexical_cast<T>( result.at(0) );
					 bounds.second = boost::lexical_cast<T>( result.at(1) );
				 } catch( boost::bad_lexical_cast const& ) {
					 BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "invalid range: " + rangeDef ) );					 
				 }
			 }
			 else
				 BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "invalid range: " + rangeDef ) );
		 }
		 else
				BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "invalid range: " + rangeDef ) );
		 if ( !( bounds.first <= bounds.second ) ) 
				BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "invalid range - min>max: " + rangeDef ) );
	 }

	 bool operator() ( T val ) const { return ( !bounds.first || val >= *bounds.first ) &&
				( !bounds.second || val <= *bounds.second ); }

	 boost::optional<T> getMin() const { return bounds.first; }
	 boost::optional<T> getMax() const { return bounds.second; }
	 void setMin( boost::optional<T> min_ ) { bounds.first = min_; }
	 void setMax( boost::optional<T> max_ ) { bounds.second = max_; }
	 
private:
	 // Field: bounds
	 // The pair of bounds of the range values
	 pair< boost::optional<T>, boost::optional<T> > bounds;

};  // class ValRange

template <typename T>
inline ValRange<T> make_val_range( T min_, T max_ ){ return ValRange<T>( min_, max_ ); }

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
	if ( valRange.getMin() ) os << valRange.getMin();
	if ( valRange.getMax() && !( valRange.getMin() == valRange.getMax() ) ) os << "-" << valRange.getMax();
	
	return os;
}

// ** Function: interpolate - linearly interpolate between two points.
template <typename TX, typename TY>
TY interpolate( TX x1, TY y1, TX x2, TY y2, TX x ) {
	assert( x1 <= x && x <= x2 );
	return y1 + ( y2 - y1 ) * ( ( x - x1 ) / ( x2 - x1 ) );
}

template <typename TXRange, typename TYRange>
typename boost::range_value<TYRange>::type
interp( const TXRange& xs, const TYRange& ys, typename boost::range_value<TXRange>::type x ) {
	assert( !boost::empty( xs ) );
	assert( boost::size( xs ) == boost::size( ys ) );
	typename boost::range_iterator<const TXRange>::type x_it = boost::lower_bound( xs, x );
	assert( x_it != boost::end( xs ) );
	if ( x_it == boost::begin( xs ) ) return *boost::begin( ys );
	typename boost::range_iterator<const TXRange>::type x_it_p = boost::prior( x_it );
	assert( *x_it_p <= x && x <= *x_it );
	typename boost::range_iterator<const TYRange>::type y_it = boost::next( boost::begin( ys ),
																																					std::distance( boost::begin( xs ), x_it ) );
	assert( y_it != boost::begin( ys ) );
	typename boost::range_iterator<const TYRange>::type y_it_p = boost::prior( y_it );
	return interpolate( *x_it_p, *y_it_p, *x_it, *y_it, x );
}

namespace tsv {


// ** Class TSVIdx - an index of a TSV file for quickly finding locations in it
class TSVIdx {
public:

   // *** Typedef index_t - the type of values in the index column.
	 typedef uint64_t index_t;
	 
	 TSVIdx( filename_t tsvFN, unsigned colNum, filename_t idxFN );

	 std::istream::streampos getStreamPos( index_t idxVal );

private:
	 boost::container::flat_map<index_t, std::istream::streampos> idx2streamPos;
};  // class TSVIdx

}  // namespace tsv

#undef ForEach

} // namespace util

}  // namespace cosi

#define cosi_using2(id1,id2)					 \
	using id1;													 \
	using id2;

#define cosi_using5(id1,id2,id3,id4,id5)  \
	using id1;													 \
	using id2;													 \
	using id3;													 \
	using id4;													 \
	using id5

#ifdef NDEBUG 
#define COSI_IF_DEBUG(x)
#define COSI_IF_NDEBUG(x) x
#define COSI_IFELSE_NDEBUG(x,y) x
#else
#define COSI_IF_DEBUG(x) x
#define COSI_IF_NDEBUG(x)
#define COSI_IFELSE_NDEBUG(x,y) y
#endif

// ** Type registration
BOOST_TYPEOF_REGISTER_TYPE(cosi::cosi_error)
BOOST_TYPEOF_REGISTER_TYPE(cosi::cosi_io_error)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::refcounted,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::SumKeeper,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::StatKeeper,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::PartialSumTree,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::ValRange,1)


#endif // #ifndef __INCLUDE_COSI_UTILS_H

