#ifndef COSI_INCLUDE_UTILDEFS_H
#define COSI_INCLUDE_UTILDEFS_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <string>
#include <boost/core/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/filesystem.hpp>

namespace cosi {


// Type: bool_t
// A Boolean value.  Normally same as bool; just, some compilers have buggy support for bool (or at least did in the past),
// so we want the option of defining this as an int.
typedef bool bool_t;

const bool_t True = true;
const bool_t False = false;

// Type: idx_t
// An index into an array.
typedef int idx_t;

typedef double frac_t;
typedef double freq_t;
typedef double prob_t;

typedef double rate_t;


#ifdef COSI_LONG_DOUBLE
typedef long double cosi_double;
#else
typedef double cosi_double;
#endif

inline bool equal_eps( cosi_double v1, cosi_double v2, cosi_double eps = 1e-14 ) { return ::fabs( v1 - v2 ) < eps; }


// Func: ToDouble
// Converts a value to the standard type 'cosi_double'.
// Helps write generic code that works without change for primitive types
// and for <typed values>.
inline cosi_double ToDouble( const cosi_double& x ) { return x; }
inline cosi_double ToDouble( int x ) { return x; }
inline cosi_double ToDouble( short x ) { return x; }
inline cosi_double ToDouble( unsigned int x ) { return x; }
inline cosi_double ToDouble( unsigned short x ) { return x; }
inline cosi_double ToDouble( long x ) { return x; }
inline cosi_double ToDouble( unsigned long x ) { return x; }
#ifdef COSI_LONG_DOUBLE
inline cosi_double ToDouble( double x ) { return x; }
#endif

inline int ToInt( const int& x ) { return x; }

typedef boost::filesystem::path filename_t;

 inline FILE *cosi_fopen( filename_t fname, const char *mode )  {
   std::string s = fname.string();
   return fopen( s.c_str(), mode );
 }

template <typename TFrom, typename TTo, typename Enable=void> struct cosi_converter: public boost::false_type { };

template <typename TTo, typename TFrom>
typename boost::enable_if< cosi_converter<TFrom,TTo>, TTo>::type cosi_cvt( TFrom x ) {
	return cosi_converter<TFrom,TTo>()( x );
}

template <typename TFrom, typename TTo>					
struct cosi_converter<TFrom,TTo,
											typename boost::enable_if< boost::is_convertible<TFrom,TTo> >::type >: public boost::true_type {
	 TTo operator()( TFrom const& x ) const { return static_cast<TTo>( x ); }
};

// use boost TTI to test for existence of this expression?  check for limitations

}  // namespace cosi

#endif // #ifndef COSI_INCLUDE_UTILDEFS_H
