//
// Header: typedval.h
//
// A simple system for dimensional analysis to catch programming errors at compile time.
// Helps you define "logical types" to use throughout the program indicating a higher-level
// purpose of each value (as opposed to simply using physical types like 'double' everywhere).
// See <defs.h> for many usage examples.
//
// The base template for constructing logical types is <TypedVal>.  In practice most logical types
// will use one of two templates derived from TypedVal: <TypedValRel> and <TypedValAbs>.
// Macros to help define logical types with less boilerplate include
// <COSI_DEFINE_TYPEDVAL> and <COSI_DEFINE_TYPEDVAL_ABSREL>.
//
// (Design decision: decided not to use Boost.Unit at this time, even though it is much more
// full-featured, because its error messages were hard to understand without deep understanding
// of template metaprogramming. Also, the formulas used in cosi are sufficiently simple and
// give rise to sufficiently simple untis, that the relevant operations can be defined manually
// without too much trouble.)
//

#ifndef __INCLUDE_COSI_TYPEDVAL_H
#define __INCLUDE_COSI_TYPEDVAL_H

#include <limits>
#include <cstdlib>
#include <cassert>
#include <boost/operators.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <cosi/general/utildefs.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {

#ifdef COSI_DEV_TYPEDVAL_DISABLE

#define COSI_DEFINE_TYPEDVAL(Type) typedef cosi_double Type
typedef cosi_double factor_t;

const factor_t ZERO_FACTOR(0.0);
const cosi_double ZERO(0.0);

#define COSI_DEFINE_TYPEDVAL_MULT(TypeA, TypeB, TypeR)
#define COSI_DEFINE_TYPEDVAL_REL(Type) typedef cosi_double Type
#define COSI_DEFINE_TYPEDVAL_ABS(Type) typedef cosi_double Type

#define COSI_DEFINE_TYPEDVAL_ABSREL(AbsT,RelT,ValT) typedef ValT AbsT; typedef ValT RelT

#define COSI_DEFINE_TYPEDVAL_ID(Type) typedef int Type

#else  // if COSI_DEV_TYPEDVAL_DISABLE not defined

// Struct: zero_t
// Represents the value zero in a way understandable at compile time.
// For example if you want to be able to initialize an object from a zero
// but not from an arbitrary floating-point value, you could use this.
struct zero_t {};

// Const: ZERO
// Represents the value zero in a way understandable at compile time.
// For example if you want to be able to initialize an object from a zero
// but not from an arbitrary floating-point value, you could use this.
//const zero_t ZERO;

inline cosi_double ToDouble( const zero_t& ) { return 0.0; }

// Struct: nan_t
// Represents the NaN value in a way understandable at compile time.
struct nan_t {};

// Const: NAN_VAL
// Represents the value NaN in a way understandable at compile time
const nan_t NAN_VAL();

//
// Class template: TypedVal
//
// A typed value; its C++ type represents its units.
// So, for example, <genid> (representing a generation) and <ploc_t>
// (representing a physical location) are both values of C++ type 'cosi_double',
// but representing them as TypedVal instantiations of different types
// lets the compiler catch nonsensical operations such as passing a location
// to a function that expects a generation.
//
// Template params:
//
//     T - the concrete C++ type, passed here using the Curiously Recurring Template Pattern
//         ( see http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern ).
//     ValT - the 'physical' (C++) value type.  Since it is almost always 'cosi_double' in cosi,
//         the parameter defaults to that.  The TypedVal must be explicitly converted
//         from and to the underlying type; there deliberately are no implicit conversions.
template <class T, typename ValT = cosi_double >
struct TypedVal  {
	 // Field: val
	 // The actual value.
	 ValT val;

	 TypedVal(): val(0.0) { }
	 explicit TypedVal( ValT val_ ): val( val_ ) { }
	 //operator ValT() const { return val; }
	 TypedVal( const zero_t& ): val( 0 ) {}
	 T& operator=( const zero_t& ) { val = 0; return (T&)*this; }
	 TypedVal( const nan_t& ): val( std::numeric_limits<ValT>::quiet_NaN() ) {}
	 T& operator=( const nan_t& ) { val = std::numeric_limits<ValT>::quiet_NaN(); return (T&)*this; }

private:	 
	 friend class boost::serialization::access;
	 template <class Archive> void serialize( Archive& ar, const unsigned int /* version */ ) {
		 ar & BOOST_SERIALIZATION_NVP( val );
	 } 
};  // struct TypedVal

template <class T, typename V> inline bool operator==( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val == f2.val; }
template <class T, typename V> inline bool operator!=( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val != f2.val; }
template <class T, typename V> inline bool operator<( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val < f2.val; }
template <class T, typename V> inline bool operator>( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val > f2.val; }
template <class T, typename V> inline bool operator<=( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val <= f2.val; }
template <class T, typename V> inline bool operator>=( const TypedVal<T,V>& f1, const TypedVal<T,V>& f2 ) { return f1.val >= f2.val; }

template <class T, typename V> inline bool operator==( const TypedVal<T,V>& f, const zero_t& ) { return f.val == 0; }
template <class T, typename V> inline bool operator!=( const TypedVal<T,V>& f, const zero_t& ) { return f.val != 0; }
template <class T, typename V> inline bool operator<( const TypedVal<T,V>& f, const zero_t& ) { return f.val < 0; }
template <class T, typename V> inline bool operator<=( const TypedVal<T,V>& f, const zero_t& ) { return f.val <= 0; }
template <class T, typename V> inline bool operator>( const TypedVal<T,V>& f, const zero_t& ) { return f.val > 0; }
template <class T, typename V> inline bool operator>=( const TypedVal<T,V>& f, const zero_t& ) { return f.val >= 0; }

template <class T, typename V> inline bool operator==( const zero_t&, const TypedVal<T,V>& f ) { return f.val == 0; }
template <class T, typename V> inline bool operator!=( const zero_t&, const TypedVal<T,V>& f ) { return f.val != 0; }
template <class T, typename V> inline bool operator<( const zero_t&, const TypedVal<T,V>& f ) { return f.val < 0; }
template <class T, typename V> inline bool operator<=( const zero_t&, const TypedVal<T,V>& f ) { return f.val <= 0; }
template <class T, typename V> inline bool operator>( const zero_t&, const TypedVal<T,V>& f ) { return f.val > 0; }
template <class T, typename V> inline bool operator>=( const zero_t&, const TypedVal<T,V>& f ) { return f.val >= 0; }

template <class T> inline cosi_double ToDouble( const TypedVal<T,cosi_double>& factor ) { return factor.val; }
template <class T, typename V> inline std::ostream& operator<<( std::ostream& s, const TypedVal<T,V>& f ) { s << f.val; return s; }
template <class T, typename V> inline std::istream& operator>>( std::istream& s, TypedVal<T,V>& f ) { s >> f.val; return s; }

template <class T> inline const TypedVal<T,cosi_double> cosi_fabs( const TypedVal<T,cosi_double>& val ) { return TypedVal<T,cosi_double>( ::fabs( ToDouble( val ) ) ); }
template <class T> inline TypedVal<T,cosi_double> cosi_fabs( TypedVal<T,cosi_double>& val ) { return TypedVal<T,cosi_double>( ::fabs( ToDouble( val ) ) ); }

inline cosi_double cosi_fabs( cosi_double x ) { return ::fabs( x ); }

// Func: equal_eps
// Test whether two values are practically equal, with a specified tolerance.
template <class T>
inline bool equal_eps(const TypedVal<T,cosi_double>& val1, const TypedVal<T,cosi_double>& val2, cosi_double eps = 1e-14) {
	return ::fabs( ToDouble( val1 ) - ToDouble( val2 ) ) < eps;
}

template <class T>
inline bool is_positive( const TypedVal<T,cosi_double>& tval ) { return ToDouble( tval ) > 0.0; }
template <class T>
inline bool is_zero( const TypedVal<T,cosi_double>& tval ) { return ToDouble( tval ) == 0.0; }

template <class T>
inline bool is_null( const TypedVal<T,cosi_double>& tval ) { return (boost::math::isnan)( ToDouble( tval ) ); }

//
// Macro: COSI_DEFINE_TYPEDVAL
//
// Define a <logical type>, for values of the physical type 'cosi_double'.
//
#define COSI_DEFINE_TYPEDVAL(Type)													\
	struct Type: public TypedVal<Type> {											\
		Type() {}																								\
		explicit Type( cosi_double val_ ): TypedVal<Type>( val_ ) {}	\
  }

//
// Logical type: factor_t
//
// A dimensionless value.  Writing factor_t(NUMBER) instead of just
// NUMBER in the code lets the system (and the reader) know that this
// is definitely a dimensionless value, rather than possibly a value
// with units which the user forgot to wrap in the correct <logical type>.
//
//COSI_DEFINE_TYPEDVAL(factor_t);
typedef double factor_t;

const factor_t ZERO_FACTOR(0.0);

//
// Macro: COSI_DEFINE_TYPEDVAL_MULT
//
// Specify a multiplication relationship between logical types: for example,
// speed multiplied by time yields distance.  The macro defines multiplication
// and division operators implied by this relationship.
//
#define COSI_DEFINE_TYPEDVAL_MULT(TypeA, TypeB, TypeR)									\
	inline TypeR operator*( const TypeA& a, const TypeB& b ) { return TypeR( ToDouble( a ) * ToDouble( b ) ); } \
	inline TypeR operator*( const TypeB& b, const TypeA& a ) { return TypeR( ToDouble( a ) * ToDouble( b ) ); } \
	inline TypeA operator/( const TypeR& r, const TypeB& b ) { return TypeA( ToDouble( r ) / ToDouble( b ) ); } \
	inline TypeB operator/( const TypeR& r, const TypeA& a ) { return TypeB( ToDouble( r ) / ToDouble( a ) ); }

// Multiplying a dimensionless <factor_t> by a <typed value> yields a value of the same type.
template <typename T, typename V> inline T operator*( const factor_t& f, const TypedVal<T,V>& v ) { return T( ToDouble( f ) * ToDouble( v ) ); }
template <typename T, typename V> inline T operator*( const TypedVal<T,V>& v, const factor_t& f ) { return T( ToDouble( f ) * ToDouble( v ) ); }
template <typename T, typename V> inline T operator/( const TypedVal<T,V>& v, const factor_t& f ) { return T( ToDouble( v ) / ToDouble( f ) ); }
template <typename T, typename V> inline factor_t operator/( const TypedVal<T,V>& v1, const TypedVal<T,V>& v2 ) { return factor_t( ToDouble( v1 ) / ToDouble( v2 ) ); }

//
// Class template: TypedValRel
//
// A "homogeneous" typed value on which all addition and subtraction operations make sense.
// For example, a length: you can add two lengths.  Or a probability: you can add two probabilities.
// Contrast with <TypedValAbs>, which represents "point" values; adding two point values does not
// make sense, and subtracting two point values yields a homogeneous value -- e.g. subtracting
// two positions gives the distance between them.
//
template <class T, typename V = cosi_double>
struct TypedValRel: public TypedVal<T,V> {
	 TypedValRel() { }
	 explicit TypedValRel( cosi_double val_ ): TypedVal<T,V>( val_ ) { }
	 explicit TypedValRel( const zero_t& val_ ): TypedVal<T,V>( val_ ) { }
	 explicit TypedValRel( const nan_t& val_ ): TypedVal<T,V>( val_ ) { }

	 T& operator=( const zero_t& ) { this->val = 0; return (T&)*this; }
	 T& operator=( const nan_t& ) { this->val = std::numeric_limits<V>::quiet_NaN(); return (T&)*this; }
	 
	 T& operator-=( const TypedValRel<T,V>& rv ) { this->val -= rv.val; return (T&)*this; }
	 T& operator+=( const TypedValRel<T,V>& rv ) { this->val += rv.val; return (T&)*this; }

	 T& operator*=( const factor_t& factor ) { this->val *= factor; return (T&)*this; }
	 T& operator/=( const factor_t& factor ) { this->val /= factor; return (T&)*this; }
	 const T operator*( const factor_t& f ) const { return T( this->val * f ); }
};

template <class T, typename V> inline
const T operator-( const TypedValRel<T,V>& rv1, const TypedValRel<T,V>& rv2 ) { return T( rv1.val - rv2.val ); }
template <class T, typename V> inline
const T operator+( const TypedValRel<T,V>& rv1, const TypedValRel<T,V>& rv2 ) { return T( rv1.val + rv2.val ); }

template <class T, typename V> inline
const T operator*( const factor_t& f, const TypedValRel<T,V>& rv ) { return T( f * rv.val  ); }

template <class T, typename V> inline
const T operator-( const TypedValRel<T,V>& rv ) { return T( -rv.val ); }


#define COSI_DEFINE_TYPEDVAL_REL(Type)													\
	struct Type: public TypedValRel<Type> {												\
		Type() {}																										\
		explicit Type( cosi_double val_ ): TypedValRel<Type>( val_ ) {}	\
  }

// 
// Class template: TypedValAbs
//
// An "absolute" typed value for which addition does not make sense.
// For example, adding two locations does not make sense; subtracting
// two locations yields not a location but a length.
template <class T, class RelT, typename V = cosi_double >
struct TypedValAbs: public TypedVal<T> {
	 typedef RelT rel_t;

	 TypedValAbs() {}
	 explicit TypedValAbs( cosi_double val_ ): TypedVal<T>( val_ ) {}
	 explicit TypedValAbs( const zero_t& val_ ): TypedVal<T,V>( val_ ) { }
	 explicit TypedValAbs( const nan_t& val_ ): TypedVal<T,V>( val_ ) { }

	 T& operator=( const zero_t& ) { this->val = 0; return (T&)*this; }
	 T& operator=( const nan_t& ) { this->val = std::numeric_limits<V>::quiet_NaN(); return (T&)*this; }
	 
	 T& operator-=( const rel_t& rv ) { this->val -= rv.val; return (T&)*this; };
	 T& operator+=( const rel_t& rv ) { this->val += rv.val; return (T&)*this; };

   T operator+( const rel_t& rv ) const { return T( this->val + rv.val ); }
   T operator-( const rel_t& rv ) const { return T( this->val - rv.val ); }

   RelT operator-( const TypedValAbs<T,RelT,V>& rv ) const { return RelT( this->val - rv.val ); }
};  // struct TypedValAbs

//
// Macro: COSI_DEFINE_TYPEDVAL_ABSREL
//
// Define a pair of related logical types, one for absolute values and
// one for relative values -- for example, location and distance.
// The two types have the same units,
// but adding two values of the absolute type does not make sense and
// subtracting two values of the absolute type yields a value of the
// relative type.
//
#define COSI_DEFINE_TYPEDVAL_ABSREL(AbsT,RelT,ValT)									\
	struct RelT: public TypedValRel<RelT,ValT> {											\
		RelT() { }																											\
		explicit RelT( ValT val_ ): TypedValRel<RelT,ValT>( val_ ) { }	\
	};																																\
																																		\
  struct AbsT: public TypedValAbs<AbsT,RelT,ValT> {									\
		AbsT() { }																											\
		explicit AbsT( ValT val_ ):																			\
			TypedValAbs<AbsT,RelT,ValT>( val_ ) { }												\
	}

//
// Macro: COSI_DEFINE_TYPEDVAL_ID
//
// Define an int <TypedVal> that would serve as an identifier.
//
#define COSI_DEFINE_TYPEDVAL_ID(Type)													\
	struct Type: public TypedVal<Type,int> {										\
		Type() {}																									\
		explicit Type( int val_ ): TypedVal<Type,int>( val_ ) {}	\
  }

template <class T>
inline int ToInt( const TypedVal<T,int>& v ) { return v.val; }

template <typename T> struct cosi_converter< TypedVal<T,int>, std::size_t >: public boost::true_type {
	 std::size_t operator()( const TypedVal<T,int>& v ) const {
		 assert( v.val >= 0 );
		 return static_cast<std::size_t>( v.val );
	 }
};

#endif  // #ifndef COSI_TYPEDVAL_DISABLE 

}  // namespace cosi

#ifndef COSI_TYPEDEVAL_DISABLE
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::TypedVal, 2 )
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::TypedValRel, 2 )
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::TypedValAbs, 2 )
BOOST_TYPEOF_REGISTER_TYPE(cosi::factor_t)
BOOST_TYPEOF_REGISTER_TYPE(cosi::zero_t)
#endif // #ifndef COSI_TYPEDVAL_DISABLE

#endif // #ifndef __INCLUDE_COSI_TYPEDVAL_H

