
// * Preamble

// ** Headers

#ifndef __INCLUDE_COSI_GENERALMATH_H
#define __INCLUDE_COSI_GENERALMATH_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include <map>
#include <functional>
#include <string>
#include <ios>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/phoenix.hpp>
#include <boost/next_prior.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/concept_check.hpp>
#include <boost/proto/functional/std/utility.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/declval.hpp>
#include <boost/fusion/include/std_pair.hpp> 
#include <boost/compressed_pair.hpp>
#include <boost/phoenix/fusion/at.hpp>
#include <boost/phoenix/function.hpp>
#include <boost/static_assert.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/mpl/int.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/exception/all.hpp>
#include <cosi/gauss_legendre.h>
#include <cosi/utils.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()


namespace cosi {
namespace math {

//
// ** Utility metafunctions
//

// Metafunction: MultType
// Returns the type of the multiplication of values of the given type.
template <typename TVal1, typename TVal2> struct MultType {
	 typedef BOOST_TYPEOF_TPL( boost::declval<TVal1>() * boost::declval<TVal2>()) type;
};

// Metafunction: DiffType
// Returns the type of the difference of values of the given type.
template <typename TVal> struct DiffType {
   typedef BOOST_TYPEOF_TPL( boost::declval<TVal>() - boost::declval<TVal>()) type;
};

// Metafunction: AreaType
// The type for an area, i.e. the result of integrating a function with given domain and range.
template <typename TDomain, typename TRange>
struct AreaType {
   typedef typename DiffType<TDomain>::type domain_diff_type;
   typedef typename DiffType<TRange>::type range_diff_type;
   typedef typename MultType<domain_diff_type, range_diff_type>::type type;
};



// * Functions
// ** Generic function concept

// *** Metafunction <<DomainType>> - Returns the domain type of a function
template <typename T> struct DomainType;

// *** Metafunction <<RangeType>> - Returns the range type of a function
template <typename T> struct RangeType;

// *** Concept <<FunctionConcept>>
//
//    Concept checking class for general math functions of one variable.
//
//    Modeled by: [[Function]] (defined below).
template <typename F>
class FunctionConcept {
public:
   //
   // Associated types: for each generic math function there is a domain and a range type.
   // Normally these would be floating-point types.
   //
   
	 typedef typename DomainType<F>::type argument_type;
   typedef typename RangeType<F>::type result_type;

   BOOST_CONCEPT_USAGE(FunctionConcept) {
     range_var = eval( const_cast<const F&>( func ), domain_var );
   }

private:
   F func;
   argument_type domain_var;
   result_type range_var;
   
};  // FunctionConcept

template <typename T> std::string getLabel( const T& );

// *** GenericClass Function - A mathematical function of one variable.
template <typename TDomain, typename TRange, typename TSpec> class Function;
        
// Template params:

//   - TDomain - the type of values in the function's domain
//   - TRange - the type of values in the function's range
//   - TSpec - the specialization defining what kind of function this is.
         
// *** Metafunction: SpecType
//    Returns the specialization type of a Function; this defines what kind of function it is
//    (constant, linear etc).
template <typename T> struct SpecType;

// Metapredicate: IsFunctionSpec
// Tests whether a given type is a valid function specialization.
template <typename T, typename Enabler = void> struct IsFunctionSpec: public boost::false_type { };

template <typename TFunc> struct DomainType { typedef typename TFunc::argument_type type; };
template <typename TFunc> struct RangeType { typedef typename TFunc::result_type type; };
template <typename TFunc> struct SpecType { typedef typename TFunc::spec_type type; };

//
// *** Generic function: eval
//
// Default implementation of eval(f,x) to evaluate a Function at a given point in its domain;
// just calls the function's operator() .
// 
template <typename TFunc, typename TArg>
typename boost::enable_if< boost::is_convertible< TArg, typename DomainType<TFunc>::type >,
                           typename RangeType<TFunc>::type >::type
eval( const TFunc& f, TArg x ) { return f( x ); }


// ** Particular function kinds
// *** Constant functions
//    We support run-time constant functions and compile-time constant functions.
            

struct RunTime;
template <typename TConstSpec = RunTime> struct Const;
template <> struct IsFunctionSpec< Const<RunTime> >: public boost::true_type {};

template <int N> struct CompileTime { static const int value = N; };

template <int N> struct IsFunctionSpec< Const< CompileTime<N> > >: public boost::true_type {};

template <typename TDomain, typename TRange, typename TSpec>
TRange evalConst( const Function< TDomain, TRange, Const< TSpec > >& f ) { return eval( f, TDomain() ); }

// template <typename TDomain, typename TRange, typename TSpec>
// std::ostream& operator<<( std::ostream& s, const Function<TDomain, TRange, TSpec>& f ) {
// 	s << "Function[ " << typeid(TSpec).name() << "]";
// 	return s;
// }

template <typename TDomain, typename TRange>
class Function<TDomain, TRange, Const< RunTime > > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef Const<RunTime> spec_type;

   Function( TRange val_ ): val( val_ ) { }

   TRange operator()( TDomain ) const { return val; }

   friend ostream& operator<<( ostream& s, const Function& f ) {
     s << "ConstRuntime[" << f.val << "]"; return s;
   }

	 bool operator==( const Function& f ) const { return val == f.val; }
	 bool operator!=( const Function& f ) const { return val != f.val; }

	 TRange& getValRef() { return val; }

private:
   TRange val;
};


template <typename TDomain, typename TRange, int val>
class Function<TDomain, TRange, Const< CompileTime<val> > > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef Const< CompileTime<val> > spec_type;
   TRange operator()( TDomain ) const { return static_cast<TRange>( val ); }

   friend ostream& operator<<( ostream& s, const Function& f ) {
     s << "ConstCompiletime[" << val << "]"; return s;
   }
   
};

// *** Monomials (functions of the form a x^k) 

template <int exponent, typename FactorType = double> struct X_To;
template <int exponent, typename FactorType>
struct IsFunctionSpec< X_To<exponent, FactorType> > : public boost::true_type {};

template <typename TDomain, typename TRange>
class Function<TDomain, TRange, X_To<0> >: public Function<TDomain, TRange,
                                                           Const< CompileTime<1> > > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef X_To<0> spec_type;
   TRange operator()( TDomain ) const { return static_cast<TRange>( 1 ); }
};


template <typename TDomain, typename TRange, typename FactorType>
class Function<TDomain, TRange,
							 X_To<1, FactorType> > {
public:
   typedef TDomain argument_type;
	 typedef typename DiffType<TDomain>::type domain_diff_type;
   //typedef typename MultType<FactorType, domain_diff_type >::type result_type;
	 typedef TRange result_type;
   typedef X_To<1,FactorType> spec_type;
   typedef FactorType factor_type;

	 BOOST_MPL_ASSERT(( boost::is_convertible< typename MultType<FactorType, domain_diff_type >::type,
											result_type> ));

   Function( ): factor( 1.0 ) {}
   Function( factor_type factor_ ): factor( factor_ ) {}
   result_type operator()( TDomain x ) const { return factor * ( x - static_cast< TDomain >( 0.0 ) ) ; }
   factor_type getFactor() const { return factor; }

   friend ostream& operator<<( ostream& s, const Function& f ) {
     s << "(" << f.factor << " * x" << ")"; return s; 
   }

   
private:
   factor_type factor;
};

template <typename TDomain, typename FactorType>
class Function<TDomain, typename MultType< FactorType,
                                           typename MultType< TDomain, TDomain >::type >::type,
               X_To<2,FactorType> > {
public:
   typedef TDomain argument_type;
   typedef typename MultType< FactorType,
                              typename MultType< TDomain, TDomain >::type >::type result_type;
   typedef X_To<2,FactorType> spec_type;
   typedef FactorType factor_type;
   Function( factor_type factor_ ): factor( factor_ ) {}
   result_type operator()( TDomain x ) const { return factor * x*x; }
   factor_type getFactor() const { return factor; }
   friend ostream& operator<<( ostream& s, const Function& f ) {
     s << "(" << f.factor << " * x^2" << ")"; return s; 
   }
private:
   const factor_type factor;
};

// ** Operations on functions

//   Construction of compound functions: e.g. given two functions we can
//   add them to create a new function that computes the sum of the two original
//   functions' results.

// *** Ops

template <typename T1, typename T2 = T1>
struct AddOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() + boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 + a2; }

   friend ostream& operator<<( ostream& s, const AddOp& ) { s << "+"; return s; }
};
template <typename T1, typename T2 = T1>
struct SubOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() - boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 - a2; }
   friend ostream& operator<<( ostream& s, const SubOp& ) { s << "-"; return s; }
};
template <typename T1, typename T2 = T1>
struct MultOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() * boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 * a2; }
   friend ostream& operator<<( ostream& s, const MultOp& ) { s << "*"; return s; }
};
template <typename T1, typename T2 = T1>
struct DivOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() / boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 / a2; }
   friend ostream& operator<<( ostream& s, const DivOp& ) { s << "/"; return s; }
};

// *** Binary Op Function

template <typename TSpec1, typename TSpec2, typename Op> struct BinOp;
template <typename TSpec1, typename TSpec2, typename Op>
struct IsFunctionSpec< BinOp< TSpec1, TSpec2, Op >,
                       typename boost::enable_if< typename boost::mpl::and_< IsFunctionSpec< TSpec1 >,
                                                                             IsFunctionSpec< TSpec2 > > >::type >: public boost::true_type {};



template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec1, typename TSpec2, template <typename,typename> class Op>
class Function<TDomain, typename Op<TRange1,TRange2>::result_type,
               BinOp<TSpec1, TSpec2, Op<TRange1,TRange2> > >: public boost::compressed_pair< Function<TDomain, TRange1, TSpec1>, Function<TDomain, TRange2, TSpec2> > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
public:
   typedef TDomain argument_type;
   typedef typename Op<TRange1,TRange2>::result_type result_type;
   typedef BinOp<TSpec1,TSpec2,Op<TRange1,TRange2> > spec_type;
   typedef TRange1 range1_type;
   typedef TRange2 range2_type;
   typedef TSpec1 spec1_type;
   typedef TSpec2 spec2_type;
   
   typedef Function<TDomain, TRange1, TSpec1> function_type_1; 
   typedef Function<TDomain, TRange2, TSpec2> function_type_2;

   Function( const function_type_1& f1_, const function_type_2& f2_ ):
     boost::compressed_pair< function_type_1, function_type_2 >( f1_, f2_ ) { }

   result_type operator()( TDomain x ) const { return binop( eval( this->first(), x ),
                                                             eval( this->second(), x ) ); }

   friend ostream& operator<<( ostream& s, const Function& f ) {
     s << "(" << f.first() << " " << f.binop << f.second() << ")"; return s;
   }
   
private:
   static const Op<TRange1,TRange2> binop;
   
};

template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec1, typename TSpec2, template<typename,typename> class Op>
const Op<TRange1,TRange2> Function<TDomain, typename Op<TRange1,TRange2>::result_type,
                                   BinOp< TSpec1, TSpec2,
                                          Op<TRange1,TRange2> > >::binop = Op<TRange1,TRange2>();

// *** Operator overloads for adding/subtracting/etc two functions
// 
// Adding two functions f and g creates a new function that for x computes f(x)+g(x)

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
Function<TDomain, typename AddOp<TRange1,TRange2>::result_type,
         BinOp<TSpec1,TSpec2, AddOp<TRange1,TRange2> > >
operator+( const Function<TDomain,TRange1,TSpec1>& f1,
           const Function<TDomain,TRange2,TSpec2>& f2 ) {
  return Function<TDomain, typename AddOp<TRange1,TRange2>::result_type,
                  BinOp<TSpec1,TSpec2,AddOp< TRange1, TRange2> > >( f1, f2 );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
Function<TDomain, typename SubOp<TRange1,TRange2>::result_type,
         BinOp<TSpec1,TSpec2, SubOp<TRange1,TRange2> > >
operator-( const Function<TDomain,TRange1,TSpec1>& f1,
           const Function<TDomain,TRange2,TSpec2>& f2 ) {
  return Function<TDomain, typename SubOp<TRange1,TRange2>::result_type,
                  BinOp<TSpec1,TSpec2,SubOp< TRange1, TRange2> > >( f1, f2 );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
Function<TDomain, typename MultOp<TRange1,TRange2>::result_type,
         BinOp<TSpec1,TSpec2, MultOp<TRange1,TRange2> > >
operator*( const Function<TDomain,TRange1,TSpec1>& f1,
           const Function<TDomain,TRange2,TSpec2>& f2 ) {
  return Function<TDomain, typename MultOp<TRange1,TRange2>::result_type,
                  BinOp<TSpec1,TSpec2,MultOp< TRange1, TRange2> > >( f1, f2 );
}
template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
Function<TDomain, typename DivOp<TRange1,TRange2>::result_type,
         BinOp<TSpec1,TSpec2, DivOp<TRange1,TRange2> > >
operator/( const Function<TDomain,TRange1,TSpec1>& f1,
           const Function<TDomain,TRange2,TSpec2>& f2 ) {
  return Function<TDomain, typename DivOp<TRange1,TRange2>::result_type,
                  BinOp<TSpec1,TSpec2,DivOp< TRange1, TRange2> > >( f1, f2 );
}

// *** Operations on const funcs
//
// Operations on constant functions: produce a new constant function, whose return value
// can be computed immediately.
//


template <typename TDomain, typename TRange1, typename TRange2>
Function<TDomain, typename AddOp<TRange1,TRange2>::result_type, Const<> >
operator+( const Function< TDomain, TRange1, Const<> >& f1,
           const Function< TDomain, TRange2, Const<> >& f2 ) {
  return Function< TDomain, typename AddOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) + evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2>
Function<TDomain, typename SubOp<TRange1,TRange2>::result_type, Const<> >
operator-( const Function< TDomain, TRange1, Const<> >& f1,
           const Function< TDomain, TRange2, Const<> >& f2 ) {
  return Function< TDomain, typename SubOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) - evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2>
Function<TDomain, typename MultOp<TRange1,TRange2>::result_type, Const<> >
operator*( const Function< TDomain, TRange1, Const<> >& f1,
           const Function< TDomain, TRange2, Const<> >& f2 ) {
  return Function< TDomain, typename MultOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) * evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2>
Function<TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >
operator/( const Function< TDomain, TRange1, Const<> >& f1,
           const Function< TDomain, TRange2, Const<> >& f2 ) {
  return Function< TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) / evalConst( f2 ) );
}

//
// *** Unary operations on functions: specificaly, negation.
//
//  For a function f, -f is a function that for value x returns -f(x).
//

template <typename TSpec, typename Op> struct UnaryOp;

template <typename TDomain, typename TRange, typename TSpec, typename Op>
class Function<TDomain, TRange, UnaryOp<TSpec, Op> > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef UnaryOp<TSpec,Op> spec_type;
   typedef TSpec unop_spec_type;
   typedef Op op_type;
   
   typedef Function<TDomain, TRange, TSpec> function_type;

   Function( const function_type& f_ ): f( f_ ) { }

   TRange operator()( TDomain x ) const { return op( f( x ) ); }

   friend ostream& operator<<( ostream& s, const Function& f ) { s << op << " " << f.f; return s; }
   
private:
   static const Op op;
   const function_type f;
   
};

template <typename TDomain, typename TRange, typename TSpec, typename Op>
const Op Function<TDomain, TRange, UnaryOp< TSpec, Op > >::op = Op();

template <typename T> 
ostream& operator<<( ostream& s, const std::negate< T >& ) { s << "-"; return s; }


template <typename TDomain, typename TRange, typename TSpec>
Function<TDomain, TRange, UnaryOp<TSpec, std::negate< TRange > > >
operator-( const Function<TDomain,TRange,TSpec>& f ) {
  return Function<TDomain, TRange, UnaryOp<TSpec, std::negate< TRange > > >( f );
}

// ** Linear functions: creating a linear function from two points through which it passes.

// *** Construct line through two points

namespace detail {
template <typename TDomain, typename TRange>
struct result_of_make_line_through {
   typedef 
   BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain, TRange, X_To<1> > >() +
                      boost::declval< Function<TDomain, TRange, Const<> > >() )) type;
};
}
template <typename TDomain, typename TRange>
typename detail::result_of_make_line_through<TDomain,TRange>::type
make_line_through( TDomain x1, TRange y1, TDomain x2, TRange y2 ) {
  assert( x1 != x2 );
  BOOST_AUTO_TPL( slope, ( y2 - y1 ) / ( x2 - x1 ) );
  BOOST_AUTO_TPL( b, y1 - ( slope * x1 ) );
  return //Function<TDomain, TRange, Const<> >( slope ) *
     Function<TDomain, TRange, X_To<1> >( slope) +
     Function<TDomain, TRange, Const<> >( b );
}


// ** Piecewise functions

//       A piecewise function is a compound function which divides the domain into a set of intervals,
//       and on each interval returns the result of a specified piece function.
//       Operations on piecewise functions can often be expressed in terms of operations on the pieces.

#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

template <typename TPieceSpec> struct Piecewise;
template <typename TPieceSpec> struct IsFunctionSpec< Piecewise<TPieceSpec>, typename boost::enable_if< IsFunctionSpec<TPieceSpec> >::type >: public boost::true_type {};

template <typename TDomain, typename TRange, typename TPieceSpec>
class Function<TDomain, TRange, Piecewise< TPieceSpec > > {
	 BOOST_MPL_ASSERT(( IsFunctionSpec<TPieceSpec> ));
public:
	 typedef TDomain argument_type;
	 typedef TRange result_type;
	 typedef Piecewise< TPieceSpec > spec_type;
	 typedef TPieceSpec piece_spec_type;
      
	 BOOST_CONCEPT_ASSERT((boost::LessThanComparable<TDomain>));
   
	 typedef Function< TDomain, TRange, TPieceSpec > piece_function_type;
	 BOOST_CONCEPT_ASSERT((FunctionConcept< piece_function_type >));
   
	 Function() { }
	 template <typename SinglePassRange>
	 Function( const SinglePassRange& r ):
		 pieces( boost::begin( r ), boost::end( r ) ) {
				
		 BOOST_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< SinglePassRange > ));
		 typedef typename boost::range_value< SinglePassRange >::type range_value_type;
		 BOOST_STATIC_ASSERT(( boost::is_convertible< range_value_type,
													 typename pieces_type::value_type >::value ));
	 }
   
	 // template <typename TDomain2>
	 // void addPiece( TDomain2 x,
	 //                 const piece_function_type& piece,
	 //                 typename boost::enable_if< boost::is_convertible< TDomain2, TDomain >, TDomain2>::type *dummy = NULL ) { (void)dummy; pieces.insert( std::make_pair( x, piece ) ); }
   
	 // Type: pieces_type
	 // The type for a map from point to the piece function starting at that point.
	 // The map is ordered from the rightmost piece to the leftmost, so that map::lower_bound()
	 // can be used to find the piece corresponding to a given domain point.
	 typedef std::map< TDomain, piece_function_type, std::greater<TDomain> > pieces_type;
   
	 // Method: getPieces
	 // Returns a map from domain point to the piece function starting at that point.
	 // The map is ordered from rightmost piece to the left.
	 // The map is mutable; change it to change the definition of this piecewise function.
	 const pieces_type& getPieces() const { return pieces; }
	 pieces_type& getPieces() { return pieces; }
   
	 TRange operator()( TDomain x ) const {
		 typename pieces_type::const_iterator it = pieces.lower_bound( x );
		 return it == pieces.end() ? std::numeric_limits<TRange>::quiet_NaN() :
				eval( it->second, x );
	 }
   
	 friend ostream& operator<<( ostream& s, const Function& f ) {
		 s << "Piecewise[";
		 for ( BOOST_AUTO_TPL( it, f.pieces.rbegin() ); it != f.pieces.rend(); it++ )
				s << "(" << it->first << ", " << it->second << ")";
		 s << "]";
		 return s;
	 }

	 bool operator==( const Function& f ) const { return pieces == f.pieces; }
	 bool operator!=( const Function& f ) const { return pieces != f.pieces; }
   
private:
	 // Field: pieces
	 // Map from domain point to the piece function starting at that point; ordered from rightmost point to leftmost.
	 pieces_type pieces;
};  // class Function<TDomain, TRange, Piecewise< TPieceSpec > >

#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic pop
#endif
#endif

// Method: loadFrom - read piecewise function from a file
template <typename TDomain, typename TRange>
void loadFrom( istream& s, Function<TDomain, TRange, Piecewise< Const<> > >& f ) {
	boost::io::ios_exception_saver save_exceptions( s );
	s.exceptions( ios::failbit | ios::badbit );
	f.getPieces().clear();
	try {
		while ( !s.eof() ) {
			TDomain d;
			TRange r;
			s >> d >> r;
			f.getPieces().insert( std::make_pair( d, r ) );
			s >> std::ws;
		}
	} catch( const std::ios_base::failure& e ) {
		throw cosi_io_error() << boost::errinfo_nested_exception( boost::current_exception() );
	}
}
	 

// *** Operations on F,P where P is a piecewise function
// Propagate to the pieces 

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
struct result_of_mult {
	 typedef BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain,TRange1,TSpec1> >() *
															boost::declval< Function<TDomain,TRange2,TPieceSpec2> >() )) piece_function_type;
	 typedef typename RangeType< piece_function_type >::type range_type;
	 typedef typename SpecType< piece_function_type >::type piece_spec_type;
	 typedef Function< TDomain, range_type, Piecewise< piece_spec_type > > type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
typename result_of_mult<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type
operator*( const Function<TDomain,TRange1,TSpec1>& f1, const Function<TDomain,TRange2,Piecewise<TPieceSpec2> >& f2 ) {
	typename result_of_mult<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type f_result;
 
	for ( typename Function< TDomain, TRange2, Piecewise<TPieceSpec2> >::pieces_type::const_iterator it =
					 f2.getPieces().begin(); it != f2.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, f1 * it->second ) );
	return f_result;
}


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
struct result_of_div {
	 typedef BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain,TRange1,TSpec1> >() /
															boost::declval< Function<TDomain,TRange2,TPieceSpec2> >() )) piece_function_type;
	 typedef typename RangeType< piece_function_type >::type range_type;
	 typedef typename SpecType< piece_function_type >::type piece_spec_type;
	 typedef Function< TDomain, range_type, Piecewise< piece_spec_type > > type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
typename result_of_div<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type
operator/( const Function<TDomain,TRange1,TSpec1>& f1, const Function<TDomain,TRange2,Piecewise<TPieceSpec2> >& f2 ) {
	typename result_of_div<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type f_result;
 
	for ( typename Function< TDomain, TRange2, Piecewise<TPieceSpec2> >::pieces_type::const_iterator it =
					 f2.getPieces().begin(); it != f2.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, f1 / it->second ) );
	return f_result;
}
   

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
struct result_of_add {
	 typedef BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain,TRange1,TSpec1> >() +
															boost::declval< Function<TDomain,TRange2,TPieceSpec2> >() )) piece_function_type;
	 typedef typename RangeType< piece_function_type >::type range_type;
	 typedef typename SpecType< piece_function_type >::type piece_spec_type;
	 typedef Function< TDomain, range_type, Piecewise< piece_spec_type > > type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
typename result_of_add<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type
operator+( const Function<TDomain,TRange1,TSpec1>& f1, const Function<TDomain,TRange2,Piecewise<TPieceSpec2> >& f2 ) {
	typename result_of_add<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type f_result;
 
	for ( typename Function< TDomain, TRange2, Piecewise<TPieceSpec2> >::pieces_type::const_iterator it =
					 f2.getPieces().begin(); it != f2.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, f1 + it->second ) );
	return f_result;
}


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
struct result_of_sub {
	 typedef BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain,TRange1,TSpec1> >() -
															boost::declval< Function<TDomain,TRange2,TPieceSpec2> >() )) piece_function_type;
	 typedef typename RangeType< piece_function_type >::type range_type;
	 typedef typename SpecType< piece_function_type >::type piece_spec_type;
	 typedef Function< TDomain, range_type, Piecewise< piece_spec_type > > type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
typename result_of_sub<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type
operator-( const Function<TDomain,TRange1,TSpec1>& f1, const Function<TDomain,TRange2,Piecewise<TPieceSpec2> >& f2 ) {
	typename result_of_sub<TDomain, TRange1, TRange2, TSpec1, TPieceSpec2>::type f_result;
 
	for ( typename Function< TDomain, TRange2, Piecewise<TPieceSpec2> >::pieces_type::const_iterator it =
					 f2.getPieces().begin(); it != f2.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, f1 - it->second ) );
	return f_result;
}


// ** Integrals


// *** Definition of indefiniteIntegral function

template <typename TDomain, typename TRange, typename TSpec> struct result_of_indefiniteIntegral {
	 BOOST_STATIC_ASSERT_MSG( sizeof( TDomain ) == 0, "don't know how to get type of indef integral of f" );
};

template <typename TDomain, typename TRange, typename TSpec>
void indefiniteIntegral( const Function< TDomain, TRange, TSpec >& f) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  BOOST_STATIC_ASSERT_MSG( sizeof( f ) == 0, "don't know how to integrate f" );
}

//
// *** Definite integrals
//

template <typename TDomain, typename TRange, typename TSpec>
typename AreaType<TDomain,TRange>::type
definiteIntegral( const Function< TDomain, TRange, TSpec >& f, TDomain a, TDomain b ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  BOOST_AUTO( f_int, indefiniteIntegral( f ) );
  return eval( f_int, b ) - eval( f_int, a );
}

template <typename TDomain, typename TRange, typename TSpec>
struct result_of_integralFromPoint {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
   typedef Function< TDomain, TRange, TSpec > function_type;
   typedef typename result_of_indefiniteIntegral< TDomain, TRange, TSpec >::type f_int;
   typedef typename f_int::result_type integralValueType;
   typedef Function< TDomain, integralValueType, Const<> > const_func_type;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f_int>() - boost::declval<const_func_type>())) type; 
};

template <typename TDomain, typename TRange, typename TSpec>
typename result_of_integralFromPoint< TDomain, TRange, TSpec >::type
integralFromPoint( const Function< TDomain, TRange, TSpec >& f, TDomain x ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  BOOST_AUTO( f_int, indefiniteIntegral( f ) );
  typedef BOOST_TYPEOF_TPL( f_int ) integralFuncType;
  typedef typename integralFuncType::result_type integral_result_type;
  integral_result_type integralAtX( eval( f_int, x ) );
  Function< TDomain, integral_result_type,
            Const<> > f_const( integralAtX );
  return f_int - f_const;
}

//
// *** Integral of a constant
//

template <typename TDomain, typename TRange>
struct result_of_indefiniteIntegral<TDomain, TRange, Const<> > {
	 typedef typename DiffType<TDomain>::type domain_diff_type;
	 typedef typename DiffType<TRange>::type range_diff_type;
	 typedef typename AreaType<TDomain,TRange>::type integral_type;
	 typedef Function< TDomain, integral_type, X_To<1, range_diff_type> > type;
};

template <typename TDomain, typename TRange>
typename result_of_indefiniteIntegral<TDomain,TRange,Const<> >::type
indefiniteIntegral( const Function< TDomain, TRange, Const<> >& f ) {
  typedef typename DiffType<TRange>::type range_diff_type;
  typedef typename AreaType<TDomain,TRange>::type integral_type;
  return Function< TDomain, integral_type,
                   X_To<1, range_diff_type> >( evalConst( f ) - static_cast<TRange>( 0 ) );
}

//
// *** Integral of a monomial
//

template <typename TDomain, typename TRange, typename TFactor>
struct result_of_indefiniteIntegral<TDomain, TRange, X_To<1,TFactor> > {

   typedef Function< TDomain,
										 typename AreaType<TDomain,TRange>::type,X_To<2,TFactor> > type;
};

template <typename TDomain, typename TRange, typename TFactor>
typename result_of_indefiniteIntegral<TDomain,TRange, X_To<1,TFactor> >::type
indefiniteIntegral( const Function< TDomain, TRange, X_To<1,TFactor> >& f ) {
	typedef typename result_of_indefiniteIntegral<TDomain,TRange,X_To<1,TFactor> >::type result_function_type;
  return result_function_type( 0.5 * f.getFactor() );
}


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
struct result_of_indefiniteIntegral<TDomain, typename MultOp<TRange1,TRange2>::result_type,
                                    BinOp< Const<>, TSpec, MultOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
   typedef Function< TDomain, TRange1, Const<> > f_const_t;
   typedef typename result_of_indefiniteIntegral<TDomain, TRange2, TSpec>::type f_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f_const_t>() * boost::declval<f_t>() )) type;
};

//
// *** Linearity of integration
//

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
typename result_of_indefiniteIntegral<TDomain, typename MultOp<TRange1, TRange2>::result_type,
                                      BinOp< Const<>, TSpec,
                                             MultOp<TRange1,TRange2> > >::type
indefiniteIntegral( const Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
                    BinOp< Const<>, TSpec,
                    MultOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  return f.first() * indefiniteIntegral( f.second() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
typename result_of_indefiniteIntegral<TDomain, typename MultOp<TRange1, TRange2>::result_type,
                                      BinOp< Const<>, TSpec,
                                             MultOp<TRange1,TRange2> > >::type
indefiniteIntegral( const Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec, Const<>,
                    MultOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  return f.second() * indefiniteIntegral( f.first() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_indefiniteIntegral<TDomain, typename AddOp<TRange1,TRange2>::result_type,
                                    BinOp<TSpec1, TSpec2, AddOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
   typedef typename result_of_indefiniteIntegral<TDomain, TRange1, TSpec1>::type f1_t;
   typedef typename result_of_indefiniteIntegral<TDomain, TRange2, TSpec2>::type f2_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f1_t>() + boost::declval<f2_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_indefiniteIntegral<TDomain, typename AddOp<TRange1, TRange2>::result_type,
                                      BinOp< TSpec1, TSpec2,
                                             AddOp<TRange1,TRange2> > >::type
indefiniteIntegral( const Function< TDomain, typename AddOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec1, TSpec2,
                    AddOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
  return indefiniteIntegral( f.first() ) + indefiniteIntegral( f.second() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_indefiniteIntegral<TDomain, typename SubOp<TRange1,TRange2>::result_type,
                                    BinOp<TSpec1, TSpec2, SubOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
   typedef typename result_of_indefiniteIntegral<TDomain, TRange1, TSpec1>::type f1_t;
   typedef typename result_of_indefiniteIntegral<TDomain, TRange2, TSpec2>::type f2_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f1_t>() - boost::declval<f2_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_indefiniteIntegral<TDomain, typename SubOp<TRange1, TRange2>::result_type,
                                      BinOp< TSpec1, TSpec2,
                                             SubOp<TRange1,TRange2> > >::type
indefiniteIntegral( const Function< TDomain, typename SubOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec1, TSpec2,
                    SubOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
  return indefiniteIntegral( f.first() ) - indefiniteIntegral( f.second() );
}

// template <typename TKey, typename TValue, typename TValueTransformer>
// MapValueTransformer<TKey, TValue, TValueTransformer>
// makeMapValueTransformer( const std::map<TKey,TValue>&, TValueTransformer transformer ) {
//  return MapValueTransformer<TKey, TValue, TValueTransformer>( transformer );
// }

//
// *** Integrating piecewise functions
//


template <typename TDomain, typename TRange, typename TPieceSpec>
struct result_of_indefiniteIntegral<TDomain, TRange, Piecewise< TPieceSpec > > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TPieceSpec> ));
   typedef typename result_of_indefiniteIntegral<TDomain,TRange,TPieceSpec>::type piece_integral_type;
   typedef typename DomainType<piece_integral_type>::type piece_integral_domain_type;
   typedef typename RangeType<piece_integral_type>::type piece_integral_range_type;
   typedef typename SpecType<piece_integral_type>::type piece_integral_spec_type;
   typedef BinOp< Const<>, piece_integral_spec_type, AddOp< piece_integral_range_type > >
   const_plus_piece_integral_spec_type;
   typedef Function<piece_integral_domain_type, piece_integral_range_type,
                    Piecewise< const_plus_piece_integral_spec_type > > type;
};

template <typename TDomain, typename TRange, typename TPieceSpec>
typename result_of_indefiniteIntegral<TDomain, TRange, Piecewise< TPieceSpec > >::type
indefiniteIntegral( const Function< TDomain, TRange, Piecewise< TPieceSpec > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TPieceSpec> ));
  
  namespace ad = boost::adaptors;
  typedef typename result_of_indefiniteIntegral<TDomain, TRange, Piecewise< TPieceSpec > >::type
     result_type;
  typedef typename DomainType<result_type>::type result_domain_type;
  typedef typename RangeType<result_type>::type result_range_type;

  result_type result;
  result_range_type definiteIntegralSoFar( 0.0 );

  for ( BOOST_AUTO_TPL( it, f.getPieces().rbegin() ); it != f.getPieces().rend(); it++ ) {
    BOOST_AUTO_TPL( pieceIndefIntegral, indefiniteIntegral( it->second ) );
		BOOST_AUTO_TPL( offsetHere, definiteIntegralSoFar - eval( pieceIndefIntegral,
																															it->first ) );
    result.getPieces().insert( 
			std::make_pair( it->first,
											Function< result_domain_type, result_range_type, Const<> >( offsetHere )
											+ pieceIndefIntegral ) );

    BOOST_AUTO_TPL( it_n, boost::next( it ) );
    if ( it_n != f.getPieces().rend() ) {
      BOOST_AUTO_TPL( at_b, eval( pieceIndefIntegral, it_n->first ) );
      BOOST_AUTO_TPL( at_a, eval( pieceIndefIntegral, it->first ) );
      BOOST_AUTO_TPL( pieceDefiniteIntegral, ( at_b - at_a ) );
      definiteIntegralSoFar += pieceDefiniteIntegral;
    }

  }
  return result;
  
  //return result_type( f.getPieces() | ad::transformed( indefiniteIntegral_pair_caller<TDomain,TRange,TPieceSpec>() ) );
                                                                   
}

template <typename TDomain, typename TRange, typename TSpec>
bool isStrictlyIncreasing( const Function<TDomain, TRange, TSpec>& f, TDomain a, TDomain b,
													 typename DiffType<TDomain>::type step = static_cast< typename DiffType<TDomain>::type >( 1e-3 ) ) {
	assert( a <= b && step >= static_cast< typename DiffType<TDomain>::type >( 0.0 ) );
	TRange prevVal = eval( f, a );
	for ( TDomain x = a + step; x < b; x += step ) {
		TRange valHere = eval( f, x );
		if ( !( prevVal < valHere ) ) {
			PRINT7( f, a, b, step, x, prevVal, valHere );
			return false;
		}
		prevVal = valHere;
	}
	return true;
}


// ** Evaluating inverse functions

//
// *** Function: evalInverse
//
// Evaluate the inverse of a given function at a given point.
//
// Params:
//
//    f - the function.   It must be strictly increasing on the interval [a,b]
//
template <typename TDomain, typename TRange, typename TSpec>
TDomain
evalInverse( const Function<TDomain, TRange, TSpec>& f, TDomain a, TDomain b, TRange targetVal,
             typename DiffType<TRange>::type eps, unsigned maxSteps = 100000 ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
	assert( isStrictlyIncreasing( f, a, b ) );
  unsigned nsteps = 0;
  //PRINT5( a, b, targetVal, eps, typeid(f).name() );
  assert( !(boost::math::isnan)( a ) );
  assert( !(boost::math::isnan)( b ) );
  assert( !(boost::math::isnan)( targetVal ) );
  while ( true ) {
    assert( a <= b );
    TDomain mid = a + ( b-a ) / 2;
    typedef typename DiffType<TRange>::type range_diff_type;
    range_diff_type curDiff = eval( f, mid ) - targetVal;
    //PRINT5( a, b, mid, targetVal, curDiff );
    if ( cosi_fabs( curDiff ) < eps ) {
      return mid;
    }
    if ( nsteps++ > maxSteps ) {
			std::ostringstream msg;
			msg.precision(16);
			msg << "findInfimum: too many steps; eps=" << eps << " curDiff=" << curDiff <<
				 " nsteps=" << nsteps << " maxSteps=" << maxSteps << " a=" << a << " b=" << b << " f=" << f;
			throw std::runtime_error( msg.str() );
		}
    assert( a < b );
    assert( eval( f, a ) < eval( f, mid ) );
    assert( eval( f, mid ) < eval( f, b ) );
    ( curDiff < static_cast<range_diff_type>( 0.0 ) ? a : b ) = mid;
  }
}

// template <typename TDomain, typename TRange, typename TPieceSpec>
// indefiniteIntegral< const Function< TDomain, TRange, Piecewise< TPieceSpec > >& f ) {
  
// }

// template <typename TDomain, typename TRange, typename TSpec>
// typename IntegrationResultType<TDomain,TRange>::type
// integrate( const TFunction<TDomain, TRange, Const<TSpec>& f, TDomain a, TDomain b ) {
//  return ( b - a ) * f.getConstValue();
// }
// Evaluating\ inverse\ functions:1 ends here

// ** Type erasure for Function: a Function that stores any Function.

template <typename TSpec = void> struct Any;
template <> struct IsFunctionSpec< Any<> >: public boost::true_type {};

template <typename TDomain, typename TRange>
class Function<TDomain, TRange, Any<> > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef Any<> spec_type;
   
private:
   struct FunctionObjectConcept {
      virtual ~FunctionObjectConcept() {}
      virtual TRange doEval( TDomain ) const = 0;
   };

   template< typename TSpec > class FunctionObjectModel : public FunctionObjectConcept {
      BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
   public:
      FunctionObjectModel( const Function<TDomain, TRange, TSpec>& f_ ) : f( f_ ) {}
      virtual ~FunctionObjectModel() {}
      virtual TRange doEval( TDomain x ) const { return eval( f, x ); }
   private:
      Function<TDomain, TRange, TSpec> f;
   };

   boost::shared_ptr<FunctionObjectConcept> object;

public:
   Function()  { }
   template< typename TSpec > Function( const Function<TDomain, TRange, TSpec>& obj ) :
     object( new FunctionObjectModel<TSpec>( obj ) ) {}

   template <typename TSpec>
   void reset( const Function<TDomain, TRange, TSpec>& obj ) {
     object.reset( new FunctionObjectModel<TSpec>( obj ) );
   }

   bool empty() const { return !object.get(); }
   operator bool() const { return !empty(); }

   TRange operator()( TDomain x ) const { assert( object.get() ); return object->doEval( x ); }
};  // end: type erasure of a Function


// * Arrival processes

//     Modeling of generic arrival processes, such as Poisson processes.


template <typename TTime, typename TSpec> class ArrivalProcess;

//
// ** Type erasure for ArrivalProcess: storing any ArrivalProcess
//

template <typename TTime, typename TRand>
class ArrivalProcess<TTime, Any< TRand > > {
   struct ArrivalProcessConcept {
      virtual ~ArrivalProcessConcept() {}
      virtual TTime doNextArrival( TTime fromTime, TTime maxTime, double rateFactor,
                                   TRand& randGen, double eps, int maxSteps ) const = 0;
			virtual ostream& print( ostream& s ) const = 0;
			virtual std::string doGetLabel() const = 0;
   };
   template <typename TSpec>
   class ArrivalProcessModel: public ArrivalProcessConcept {
   public:
      ArrivalProcessModel( const ArrivalProcess<TTime, TSpec>& proc_ ): proc( proc_ ) {}
      virtual ~ArrivalProcessModel() {}
      virtual TTime doNextArrival( TTime fromTime, TTime maxTime, double rateFactor,
                                   TRand& randGen, double eps, int maxSteps ) const {
        return proc.nextArrivalTime( fromTime, maxTime, rateFactor, randGen, eps, maxSteps );
      }
			virtual ostream& print( ostream& s ) const { s << proc; return s; }
			virtual std::string doGetLabel() const { return getLabel( proc ); }
   private:
      ArrivalProcess<TTime, TSpec> proc;
   };

   boost::shared_ptr<ArrivalProcessConcept> object;
public:
   ArrivalProcess() { }
   
   template< typename TSpec >
   ArrivalProcess( const ArrivalProcess< TTime, TSpec >& proc ):
     object( new ArrivalProcessModel<TSpec>( proc ) ) { }
   
   template< typename TSpec >
   void reset( const ArrivalProcess< TTime, TSpec >& proc ) {
     object.reset( new ArrivalProcessModel<TSpec>( proc ) );
   }
   bool empty() const { return !object.get(); }
   operator bool() const { return !empty(); }
   void reset() { object.reset(); }
   
   TTime nextArrivalTime( TTime fromTime, TTime maxTime, double rateFactor,
                          TRand& randGen, double eps = 1e-5, int maxSteps = 100000 ) {
     return object->doNextArrival( fromTime, maxTime, rateFactor, randGen, eps, maxSteps );
   }
	 ostream& print( ostream& s ) const { return object->print( s ); }
   
   friend ostream& operator<<( ostream& s, const ArrivalProcess& f ) {
		 return f.print( s );
	 }
	 friend std::string getLabel( const ArrivalProcess& p ) { return p.doGetLabel(); }
}; // class ArrivalProcess<TTime, Any< TRand > >


//
// ** Poisson processes
//


template <typename TRate, typename TSpec> class Poisson;

struct CinlarsMethod;
template <typename TRateFunctionSpec, typename TMethodSpec = CinlarsMethod> class NonHomogeneous;

template <typename TTime, typename TRate, typename TRateFunctionSpec>
class ArrivalProcess<TTime, Poisson< TRate, NonHomogeneous< TRateFunctionSpec, CinlarsMethod > > > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TRateFunctionSpec> ));
public:
   typedef TTime result_type;
   typedef Function< TTime, TRate, TRateFunctionSpec > rate_function_type;
   typedef BOOST_TYPEOF_TPL(( integralFromPoint( boost::declval<rate_function_type>(),
                                                 boost::declval<TTime>() ) ))
	 rate_function_integral_type;

   typedef typename rate_function_integral_type::result_type integral_val_type;

   typedef typename DiffType< integral_val_type >::type
   integral_diff_type;

   //typedef typename DiffType<TTime> time_diff_type;
   
   ArrivalProcess( const rate_function_type& rateFunction_, TTime startTime_, std::string label_ = "" ):
     rateFunctionIntegral( integralFromPoint( rateFunction_, startTime_ ) ), label( label_ ) {

		 PRINT4( label_, rateFunction_, startTime_, rateFunctionIntegral );
   }

   // Method: function call operator
   // Generate next event from this process, starting from time 'fromTime' but no further than
   // time 'maxTime'.
   //
   // Template params:
   //
   //     TRand - type of the random numbers source
   //
   // Returns:
   //
   //   time of next event, or maxTime if no events happen before maxTime.
   //
   template <typename TFactor, typename TRand>
   TTime nextArrivalTime( TTime fromTime, TTime maxTime, TFactor rateFactor, TRand& randGen,
													double eps = 1e-5, int maxSteps = 100000  ) const {
		 assert( rateFactor >= 0 );
     assert( fromTime >= TTime(0.0) );
     assert( maxTime >= fromTime );

		 if ( rateFactor == 0  ||  fromTime == maxTime ) return maxTime;
     BOOST_AUTO_TPL(integralAtFromTime, ( eval( rateFunctionIntegral, fromTime )  ));
     BOOST_AUTO_TPL(integralAtMaxTime, ( eval( rateFunctionIntegral, maxTime ) ));
     assert( integralAtMaxTime >= integralAtFromTime );
		 if ( integralAtMaxTime == integralAtFromTime ) return maxTime;
#ifndef NDEBUG
		 PRINT7( label, fromTime, maxTime, rateFactor, eps, integralAtFromTime, integralAtMaxTime );
#ifdef COSI_DEV_PRINT_GENERALMATH
		 PRINT( rateFunctionIntegral );
#endif		 
     BOOST_TYPEOF_TPL( integralAtFromTime ) zeroIntegral(0.0);
     assert( integralAtFromTime >= zeroIntegral );
     assert( integralAtMaxTime >= zeroIntegral );
#endif
		 
     typename boost::random::uniform_01<> uniformDistr;
     BOOST_AUTO_TPL( u,-std::log( uniformDistr( randGen ) ) );
     BOOST_AUTO_TPL(targetIntegralValue, integralAtFromTime +
                    static_cast< integral_diff_type >( u / rateFactor ) );
     if ( integralAtMaxTime < targetIntegralValue ) return maxTime;
     TTime  targetTime =  evalInverse( rateFunctionIntegral, fromTime, maxTime, targetIntegralValue,
                                       static_cast<integral_diff_type>( eps / ( 10 * rateFactor ) ),
																			 maxSteps );

     // if ( !(std::abs( targetTime - targetIntegralValue ) < eps ) ) {
     //     PRINT2( std::abs( targetTime - targetIntegralValue ), eps );
     // }
     assert(( equal_eps( rateFunctionIntegral( targetTime ), targetIntegralValue ), eps ));
#ifndef NDEBUG
		 if ( !( equal_eps( ( rateFunctionIntegral( targetTime ) - rateFunctionIntegral( fromTime ) )
												* rateFactor, static_cast< integral_diff_type >(u),  eps ) ) ) {
			 PRINT12( fromTime, targetTime, targetTime - fromTime, u, rateFunctionIntegral( targetTime ), rateFunctionIntegral( fromTime ),
								targetIntegralValue, rateFunctionIntegral( targetTime ) - targetIntegralValue,
								rateFactor, rateFunctionIntegral( targetTime ) - rateFunctionIntegral( fromTime ),
								( rateFunctionIntegral( targetTime ) - rateFunctionIntegral( fromTime ) ) * rateFactor,
								ToDouble( ( rateFunctionIntegral( targetTime ) - rateFunctionIntegral( fromTime ) ) * rateFactor ) - u );
							
		 }
#endif		 
     assert(( equal_eps( ( rateFunctionIntegral( targetTime ) - rateFunctionIntegral( fromTime ) )
                         * rateFactor, static_cast< integral_diff_type >(u), eps ) ));
     return targetTime;
//     return std::min( maxTime, findIntegralInverse( rateFunctionIntegral( fromTime ) ;
   }
   const rate_function_integral_type& getRateFunctionIntegral() const { return rateFunctionIntegral; }
   
   friend ostream& operator<<( ostream& s, const ArrivalProcess& f ) {
		 s << "ArrivalProcess[label=" << f.label << ", rateFunctionIntegral=" << f.rateFunctionIntegral << "]";
		 return s;
	 }

	 friend std::string getLabel( const ArrivalProcess& p ) { return p.label; }

private:
   rate_function_integral_type rateFunctionIntegral;
	 std::string label;
   //integral_val_type integralAtZeroTime;
};  // class ArrivalProcess<TTime, Poisson< TRate, NonHomogeneous< TRateFunctionSpec > > >

template <typename TTime, typename TRate, typename TRateFunctionSpec >
struct result_of_makeNonHomogeneousPoissonProcess {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TRateFunctionSpec> ));
   typedef ArrivalProcess<TTime, Poisson< TRate, NonHomogeneous< TRateFunctionSpec > > > type;
};

template <typename TTime, typename TRate, typename TRateFunctionSpec >
typename result_of_makeNonHomogeneousPoissonProcess<TTime, TRate, TRateFunctionSpec>::type
makeNonHomogeneousPoissonProcess( const Function< TTime, TRate, TRateFunctionSpec >& rateFunction,
                                  TTime fromTime, std::string label ) {
  typedef typename result_of_makeNonHomogeneousPoissonProcess<TTime, TRate, TRateFunctionSpec>::type
     result_type;
  return result_type( rateFunction, fromTime, label );
}


// * Utils
// ** Numerical integration

template <typename TDomain, typename TRange, typename TSpec>
double FuncCaller( double x, void *data ) {
  const Function<TDomain,TRange,TSpec> *f =  (const Function<TDomain,TRange,TSpec> *)data;
  assert(data);
  return ToDouble( eval( *f, static_cast<TDomain>( x ) ) );
}


template <typename TDomain, typename TRange, typename TSpec>
typename AreaType<TDomain,TRange>::type
inline
integrateNumerically( const Function<TDomain,TRange,TSpec>& f,
                      TDomain a, TDomain b, int n ) {
  Function<TDomain,TRange,TSpec>& f_nonconst = const_cast<Function<TDomain,TRange,TSpec>&>( f );
  typedef double (*f_t)(double, void *);
  f_t fptr = &FuncCaller<TDomain, TRange, TSpec>;
  typedef typename AreaType<TDomain,TRange>::type result_type;
  return static_cast<result_type>( gauss_legendre( n, fptr, &f_nonconst, ToDouble( a ), ToDouble( b ) ) );
  
}

// ** Interpolated functions

// *** Class: InterpFn
//
// An interpolated function: stores a vector of (x,f(x)) pairs and
// can compute the values of f at new points by simple linear interpolation.
//
template <typename TDomain = double, typename TRange = double>
class InterpFn /*: public AbstractFunction<TDomain, TRange>*/ {
public:
   // Method: clear
   // Remove all (x,f(x)) pairs from this function.
   void clear() { x2y.clear(); }

   // Method: addPt
   // Store an (x,f(x)) pair specifying the value of the function at a given point 'x'.
   // Points must be added in order of increasing 'x'.
   void addPt( TDomain x, TRange y ) {
     assert( x2y.empty() || x > boost::prior( x2y.end() )->first );
     x2y.insert( x2y.end(), make_pair( x, y ) );
   }

   // Method: eval
   // Evaluate the function at the specified point 'x', interpolating as necessary
   // between the neighboring values.
   TRange operator()( TDomain x ) const {
     assert( !x2y.empty() );
     assert( x >= x2y.begin()->first );
     assert( x <= x2y.rbegin()->first );
		 BOOST_AUTO_TPL(it, x2y.lower_bound( x ) );
     assert( it != x2y.end() );
     if ( it->first == x ) return it->second;
     assert( it != x2y.begin() );
     BOOST_AUTO_TPL( it_p, boost::prior( it ) );
     double frac = ( ToDouble( x ) - ToDouble( it_p->first ) ) / ( ToDouble( it->first ) - ToDouble( it_p->first ) );
     double result = ( ToDouble( it_p->second ) * (1.0-frac) + ToDouble( it->second ) * frac );
     return TRange( result );
   }

   virtual TRange eval( TDomain x ) const { return (*this)( x ); }

   size_t getNumPts() const { return x2y.size(); }
   
#if 0
   struct Pt {
      TDomain x;
      TRange y;
      typename map<TDomain,TRange>::const_iterator lwrBnd, uprBnd;
   };

   Pt evalPt( TDomain x ) const {
		 BOOST_AUTO_TPL( it, x2y.lower_bound( x ) );
     assert( it != x2y.end() );
     Pt p;
     p.x = x;
     p.lwrBnd = it;
     if ( it->first == x ) {
       p.uprBnd = it;
       p.y = it->second;
     } else {
       assert( it != x2y.begin() );
       BOOST_AUTO_TPL( it_p, boost::prior( it ) );
       p.uprBnd = it_p;
       double frac = ( x - it_p->first ) / ( it->first - it_p->first );
       double result = ( it_p->second * (1.0-frac) + it->second * frac );
       p.y = result;
     }
     
     return p;
   }
   
   Pt evalPt( Pt from, Pt to, TDomain x ) const {
		 BOOST_AUTO_TPL( it, lower_bound( from.uprBnd, to.lwrBnd, x ) );

     assert( it != x2y.end() );
     Pt p;
     p.x = x;
     p.lwrBnd = it;
     if ( it->first == x ) {
       p.uprBnd = it;
       p.y = it->second;
     } else {
       assert( it != x2y.begin() );
       BOOST_AUTO_TPL( it_p, boost::prior( it ) );
       p.uprBnd = it_p;
       double frac = ( x - it_p->first ) / ( it->first - it_p->first );
       double result = ( it_p->second * (1.0-frac) + it->second * frac );
       p.y = result;
     }
     
     return p;
   }
#endif

   void save( filename_t filename ) const {
     boost::filesystem::ofstream out( filename );
     out.exceptions( std::ios_base::failbit | std::ios_base::badbit );
     out.precision( 20 );
     out << "x\ty\n";
     for ( x2y_const_iterator it = x2y.begin(); it != x2y.end(); it++ )
        out << it->first << "\t" << it->second << "\n";
   }
   
   void load( filename_t filename ) {
     boost::filesystem::ifstream inFile( filename );
     inFile.exceptions( std::ios_base::failbit | std::ios_base::badbit );
     clear();
     size_t lineNum = 0;
     while ( inFile ) {
       std::string line;
       std::getline( inFile, line );
       if ( lineNum == 0 ) continue;
       std::istringstream lineStr( line );
       lineStr.exceptions( std::ios_base::failbit | std::ios_base::badbit );
       TDomain x;
       TRange y;
       lineStr >> x >> y;
       addPt( x, y );
       lineNum++;
     }
   }

private:

   typedef std::map<TDomain,TRange> x2y_type;
   typedef typename x2y_type::const_iterator x2y_const_iterator;
   
   // Field: x2y
   // Map from point x to the value f(x), for the points at which the function value
   // is explicitly specified.
   x2y_type x2y;

	 friend class boost::serialization::access;
   template <class Archive> void serialize( Archive& ar, const unsigned int /* version */ ) {
     ar & BOOST_SERIALIZATION_NVP( x2y );
   }
   
};  // class InterpFn


//
// *** Class: InterpBiFun
//
// An bijective interpolated function, which can map values both from
// domain to range and from range to domain.
//
template <typename TDomain = double, typename TRange = double >
class InterpBiFun: public InterpFn<TDomain, TRange> {
public:
   typedef InterpFn<TDomain,TRange> PARENT;
   
   void clear() { PARENT::clear(); inverse.clear(); }
   void addPt( TDomain x, TRange y ) { PARENT::addPt( x, y ); inverse.addPt( y, x ); }

   const InterpFn<TRange,TDomain>& getInverse() const { return inverse; }
   
private:
   InterpFn<TRange,TDomain> inverse;
};

template <typename TRange>
inline TRange interpolate( const std::vector<TRange>& f, double i ) {
  int int_i = static_cast<int>( floor( i ) );
  double frac_i = i - int_i;
  return static_cast<TRange>( ( fabs( int_i - i ) < 1e-5 ) ? f[ int_i ] :
                              (1-frac_i) * f[int_i] + frac_i*f[int_i+1] ); 
}

// ** Numerical utils: integrate, findWhere

template <typename TDomain, typename TRange> inline
typename AreaType<TDomain,TRange>::type integrate( const std::vector<TRange>& f, const std::vector<TDomain>& x, double from, double to ) {

  typedef typename AreaType<TDomain,TRange>::type result_type;

  assert( 0 <= from && from <= f.size()-1 );
  assert( 0 <= to && to <= f.size()-1 );

  TRange firstF = interpolate( f, from );
  TDomain firstX = interpolate( x, from );
  TRange lastF = interpolate( f, to );
  TDomain lastX = interpolate( x, to );

  if ( static_cast<int>(from) == static_cast<int>(to) )
     // The integral doesn't cross any of the f samples,
     // handle this case separately
     return static_cast<result_type>( 0.5*(firstF+lastF)*(lastX-firstX) );

  int ceilFrom = static_cast<int>( ceil( from ) );
  int floorTo = static_cast<int>( floor( to ) );

  result_type I = 0.5*(f[ceilFrom]+firstF)*(x[ceilFrom]-firstX);
  for ( int i = ceilFrom; i < floorTo; i++ ) {
    result_type incr = 0.5*(f[i]+f[i+1])*(x[i+1]-x[i]);
    I += incr;
  }

  result_type finalAdd = 0.5*(f[floorTo]+lastF)*(lastX-x[floorTo]);
  I += finalAdd;
  return I;
}

template <typename T> inline
T findWhere( const std::vector<T>& f, T c ) {

  for ( int i = 0; i < static_cast<int>( f.size() )-1; i++ ) {
    if( (f[i] <= c && f[i+1] >= c)
        || (f[i] >= c && f[i+1] <= c) ) {
      // Interpolate between i and i+1
      T frac = fabs( (f[i] - c)/(f[i+1]-f[i]) );
      return (1-frac)*i + frac*(i+1);
    }
  }
  // f never crosses c
  return -1;
}

// * Postamble

}  // namespace math
}  // namespace cosi

//////////////////////////////////////////////////////////////////////
// Boost.TypeOf registration
//////////////////////////////////////////////////////////////////////


//BOOST_TYPE_ERASURE_FREE((has_eval), cosi::math::eval, 2)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::MultType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::DiffType,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::DomainType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::RangeType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::SpecType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::IsFunctionSpec,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::FunctionConcept,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Function,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Const,1)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::RunTime)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::CompileTime,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::X_To,(int)(typename))

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::BinOp,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::AddOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::SubOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::MultOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::DivOp,2)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Piecewise,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::AreaType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::detail::result_of_make_line_through,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::UnaryOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_indefiniteIntegral,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_integralFromPoint,3)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::Any)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::ArrivalProcess,2);
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Poisson,2)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::CinlarsMethod)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::NonHomogeneous,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_makeNonHomogeneousPoissonProcess,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::InterpFn,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::InterpBiFun,2)

#endif  // #ifndef __INCLUDE_COSI_GENERALMATH_H
// Postamble:1 ends here
