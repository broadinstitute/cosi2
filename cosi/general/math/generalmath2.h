
// * Preamble

// ** Headers

#ifndef __INCLUDE_COSI_GENERALMATH2_H
#define __INCLUDE_COSI_GENERALMATH2_H

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
#include <boost/mpl/int.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <cosi/gauss_legendre.h>
#include <cosi/utils.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()


namespace cosi {
namespace math2 {

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
template <typename TDom, typename TRan>
struct AreaType {
   typedef typename DiffType<TDom>::type domain_diff_type;
   typedef typename DiffType<TRan>::type range_diff_type;
   typedef typename MultType<domain_diff_type, range_diff_type>::type type;
};



// * Functions
// ** Generic function concept

// *** Concept <<FuncConcept>>
//
//    Concept checking class for general math functions of one variable.
//
//    Modeled by: [[Func]] (defined below).
template <typename F>
class FuncConcept {
public:
   //
   // Associated types: for each generic math function there is a domain and a range type.
   // Normally these would be floating-point types.
   //
   
   BOOST_CONCEPT_USAGE(FuncConcept) {
		 double x;
     (void)eval( const_cast<const F&>( func ), x );
   }

private:
   F func;
   
};  // FuncConcept

// *** GenericClass Func - A mathematical function.
template <typename TSpec> class Func;

        
// Template params:

//   - TSpec - the specialization defining what kind of function this is.
         
// *** Metafunction: SpecType
//    Returns the specialization type of a Func; this defines what kind of function it is
//    (constant, linear etc).
template <typename T> struct SpecType;

// Metapredicate: IsFuncSpec
// Tests whether a given type is a valid function specialization.
template <typename T, typename Enabler = void> struct IsFuncSpec: public boost::false_type { };

struct SpecBase {};

template <typename TFunc> struct DomainType { typedef typename TFunc::argument_type type; };
template <typename TFunc> struct RangeType { typedef typename TFunc::result_type type; };
template <typename TFunc> struct SpecType { typedef typename TFunc::spec_type type; };

//
// *** Generic function: eval
//
// Default implementation of eval(f,x) to evaluate a Func at a given point in its domain;
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
template <typename TRan, typename TConstSpec = RunTime> struct Const;
template <typename TRan> struct IsFuncSpec< Const<TRan, RunTime> >: public boost::true_type {};

template <int N> struct CompileTime { static const int value = N; };

template <typename TRan, int N> struct IsFuncSpec< Const< TRan, CompileTime<N> > >: public boost::true_type {};

template <typename TDom, typename TRan, typename TSpec>
TRan evalConst( const Func< Const< TRan, TSpec > >& f ) { return eval( f, TDom() ); }

template <typename TRan>
class Func<Const< TRan, RunTime > > {
public:
   typedef Const<TRan, RunTime> spec_type;

   Func( TRan val_ ): val( val_ ) { }

	 template <typename TDom>
   TRan operator()( TDom ) const { return val; }

   friend ostream& operator<<( ostream& s, const Func& f ) {
     s << "ConstRuntime[" << f.val << "]"; return s;
   }

private:
	 TRan val;
};

template <typename TRan>
Func< Const<TRan> > makeConst( TRan x ) { return Func< Const<TRan> >( x ); }

template <typename TRan, int val>
class Func<Const< TRan, CompileTime<val> > > {
public:
   typedef Const< TRan, CompileTime<val> > spec_type;
	 template <typename TDom>
   TRan operator()( TDom ) const { return static_cast<TRan>( val ); }

   friend ostream& operator<<( ostream& s, const Func& f ) {
     s << "ConstCompiletime[" << val << "]"; return s;
   }
   
};

// *** Monomials (functions of the form a x^k) 

template <int exponent, typename FactorType = double> struct X_To;
template <int exponent, typename FactorType>
   struct IsFuncSpec< X_To<exponent, FactorType> > : public boost::true_type {};

template <typename FactorType>
class Func<X_To<1, FactorType> > {
public:
   typedef X_To<1,FactorType> spec_type;
   typedef FactorType factor_type;

   Func( ): factor( 1.0 ) {}
   Func( factor_type factor_ ): factor( factor_ ) {}
	 template <typename TDom>
	 BOOST_TYPEOF_TPL( boost::declval<FactorType>() * ( boost::declval<TDom>() - boost::declval<TDom>() )  )
    operator()( TDom x ) const { return factor * ( x - static_cast< TDom >( 0.0 ) ) ; }
   factor_type getFactor() const { return factor; }

   friend ostream& operator<<( ostream& s, const Func& f ) {
     s << "(" << f.factor << " * x" << ")"; return s; 
   }

   
private:
   factor_type factor;
};

template <typename TFactor>
Func< X_To<1,TFactor> > makeX( TFactor factor ) { return Func< X_To<1,TFactor> >( factor ); }

template <typename FactorType>
class Func<X_To<2,FactorType> > {
public:
   typedef X_To<2,FactorType> spec_type;
   typedef FactorType factor_type;
   Func( factor_type factor_ ): factor( factor_ ) {}
	 template <typename TDom>
   BOOST_TYPEOF_TPL( boost::declval<FactorType>() * ( boost::declval<TDom>() - boost::declval<TDom>() ) *  ( boost::declval<TDom>() - boost::declval<TDom>() )  )
		 operator()( TDom x ) const { return factor * x*x; }
   factor_type getFactor() const { return factor; }
   friend ostream& operator<<( ostream& s, const Func& f ) {
     s << "(" << f.factor << " * x^2" << ")"; return s; 
   }
private:
   const factor_type factor;
};

template <typename TFactor>
Func< X_To<2,TFactor> > makeX2( TFactor factor ) { return Func< X_To<2,TFactor> >( factor ); }

// ** Operations on functions

//   Construction of compound functions: e.g. given two functions we can
//   add them to create a new function that computes the sum of the two original
//   functions' results.

// *** Ops

struct AddOp {
	 template <typename T1, typename T2>
   BOOST_TYPEOF_TPL( boost::declval<T1>() + boost::declval<T2>() )
		 operator()( T1 a1, T2 a2 ) const { return a1 + a2; }
   friend ostream& operator<<( ostream& s, const AddOp& ) { s << "+"; return s; }
};
struct SubOp {
	 template <typename T1, typename T2>
   BOOST_TYPEOF_TPL( boost::declval<T1>() - boost::declval<T2>() )
		 operator()( T1 a1, T2 a2 ) const { return a1 - a2; }
   friend ostream& operator<<( ostream& s, const SubOp& ) { s << "-"; return s; }
};
struct MultOp {
	 template <typename T1, typename T2>
   BOOST_TYPEOF_TPL( boost::declval<T1>() * boost::declval<T2>() ) operator()( T1 a1, T2 a2 ) const { return a1 * a2; }
   friend ostream& operator<<( ostream& s, const MultOp& ) { s << "*"; return s; }
};
struct DivOp {
	 template <typename T1, typename T2>
   BOOST_TYPEOF_TPL( boost::declval<T1>() / boost::declval<T2>() ) operator()( T1 a1, T2 a2 ) const { return a1 / a2; }
   friend ostream& operator<<( ostream& s, const DivOp& ) { s << "/"; return s; }
};

// *** Binary Op Func

template <typename TSpec1, typename TSpec2, typename Op> struct BinOp;
template <typename TSpec1, typename TSpec2, typename Op>
struct IsFuncSpec< BinOp< TSpec1, TSpec2, Op >,
                       typename boost::enable_if< typename boost::mpl::and_< IsFuncSpec< TSpec1 >,
                                                                             IsFuncSpec< TSpec2 > > >::type >: public boost::true_type {};



template <typename TSpec1, typename TSpec2, typename Op>
class Func<BinOp<TSpec1, TSpec2, Op > >:
		 public boost::compressed_pair< Func<TSpec1>, Func<TSpec2> > {
   BOOST_MPL_ASSERT(( IsFuncSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFuncSpec<TSpec2> ));
public:
   typedef BinOp<TSpec1,TSpec2,Op > spec_type;
   typedef TSpec1 spec1_type;
   typedef TSpec2 spec2_type;
   
   typedef Func<TSpec1> function_type_1; 
   typedef Func<TSpec2> function_type_2;

   Func( const function_type_1& f1_, const function_type_2& f2_ ):
     boost::compressed_pair< function_type_1, function_type_2 >( f1_, f2_ ) { }

	 template <typename TDom>
	 BOOST_TYPEOF_TPL( binop( eval( boost::declval<function_type_1>(), boost::declval<TDom>() ),
														eval( boost::declval<function_type_2>(), boost::declval<TDom>() ) ) )
	 operator()( TDom x ) const { return binop( eval( this->first(), x ),
																								 eval( this->second(), x ) ); }

   friend ostream& operator<<( ostream& s, const Func& f ) {
     s << "(" << f.first() << " " << f.binop << f.second() << ")"; return s;
   }
   
private:
   static const Op binop;
   
};

template <typename TSpec1, typename TSpec2, typename Op>
const Op Func< BinOp< TSpec1, TSpec2, Op > >::binop = Op();

// *** Operator overloads for adding/subtracting/etc two functions
// 
// Adding two functions f and g creates a new function that for x computes f(x)+g(x)


template <typename TSpec1, typename TSpec2, typename Op>
Func< BinOp<TSpec1,TSpec2, Op > >
ApplyOp( const Func<TSpec1>& f1, const Func<TSpec2>& f2, Op ) {
  return Func< BinOp<TSpec1,TSpec2,Op > >( f1, f2 );
}

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( ApplyOp( boost::declval<TSpec1>(), boost::declval<TSpec2>(), boost::declval<AddOp>() ) )
operator+( const Func<TSpec1>& f1, const Func<TSpec2>& f2 ) { return ApplyOp( f1, f2, AddOp() ); }

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( ApplyOp( boost::declval<TSpec1>(), boost::declval<TSpec2>(), boost::declval<SubOp>() ) )
operator-( const Func<TSpec1>& f1, const Func<TSpec2>& f2 ) { return ApplyOp( f1, f2, SubOp() ); }

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( ApplyOp( boost::declval<TSpec1>(), boost::declval<TSpec2>(), boost::declval<MultOp>() ) )
operator*( const Func<TSpec1>& f1, const Func<TSpec2>& f2 ) { return ApplyOp( f1, f2, MultOp() ); }

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( ApplyOp( boost::declval<TSpec1>(), boost::declval<TSpec2>(), boost::declval<DivOp>() ) )
operator/( const Func<TSpec1>& f1, const Func<TSpec2>& f2 ) { return ApplyOp( f1, f2, DivOp() ); }


// *** Operations on const funcs
//
// Operations on constant functions: produce a new constant function, whose return value
// can be computed immediately.
//

template <typename TRan1, typename TRan2, typename Op>
Func< Const< BOOST_TYPEOF_TPL( Op( boost::declval<TRan1>(), boost::declval<TRan2>() ) )  > >
ApplyOp( const Func< Const< TRan1 > >& f1, const Func< Const< TRan2 > >& f2 ) {
  return
		 Func< Const< BOOST_TYPEOF_TPL( Op( boost::declval<TRan1>(),
																						boost::declval<TRan2>() ) ) > >( Op( evalConst( f1 ),
																																									 evalConst( f2 ) ) );
}

//
// *** Unary operations on functions: specificaly, negation.
//
//  For a function f, -f is a function that for value x returns -f(x).
//

template <typename TSpec, typename Op> struct UnaryOp;

template <typename TSpec, typename Op>
class Func<UnaryOp<TSpec, Op> > {
   BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
public:
   typedef UnaryOp<TSpec,Op> spec_type;
   typedef TSpec unop_spec_type;
   typedef Op op_type;
   typedef Func<TSpec> function_type;


   Func( const function_type& f_ ): f( f_ ) { }

	 template <typename TDom>
	 BOOST_TYPEOF_TPL( op( eval( boost::declval<function_type>(), boost::declval<TDom>() ) ) )
   operator()( TDom x ) const { return op( eval( f, x ) ); }

   friend ostream& operator<<( ostream& s, const Func& f ) { s << op << " " << f.f; return s; }
   
private:
   static const Op op;
   function_type f;
};

template <typename TSpec, typename Op>
  const Op Func<UnaryOp< TSpec, Op > >::op = Op();


template <typename TSpec, typename Op>
Func< UnaryOp<TSpec, Op > >
ApplyOp( const Func<TSpec>& f, Op ) {
  return Func< UnaryOp<TSpec,Op > >( f );
}

struct NegateOp {
	 template <typename T> BOOST_TYPEOF_TPL( -boost::declval<T>() ) operator()( T a ) const { return -a; }
   friend ostream& operator<<( ostream& s, const NegateOp& ) { s << "-"; return s; }
};


template <typename TSpec>
BOOST_TYPEOF_TPL( ApplyOf( boost::declval< Func<TSpec> >(), boost::declval<NegateOp>() ) )
operator-( const Func<TSpec>& f ) {
  return ApplyOp( f, NegateOp() );
}

// ** Linear functions: creating a linear function from two points through which it passes.



// *** Construct line through two points

template <typename TDom, typename TRan>
BOOST_TYPEOF_TPL(( boost::declval< Func<X_To<1> > >() +
									 boost::declval< Func<Const<TDom> > >() ))
make_line_through( TDom x1, TRan y1, TDom x2, TRan y2 ) {
  assert( x1 != x2 );
  BOOST_AUTO_TPL( slope, ( y2 - y1 ) / ( x2 - x1 ) );
  BOOST_AUTO_TPL( b, y1 - ( slope * x1 ) );
  return Func<X_To<1> >( slope ) + makeConst( b );
}


// ** Piecewise functions

//       A piecewise function is a compound function which divides the domain into a set of intervals,
//       and on each interval returns the result of a specified piece function.
//       Operations on piecewise functions can often be expressed in terms of operations on the pieces.


template <typename TDom, typename TPieceSpec> struct Piecewise;
template <typename TDom, typename TPieceSpec>
struct IsFuncSpec< Piecewise< TDom, TPieceSpec >,
									 typename boost::enable_if< IsFuncSpec< TPieceSpec > > >: public boost::true_type {};

template <typename TDom, typename TPieceSpec>
class Func<Piecewise< TDom, TPieceSpec > > {
//	 BOOST_MPL_ASSERT(( IsFuncSpec<TPieceSpec> ));
public:
	 typedef Piecewise< TDom, TPieceSpec > spec_type;
	 typedef TPieceSpec piece_spec_type;
   
	 typedef Func< TPieceSpec > piece_function_type;
	 BOOST_CONCEPT_ASSERT((FuncConcept< piece_function_type >));
   
	 Func() { }
	 template <typename SinglePassRange>
	 Func( const SinglePassRange& r ):
		 pieces( boost::begin( r ), boost::end( r ) ) {
		 BOOST_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< SinglePassRange > ));
		 typedef typename boost::range_value< SinglePassRange >::type range_value_type;
		 BOOST_STATIC_ASSERT(( boost::is_convertible< range_value_type,
													 typename pieces_type::value_type >::value ));
	 }
   
	 // template <typename TDom2>
	 // void addPiece( TDom2 x,
	 //                 const piece_function_type& piece,
	 //                 typename boost::enable_if< boost::is_convertible< TDom2, TDom >, TDom2>::type *dummy = NULL ) { (void)dummy; pieces.insert( std::make_pair( x, piece ) ); }
   
	 // Type: pieces_type
	 // The type for a map from point to the piece function starting at that point.
	 // The map is ordered from the rightmost piece to the leftmost, so that map::lower_bound()
	 // can be used to find the piece corresponding to a given domain point.
	 typedef std::map< TDom, piece_function_type, std::greater<TDom> > pieces_type;
   
	 // Method: getPieces
	 // Returns a map from domain point to the piece function starting at that point.
	 // The map is ordered from rightmost piece to the left.
	 // The map is mutable; change it to change the definition of this piecewise function.
	 const pieces_type& getPieces() const { return pieces; }
	 pieces_type& getPieces() { return pieces; }
	 
	 typedef BOOST_TYPEOF_TPL(( eval( boost::declval< piece_function_type >(), boost::declval< TDom >() ) ))
	 result_type;
	 
   result_type operator()( TDom x ) const {
		 BOOST_CONCEPT_ASSERT((boost::LessThanComparable<TDom>));
		 
		 typename pieces_type::const_iterator it = pieces.lower_bound( x );
		 return it == pieces.end() ? std::numeric_limits<result_type>::quiet_NaN() :
				eval( it->second, x );
	 }
   
	 friend ostream& operator<<( ostream& s, const Func& f ) {
		 s << "Piecewise[";
		 for ( BOOST_AUTO_TPL( it, f.pieces.rbegin() ); it != f.pieces.rend(); it++ )
				s << "(" << it->first << ", " << it->second << ")";
		 s << "]";
		 return s;
	 }
   
   
private:
	 // Field: pieces
	 // Map from domain point to the piece function starting at that point; ordered from rightmost point to leftmost.
	 pieces_type pieces;
};  // class Func<TDom, TRan, Piecewise< TPieceSpec > >

#if 1

// *** Operations on F,P where P is a piecewise function
// Propagate to the pieces 

template <typename TDom, typename TSpec1, typename TPieceSpec2, typename Op>
struct result_of_ApplyOp {
  typedef BOOST_TYPEOF_TPL(
		ApplyOp( boost::declval< Func<TSpec1> >(),
						 boost::declval< Func< Piecewise< TDom, TPieceSpec2 > > >(),
						 boost::declval<Op>() ) ) type;
};

template <typename TDom, typename TSpec1, typename TPieceSpec2, typename Op>
typename result_of_ApplyOp<TDom,TSpec1,TPieceSpec2,Op>::type
ApplyOp( const Func<TSpec1>& f1, const Func<Piecewise<TDom,TPieceSpec2> >& f2, Op ) {
	typename result_of_ApplyOp<TDom,TSpec1,TPieceSpec2,Op>::type f_result;
	for ( BOOST_AUTO_TPL( it, f2.getPieces().begin() ); it != f2.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, Op( it->second ) ) );
	return f_result;
}

// ** Integrals


// *** Definition of indefiniteIntegral function

template <typename TSpec> struct result_of_indefiniteIntegral {
	 BOOST_STATIC_ASSERT_MSG( sizeof( TSpec ) + 1 == 0, "don't know how to get type of indef integral of f" );
};

template <typename TSpec>
void indefiniteIntegral( const Func<TSpec >& f) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  BOOST_STATIC_ASSERT_MSG( sizeof( TSpec ) + 1 == 0, "don't know how to integrate f" );
}

//
// *** Definite integrals
//



template <typename TDom, typename TSpec>
typename AreaType< TDom, BOOST_TYPEOF_TPL( eval( boost::declval< Func<TSpec> >(),
																								 boost::declval< TDom > ) ) >::type
definiteIntegral( const Func< TSpec >& f, TDom a, TDom b ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  BOOST_AUTO( f_int, indefiniteIntegral( f ) );
  return eval( f_int, b ) - eval( f_int, a );
}

template <typename TDom, typename TSpec>
BOOST_TYPEOF_TPL( indefiniteIntegral( boost::declval< Func<TSpec> >() ) -
									makeConst( eval( indefiniteIntegral( boost::declval< Func<TSpec> >() ) ),
														 boost::declval<TDom>() ) )
integralFromPoint( const Func< TSpec >& f, TDom x ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  BOOST_AUTO( f_int, indefiniteIntegral( f ) );
  return f_int - makeConst( eval( f_int, x ) );
}

//
// *** Integral of a constant
//

template <typename TRan>
BOOST_TYPEOF_TPL( makeX( evalConst( boost::declval< Func< Const< TRan > > >() ) - boost::declval<TRan>() ) )
indefiniteIntegral( const Func< Const< TRan > >& f ) {
  return makeX( evalConst( f ) - static_cast<TRan>( 0 ) );
}

//
// *** Integral of a monomial
//

template <typename TFactor>
BOOST_TYPEOF_TPL( makeX2( 0.5 * boost::declval<TFactor>() ) )
indefiniteIntegral( const Func< X_To<1,TFactor> >& f ) {
  return makeX2( 0.5 * f.getFactor() );
}


//
// *** Linearity of integration
//

template <typename TRan, typename TSpec>
BOOST_TYPEOF_TPL( boost::declval< Func< Const<TRan> > >() * indefiniteIntegral( boost::declval< Func< TSpec > >() ) )
indefiniteIntegral( const Func< BinOp< Const<TRan>, TSpec, MultOp > >& f ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  return f.first() * indefiniteIntegral( f.second() );
}

template <typename TRan, typename TSpec>
BOOST_TYPEOF_TPL( boost::declval< Func< Const<TRan> > >() * indefiniteIntegral( boost::declval< Func< TSpec > >() ) )
indefiniteIntegral( const Func< BinOp< TSpec, Const<TRan>, MultOp > >& f ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  return f.second() * indefiniteIntegral( f.first() );
}

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( indefiniteIntegral( boost::declval< Func<TSpec1> >() ) +
									indefiniteIntegral( boost::declval< Func<TSpec2> >() ) )
indefiniteIntegral( const Func< BinOp< TSpec1, TSpec2, AddOp > >& f ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec2> ));
  return indefiniteIntegral( f.first() ) + indefiniteIntegral( f.second() );
}

template <typename TSpec1, typename TSpec2>
BOOST_TYPEOF_TPL( indefiniteIntegral( boost::declval< Func<TSpec1> >() ) -
									indefiniteIntegral( boost::declval< Func<TSpec2> >() ) )
indefiniteIntegral( const Func< BinOp< TSpec1, TSpec2, SubOp > >& f ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec2> ));
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

template <typename TDom, typename TPieceSpec>
Func< Piecewise< TDom, typename SpecType< BOOST_TYPEOF_TPL( indefiniteIntegral( boost::declval< Func< TPieceSpec > >() ) ) >::type > >
indefiniteIntegral( const Func< Piecewise< TDom, TPieceSpec > >& f ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TPieceSpec> ));
  
  namespace ad = boost::adaptors;
  typedef BOOST_TYPEOF_TPL( indefiniteIntegral( boost::declval< Func< TPieceSpec > >() ) )
     piece_integral_type;
	typedef Func< Piecewise< TDom, typename piece_integral_type::spec_type > > result_type;
  typedef BOOST_TYPEOF_TPL( eval( boost::declval< piece_integral_type >(), boost::declval< TDom >() ) ) result_range_type;

  result_type result;
  result_range_type definiteIntegralSoFar( 0.0 );

  for ( BOOST_AUTO_TPL( it, f.getPieces().rbegin() ); it != f.getPieces().rend(); it++ ) {
    BOOST_AUTO_TPL( pieceIndefIntegral, indefiniteIntegral( it->second ) );
		BOOST_AUTO_TPL( offsetHere, definiteIntegralSoFar - eval( pieceIndefIntegral,
																															it->first ) );
    result.getPieces().insert( 
			std::make_pair( it->first,
											makeConst( offsetHere )
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
  //return result_type( f.getPieces() | ad::transformed( indefiniteIntegral_pair_caller<TDom,TRan,TPieceSpec>() ) );
}


// ** Evaluating inverse functions

//
// *** Func: evalInverse
//
// Evaluate the inverse of a given function at a given point.
//
// Params:
//
//    f - the function.   It must be strictly increasing on the interval [a,b]
//
template <typename TDom, typename TRan, typename TSpec>
TDom
evalInverse( const Func<TSpec>& f, TDom a, TDom b, TRan targetVal,
             typename DiffType<TRan>::type eps, unsigned maxSteps = 10000 ) {
  BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
  unsigned nsteps = 0;
  //PRINT5( a, b, targetVal, eps, typeid(f).name() );
  assert( !(boost::math::isnan)( a ) );
  assert( !(boost::math::isnan)( b ) );
  assert( !(boost::math::isnan)( targetVal ) );
  while ( true ) {
    assert( a <= b );
    TDom mid = a + ( b-a ) / 2;
    typedef typename DiffType<TRan>::type range_diff_type;
    range_diff_type curDiff = eval( f, mid ) - targetVal;
    //PRINT5( a, b, mid, targetVal, curDiff );
    if ( cosi_fabs( curDiff ) < eps ) {
      return mid;
    }
    if ( nsteps++ > maxSteps ) throw std::runtime_error( "findInfimum: too many steps" );
    assert( a < b );
    assert( eval( f, a ) < eval( f, mid ) );
    assert( eval( f, mid ) < eval( f, b ) );
    ( curDiff < static_cast<range_diff_type>( 0.0 ) ? a : b ) = mid;
  }
}

// template <typename TDom, typename TRan, typename TPieceSpec>
// indefiniteIntegral< const Func< TDom, TRan, Piecewise< TPieceSpec > >& f ) {
  
// }

// template <typename TDom, typename TRan, typename TSpec>
// typename IntegrationResultType<TDom,TRan>::type
// integrate( const TFunction<TDom, TRan, Const<TSpec>& f, TDom a, TDom b ) {
//  return ( b - a ) * f.getConstValue();
// }
// Evaluating\ inverse\ functions:1 ends here

// ** Type erasure for Func: a Func that stores any Func.

template <typename TDom, typename TRan, typename TSpec = void> struct AnyFunc;
template <typename TDom, typename TRan> struct IsFuncSpec< AnyFunc< TDom, TRan > >: public boost::true_type {};

template <typename TDom, typename TRan>
class Func< AnyFunc< TDom, TRan > > {
public:
   typedef TDom argument_type;
   typedef TRan result_type;
   typedef AnyFunc< TDom, TRan > spec_type;
   
private:
   struct FunctionObjectConcept {
      virtual ~FunctionObjectConcept() {}
      virtual TRan doEval( TDom ) const = 0;
   };

   template< typename TSpec > class FunctionObjectModel : public FunctionObjectConcept {
      BOOST_MPL_ASSERT(( IsFuncSpec<TSpec> ));
   public:
      FunctionObjectModel( const Func<TSpec>& f_ ) : f( f_ ) {}
      virtual ~FunctionObjectModel() {}
      virtual TRan doEval( TDom x ) const { return eval( f, x ); }
   private:
      Func<TSpec> f;
   };

   boost::shared_ptr<FunctionObjectConcept> object;

public:
   Func()  { }
   template< typename TSpec > Func( const Func<TSpec>& obj ) :
     object( new FunctionObjectModel<TSpec>( obj ) ) {}

   template <typename TSpec>
   void reset( const Func<TSpec>& obj ) {
     object.reset( new FunctionObjectModel<TSpec>( obj ) );
   }

   bool empty() const { return !object.get(); }
   operator bool() const { return !empty(); }

   TRan operator()( TDom x ) const { assert( object.get() ); return object->doEval( x ); }
};  // end: type erasure of a Func


// * Arrival processes

//     Modeling of generic arrival processes, such as Poisson processes.


template <typename TTime, typename TSpec> class ArrivalProcess;

//
// ** Type erasure for ArrivalProcess: storing any ArrivalProcess
//
template <typename TSpec = void> struct AnyProc;

template <typename TTime, typename TRand>
class ArrivalProcess<TTime, AnyProc< TRand > > {
   struct ArrivalProcessConcept {
      virtual ~ArrivalProcessConcept() {}
      virtual TTime doNextArrival( TTime fromTime, TTime maxTime, double rateFactor,
                                   TRand& randGen, double eps ) const = 0;
   };
   template <typename TSpec>
   class ArrivalProcessModel: public ArrivalProcessConcept {
   public:
      ArrivalProcessModel( const ArrivalProcess<TTime, TSpec>& proc_ ): proc( proc_ ) {}
      virtual ~ArrivalProcessModel() {}
      virtual TTime doNextArrival( TTime fromTime, TTime maxTime, double rateFactor,
                                   TRand& randGen, double eps ) const {
        return proc.nextArrivalTime( fromTime, maxTime, rateFactor, randGen, eps );
      }
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
                          TRand& randGen, double eps = 1e-5 ) {
     return object->doNextArrival( fromTime, maxTime, rateFactor, randGen, eps );
   }
   
}; // class ArrivalProcess<TTime, AnyProc< TRand > >


//
// ** Poisson processes
//


template <typename TSpec> class Poisson;

struct CinlarsMethod;
template <typename TRateFunctionSpec, typename TMethodSpec = CinlarsMethod> class NonHomogeneous;

template <typename TTime, typename TRateFunctionSpec>
class ArrivalProcess<TTime, Poisson< NonHomogeneous< TRateFunctionSpec, CinlarsMethod > > > {
   BOOST_MPL_ASSERT(( IsFuncSpec<TRateFunctionSpec> ));
public:
   typedef TTime result_type;
   typedef Func< TRateFunctionSpec > rate_function_type;
   typedef BOOST_TYPEOF_TPL(( integralFromPoint( boost::declval<rate_function_type>(),
                                                 boost::declval<TTime>() ) ))
    rate_function_integral_type;

   typedef typename rate_function_integral_type::result_type integral_val_type;

   typedef typename DiffType< integral_val_type >::type
   integral_diff_type;

   //typedef typename DiffType<TTime> time_diff_type;
   
   ArrivalProcess( const rate_function_type& rateFunction_, TTime startTime_ ):
     rateFunctionIntegral( integralFromPoint( rateFunction_, startTime_ ) ) {
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
                           double eps = 1e-5 ) const {
     BOOST_AUTO_TPL(integralAtFromTime, ( eval( rateFunctionIntegral, fromTime )  ));
     BOOST_AUTO_TPL(integralAtMaxTime, ( eval( rateFunctionIntegral, maxTime ) ));
#ifndef NDEBUG     
     BOOST_TYPEOF_TPL( integralAtFromTime ) zeroIntegral(0.0);
     assert( integralAtFromTime >= zeroIntegral );
     assert( integralAtMaxTime > zeroIntegral );
     assert( fromTime >= TTime(0.0) );
     assert( maxTime > fromTime );
     assert( integralAtMaxTime > integralAtFromTime );
#endif     
     typename boost::random::uniform_01<> uniformDistr;
     BOOST_AUTO_TPL( u,-std::log( uniformDistr( randGen ) ) );
     BOOST_AUTO_TPL(targetIntegralValue, integralAtFromTime +
                    static_cast< integral_diff_type >( u / rateFactor ) );
     if ( integralAtMaxTime < targetIntegralValue ) return maxTime;
     TTime  targetTime =  evalInverse( rateFunctionIntegral, fromTime, maxTime, targetIntegralValue,
                                       static_cast<integral_diff_type>( eps / ( rateFactor ) ) );

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
   
private:
   rate_function_integral_type rateFunctionIntegral;
   //integral_val_type integralAtZeroTime;
};  // class ArrivalProcess<TTime, Poisson< NonHomogeneous< TRateFunctionSpec > > >

template <typename TTime, typename TRateFunctionSpec >
struct result_of_makeNonHomogeneousPoissonProcess {
   BOOST_MPL_ASSERT(( IsFuncSpec<TRateFunctionSpec> ));
   typedef ArrivalProcess<TTime, Poisson< NonHomogeneous< TRateFunctionSpec > > > type;
};

template <typename TTime, typename TRateFunctionSpec >
typename result_of_makeNonHomogeneousPoissonProcess< TTime, TRateFunctionSpec >::type
makeNonHomogeneousPoissonProcess( const Func< TRateFunctionSpec >& rateFunction,
                                  TTime fromTime ) {
  typedef typename result_of_makeNonHomogeneousPoissonProcess<TTime, TRateFunctionSpec>::type
     result_type;
  return result_type( rateFunction, fromTime );
}


// * Utils
// ** Numerical integration

template <typename TSpec>
double FuncCaller( double x, void *data ) {
  const Func<TSpec> *f =  (const Func<TSpec> *)data;
  assert(data);
  return ToDouble( eval( *f, x ) );
}

template <typename TDom, typename TSpec>
typename AreaType<TDom, BOOST_TYPEOF_TPL( eval( boost::declval< Func<TSpec> >(), boost::declval< TDom >() ) )  >::type
inline
integrateNumerically( const Func<TSpec>& f,
                      TDom a, TDom b, int n ) {
  Func<TSpec>& f_nonconst = const_cast<Func<TSpec>&>( f );
  typedef double (*f_t)(double, void *);
  f_t fptr = &FuncCaller<TSpec>;
	typedef BOOST_TYPEOF_TPL( eval( boost::declval< Func<TSpec> >(), boost::declval< TDom >() ) ) TRan;
  typedef typename AreaType<TDom, TRan >::type result_type;
  return static_cast<result_type>( gauss_legendre( n, fptr, &f_nonconst, ToDouble( a ), ToDouble( b ) ) );
  
}

// ** Interpolated functions

// *** Class: InterpFn
//
// An interpolated function: stores a vector of (x,f(x)) pairs and
// can compute the values of f at new points by simple linear interpolation.
//
template <typename TDom = double, typename TRan = double>
class InterpFn /*: public AbstractFunction<TDom, TRan>*/ {
public:
   // Method: clear
   // Remove all (x,f(x)) pairs from this function.
   void clear() { x2y.clear(); }

   // Method: addPt
   // Store an (x,f(x)) pair specifying the value of the function at a given point 'x'.
   // Points must be added in order of increasing 'x'.
   void addPt( TDom x, TRan y ) {
     assert( x2y.empty() || x > boost::prior( x2y.end() )->first );
     x2y.insert( x2y.end(), make_pair( x, y ) );
   }

   // Method: eval
   // Evaluate the function at the specified point 'x', interpolating as necessary
   // between the neighboring values.
   TRan operator()( TDom x ) const {
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
     return TRan( result );
   }

   virtual TRan eval( TDom x ) const { return (*this)( x ); }

   size_t getNumPts() const { return x2y.size(); }
   
#if 0
   struct Pt {
      TDom x;
      TRan y;
      typename map<TDom,TRan>::const_iterator lwrBnd, uprBnd;
   };

   Pt evalPt( TDom x ) const {
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
   
   Pt evalPt( Pt from, Pt to, TDom x ) const {
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
     std::ofstream out( filename );
     out.exceptions( std::ios_base::failbit | std::ios_base::badbit );
     out.precision( 20 );
     out << "x\ty\n";
     for ( x2y_const_iterator it = x2y.begin(); it != x2y.end(); it++ )
        out << it->first << "\t" << it->second << "\n";
   }
   
   void load( filename_t filename ) {
     std::ifstream inFile( filename );
     inFile.exceptions( std::ios_base::failbit | std::ios_base::badbit );
     clear();
     size_t lineNum = 0;
     while ( inFile ) {
       std::string line;
       std::getline( inFile, line );
       if ( lineNum == 0 ) continue;
       std::istringstream lineStr( line );
       lineStr.exceptions( std::ios_base::failbit | std::ios_base::badbit );
       TDom x;
       TRan y;
       lineStr >> x >> y;
       addPt( x, y );
       lineNum++;
     }
   }

private:

   typedef std::map<TDom,TRan> x2y_type;
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
// *** Class: InterpBiFunc
//
// An bijective interpolated function, which can map values both from
// domain to range and from range to domain.
//
template <typename TDom = double, typename TRan = double >
class InterpBiFunc: public InterpFn<TDom, TRan> {
public:
   typedef InterpFn<TDom,TRan> PARENT;
   
   void clear() { PARENT::clear(); inverse.clear(); }
   void addPt( TDom x, TRan y ) { PARENT::addPt( x, y ); inverse.addPt( y, x ); }

   const InterpFn<TRan,TDom>& getInverse() const { return inverse; }
   
private:
   InterpFn<TRan,TDom> inverse;
};

template <typename TRan>
inline TRan interpolate( const std::vector<TRan>& f, double i ) {
  int int_i = static_cast<int>( floor( i ) );
  double frac_i = i - int_i;
  return static_cast<TRan>( ( fabs( int_i - i ) < 1e-5 ) ? f[ int_i ] :
                              (1-frac_i) * f[int_i] + frac_i*f[int_i+1] ); 
}

// ** Numerical utils: integrate, findWhere

template <typename TDom, typename TRan> inline
typename AreaType<TDom,TRan>::type integrate( const std::vector<TRan>& f, const std::vector<TDom>& x, double from, double to ) {

  typedef typename AreaType<TDom,TRan>::type result_type;

  assert( 0 <= from && from <= f.size()-1 );
  assert( 0 <= to && to <= f.size()-1 );

  TRan firstF = interpolate( f, from );
  TDom firstX = interpolate( x, from );
  TRan lastF = interpolate( f, to );
  TDom lastX = interpolate( x, to );

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
#endif
}  // namespace math
}  // namespace cosi

//////////////////////////////////////////////////////////////////////
// Boost.TypeOf registration
//////////////////////////////////////////////////////////////////////


//BOOST_TYPE_ERASURE_FREE((has_eval), cosi::math::eval, 2)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::MultType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::DiffType,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::DomainType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::RangeType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::SpecType,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::IsFuncSpec,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::FuncConcept,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::Func,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::Const,1)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math2::RunTime)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::CompileTime,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::X_To,(int)(typename))

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::BinOp,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::AddOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::SubOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::MultOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::DivOp,2)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::Piecewise,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::AreaType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::detail::result_of_make_line_through,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::UnaryOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::result_of_indefiniteIntegral,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::result_of_integralFromPoint,3)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math2::Any)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::ArrivalProcess,2);
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::Poisson,2)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math2::CinlarsMethod)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::NonHomogeneous,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::result_of_makeNonHomogeneousPoissonProcess,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::InterpFn,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math2::InterpBiFunc,2)

#endif  // #ifndef __INCLUDE_COSI_GENERALMATH2_H
// Postamble:1 ends here
