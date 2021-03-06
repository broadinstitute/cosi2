#ifndef __INCLUDE_COSI_GENERALMATH_H
#define __INCLUDE_COSI_GENERALMATH_H

// #+TITLE: Module generalmath
//
// A general representation of functions of one variable.  Functions are represented as trees, and
// various generic operations on functions are implemented: composition, inversion, differentiation,
// integration etc.  Functions are represnted analytically, as opposed to just black boxes,
// which allows exact analytic computation of many operations on functions, and permits detailed
// type-checking.

// * Interface

#include <boost/core/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>

namespace cosi {
namespace math {

// ** Template: Function - a function of one variable, from a specific domain to a specific range
// Template params:
//    TDomain - the domain type of the function.  Even if the function is constant, it still has a
//      specific domain type.
//    TRange - the range type of the function.
//    TSpec - definition of the function: what the function is, i.e. how it computes its result from its argument.
template <typename TDomain, typename TRange, typename TSpec> class Function;

template <typename TFunc> struct DomainType;
template <typename TFunc> struct RangeType;
template <typename TFunc> struct SpecType;

// ** Function: eval - evaluate a function at a given point.
template <typename TDomain, typename TRange, typename TSpec, typename TArg>
typename boost::enable_if< boost::is_convertible< TArg, TDomain >,
													 TRange >::type
eval( Function<TDomain, TRange, TSpec> const&, TArg val );

template <typename TDomain, typename TRange, typename TSpec, typename Enable=void>
struct result_of_indefiniteIntegral;

template <typename TDomain, typename TRange, typename TSpec>
typename result_of_indefiniteIntegral<TDomain, TRange, TSpec>::type
indefiniteIntegral( const Function< TDomain, TRange, TSpec >& f);

template <typename TDomain, typename TRange, typename TSpec> struct result_of_differentiate;

template <typename TDomain, typename TRange, typename TSpec>
typename result_of_differentiate<TDomain, TRange, TSpec>::type
differentiate( const Function< TDomain, TRange, TSpec >& f);

template <typename TDomain, typename TRange, typename TSpec> struct result_of_inverse;

template <typename TDomain, typename TRange, typename TSpec>
typename result_of_inverse<TDomain, TRange, TSpec>::type
inverse( const Function< TDomain, TRange, TSpec >& f);


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_compose;

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_compose<TDomain,TRange1,TRange2,TSpec1,TSpec2>::type
compose( const Function<TRange1,TRange2,TSpec1>& f1, const Function<TDomain,TRange1,TSpec2>& f2 );

} // namespace math
} // namespace cosi

// * Implementation

// * Preamble

// ** Headers

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
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/has_multiplies.hpp>
#include <boost/type_traits/has_divides.hpp>
#include <boost/type_traits/has_minus.hpp>
#include <boost/type_traits/has_plus.hpp>
#include <boost/core/demangle.hpp>
#include <boost/tti/has_type.hpp>
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
#include <boost/cstdint.hpp>
#include <boost/type_index.hpp>
#include <cosi/general/utils.h>
#include <cosi/general/math/gauss_legendre.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {
namespace math {


template <typename Op, typename T1, typename T2=void, typename T3=void, typename T4=void> struct rez;

namespace ops {
struct Add;
struct Sub;
struct Div;
struct Mult;
}

template <typename TDomain1, typename TRange1, typename TSpec1,
					typename TDomain2, typename TRange2, typename TSpec2> struct result_of_multiplies;

BOOST_TTI_HAS_TYPE(type)

//
// ** Utility metafunctions
//

// Metafunction: MultType
// Returns the type of the multiplication of values of the given type.
template <typename TVal1, typename TVal2, typename Enable=void> struct MultType;
template <typename TVal1, typename TVal2>
struct MultType<TVal1, TVal2, typename boost::enable_if< boost::has_multiplies< TVal1, TVal2 > >::type > {
	 typedef BOOST_TYPEOF_TPL( boost::declval<TVal1>() * boost::declval<TVal2>()) type;
};

// Metafunction: AddType
// Returns the type of the multiplication of values of the given type.
template <typename TVal1, typename TVal2=TVal1> struct AddType {
	 typedef BOOST_TYPEOF_TPL( boost::declval<TVal1>() + boost::declval<TVal2>()) type;
};


// Metafunction: DiffType
// Returns the type of the difference of values of the given type.
template <typename TVal1, typename TVal2 = TVal1> struct DiffType {
   typedef BOOST_TYPEOF_TPL( boost::declval<TVal1>() - boost::declval<TVal2>()) type;
};

// Metafunction: DivType
// Returns the type of the division of values of the given type.
template <typename TVal1, typename TVal2> struct DivType {
	 typedef BOOST_TYPEOF_TPL( boost::declval<TVal1>() / boost::declval<TVal2>()) type;
};

template <typename T> struct InvType: public DivType< double, T> {
};

// Metafunction: DerivType
// The type for a derivative, i.e. the result of differentiating a function with given domain and range.
template <typename TDomain, typename TRange>
struct DerivType {
   typedef typename DiffType<TDomain>::type domain_diff_type;
   typedef typename DiffType<TRange>::type range_diff_type;
   typedef typename DivType<range_diff_type, domain_diff_type>::type type;
};

// Metafunction: AreaType
// The type for an area, i.e. the result of integrating a function with given domain and range.


template <typename TDomain, typename TRange>
struct has_area: public boost::mpl::and_< boost::has_minus< TDomain >,
																					boost::has_minus< TRange >,
																					boost::has_multiplies< typename DiffType<TDomain>::type,
																																 typename DiffType<TRange>::type > > { };

template <typename TDomain, typename TRange, typename TSpec>
struct can_integrate: public boost::mpl::and_< has_area<TDomain, TRange>,
																							 has_type_type< result_of_indefiniteIntegral<TDomain,TRange,TSpec> > > {
};


template <typename TDomain, typename TRange, typename Enable=void>
struct AreaType;

template <typename TDomain, typename TRange>
struct AreaType<TDomain, TRange, typename boost::enable_if< has_area< TDomain, TRange > >::type > {
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

//template <typename T> std::string getLabel( const T& );

// *** GenericClass Function - A mathematical function of one variable.
//template <typename TDomain, typename TRange, typename TSpec> class Function;
        
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

template <typename TDomain, typename TRange, typename TSpec>
struct DomainType< Function<TDomain, TRange, TSpec> > { typedef TDomain type; };
template <typename TDomain, typename TRange, typename TSpec>
struct RangeType< Function<TDomain, TRange, TSpec> > { typedef TRange type; };
template <typename TDomain, typename TRange, typename TSpec>
struct SpecType< Function<TDomain, TRange, TSpec> > { typedef TSpec type;
	    BOOST_MPL_ASSERT(( IsFunctionSpec<type> ));
};

// template <typename TFunc> struct RangeType { typedef typename TFunc::result_type type; };
// template <typename TFunc> struct SpecType {
// 	 typedef typename TFunc::spec_type type;
//    BOOST_MPL_ASSERT(( IsFunctionSpec<type> ));
// };

//
// *** Generic function: eval
//
// Default implementation of eval(f,x) to evaluate a Function at a given point in its domain;
// just calls the function's operator() .
// 
// template <typename TFunc, typename TArg>
// typename boost::enable_if< boost::is_convertible< TArg, typename DomainType<TFunc>::type >,
//                            typename RangeType<TFunc>::type >::type
// eval( const TFunc& f, TArg x ) { return f( x ); }

// template <typename TDomain, typename TRange, typename TSpec, typename TArg>
// typename boost::enable_if< boost::is_convertible< TArg, TDomain >,
//                            TRange >::type
// eval( Function<TDomain, TRange, TSpec> const& f, TArg x ) { return f( x ); }

// ** Particular function kinds
// *** Constant functions
//    We support run-time constant functions and compile-time constant functions.

template <typename TDomain, typename TRange, typename TSpec, typename TArg>
typename boost::enable_if< boost::is_convertible< TArg, TDomain >,
													 TRange >::type
eval( Function<TDomain, TRange, TSpec> const& f, TArg x ) { return f(x); }

struct RunTime;
template <typename TConstSpec = RunTime> struct Const;
template <> struct IsFunctionSpec< Const<RunTime> >: public boost::true_type {};

template <int N> struct CompileTime { static const int value = N; };

template <int N> struct IsFunctionSpec< Const< CompileTime<N> > >: public boost::true_type {};

template <typename TDomain, typename TRange, typename TSpec>
TRange evalConst( const Function< TDomain, TRange, Const< TSpec > >& f ) { return eval<TDomain, TRange, Const<TSpec>,
																																											 TDomain>( f, TDomain() ); }

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

	 Function( ): val( std::numeric_limits<TRange>::quiet_NaN() ) { }
   Function( TRange val_ ): val( val_ ) { }

   TRange operator()( TDomain ) const { return val; }

   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
     s << f.val; return s;
   }

	 bool operator==( const Function& f ) const { return val == f.val; }
	 bool operator!=( const Function& f ) const { return val != f.val; }

	 TRange& getValRef() { return val; }

private:
   TRange val;
};

template <typename TDomain, typename TRange> inline
Function< TDomain, TRange, Const<> >
fn_const( TRange v ) { return Function< TDomain, TRange, Const<> >( v ); }


template <typename TDomain, typename TRange, int val>
class Function<TDomain, TRange, Const< CompileTime<val> > > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef Const< CompileTime<val> > spec_type;
   TRange operator()( TDomain ) const { return static_cast<TRange>( val ); }

   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
     s << val; return s;
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

	 // BOOST_MPL_ASSERT(( boost::is_convertible< typename MultType<FactorType, TDomain >::type,
	 // 										result_type >
	 // 										));

	 
	 // BOOST_MPL_ASSERT(( boost::is_convertible< typename MultType<FactorType, domain_diff_type >::type,
	 // 										result_type >
	 // 										));

   Function( ): factor( 1.0 ) {}
   Function( factor_type factor_ ): factor( factor_ ) {}
   //result_type operator()( TDomain x ) const { return factor * ( x - static_cast< TDomain >( 0.0 ) ) ; }
	 result_type operator()( TDomain x ) const { return result_type( ToDouble( factor ) * ToDouble( x ) ); }
   factor_type getFactor() const { return factor; }

   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
		 typedef BOOST_TYPEOF_TPL( f.factor ) f_t;
		 if ( f.factor == f_t(1.0) ) s << "x";
		 else
				s << "(" << f.factor << " * x" << ")";
		 return s; 
   }

   
private:
   factor_type factor;
};

template <typename TDomain, typename TRange>
Function< TDomain, TRange, X_To<1> > fn_x() { return Function< TDomain, TRange, X_To<1> >(); }

template <typename TDomain, typename TFactor>
Function< TDomain, typename MultType<TFactor, typename DiffType<TDomain>::type >::type,
						X_To<1,TFactor> >
fn_x( TFactor c ) {
	return Function< TDomain, typename MultType<TFactor, typename DiffType<TDomain>::type >::type,
									 X_To<1,TFactor> >( c );
}


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
   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
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

   friend std::ostream& operator<<( std::ostream& s, const AddOp& ) { s << "+"; return s; }
};
template <typename T1, typename T2 = T1>
struct SubOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() - boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 - a2; }
   friend std::ostream& operator<<( std::ostream& s, const SubOp& ) { s << "-"; return s; }
};
template <typename T1, typename T2 = T1>
struct MultOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() * boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 * a2; }
   friend std::ostream& operator<<( std::ostream& s, const MultOp& ) { s << "*"; return s; }
};
template <typename T1, typename T2 = T1>
struct DivOp {
   typedef BOOST_TYPEOF_TPL( boost::declval<T1>() / boost::declval<T2>() ) result_type;
   result_type operator()( T1 a1, T2 a2 ) const { return a1 / a2; }
   friend std::ostream& operator<<( std::ostream& s, const DivOp& ) { s << "/"; return s; }
};

// *** Binary Op Function

template <typename TSpec1, typename TSpec2, typename Op> struct BinOp;
template <typename TSpec1, typename TSpec2, typename Op>
struct IsFunctionSpec< BinOp< TSpec1, TSpec2, Op >,
                       typename boost::enable_if< boost::mpl::and_< IsFunctionSpec< TSpec1 >,
																																		IsFunctionSpec< TSpec2 > > >::type >:
		 public boost::true_type {};

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

   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
     s << "(" << f.first() << " " << f.binop << " " << f.second() << ")"; return s;
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


template <typename TDomain, typename TRange1, typename TRange2, typename TConstSpec1, typename TConstSpec2>
Function<TDomain, typename AddOp<TRange1,TRange2>::result_type, Const<> >
operator+( const Function< TDomain, TRange1, Const<TConstSpec1> >& f1,
           const Function< TDomain, TRange2, Const<TConstSpec2> >& f2 ) {
  return Function< TDomain, typename AddOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) + evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TConstSpec1, typename TConstSpec2>
Function<TDomain, typename SubOp<TRange1,TRange2>::result_type, Const<> >
operator-( const Function< TDomain, TRange1, Const<TConstSpec1> >& f1,
           const Function< TDomain, TRange2, Const<TConstSpec2> >& f2 ) {
  return Function< TDomain, typename SubOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) - evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TConstSpec1, typename TConstSpec2>
Function<TDomain, typename MultOp<TRange1,TRange2>::result_type, Const<> >
operator*( const Function< TDomain, TRange1, Const<TConstSpec1> >& f1,
           const Function< TDomain, TRange2, Const<TConstSpec2> >& f2 ) {
  return Function< TDomain, typename MultOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) * evalConst( f2 ) );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TConstSpec1, typename TConstSpec2>
Function<TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >
operator/( const Function< TDomain, TRange1, Const<TConstSpec1> >& f1,
           const Function< TDomain, TRange2, Const<TConstSpec2> >& f2 ) {
  return Function< TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) / evalConst( f2 ) );
}

// template <typename TDomain1, typename TRange1, typename TSpec1, typename TDomain2, typename TRange2, typename TSpec2>
// struct rez<ops::Mult, Function<TDomain1,TRange1,Const<> >, Function<TDomain2,TRange2,X_To<1> > > {
// 	 typedef Function< TDomain, TRange, Const<> > f1_t;
// 	 typedef Function<TDomain, TRange, X_To<1> > f2_t;
// 	 typedef BOOST_TYPEOF_TPL(( 
// };

//typename rez<ops::Mult, Function< TDomain, TRange, Const<> >, Function<TDomain, TRange, X_To<1> > >::type
template <typename TDomain, typename TRange1, typename TRange2>
Function<TDomain, typename MultType<TRange1,TRange2>::type, X_To<1> >
operator*( const Function< TDomain, TRange1, Const<> >& f1, const Function<TDomain, TRange2, X_To<1> >& f2 ) {
	return Function<TDomain, typename MultType<TRange1,TRange2>::type, X_To<1> >( evalConst( f1 ) * f2.getFactor() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TRange3, typename TSpec>
struct multAresult {
	 typedef Function< TDomain, TRange1, Const<> > f1_t;
	 typedef Function< TDomain, typename MultOp<TRange2,TRange3>::result_type,
										 BinOp< Const<>,
														TSpec,
														MultOp< TRange2, TRange3 >
														>
										 > f2_t;
	 typedef BOOST_TYPEOF_TPL(( ( boost::declval<f1_t>() * boost::declval<f2_t>().first() ) *
															boost::declval<f2_t>().second() )) type;
};

// template <typename TDomain, typename TRange1, typename TRange2, typename TRange3, typename TSpec>
// // Function< TDomain, typename MultOp< typename MultOp<TRange1,TRange2>::result_type, TRange3 >::result_type,
// // 					BinOp< Const<>, TSpec, MultOp< typename MultOp< TRange1, TRange2 >::result_type, TRange3 >

template <typename TDomain, typename TRange1, typename TRange2, typename TRange3, typename TSpec>
typename multAresult<TDomain, TRange1, TRange2, TRange3, TSpec>::type
operator*( Function< TDomain, TRange1, Const<> > const& f1,
					 Function< TDomain, typename MultOp<TRange2,TRange3>::result_type,
										 BinOp< Const<>,
														TSpec,
														MultOp< TRange2, TRange3 >
														>
										 > const& f2
	) {
	return ( f1 * f2.first() ) * f2.second();
}



// template <typename TDomain, typename TRange1, typename TRange2, typename TConstSpec1, typename TConstSpec2>
// Function<TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >
// operator/( const Function< TDomain, TRange1, Const<TConstSpec1> >& f1,
//            const Function< TDomain, TRange2, BinConst<TConstSpec2> >& f2 ) {
//   return Function< TDomain, typename DivOp<TRange1,TRange2>::result_type, Const<> >( evalConst( f1 ) / evalConst( f2 ) );
// }



// *** cval: convenience wrapper for const values
//
// We define a convenience wrapper cval(x), which saves the value and its type.
// The result can then be used in binary operations with functions, without worrying about specifying
// the domain and range: the domain is assumed the same as the function with which we're operating,
// and the range is the saved type of the value.

template <typename TVal>
struct CVal {
	 TVal val;

	 explicit CVal( TVal val_ ): val( val_ ) { }

	 template <typename TDomain>
	 operator Function< TDomain, TVal, Const<> >() const {
		 return Function< TDomain, TVal, Const<> >( val );
	 }
};

template <typename TVal>
 inline CVal<TVal> cval( TVal x ) { return CVal<TVal>( x ); }

namespace detail {

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_mult {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TVal, Const<> > >() *
											boost::declval< Function< TDomain, TRange, TSpec> >() )) type;
};

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_sub {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TVal, Const<> > >() -
											boost::declval< Function< TDomain, TRange, TSpec> >() )) type;
};

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_sub_r {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TRange, TSpec> >() -
											boost::declval< Function< TDomain, TVal, Const<> > >() )) type;
};

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_div {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TVal, Const<> > >() /
											boost::declval< Function< TDomain, TRange, TSpec> >() )) type;
};

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_div_r {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TRange, TSpec> >() /
											boost::declval< Function< TDomain, TVal, Const<> > >() )) type;
};

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
struct result_of_factor_add {
	 typedef
	 BOOST_TYPEOF_TPL(( boost::declval< Function< TDomain, TVal, Const<> > >() +
											boost::declval< Function< TDomain, TRange, TSpec> >() )) type;
};

}  // namespace detail

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_multiplies< TVal, TRange>,
													 typename detail::result_of_factor_mult<TVal,TDomain,TRange,TSpec>::type >::type
operator*( CVal<TVal> v,
           const Function< TDomain, TRange, TSpec >& f ) {
  return Function< TDomain, TVal, Const<> >( v.val ) * f;
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_multiplies< TVal, TRange>,
													 typename detail::result_of_factor_mult<TVal,TDomain,TRange,TSpec>::type >::type
operator*( const Function< TDomain, TRange, TSpec >& f, CVal<TVal> v ) {
  return Function< TDomain, TVal, Const<> >( v.val ) * f;
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_minus< TVal, TRange>,
													 typename detail::result_of_factor_sub<TVal,TDomain,TRange,TSpec>::type >::type
operator-( CVal<TVal> v,
           Function< TDomain, TRange, TSpec > const& f ) {
  return Function< TDomain, TVal, Const<> >( v.val ) - f;
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_minus< TRange, TVal >,
													 typename detail::result_of_factor_sub_r<TVal,TDomain,TRange,TSpec>::type >::type
operator-( const Function< TDomain, TRange, TSpec >& f, CVal<TVal> v ) {
  return f - Function< TDomain, TVal, Const<> >( v.val );
}


template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_divides< TVal, TRange>,
													 typename detail::result_of_factor_div<TVal,TDomain,TRange,TSpec>::type >::type
operator/( CVal<TVal> v,
           const Function< TDomain, TRange, TSpec >& f ) {
  return Function< TDomain, TVal, Const<> >( v.val ) / f;
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_divides< TRange, TVal>,
													 typename detail::result_of_factor_div_r<TVal,TDomain,TRange,TSpec>::type >::type
operator/( const Function< TDomain, TRange, TSpec >& f, CVal<TVal> v ) {
  return f / Function< TDomain, TVal, Const<> >( v.val );
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_plus< TVal, TRange>,
													 typename detail::result_of_factor_add<TVal,TDomain,TRange,TSpec>::type >::type
operator+( CVal<TVal> v,
           const Function< TDomain, TRange, TSpec >& f ) {
  return Function< TDomain, TVal, Const<> >( v.val ) + f;
}

template <typename TVal, typename TDomain, typename TRange, typename TSpec>
typename boost::enable_if< boost::has_plus< TRange, TVal >,
													 typename detail::result_of_factor_add<TVal,TDomain,TRange,TSpec>::type >::type
operator+( const Function< TDomain, TRange, TSpec >& f, CVal<TVal> v ) {
  return Function< TDomain, TVal, Const<> >( v.val ) + f;
}



//
// *** Unary operations on functions
//
//  For a function f, -f is a function that for value x returns -f(x).
//

template <typename TSpec, typename Op> struct UnaryOp;
template <typename TSpec, typename Op>
struct IsFunctionSpec< UnaryOp<TSpec, Op>,
											 typename boost::enable_if< IsFunctionSpec< TSpec > >::type >: public boost::true_type {};
																						 

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

//   friend std::ostream& operator<<( std::ostream& s, const Function& f ) { s << op << " " << f.f; return s; }

	 const function_type& getFunction() const { return f; }
   
private:
   static const Op op;
   const function_type f;
   
};

template <typename TDomain, typename TRange, typename TSpec, typename Op>
const Op Function<TDomain, TRange, UnaryOp< TSpec, Op > >::op = Op();

template <typename T> 
std::ostream& operator<<( std::ostream& s, const std::negate< T >& ) { s << "-"; return s; }


template <typename TDomain, typename TRange, typename TSpec>
Function<TDomain, TRange, UnaryOp<TSpec, std::negate< TRange > > >
operator-( const Function<TDomain,TRange,TSpec>& f ) {
  return Function<TDomain, TRange, UnaryOp<TSpec, std::negate< TRange > > >( f );
}

template <typename TDomain, typename TRange, typename TSpec>
std::ostream& operator<<( std::ostream& s, const Function<TDomain,TRange,UnaryOp<TSpec, std::negate<TRange> > >& f ) {
	s << "-(" << f.getFunction() << ")"; return s;
}


struct Exp
{
	 template <typename T>
	 T operator()(T x) const { return exp(x);}
};

template <typename TDomain, typename TSpec>
std::ostream& operator<<( std::ostream& s, const Function<TDomain,double,UnaryOp<TSpec, Exp> >& f ) {
	s << "e^(" << f.getFunction() << ")"; return s;
}

struct Log
{
	 template <typename T>
	 T operator()(T x) const { return log(x);}
};
template <typename TDomain, typename TSpec>
std::ostream& operator<<( std::ostream& s, const Function<TDomain,double,UnaryOp<TSpec, Log> >& f ) {
	s << "ln(" << f.getFunction() << ")"; return s;
}


template <typename TDomain, typename TSpec>
Function<TDomain, double, UnaryOp<TSpec,Exp> >
exp_( const Function<TDomain,double,TSpec>& f ) {
	return Function<TDomain, double, UnaryOp<TSpec,Exp> >( f );
}

template <typename TDomain, typename ConstKind>
Function< TDomain, double, Const<> >
exp_( const Function<TDomain, double, Const<ConstKind> >& f ) {
	return Function< TDomain, double, Const<> >( exp( evalConst( f ) ) );
}

template <typename TDomain, typename TSpec>
Function<TDomain, double, UnaryOp<TSpec,Log> >
log_( const Function<TDomain,double,TSpec>& f ) {
	return Function<TDomain, double, UnaryOp<TSpec,Log> >( f );
}

template <typename TDomain, typename TSpec>
Function<TDomain, double, TSpec>
exp_( const Function<TDomain, double, UnaryOp<TSpec,Log> >& f ) {
	return f.getFunction();
}

template <typename TDomain, typename TSpec>
Function<TDomain, double, TSpec>
log_( const Function<TDomain, double, UnaryOp<TSpec,Exp> >& f ) {
	return f.getFunction();
}

template <typename TDomain, typename TRange, typename TSpec>
struct result_of_exp;
#if 0
template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_exp<TDomain, double,
										 BinOp<TSpec1, TSpec2, MultOp<TRange1,TRange2> > > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
	 typedef Function< TDomain, TRange1, TSpec1 > f1_t;
	 typedef Function< TDomain, TRange2, TSpec2 > f2_t;

	 typedef BOOST_TYPEOF_TPL(( exp_( boost::declval<f1_t>() ) * exp_( boost::declval<f2_t>() ) ) ) type;

   // typedef BOOST_TYPEOF_TPL(( exp_( boost::declval< typename result_of_indefiniteIntegral<TDomain, TRange1, TSpec1>::type f1_t;
   // typedef typename result_of_indefiniteIntegral<TDomain, TRange2, TSpec2>::type f2_t;
   // typedef BOOST_TYPEOF_TPL(( boost::declval<f1_t>() + boost::declval<f2_t>() )) type;
};

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_exp< TDomain, double, 
												BinOp< TSpec1, TSpec2, MultOp<TRange1,TRange2> > >::type
exp_( const Function< TDomain, double,
			BinOp< TSpec1, TSpec2, MultOp<TRange1, TRange2> > >& f ) {
	return exp_( f.first() ) * exp_( f.second() );
}
#endif

template <typename TDomain, typename TRange1, typename TSpec1,
					typename TRange2, typename TSpec2> struct result_of_pow {
	 //typedef Function< TDomain, double, UnaryOp< BinOp< TSpec2, UnaryOp< Log, TSpec1 >, MultOp
	 typedef Function< TDomain, TRange1, TSpec1 > f_base_t;
	 typedef Function< TDomain, TRange2, TSpec2 > f_exponent_t;
	 typedef BOOST_TYPEOF_TPL( exp_( boost::declval<f_exponent_t>() * log_( boost::declval<f_base_t>() ) ) ) type;
};

template< typename TDomain, typename TRange1, typename TSpec1, typename TRange2, typename TSpec2 >
typename result_of_pow< TDomain, TRange1, TSpec1, TRange2, TSpec2>::type
pow( Function< TDomain, TRange1, TSpec1 > const& f_base,
		 Function< TDomain, TRange2, TSpec2 > const& f_exponent ) {
	return exp_( f_exponent * log_( f_base ) );
}

template <typename TDomain, typename TRange1, typename TSpec1, typename TRange2, typename TSpec2>
struct result_of_pow< TDomain, typename MultOp<TRange1,TRange2>::result_type,
											BinOp< TSpec1, TSpec2, MultOp<TRange1,TRange2> >,
											int, Const< CompileTime< -1 > > > {
	typedef Function< TDomain, TRange1, TSpec1 > f1_t;
	typedef Function< TDomain, TRange2, TSpec2 > f2_t;
	typedef Function< TDomain, double, Const< CompileTime< -1 > > > f_exponent_t;
	typedef BOOST_TYPEOF_TPL(( pow( boost::declval<f1_t>(), boost::declval<f_exponent_t>() ) *
														 pow( boost::declval<f2_t>(), boost::declval<f_exponent_t>() ) )) type;
};

// template <typename TDomain, typename TRange1, typename TSpec1, typename TRange2, typename TSpec2>
// typename result_of_pow< TDomain, typename MultOp<TRange1,TRange2>::result_type,
// 												BinOp< TSpec1, TSpec2, MultOp<TRange1,TRange2> >,
// 												int, Const< CompileTime< -1 > > >::type
// pow( Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
// 		 BinOp< TSpec1, TSpec2, MultOp<TRange1,TRange2> > > const& f_base,
// 		 Function< TDomain, double, Const< CompileTime< -1 > > > const& f_exponent ) {
// 	return pow( f_base.first(), f_exponent ) * pow( f_base.second(), f_exponent );
// }

// template <typename TDomain, typename TSpec>
// Function< TDomain, double, UnaryOp< UnaryOp<TSpec,Exp>, std::negate<double> >
// pow( Function< TDomain, double, UnaryOp<TSpec,Exp> > const& f_base,
// 		 Function< TDomain, int, Const< CompileTime< -1 > > > const& f_exponent ) {
// 	return exp_( -f_base.getFunction() );
// }




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
// #if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 6 ) )
// #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

template <typename TPieceSpec> struct Piecewise;
template <typename TPieceSpec> struct IsFunctionSpec< Piecewise<TPieceSpec>,
																											typename boost::enable_if< IsFunctionSpec<TPieceSpec> >::type >:
 public boost::true_type {};

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
	 explicit Function( const SinglePassRange& r ):
		 pieces( boost::begin( r ), boost::end( r ) ) {
				
		 BOOST_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< SinglePassRange > ));
		 typedef typename boost::range_value< SinglePassRange >::type range_value_type;
		 BOOST_STATIC_ASSERT(( boost::is_convertible< range_value_type,
													 typename pieces_type::value_type >::value ));
	 }
	 Function( Function const& f ):
		 pieces( f.pieces ) { }
   
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
   
	 friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
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

template <typename TDomain, typename TRange>
void set( Function< TDomain, TRange, Const<> >& f, TDomain x, TRange val ) {
	f.getPieces.insert( std::make_pair( x, fn_const<TDomain>( val ) ) );
}

template <typename TDomain, typename TRange, typename TPieceSpec>
typename boost::enable_if< boost::is_convertible< Function< TDomain, TRange, Const<> >,
																									Function< TDomain, TRange, TPieceSpec > > >::type
set( Function< TDomain, TRange, Piecewise< TPieceSpec > >& f, TDomain x, TRange val ) {
	f.getPieces().insert( std::make_pair( x, fn_const<TDomain>( val ) ) );
}


#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic pop
#endif
#endif

// Method: loadFrom - read piecewise function from a file
template <typename TDomain, typename TRange>
void loadFrom( std::istream& s, Function<TDomain, TRange, Piecewise< Const<> > >& f ) {
	boost::io::ios_exception_saver save_exceptions( s );
	s.exceptions( std::ios::failbit | std::ios::badbit );
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

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TPieceSpec2>
struct result_of_log {
	 typedef BOOST_TYPEOF_TPL(( boost::declval< Function<TDomain,TRange1,TSpec1> >() *
															boost::declval< Function<TDomain,TRange2,TPieceSpec2> >() )) piece_function_type;
	 typedef typename RangeType< piece_function_type >::type range_type;
	 typedef typename SpecType< piece_function_type >::type piece_spec_type;
	 typedef Function< TDomain, range_type, Piecewise< piece_spec_type > > type;
};


template <typename TDomain, typename TPieceSpec>
Function< TDomain, double, Piecewise< UnaryOp< TPieceSpec, Log > > >
log_( const Function<TDomain,double,Piecewise<TPieceSpec> >& f ) {
	Function< TDomain, double, Piecewise< UnaryOp< TPieceSpec, Log > > > f_result;
 
	for ( typename Function< TDomain, double, Piecewise<TPieceSpec> >::pieces_type::const_iterator it =
					 f.getPieces().begin(); it != f.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, log_( it->second ) ) );
	return f_result;
}


template <typename TDomain, typename TPieceSpec>
Function< TDomain, double, Piecewise< UnaryOp< TPieceSpec, Exp > > >
exp_( const Function<TDomain,double,Piecewise<TPieceSpec> >& f ) {
	Function< TDomain, double, Piecewise< UnaryOp< TPieceSpec, Exp > > > f_result;
 
	for ( typename Function< TDomain, double, Piecewise<TPieceSpec> >::pieces_type::const_iterator it =
					 f.getPieces().begin(); it != f.getPieces().end(); it++ )
		 f_result.getPieces().insert( std::make_pair( it->first, exp_( it->second ) ) );
	return f_result;
}

namespace detail {
template <typename TDomain, typename TRange, typename TPieceSpec>
struct transformPiece {
	 typedef Function<TDomain,TRange,TPieceSpec> piece_type;
	 typedef typename result_of_inverse<TDomain, TRange, TPieceSpec>::type piece_inv_type;

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<piece_inv_type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<piece_inv_type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<piece_inv_type>::type > ));
	 
	 typedef typename SpecType<piece_inv_type>::type piece_inv_spec_type;

	 typedef std::pair<const TRange,piece_inv_type> result_type;
	 result_type operator()( const std::pair<const TDomain, piece_type>& p  ) const {
		 const TDomain& a = p.first;
		 const piece_type& piece = p.second;
		 TRange piece_at_a = eval( piece, a );
		 piece_inv_type piece_inv = inverse( piece );
		 return std::make_pair( piece_at_a, piece_inv );
	 }
};
}

template <typename TDomain, typename TRange, typename TPieceSpec>
struct result_of_inverse<TDomain, TRange, Piecewise<TPieceSpec> > {
	 typedef Function<TRange,TDomain,
										Piecewise< typename detail::transformPiece<TDomain, TRange, TPieceSpec>::piece_inv_spec_type > >
		type;
	 
	 type operator()( const Function<TDomain, TRange, Piecewise<TPieceSpec> >& f ) {
		 return
				type( boost::adaptors::transform( f.getPieces(),
																					detail::transformPiece<TDomain,TRange,TPieceSpec>() ) );
	 }

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<type>::type > ));
};

// struct invert_piece {
// 	 operator
// };


// template <typename TDomain, typename TPieceSpec>
// exp_( Function< TDomain, double, UnaryOp< Piecewise<TPieceSpec>, Exp> const& f ) {
// 	typename result_of_exp< TDomain, TPieceSpec >::type f_result;
	
// }




// ** Integrals


// *** Definition of indefiniteIntegral function

// template <typename TDomain, typename TRange, typename TSpec> struct result_of_indefiniteIntegral {
// 	 BOOST_STATIC_ASSERT_MSG( sizeof( TDomain ) == 0, "don't know how to get type of indef integral of f" );
// };

template <typename TDomain, typename TRange, typename TSpec>
void indefiniteIntegral( const Function< TDomain, TRange, TSpec >& f) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  BOOST_STATIC_ASSERT_MSG( sizeof( f ) == 0, "don't know how to integrate f" );
}

template <typename TDomain, typename TRange, typename TSpec> struct result_of_differentiate {
	 BOOST_STATIC_ASSERT_MSG( sizeof( TDomain ) == 0, "don't know how to get type of differentiation of f" );
};

template <typename TDomain, typename TRange, typename TSpec>
void differentiate( const Function< TDomain, TRange, TSpec >& f) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  BOOST_STATIC_ASSERT_MSG( sizeof( f ) == 0, "don't know how to differentiate f" );
}

template <typename TDomain, typename TRange, typename TConstSpec>
struct result_of_differentiate<TDomain, TRange, Const<TConstSpec> > {
	 typedef typename DerivType<TDomain,TRange>::type deriv_type;
	 typedef Function< TDomain, deriv_type, Const< CompileTime<0> > > type;
};

template <typename TDomain, typename TRange, typename TConstSpec >
typename result_of_differentiate<TDomain,TRange,Const<TConstSpec> >::type
differentiate( const Function< TDomain, TRange, Const<TConstSpec> >& ) {
	return typename result_of_differentiate<TDomain,TRange,Const<TConstSpec> >::type();
}

template <typename TDomain, typename TRange, typename TFactor>
struct result_of_differentiate<TDomain, TRange, X_To<1,TFactor> > {
	 typedef typename DerivType<TDomain,TRange>::type deriv_type;
   typedef Function< TDomain, deriv_type, Const<> > type;
};

template <typename TDomain, typename TRange, typename TFactor>
typename result_of_differentiate<TDomain,TRange, X_To<1,TFactor> >::type
differentiate( const Function< TDomain, TRange, X_To<1,TFactor> >& f ) {
	return typename result_of_differentiate<TDomain,TRange,X_To<1,TFactor> >::type( f.getFactor() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
struct result_of_differentiate<TDomain, typename MultOp<TRange1,TRange2>::result_type,
															 BinOp< Const<>, TSpec, MultOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
   typedef Function< TDomain, TRange1, Const<> > f_const_t;
   typedef typename result_of_differentiate<TDomain, TRange2, TSpec>::type f_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f_const_t>() * boost::declval<f_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
typename result_of_differentiate<TDomain, typename MultOp<TRange1, TRange2>::result_type,
																 BinOp< Const<>, TSpec,
																				MultOp<TRange1,TRange2> > >::type
differentiate( const Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
                    BinOp< Const<>, TSpec,
                    MultOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  return f.first() * differentiate( f.second() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
typename result_of_differentiate<TDomain, typename MultOp<TRange1, TRange2>::result_type,
                                      BinOp< Const<>, TSpec,
                                             MultOp<TRange1,TRange2> > >::type
differentiate( const Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec, Const<>,
                    MultOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  return f.second() * differentiate( f.first() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_differentiate<TDomain, typename AddOp<TRange1,TRange2>::result_type,
                                    BinOp<TSpec1, TSpec2, AddOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
   typedef typename result_of_differentiate<TDomain, TRange1, TSpec1>::type f1_t;
   typedef typename result_of_differentiate<TDomain, TRange2, TSpec2>::type f2_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f1_t>() + boost::declval<f2_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_differentiate<TDomain, typename AddOp<TRange1, TRange2>::result_type,
                                      BinOp< TSpec1, TSpec2,
                                             AddOp<TRange1,TRange2> > >::type
differentiate( const Function< TDomain, typename AddOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec1, TSpec2,
                    AddOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
  return differentiate( f.first() ) + differentiate( f.second() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_differentiate<TDomain, typename SubOp<TRange1,TRange2>::result_type,
                                    BinOp<TSpec1, TSpec2, SubOp<TRange1,TRange2> > >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
   typedef typename result_of_differentiate<TDomain, TRange1, TSpec1>::type f1_t;
   typedef typename result_of_differentiate<TDomain, TRange2, TSpec2>::type f2_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f1_t>() - boost::declval<f2_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
typename result_of_differentiate<TDomain, typename SubOp<TRange1, TRange2>::result_type,
                                      BinOp< TSpec1, TSpec2,
                                             SubOp<TRange1,TRange2> > >::type
differentiate( const Function< TDomain, typename SubOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec1, TSpec2,
                    SubOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
  return differentiate( f.first() ) - differentiate( f.second() );
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

template <typename TDomain, typename TRange, typename TSpec, typename Enable = void>
struct result_of_integralFromPoint;

template <typename TDomain, typename TRange, typename TSpec>
struct result_of_integralFromPoint< TDomain, TRange, TSpec,
																		typename boost::enable_if< can_integrate< TDomain, TRange, TSpec > >::type > {
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
struct result_of_indefiniteIntegral<TDomain, TRange, Const<>,
																		typename boost::enable_if< has_area< TDomain, TRange > >::type > {
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
struct result_of_indefiniteIntegral<TDomain, TRange, X_To<1,TFactor>,
																		typename boost::enable_if< has_area< TDomain, TRange > >::type > {

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
                                    BinOp< Const<>, TSpec, MultOp<TRange1,TRange2> >,
																		typename boost::enable_if< can_integrate< TDomain, TRange2, TSpec > >::type >  {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
   typedef Function< TDomain, TRange1, Const<> > f_const_t;
   typedef typename result_of_indefiniteIntegral<TDomain, TRange2, TSpec>::type f_t;
   typedef BOOST_TYPEOF_TPL(( boost::declval<f_const_t>() * boost::declval<f_t>() )) type;
};


template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
struct result_of_indefiniteIntegral<TDomain, typename MultOp<TRange1,TRange2>::result_type,
                                    BinOp< TSpec, Const<>, MultOp<TRange1,TRange2> >,
																		typename boost::enable_if< can_integrate< TDomain, TRange2, TSpec > >::type >  {
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
                                      BinOp< TSpec, Const<>,
                                             MultOp<TRange1,TRange2> > >::type
indefiniteIntegral( const Function< TDomain, typename MultOp<TRange1,TRange2>::result_type,
                    BinOp< TSpec, Const<>,
                    MultOp<TRange1,TRange2> > >& f ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
  return f.second() * indefiniteIntegral( f.first() );
}

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec1, typename TSpec2>
struct result_of_indefiniteIntegral<TDomain, typename AddOp<TRange1,TRange2>::result_type,
                                    BinOp<TSpec1, TSpec2, AddOp<TRange1,TRange2> >,
																		typename boost::enable_if< boost::mpl::and_< can_integrate< TDomain, TRange1, TSpec1 >,
																																								 can_integrate< TDomain, TRange2, TSpec2 > > >::type >  {
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
                                    BinOp<TSpec1, TSpec2, SubOp<TRange1,TRange2> >,
																		typename boost::enable_if< boost::mpl::and_< can_integrate< TDomain, TRange1, TSpec1 >,
																																								 can_integrate< TDomain, TRange2, TSpec2 > > >::type >  {
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


template <typename TDomain, typename TSpec>
struct result_of_indefiniteIntegral<TDomain,double, UnaryOp<TSpec,Exp>,
																		typename boost::enable_if< boost::is_same<
														 typename SpecType< typename result_of_differentiate<TDomain,double,TSpec>::type >::type,
														 Const<>
																																 > >::type
																		> {
	 typedef typename result_of_differentiate<TDomain, double, TSpec>::type deriv_type;
	 typedef typename DivType< Function< TDomain, double, Const< CompileTime< 1 > > >, deriv_type >::type fctr_type;
	 typedef Function<TDomain,double,TSpec> exponent_type;
	 typedef Function<TDomain,double, UnaryOp<TSpec,Exp> > orig_func_type;
	 typedef typename MultType< fctr_type, orig_func_type >::type type;
};


template <typename TDomain, typename TSpec>
typename result_of_indefiniteIntegral<TDomain,double, UnaryOp< TSpec, Exp > >::type
//typename result_of_indefiniteIntegral<TDomain,double, UnaryOp< TSpec, Exp > >::type
indefiniteIntegral( const Function< TDomain, double, UnaryOp< TSpec, Exp > >& f ) {
	BOOST_AUTO_TPL( d, differentiate( f.getFunction() ) );
//	typedef BOOST_TYPEOF_TPL( d ) d_type;
//	typedef Const<> const_t;
//	BOOST_MPL_ASSERT(( boost::is_same< typename SpecType<d_type>::type, const_t ));
	return ( Function< TDomain, double, Const< CompileTime<1> > >() / d ) * f;
}


//
// *** Integrating piecewise functions
//


template <typename TDomain, typename TRange, typename TPieceSpec >
struct result_of_indefiniteIntegral<TDomain, TRange, Piecewise< TPieceSpec >,
																		typename boost::enable_if< can_integrate< TDomain, TRange, TPieceSpec > >::type
> {
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

template <typename TSpec1, typename TSpec2, typename TRange> struct Compose;
template <typename TSpec1, typename TSpec2, typename TRange>
struct IsFunctionSpec< Compose< TSpec1, TSpec2, TRange >,
											 typename boost::enable_if< boost::mpl::and_< IsFunctionSpec< TSpec1 >,
																																		IsFunctionSpec< TSpec2 > > >::type >:
		 public boost::true_type {};

template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec1, typename TSpec2>
class Function<TDomain, TRange2, Compose<TSpec1, TSpec2, TRange1> >:
		 public boost::compressed_pair< Function<TRange1, TRange2, TSpec1>, Function<TDomain, TRange1, TSpec2> > {
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec1> ));
   BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec2> ));
public:
   typedef TDomain argument_type;
   typedef TRange2 result_type;
   typedef Compose<TSpec1,TSpec2,TRange1> spec_type;
   typedef TRange1 range1_type;
   typedef TRange2 range2_type;
   typedef TSpec1 spec1_type;
   typedef TSpec2 spec2_type;
   
   typedef Function<TRange1, TRange2, TSpec1> function_type_1;
   typedef Function<TDomain, TRange1, TSpec2> function_type_2; 

   Function( const function_type_1& f1_, const function_type_2& f2_ ):
     boost::compressed_pair< function_type_1, function_type_2 >( f1_, f2_ ) { }

   result_type operator()( TDomain x ) const { return eval( this->first(),
																														eval( this->second(), x ) ); }
	 
   friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
     s << "(" << f.first() << " o " << f.second() << ")"; return s;
   }
};

template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec1, typename TSpec2>
struct result_of_compose {
	 typedef Function<TDomain,TRange2,Compose<TSpec1,TSpec2,TRange1> > type;
	 type operator()( const Function<TRange1,TRange2,TSpec1>& f1, const Function<TDomain,TRange1,TSpec2>& f2 ) {
		 return type( f1, f2 );
	 }
};

template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec1, typename TSpec2>
typename result_of_compose<TDomain,TRange1,TRange2,TSpec1,TSpec2>::type
compose( const Function<TRange1,TRange2,TSpec1>& f1, const Function<TDomain,TRange1,TSpec2>& f2 ) {
	return result_of_compose<TDomain,TRange1,TRange2,TSpec1,TSpec2>()( f1, f2 );
}

#if 0
template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec2>
struct result_of_compose<TDomain,TRange1,TRange2,X_To<1>,TSpec2> {
	 typedef Function<TRange1,TRange2,X_To<1> > f1_t;
	 typedef Function<TDomain,TRange1,TSpec2> f2_t;
	 typedef BOOST_TYPEOF_TPL(( fn_const<TDomain>( boost::declval<f1_t>().getFactor() ) * boost::declval<f2_t>() )) type;
};

template <typename TDomain, typename TRange1, typename TRange2,
          typename TSpec2>
typename result_of_compose<TDomain,TRange1,TRange2, X_To<1>,TSpec2>::type
compose( const Function<TRange1,TRange2,X_To<1> >& f1, const Function<TDomain,TRange1,TSpec2>& f2 ) {
	return fn_const<TDomain>( f1.getFactor() ) * f2;
}
#endif

#if 0
template <typename TDomain, typename TFactor, typename TRange2, typename TSpec2>
struct result_of_compose<TDomain, typename MultType<TFactor,TDomain>::type, TRange2,
												 X_To<1,TFactor>, TSpec2 > {
	 typedef Function<TDomain, typename MultType<TDomain,TFactor>::type,X_To<1,TFactor> > f1_t;
	 typedef Function< typename MultType<TDomain,TFactor>::type,TRange2,TSpec2>	f2_t;
	 typedef BOOST_TYPEOF_TPL( boost::declval<f1_t>().getFactor() * boost::declval<f2_t>() );
};


template <typename TDomain, typename TRange1, typename TFactor, typename TRange2, typename TSpec2>
typename result_of_compose<TDomain, TRange1, TRange2,
													 X_To<1,TFactor>, TSpec2>::type
compose( const Function<TDomain, TRange1,X_To<1,TFactor> >& f1,
				 const Function< typename MultType<TDomain,TFactor>::type,TRange2,TSpec2>& f2 ) {
	return f1.getFactor() * f2;
}
#endif

// ** Evaluating inverse functions

template <typename TDomain, typename TRange, typename TFactor>
struct result_of_inverse<TDomain, TRange, X_To<1, TFactor> > {
	 typedef typename InvType<TFactor>::type TFactorInv;
	 typedef Function< TRange, TDomain, X_To<1, TFactorInv> >  type;
	 type operator()( const Function< TDomain, TRange, X_To<1, TFactor> >& f ) {
		 return type( 1.0 / f.getFactor() );
	 }
	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<type>::type > ));	 
};

namespace detail {
template < template <typename,typename> class Op> struct slv;
struct Arg1 {};
struct Arg2 {};
template <typename T1, typename T2> const T1& aux_get( Arg1, const boost::compressed_pair<T1,T2>& p )
{ return p.first(); }
template <typename T1, typename T2> const T2& aux_get( Arg2, const boost::compressed_pair<T1,T2>& p )
{ return p.second(); }

template <> struct slv< AddOp > {
	 template <typename T2, typename TR>
	 static typename SubOp<TR,T2>::result_type slv_for(Arg1, T2 v2, TR vr  ) { return vr - v2; }
	 template <typename T1, typename TR>
	 static typename SubOp<TR,T1>::result_type slv_for(Arg2, T1 v1, TR vr  ) { return vr - v1; }
};
template <> struct slv< MultOp > {
	 template <typename T2, typename TR>
	 static typename DivOp<TR,T2>::result_type slv_for(Arg1, T2 v2, TR vr  ) { return vr / v2; }
	 template <typename T1, typename TR>
	 static typename DivOp<TR,T1>::result_type slv_for(Arg2, T1 v1, TR vr  ) { return vr / v1; }
};
template <> struct slv< SubOp > {
	 // v1 - v2 = vr
	 template <typename T2, typename TR>
	 static typename AddOp<T2,TR>::result_type slv_for(Arg1, T2 v2, TR vr  ) { return v2 + vr; }
	 template <typename T1, typename TR>
	 static typename SubOp<T1,TR>::result_type slv_for(Arg2, T1 v1, TR vr  ) { return v1 - vr; }
};
template <> struct slv< DivOp > {
	 // v1 / v2 = vr
	 template <typename T2, typename TR>
	 static typename MultOp<T2,TR>::result_type slv_for(Arg1, T2 v2, TR vr  ) { return v2 * vr; }
	 template <typename T1, typename TR>
	 static typename DivOp<T1,TR>::result_type slv_for(Arg2, T1 v1, TR vr  ) { return v1 / vr; }
};

template <typename TDomain, typename TRange, typename TRangeC, typename TRangeF, typename TSpec,
					template <typename,typename> class Op, typename ArgC, typename ArgF, typename bin_op_t>
struct result_of_inverse_aux {
	 typedef Function< TDomain, TRange, bin_op_t > f_t;
	 typedef typename result_of_inverse<TDomain, TRangeF, TSpec>::type f_inv_t;

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<f_inv_t>::type, TRangeF > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<f_inv_t>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<f_inv_t>::type > ));
	 
	 typedef Function<TRange,TRange,X_To<1> > x_t;
	 typedef Function<TRange,TRangeC,Const<> >	 c_t;
	 typedef BOOST_TYPEOF_TPL(( slv<Op>::slv_for( boost::declval<ArgF>(),
																								boost::declval<c_t>(), boost::declval<x_t>() ) ))
	 r_t;

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<r_t>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<r_t>::type, TRangeF > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<r_t>::type > ));

	 typedef BOOST_TYPEOF_TPL(( compose( boost::declval< f_inv_t >(), boost::declval< r_t >()
																))) type;

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<type>::type > ));
	 
	 type operator()( const f_t& f ) {
		 return compose( inverse( aux_get( ArgF(), f) ),
										 slv<Op>::slv_for( ArgF(), c_t( evalConst( aux_get( ArgC(), f ) ) ), x_t() ) );
	 }
};
}  // namespace detail
template <typename TDomain, typename TRange, typename TRangeC, typename TRangeF, typename TSpec,
					template <typename,typename> class Op>
struct result_of_inverse<TDomain, TRange, BinOp< Const<>, TSpec, Op< TRangeC, TRangeF > > >:
		 public detail::result_of_inverse_aux<TDomain,TRange, TRangeC, TRangeF, TSpec, Op, detail::Arg1, detail::Arg2,
																					BinOp< Const<>, TSpec, Op< TRangeC, TRangeF > > > {
};

template <typename TDomain, typename TRange, typename TRangeC, typename TRangeF, typename TSpec,
					template <typename,typename> class Op>
struct result_of_inverse<TDomain, TRange, BinOp< TSpec, Const<>, Op< TRangeF, TRangeC > > >:
		 public detail::result_of_inverse_aux<TDomain,TRange, TRangeC, TRangeF, TSpec, Op, detail::Arg2, detail::Arg1,
																					BinOp< TSpec, Const<>, Op< TRangeF, TRangeC > > > {
};

template <typename TDomain, typename TRange, typename TSpec>
typename result_of_inverse<TDomain, TRange, TSpec>::type
inverse( const Function<TDomain, TRange, TSpec>& f ) {
	typedef typename result_of_inverse<TDomain, TRange, TSpec>::type type;
	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<type>::type > ));
	
	return result_of_inverse<TDomain,TRange,TSpec>()( f );
}

template <typename TDomain, typename TRange, typename TSpec>
struct result_of_inverse<TDomain, TRange, UnaryOp< TSpec, Exp > > {
	 BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
	 BOOST_MPL_ASSERT(( boost::is_same< TRange, double > ));
	 typedef Function< TDomain, TRange, UnaryOp< TSpec, Exp > > f_t;
	 typedef Function< TDomain, TRange, TSpec > f_exponent_t;

	 typedef BOOST_TYPEOF_TPL(( boost::declval<f_t>().getFunction() )) f_exponent_t2;
	 
	 BOOST_MPL_ASSERT(( boost::is_same< f_exponent_t, f_exponent_t2 > ));
										 
	 typedef BOOST_TYPEOF_TPL( inverse( boost::declval<f_exponent_t>() ) ) f_inv_t;

	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<f_inv_t>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<f_inv_t>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<f_inv_t>::type > ));
	 
	 typedef BOOST_TYPEOF_TPL( compose( boost::declval<f_inv_t>(),
																			log_( fn_x<TRange,TRange>() ) ) ) type;
	 
	 BOOST_MPL_ASSERT(( boost::is_same< typename DomainType<type>::type, TRange > ));
	 BOOST_MPL_ASSERT(( boost::is_same< typename RangeType<type>::type, TDomain > ));
	 BOOST_MPL_ASSERT(( IsFunctionSpec< typename SpecType<type>::type > ));


	 type operator()( const f_t& f ) { return compose( inverse( f.getFunction() ), log_( fn_x<TRange,TRange>() ) ); }
};

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
             typename DiffType<TRange>::type eps, boost::uint32_t maxSteps = 100000 ) {
  BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
	assert( isStrictlyIncreasing( f, a, b ) );
  boost::uint32_t nsteps = 0;
  //PRINT5( a, b, targetVal, eps, typeid(f).name() );
	using util::chkCond;
  chkCond( !(boost::math::isnan)( targetVal ), "evalInverse: targetVal is NaN" );
  while ( true ) {
		chkCond( !(boost::math::isnan)( a ), "evalInverse: a is NaN" );
		chkCond( !(boost::math::isnan)( b ), "evalInverse: b is NaN" );
		chkCond( a <= b, "evalInverse: a <= b is false" );

    TDomain mid = a + ( b-a ) / 2;
    typedef typename DiffType<TRange>::type range_diff_type;
		TRange curVal = eval( f, mid );
		chkCond( !(boost::math::isnan)( curVal ), "evalInverse: eval returned NaN" );
		BOOST_AUTO_TPL( curDiff, curVal - targetVal );
		chkCond( !(boost::math::isnan)( curDiff ), "evalInverse: curDiff is NaN" );
    BOOST_AUTO_TPL( curDiff_abs, cosi_fabs( curVal - targetVal ) );
		chkCond( !(boost::math::isnan)( curDiff_abs ), "evalInverse: curDiff_abs is NaN" );
		typedef BOOST_TYPEOF_TPL( curDiff_abs ) curDiff_abs_t;
		if ( !( curDiff_abs >= curDiff_abs_t(0.) ) ) {
			std::cerr << "findInfimum: curDiff_abs negative; eps=" << eps << " curDiff=" << curDiff <<
				 " nsteps=" << nsteps << " maxSteps=" << maxSteps << " a=" << a << " b=" << b <<
				 " targetVal=" << targetVal << " curVal=" << curVal <<
				 " curDiff_abs >= curDiff_abs_t(0.) = " << ( curDiff_abs >= curDiff_abs_t(0.) ) << " f=" << f;
		}
		chkCond( curDiff_abs >= curDiff_abs_t(0.), "evalInverse: curDiff is negative" );
		// check that it is not NAN
		
    //PRINT5( a, b, mid, targetVal, curDiff );
    if ( curDiff_abs < eps ) {
      return mid;
    }
    if ( nsteps++ > maxSteps ) {
			std::ostringstream msg;
			msg.precision(16);
			msg << "findInfimum: too many steps; eps=" << eps << " curDiff=" << curDiff <<
				 " nsteps=" << nsteps << " maxSteps=" << maxSteps << " a=" << a << " b=" << b <<
				 " targetVal=" << targetVal << " f=" << f;
			throw std::runtime_error( msg.str() );
		}
    assert( a < b );
    assert( eval( f, a ) < eval( f, mid ) );
    assert( eval( f, mid ) < eval( f, b ) );
    ( curDiff < static_cast<range_diff_type>( 0.0 ) ? a : b ) = mid;
  }
}


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
struct result_of_indefiniteIntegral< TDomain, TRange, Any<>,
																		 typename boost::enable_if< has_area< TDomain, TRange > >::type > {
	 typedef typename AreaType<TDomain,TRange>::type integral_type;
	 typedef Function< TDomain, integral_type, Any<> > type;
};

namespace detail_fn_any {

template <typename TDomain, typename TRange, typename Enable = void>
struct FunctionObjectConcept {
	 virtual ~FunctionObjectConcept() {}
	 virtual TRange doEval( TDomain ) const = 0;
	 virtual Function< TRange, TDomain, Any<> > doInverse() const = 0;
	 virtual void doOutput( std::ostream& ) const = 0;
};

template <typename TDomain, typename TRange>
struct FunctionObjectConcept<TDomain, TRange, typename boost::enable_if< has_area< TDomain, TRange > >::type > {
	 virtual ~FunctionObjectConcept() {}
	 virtual TRange doEval( TDomain ) const = 0;
	 virtual Function< TRange, TDomain, Any<> > doInverse() const = 0;
	 virtual void doOutput( std::ostream& ) const = 0;	 
	 virtual Function< TDomain, typename AreaType< TDomain, TRange >::type, Any<> >
	 do_indefiniteIntegral() const = 0;
};

template< typename TDomain, typename TRange, typename TSpec, typename Enable=void>
class FunctionObjectModel:
		 public FunctionObjectConcept< TDomain, TRange > {
	 BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
public:
	 FunctionObjectModel( const Function<TDomain, TRange, TSpec>& f_ ) : f( f_ ) {}
	 virtual ~FunctionObjectModel() {}
	 virtual TRange doEval( TDomain x ) const { return eval( f, x ); }
	 virtual Function< TRange, TDomain, Any<> > doInverse() const { return inverse( f ); }
	 virtual void doOutput( std::ostream& s ) const { s << f; }
	 
protected:
	 Function<TDomain, TRange, TSpec> f;
};

template <typename TDomain, typename TRange, typename TSpec, typename Enable=void> struct do_intgrl {
	 static Function< TDomain, typename AreaType< TDomain, TRange >::type, Any<> > get_integral( const Function< TDomain, TRange, TSpec >& ) {
		 return Function< TDomain, typename AreaType< TDomain, TRange >::type, Any<> >();
	 }
};

template <typename TDomain, typename TRange, typename TSpec>
struct do_intgrl< TDomain, TRange, TSpec, typename boost::enable_if< can_integrate< TDomain, TRange, TSpec > >::type > {
	 static Function< TDomain, typename AreaType< TDomain, TRange >::type, Any<> > get_integral( const Function< TDomain, TRange, TSpec >& f ) {
		 return indefiniteIntegral( f );
	 }
};


template< typename TDomain, typename TRange, typename TSpec>
class FunctionObjectModel< TDomain, TRange, TSpec, typename boost::enable_if< has_area< TDomain, TRange > >::type  >:
		 public FunctionObjectConcept< TDomain, TRange > {
	 BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
public:
	 FunctionObjectModel( const Function<TDomain, TRange, TSpec>& f_ ) : f( f_ ) {}
	 virtual ~FunctionObjectModel() {}
	 virtual TRange doEval( TDomain x ) const { return eval( f, x ); }
	 virtual Function< TRange, TDomain, Any<> > doInverse() const { return inverse( f ); }
	 virtual void doOutput( std::ostream& s ) const { s << f; }

	 
	 virtual Function< TDomain, typename AreaType< TDomain, TRange >::type, Any<> >
	 do_indefiniteIntegral() const { return do_intgrl<TDomain,TRange,TSpec>::get_integral( f ); }

protected:
	 Function<TDomain, TRange, TSpec> f;
};

}  // namespace detail_fn_any

template <typename TDomain, typename TRange>
class Function<TDomain, TRange, Any<> > {
public:
   typedef TDomain argument_type;
   typedef TRange result_type;
   typedef Any<> spec_type;

public:
   boost::shared_ptr< detail_fn_any::FunctionObjectConcept<TDomain, TRange> > object;

public:
   Function()  { }
   template< typename TSpec > Function( const Function<TDomain, TRange, TSpec>&  obj ) :
     object( new detail_fn_any::FunctionObjectModel<TDomain, TRange, TSpec>( obj ) ) {}

   template <typename TSpec>
   Function& operator=( const Function<TDomain, TRange, TSpec>& obj ) {
     reset( obj );
		 return *this;
   }

   template <typename TSpec>
   void reset( const Function<TDomain, TRange, TSpec>& obj ) {
     object.reset( new detail_fn_any::FunctionObjectModel<TDomain, TRange, TSpec>( obj ) );
   }

   bool empty() const { return !object.get(); }
   operator bool() const { return !empty(); }

   TRange operator()( TDomain x ) const { assert( object.get() ); return object->doEval( x ); }

	 friend std::ostream& operator<<( std::ostream& s, const Function& f ) {
		 assert( f.object.get() );
		 f.object->doOutput( s );
		 return s;
	 }
	 
};  // end: type erasure of a Function

template <typename TDomain, typename TRange>
	 typename boost::lazy_enable_if< has_area< TDomain, TRange  >,
																	 result_of_indefiniteIntegral< TDomain, TRange, Any<> > >::type
indefiniteIntegral( const Function< TDomain, TRange, Any<> >& f ) {
		 return f.object->do_indefiniteIntegral();
}

template <typename TDomain, typename TRange>
struct result_of_inverse<TDomain, TRange, Any<> > {
	 typedef Function<TRange, TDomain, Any<> > type;
	 type operator()( const Function< TDomain, TRange, Any<> >& f ) { return f.object->doInverse(); }
};

template <typename TDomain, typename TRange, typename TSpec> inline
Function< TDomain, TRange, Any<> >
fn_any( Function< TDomain, TRange, TSpec > const& f ) { return Function< TDomain, TRange, Any<> >( f ); }


struct Bad;
template <> struct IsFunctionSpec< Bad >: public boost::true_type {};
template <typename TDomain, typename TRange>
class Function<TDomain, TRange, Bad > {
public:
	 typedef TDomain argument_type;
	 typedef TRange result_type;
	 typedef Bad spec_type;

	 Function() { }
	 TRange operator()( TDomain ) const { return std::numeric_limits<TRange>::quiet_NaN(); }
	 friend std::ostream& operator<<( std::ostream& s, const Function& ) {
		 s << "Bad()";
		 return s;
	 }
};

template <typename TDomain, typename TRange, typename TSpec> struct result_of_inverse {
 	 typedef Function<TRange,TDomain,Bad> type;
	 type operator()( const Function< TDomain, TRange, TSpec >& f ) {
		 BOOST_MPL_ASSERT(( IsFunctionSpec<TSpec> ));
		 //BOOST_STATIC_ASSERT_MSG( sizeof( f ) == 0, "don't know how to invert f" );
		 throw std::invalid_argument( "cannot compute inverse for " +
																	boost::core::demangle( typeid( f ).name() ) );
		 return type();
	 }
};


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
			virtual std::ostream& print( std::ostream& s ) const = 0;
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
			virtual std::ostream& print( std::ostream& s ) const { s << proc; return s; }
			virtual std::string doGetLabel() const { return proc.getLabel(); }
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
	 std::ostream& print( std::ostream& s ) const { return object->print( s ); }

	 std::string getLabel() const { return object->doGetLabel(); }
   
   friend std::ostream& operator<<( std::ostream& s, const ArrivalProcess& f ) {
		 return f.print( s );
	 }
//	 friend std::string getLabel<>( const ArrivalProcess& p ); // { return p.object->doGetLabel(); }
}; // class ArrivalProcess<TTime, Any< TRand > >

// template <typename TTime, typename TRand>
// inline std::string getLabel( const ArrivalProcess<TTime, Any<TRand> >& p ) {
//   return p.object->doGetLabel();
//  }


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
		 PRINT7( label, fromTime, maxTime, rateFactor, eps, integralAtFromTime, integralAtMaxTime );
		 if ( integralAtMaxTime == integralAtFromTime ) return maxTime;
#ifndef NDEBUG
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
   
   friend std::ostream& operator<<( std::ostream& s, const ArrivalProcess& f ) {
		 s << "ArrivalProcess[label=" << f.label << ", rateFunctionIntegral=" << f.rateFunctionIntegral << "]";
		 return s;
	 }

	 std::string getLabel() const { return label; }

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
     x2y.insert( x2y.end(), std::make_pair( x, y ) );
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

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
struct divAresult {
	 typedef Function< TDomain, TRange1, Const<> > f1_t;
	 typedef Function< TDomain, typename MultOp<TRange2,double>::result_type,
										 BinOp< Const<>,
														UnaryOp< TSpec, Exp >,
														MultOp< TRange2, double >
														>
										 > f2_t;
	 typedef BOOST_TYPEOF_TPL(( ( boost::declval<f1_t>() / boost::declval<f2_t>().first() ) *
															exp_( Function< TDomain, double, Const<> >( -1.0 ) *
																		boost::declval<f2_t>().second().getFunction() )) ) type;
};

template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
typename divAresult<TDomain,TRange1,TRange2,TSpec>::type
operator/( Function< TDomain, TRange1, Const<> > const& f1,
					 Function< TDomain, typename MultOp<TRange2,double>::result_type,
										 BinOp< Const<>,
														UnaryOp< TSpec, Exp >,
														MultOp< TRange2, double >
														>
					 > const& f2 ) {
	return ( f1 / f2.first() ) * exp_( Function< TDomain, double, Const<> >( -1.0 ) * f2.second().getFunction() );
}

// template <typename TDomain, typename TRange1, typename TRange2, typename TSpec>
// typename divAresult<TDomain, TRange1, TRange2, TSpec>::type
// operator/( Function< TDomain, TRange1, Const<> > const& f1,
// 					 Function< TDomain, typename MultOp<TRange2,double>::result_type,
// 										 BinOp< Const<>,
// 					                  UnaryOp< TSpec, Exp >,
// 														MultOp< TRange2, double >
// 														>
// 										 > const& f2
// 	) {
// 	return ( f1 / f2.first() ) * exp_( -f2.second().getFunction() ) ;
// }



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
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::detail::Arg1);
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::detail::Arg2);
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::detail::slv,BOOST_TYPEOF_TEMPLATE(2))
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::detail::result_of_inverse_aux,
															 (typename)(typename)(typename)(typename)(typename)BOOST_TYPEOF_TEMPLATE(2)
															 (typename)(typename)(typename))

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::BinOp,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::AddOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::SubOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::MultOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::DivOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Compose,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Piecewise,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::AreaType,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::detail::result_of_make_line_through,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::UnaryOp,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_indefiniteIntegral,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_inverse,3)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_integralFromPoint,3)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::Any)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::ArrivalProcess,2);
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Poisson,2)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::CinlarsMethod)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::NonHomogeneous,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::result_of_makeNonHomogeneousPoissonProcess,3)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::InterpFn,2)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::InterpBiFun,2)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Exp,1)
BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::Log,1)

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::math::CVal,1)
BOOST_TYPEOF_REGISTER_TYPE(cosi::math::Bad)

#endif  // #ifndef __INCLUDE_COSI_GENERALMATH_H
// Postamble:1 ends here
