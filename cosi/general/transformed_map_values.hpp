// Boost.Range library
//
//  Copyright Thorsten Ottosen, Neil Groves 2006 - 2008. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// For more information, see http://www.boost.org/libs/range/
//

#ifndef COSI_RANGE_ADAPTOR_TRANSFORMED_MAP_VALUES_HPP
#define COSI_RANGE_ADAPTOR_TRANSFORMED_MAP_VALUES_HPP

#include <utility>
#include <boost/range/adaptor/argument_fwd.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/static_assert.hpp>
#include <boost/range/concepts.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost
{

template <typename MapValueType, typename TValueTransformer>
struct MapValueTransformer {


typedef MapValueType pair_t;
typedef typename MapValueType::first_type first_type;
typedef typename MapValueType::second_type second_type;
BOOST_STATIC_ASSERT(( boost::is_same< pair_t, std::pair< first_type, second_type > >::value ));

typedef typename boost::result_of<TValueTransformer(second_type)>::type result_type;
BOOST_CONCEPT_ASSERT((boost::UnaryFunction<TValueTransformer, result_type, pair_t>));

	 MapValueTransformer( TValueTransformer transformer_ ): transformer( transformer_ ) { }
	 
	 // result_type operator()( pair_t& p ) { return std::make_pair( p.first(),
	 // 																															transformer( p.second() ) ); }
	 result_type operator() ( const pair_t& p ) const
			{ return std::make_pair( p.first(),
			 												 transformer( p.second() ) ); }

	 TValueTransformer transformer;
};

template <typename UnaryFunction, typename InputRange> struct TransformedRangeType {
typedef typename boost::range_value<InputRange>::type range_value_type;
typedef MapValueTransformer< range_value_type, UnaryFunction > map_value_transformer_type;
typedef transformed_range<map_value_transformer_type, InputRange> type;
};

template <typename UnaryFunction, typename InputRange> struct TransformedRangeTypeConst {
typedef typename boost::range_value<const InputRange>::type range_value_type;
typedef MapValueTransformer< range_value_type, UnaryFunction > map_value_transformer_type;
typedef transformed_range<map_value_transformer_type, const InputRange> type;
};



namespace range_detail {
        template< class T >
        struct transform_map_values_holder : holder<T>
        {
            transform_map_values_holder( T r ) : holder<T>(r)
            { }
        };

        template< class InputRng, class UnaryFunctionType >
        inline typename TransformedRangeType<UnaryFunctionType,InputRng>::type
        operator|( InputRng& r,
                   const transform_map_values_holder<UnaryFunctionType>& f )
        {
typedef typename TransformedRangeType<UnaryFunctionType,InputRng>::type result_type;
return result_type(f, r);
        }

        template< class InputRng, class UnaryFunctionType >
        inline typename TransformedRangeTypeConst<UnaryFunctionType,InputRng>::type
        operator|( const InputRng& r,
                   const transform_map_values_holder<UnaryFunctionType>& f )
        {
typedef typename TransformedRangeTypeConst<UnaryFunctionType,InputRng>::type result_type;
return result_type(f, r);
        }


}

    namespace adaptors
    {



        namespace
        {
            const range_detail::forwarder<range_detail::transform_map_values_holder>
                    transformed_map_values =
                      range_detail::forwarder<range_detail::transform_map_values_holder>();
        }

        template<class UnaryFunctionType, class InputRange>
        inline typename TransformedRangeType<UnaryFunctionType,InputRange>::type
        transform_map_values(InputRange& rng, UnaryFunctionType fn)
        {
//BOOST_CONCEPT_ASSERT((boost::UnaryFunction<UnaryFunctionType>));
BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<InputRange>));
//STATIC_ASSERT( boost::is_convertible< typename boost::range_value<rng>::type
typedef typename TransformedRangeType<UnaryFunctionType,InputRange>::type result_type;
return result_type(fn, rng);
        }

        template<class UnaryFunctionType, class InputRange>
        inline
typename TransformedRangeTypeConst<UnaryFunctionType,InputRange>::type				
        transform_map_values(const InputRange& rng, UnaryFunctionType fn)
        {
//BOOST_CONCEPT_ASSERT((boost::UnaryFunction<UnaryFunctionType>));
BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<const InputRange>));
//STATIC_ASSERT( boost::is_convertible< typename boost::range_value<rng>::type
typedef typename TransformedRangeTypeConst<UnaryFunctionType,InputRange>::type result_type;
return result_type(fn, rng);

        }
    } // 'adaptors'

}

#endif
