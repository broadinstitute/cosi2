#ifndef COSI_INCLUDE_TYPEUTIL_H
#define COSI_INCLUDE_TYPEUTIL_H

#include <boost/utility/declval.hpp>

namespace cosi {

// Metafunction: DiffType
// Returns the type of the difference of values of the given type.
template <typename TVal> struct DiffType {
   typedef BOOST_TYPEOF_TPL( boost::declval<TVal>() - boost::declval<TVal>()) type;
};


}  // namespace cosi

#endif // #ifndef COSI_INCLUDE_TYPEUTIL_H
