#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cosi/general/math/generalmath.h>
#include <cosi/ehh.h>

namespace cosi {

namespace ehh {
  
namespace lam = ::boost::lambda;
namespace ad = ::boost::adaptors;
namespace ran = ::boost::range;
using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using ad::transformed;
using ad::sliced;
using std::make_pair;
using ran::push_back;



}  // namespace ehh

}  // namespace cosi
