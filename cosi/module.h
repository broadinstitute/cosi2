#ifndef COSI_INCLUDE_MODULE_H
#define COSI_INCLUDE_MODULE_H

#include <string>
#include <map>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/program_options.hpp>
#include <cosi/decls.h>
#include <cosi/general/math/cosirand.h>

namespace cosi {

typedef boost::mpl::vector< Demography, MSweep > cosi_modules;

class Modules;


// * Class: Module - abstract base class for the different cosi modules
class Module: HasRandGen {
	 // Fn: defineParams - defines command-line params relevant to this module
	 virtual void defineParams( boost::program_options::options_description& ) { }

	 // Fn: init() - initialize the module
	 virtual void init() { }

	 // Fn: processParamFileLine - process one line of the parameter file
	 virtual void processParamFileLine( std::string ) { }

	 template <typename M>
	 boost::shared_ptr<M> getModule( Modules * );

	 Modules *modules;

};  // class Module

// namespace dtl {

// namespace fusion = boost::fusion;
// namespace mpl = boost::mpl;

// template<typename T>
// struct make_out_signature_pair
// {
// 	 typedef typename fusion::result_of::make_pair
//     <
// 		 T, boost::shared_ptr<T>
// 		 >::type type;
// };

// template<typename OutSignatures>
// struct module_map_impl
// {
// 	 typedef typename fusion::result_of::as_map
//     <
// 		 typename mpl::transform
//         <
// 		 OutSignatures,
// 			 make_out_signature_pair<mpl::_1>
// 			 >::type
// 		 >::type OutMap;
// };}  // namespace dtl



// // * Class: Modules - a collection of modules.
// class Modules {
// public:
// 	 typename dtl::module_map_impl< cosi_modules >::OutMap moduleMap;

// 	 template <typename M>
// 	 boost::shared_ptr<M> addModule( boost::shared_ptr<M> mptr ) {
// 		 boost::fusion::map<M>::at_key( moduleMap ) = mptr;
// 		 return mptr;
// 	 }

// };

// // template <typename M>
// // boost::shared_ptr<M> Module::getModule( Modules *mdls ) {
// // 	return boost::fusion::map<M>::at_key( mdls->moduleMap );
// // }


}  // namespace cosi

#endif // #ifndef COSI_INCLUDE_MODULE_H
