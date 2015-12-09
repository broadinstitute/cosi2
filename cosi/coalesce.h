//
// Header: coalesce.h
//
// Keeping track of coalescence probabilities, and choosing coalescent events.

/* $Id: coalesce.h,v 1.2 2011/05/03 18:50:54 sfs Exp $ */
#ifndef __COSI_INCLUDE_COALESCE_H
#define __COSI_INCLUDE_COALESCE_H

#include <functional>
#include <limits>
#include <vector>
#include <boost/random/discrete_distribution.hpp>
#include <cosi/decls.h>
#include <cosi/cosirand.h>
#include <cosi/arrproc2.h>
#include <cosi/basemodel.h>
#include <cosi/generalmath.h>

namespace cosi {

namespace coal {

//
// Class: Coalesce
//
// Keeps track of coalescence probabilities, and choosese coalescent events.
//
class Coalesce: public HasRandGen {
public:
	 Coalesce( DemographyP demography_ );

	 
	 double coalesce_get_rate (void) const;
	 
	 pop_idx_t coalesce_pick_popindex () const;

	 gens_t coalesce_get_wait_time_nonhomog( genid, gens_t /* maxWaitTime */ ) const;
	 pop_idx_t coalesce_pick_popindex_nonhomog() const;

	 void setBaseModel( BaseModelP baseModel_ ) { baseModel = baseModel_; }

	 typedef
	 arrival2::ArrivalProcess< genid,
														 arrival2::Stoch< RandGen,
																							arrival2::Compound<
																								arrival2::Poisson<
																									math::Piecewise< math::Any<> >, nchromPairs_float_t > > > >
	 coal_processes_type;

	 boost::shared_ptr< coal_processes_type >
	 createCoalProcesses();
	 
	 
private:
	 DemographyP demography;
	 BaseModelP baseModel;
	 mutable double lastrate;
	 mutable std::vector<double> popRates;
	 mutable pop_idx_t nonhomogPopIdx;
	 boost::random::discrete_distribution<> dist;	 

};  // class Coalesce

} // namespace coal

}  // namespace cosi

#endif
// #ifndef __COSI_INCLUDE_COALESCE_H
