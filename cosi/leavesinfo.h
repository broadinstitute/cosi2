#ifndef COSI_INCLUDE_LEAVESINFO_H
#define COSI_INCLUDE_LEAVESINFO_H

#include <vector>
#include <iostream>
#include <boost/typeof/typeof.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace cosi {

struct LeavesInfo {
	std::vector< leaf_id_t > leafOrder;
	std::vector< nchroms_t > sampleSizes;
	std::vector< popid > popNames;
};

LeavesInfoP computeStdLeavesInfo( std::vector<nchroms_t> const& sampleSizes, 
																	std::vector<popid> const& popNames );

	
std::ostream& operator<<( std::ostream&, LeavesInfo const& );


}  // namespace cosi

BOOST_TYPEOF_REGISTER_TYPE(cosi::LeavesInfo)

#endif
// #ifndef COSI_INCLUDE_LEAVESINFO_H
