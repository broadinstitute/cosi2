#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

#include <vector>
#include <cosi/decls.h>
#include <cosi/basemodel.h>
#include <cosi/leafset.h>

namespace cosi {
MSweepP make_MSweep( DemographyP, BaseModelP, RandGenP );
BaseModelP getSweepModel( MSweepP );
boost::shared_ptr< std::vector< leaf_id_t > >
computeLeafOrder( MSweepP );
void addSelMut( MSweepP, MutlistP );
}  // namespace cosi


#endif // #ifndef COSI_INCLUDE_MSWEEP_H
