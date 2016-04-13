#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

#include <vector>
#include <cosi/general/utils.h>
#include <cosi/decls.h>
#include <cosi/basemodel.h>
#include <cosi/leafset.h>

namespace cosi {

// * Class: MSweepP - implements functionality related to selected sweeps.

// ** Ctor: make_MSweep - construct an MSweep object.
MSweepP make_MSweep( DemographyP, GenMapP, BaseModelP, RandGenP );

// ** Fn: as_Module - cast this MSweep object to a Module
ModuleP as_Module( MSweepP );

// ** Fn: getSweepModel - returns the [[sweep model]].
BaseModelP getSweepModel( MSweepP );

// ** Fn: computeLeafOrder - compute the proper order of [[leaves]] in the output.
//
boost::shared_ptr< std::vector< leaf_id_t > >
computeLeafOrder( MSweepP );

// ** Fn: addSelMut - add to the list of mutations the selected mutation.
void addSelMut( MSweepP, MutlistP );

void setSweepInfo( BaseModel& baseModel,
									 genid selGen, double selCoeff, loc_t selPos, popid selPop,
									 util::ValRange<freq_t> final_sel_freq );

}  // namespace cosi


#endif // #ifndef COSI_INCLUDE_MSWEEP_H
