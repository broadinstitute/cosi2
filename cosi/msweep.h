#ifndef COSI_INCLUDE_MSWEEP_H
#define COSI_INCLUDE_MSWEEP_H

#include <vector>
#include <cosi/general/utils.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>

namespace cosi {

// * Class: MSweepP - implements functionality related to selected sweeps.

// ** Ctor: make_MSweep - construct an MSweep object.
MSweepP make_MSweep( DemographyP, BaseModelP, RandGenP );

// ** Fn: as_Module - cast this MSweep object to a Module
ModuleP as_Module( MSweepP );

// ** Fn: getSweepModel - returns the [[sweep model]].
BaseModelP getSweepModel( MSweepP );

// ** Fn: computeLeavesInfo - compute the proper order of [[leaves]] in the output.
//
LeavesInfoP
computeLeavesInfo( MSweepP );

// ** Fn: addSelMut - add to the list of mutations the selected mutation.
void addSelMut( MSweepP, MutlistP );

void setSweepInfo( BaseModel& baseModel,
									 genid selGen, double selCoeff, loc_t selPos, popid selPop,
									 util::ValRange<freq_t> final_sel_freq, popid selBegPop, genid selBegGen );

}  // namespace cosi


#endif // #ifndef COSI_INCLUDE_MSWEEP_H
