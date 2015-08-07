#ifndef INCLUDE_OUTPUT_COAL_DATA_COSI_H
#define INCLUDE_OUTPUT_COAL_DATA_COSI_H

#include <cosi/defs.h>
#include <cosi/decls.h>

namespace cosi {

namespace customstats {

void init( DemographyP demography, size_t nsims_, int seqlen_ );

void record_sim(DemographyP demography, GenMapP genMap, len_bp_int_t length, MutlistP mutlist,
								bool_t inf_sites );

void finish();

}  // namespace customstats

}  // namespace cosi

#endif // #ifndef INCLUDE_OUTPUT_COAL_DATA_COSI_H
