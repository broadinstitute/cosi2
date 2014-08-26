//
// Header: sweep1.h
//
// Simulation of selective sweep in a single population, by reusing most of the
// neutral-coalescent-simulation machinery.
//

#ifndef COSI_INCLUDE_SWEEP2_H
#define COSI_INCLUDE_SWEEP2_H

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <cosi/defs.h>
#include <cosi/decls.h>

namespace cosi {

namespace sweep2 {

using std::vector;
using std::pair;
using std::istream;

void sweep2_setTrajFN( filename_t fname );

void sweep2_set_sweepFracSample( bool sweepFracSample_ );

// Function: Event_SweepOnePop_typeStr
// Returns the string which identifies a single-population sweep event in a cosi parameter file;
// the string that comes after pop_event.
const char *Event_SweepOnePop2_typeStr();

}  // namespace sweep2

}  // namespace cosi

#endif  // #ifndef COSI_INCLUDE_SWEEP2_H
