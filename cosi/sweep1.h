//
// Header: sweep1.h
//
// Simulation of selective sweep in a single population, by reusing most of the
// neutral-coalescent-simulation machinery.
//

#ifndef __COSI_INCLUDE_SWEEP1_H
#define __COSI_INCLUDE_SWEEP1_H

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <cosi/defs.h>
#include <cosi/decls.h>

namespace cosi {

namespace sweep1 {

using std::vector;
using std::pair;
using std::istream;

void sweep1_setTrajFN( filename_t fname );
void sweep1_set_deltaTfactor( double deltaTfactor_ );

// Function: Event_SweepOnePop_typeStr
// Returns the string which identifies a single-population sweep event in a cosi parameter file;
// the string that comes after pop_event.
const char *Event_SweepOnePop_typeStr();

}  // namespace sweep1

}  // namespace cosi

#endif  // #ifndef __COSI_INCLUDE_SWEEP1_H
