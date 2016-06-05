#ifndef INCLUDE_COSI_NTHEVT_H
#define INCLUDE_COSI_NTHEVT_H

#include  <cosi/defs.h>

namespace cosi {
namespace nthevt {

enum evt_kind_t { EVT_NONE, EVT_COAL, EVT_RECOMB, EVT_MIGR };

void record_evt( evt_kind_t, genid );
bool skipRest();

}  // namespace nthevt
}  // namespace cosi


#endif
//#ifndef INCLUDE_COSI_NTHEVT_H

