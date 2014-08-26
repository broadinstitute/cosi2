#ifndef __COSI_INCLUDE_SEGLISTFWD_H
#define __COSI_INCLUDE_SEGLISTFWD_H

namespace cosi {
namespace seglist2 {
namespace seglist_detail {
struct seglist;
typedef struct seglist Seglist;

struct seg;
typedef struct seg Seg;

struct segptr;
typedef struct segptr Segptr;

typedef Segptr Finger_elt;
typedef Finger_elt Finger;

}  // namespace seglist_detail
using seglist_detail::Seglist;
using seglist_detail::Seg;
using seglist_detail::Finger;
}  // namespace seglist

}  // namespace cosi

#endif  // #ifndef __COSI_INCLUDE_SEGLISTFWD_H
