#ifndef __COSI_INCLUDE_NODEUTIL_H
#define __COSI_INCLUDE_NODEUTIL_H

#include <cosi/nodefwd.h>

namespace cosi {
namespace node {

int nodelist_add(NodeList *, Node *);
int nodelist_remove(NodeList *, Node *);
void nodelist_remove_idx(NodeList *, int);

inline Node * nodelist_get_node (int idx, const NodeList *nodes) { return (*nodes)[idx]; }

} // namespace node
} // namespace cosi

#endif // #ifndef __COSI_INCLUDE_NODEUTIL_H
