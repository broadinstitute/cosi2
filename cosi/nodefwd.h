#ifndef __COSI_INCLUDE_NODEFWD_H
#define __COSI_INCLUDE_NODEFWD_H

#include <vector>
#include <boost/shared_ptr.hpp>

namespace cosi {
namespace node {
class Node;

//
// Type: NodeList
//
// A vector of <Nodes>.
//
typedef std::vector<Node *> NodeList;
class NodePool;
typedef boost::shared_ptr<node::NodePool> NodePoolP;

}
}

#endif // #ifndef __COSI_INCLUDE_NODEFWD_H 
