//
// Header: mutcontext.h
//
// Code for computing mutcontext neighborhood on which a mutation appears.
//

#ifndef __INCLUDE_COSI_MUTCONTEXT_H
#define __INCLUDE_COSI_MUTCONTEXT_H

#include <vector>
#include <map>
#include <utility>
#include <boost/multi_array.hpp>
#include <cosi/defs.h>
#include <cosi/seglist.h>

namespace cosi {
namespace mutcontext {

using std::vector;
using std::map;
using seglist::Seglist;

//
// Function: computeMutContexts
//
// For a mutation placed at location 'loc' of 'seglist', for leaves in the leafset of the seg
// of 'seglist' containing 'loc', return the boundaries of the mutcontext neighborhood
// around that mutation at that leaf.  In other words: the mutation was born on a specific
// ancestral chromosome; this shows the part of that ancestral chromosome 
// 
//
vector< BasicSeg_loc > computeMutContexts( const Seglist *seglist, loc_t loc );

void saveMutContexts( const Seglist *seglist, loc_t loc );

typedef map< loc_t, vector< BasicSeg_loc > > mutContexts_t;

const mutContexts_t& getSavedMutContexts();

}  // namespace mutcontext
}  // namespace cosi

#endif // #ifndef __INCLUDE_COSI_MUTCONTEXT_H
