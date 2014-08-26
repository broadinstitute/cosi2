/*
 * Header: leafset.h
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 *
 * This file chooses one of several leafset implementation based on preprocessor directives:
 *
 *    COSI_LEAFSET_TREE - representing leafsets as trees of nodes
 *    COSI_LEAFSET_BITSET - representing leafsets as bitsets
 *    COSI_LEAFSET_SIZEONLY - keeping track of only the size of each leafset
 */

#ifndef __INCLUDE_COSI_LEAFSET_H
#define __INCLUDE_COSI_LEAFSET_H

#if ( defined(COSI_LEAFSET_TREE) + defined(COSI_LEAFSET_BITSET) + defined( COSI_LEAFSET_SIZEONLY ) ) > 1
#error "Must choose at most one leafset implementation"
#endif

#if ( defined(COSI_LEAFSET_TREE) + defined(COSI_LEAFSET_BITSET) + defined( COSI_LEAFSET_SIZEONLY ) ) == 0
#define COSI_LEAFSET_TREE
#endif

#if defined(COSI_LEAFSET_BITSET)
#include <cosi/leafset-bitset.h>
#elif defined(COSI_LEAFSET_SIZEONLY)
#include <cosi/leafset-sizeonly.h>
#elif defined(COSI_LEAFSET_TREE)
#include <cosi/leafset-tree.h>
#else
#error "Must choose a leafset implementation"
#endif
  
#endif
// #ifndef __INCLUDE_COSI_LEAFSET_H
