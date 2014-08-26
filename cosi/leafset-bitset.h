
/*
 * Header: leafset.h
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 */

#ifndef __INCLUDE_COSI_LEAFSET_BITSET_H
#define __INCLUDE_COSI_LEAFSET_BITSET_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/scoped_array.hpp>
#include <boost/intrusive_ptr.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/mempool.h>

namespace cosi {

/* Type: leaf_id_t */
/* Identifier of one leaf node of the ARG, representing a present-day <chrom>. */
typedef nodeid leaf_id_t;

#define NULL_LEAF_ID -1

/* Type: leafset_word_t */
/* Used for bitset representation of leafsets. */
typedef unsigned long leafset_word_t;

#define BITS_PER_LEAFSET_WORD ((sizeof(leafset_word_t) * BITS_PER_BYTE))

/* Struct: leafset_t */
/* A set of leaves, i.e. ARG nodes representing present-day chroms. */
/**/
/* One common use of leafsets: For each <seg> of each <Node>'s <Seglist>, */
/*	we store the set of leaves */
/* that inherit this seg -- represented as a leafset. */
struct leafset_struct: public refcounted<leafset_struct> {
private:
	 friend boost::intrusive_ptr<leafset_struct> leafset_alloc(void);

public:
	 /* Field: bits */
	 /* A bitset representation of the leaves in the leafset. */
	 leafset_word_t bits[];

	 void* operator new (size_t size);
	 void operator delete (void *p);

	 leafset_struct();
};

typedef leafset_struct leafset_t;

#define LEAFSET_NULL ((leafset_t *)NULL)

typedef boost::intrusive_ptr<leafset_t> leafset_p;

extern size_t leafset_nwords;
extern boost::scoped_array<leafset_word_t> leaf_set_mask, leaf_clear_mask;

/* FuncP: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
/* Must be called once, prior to the creation of any leafsets. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ );
leaf_id_t leafset_get_max_leaf_id(void);  

leafset_p make_empty_leafset(void);

/* FuncP: make_singleton_leafset */
/* Returns the leafset containing only the specified leaf. */
leafset_p make_singleton_leafset( leaf_id_t leaf );

/* FuncP: make_range_leafset */
/* Returns the leafset containing the given range of leaves. */
leafset_p make_range_leafset( leaf_id_t fromLeaf, leaf_id_t toLeaf );

void leafset_add( leafset_p leafset, leaf_id_t leaf );

/* FuncP: leafset_union */
/* Returns the leafset containing the union of the given leafsets. */
leafset_p leafset_union( leafset_p , leafset_p  );

/* FuncP: leafset_intersection */
/* Returns the leafset containing the intersection of the given leafsets. */
leafset_p leafset_intersection( leafset_p , leafset_p  );
bool_t leafset_intersects( leafset_p , leafset_p  );
inline bool_t leafset_disjoint( leafset_p s1, leafset_p s2 ) { return !leafset_intersects( s1, s2 ); }


leafset_p leafset_difference( leafset_p s1, leafset_p s2);

/* FuncP: leafset_complement */
/* Returns the leafset containing the complement of the given leafset. */
leafset_p leafset_complement( leafset_p  );

bool_t leafset_contains( leafset_p , leaf_id_t leaf );

nchroms_t leafset_cardinality( leafset_p  );
bool_t leafset_is_empty( leafset_p leafset );
std::string leafset_str( leafset_p  );

bool_t leafset_equal( leafset_p leafset1, leafset_p leafset2);
bool_t leafset_is_singleton( leafset_p leafset, leaf_id_t leaf );

void leafset_print( leafset_p leafset );

int leafset_num_alloced(void);

bool_t leafset_is_subset( leafset_p L1, leafset_p L2 );

leafset_p leafset_from_str( const char *s );

void leafset_write_bin( FILE *f, leafset_p  );
leafset_p leafset_read_bin( FILE *f );

#define FOR_LEAFSET(leafset,leaf_var,body) do {		\
		const leafset_word_t *bits_p = leafset->bits;	\
		leaf_id_t words_left = leafset_nwords;				\
		leaf_id_t leaf_var = 0;												\
		while( words_left-- ) {												\
			if ( *bits_p ) {														\
				leafset_word_t mask = 1;									\
				while( mask ) {														\
					if ( *bits_p & mask ) {									\
						body;																	\
					}																				\
					mask <<= 1;															\
					leaf_var++;															\
				}																					\
			}																						\
			else  leaf_var += BITS_PER_LEAFSET_WORD;		\
      bits_p++;																		\
		}																							\
  } while(0)

} // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_LEAFSET_BITSET_H
