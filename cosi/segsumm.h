/*
 * Header: segsumm.h
 *
 * A coarse approximation to a <Seglist>, kept with each seglist; lets us quickly detect when two lists do not intersect.
 * We divide the interval [0,1] into a given number of sub-intervals( "seglets" ).  A segsumm for a seglist is then a bitmask
 * with one bit per seglet, representing whether any part of the seglist might fall inside the seglet.  The approximation
 * is conservative in that a segsumm might show a 1 bit in a seglet which the seglist does not actually overlap.
 * However, a zero bit in a seglet guarantees that the seglist does not intersect that seglet.
 * So, if two segsumms do not intersect, then the corresponding seglists definitely do not intersect.
 * Segsumm intersection can be computed quickly using bitwise operations.
 *
 * To facilitate computing the segsumm of a seglist's complement given the seglist's segsumm, we also keep for each seglist
 * a summary of its complement.
 */

#ifndef __INCLUDE_COSI_SEGSUMM_H
#define __INCLUDE_COSI_SEGSUMM_H

#include <iostream>
#include <cosi/defs.h>
#include <cosi/utils.h>

namespace cosi {
/* FuncP: segsumm_init_module */
/* Initialize the segsumm module.  Must be called once, before calling any other routines in this module. */ 
void segsumm_init_module(void);

/* Const: NUM_SEGLETS */
/* Number of intervals into which we divide the [0,1] interval for summary purposes. */
/* Current code assumes this value is 64, so that seglet operations can be done with 64-bit integers. */
/* But it would not be hard to change the code to work with multiple machine words. */
const int NUM_SEGLETS = 64;

/* Type: seglet_idx_t */
/* Id of a seglet */
typedef int seglet_idx_t;

/* Type: seglet_word_t */
/* One machine word, as used in the bitset representation of <segsumm>. */
typedef unsigned long long segsumm_word_t;

/* Struct: segsumm */
/* A coarse summary of what <seglets> a <Seglist> enters. */
typedef struct segsumm {
	 /* Field: summ */
	 /* For each <seglet>, whether the given <Seglist> reaches inside that seglet. */
	 /* A zero bit means the corresponding seglet is definitely clear; a one bit says nothing. */
	 segsumm_word_t summ;

	 /* Field: inv_summ */
	 /* A summary of the inverse: a zero bit means the corresponding seglet is definitely completely full in <summ>. */
	 /* A one bit says nothing. */
	 segsumm_word_t inv_summ;
} segsumm_t;

extern segsumm_t empty_segsumm, full_segsumm;

bool_t segsumm_split_on_seg( const segsumm_t *s, loc_t loc1, loc_t loc2,
														 segsumm_t *s_inside, segsumm_t *s_outside );
std::ostream& operator<<( ostream&, const segsumm_t& );
void segsumm_split( const segsumm_t *s, loc_t loc, segsumm_t *s_left, segsumm_t *s_right );
bool_t segsumm_is_seglet_definitely_full( const segsumm_t *s, seglet_idx_t i );
bool_t segsumm_is_seglet_definitely_empty( const segsumm_t *s, seglet_idx_t i );
seglet_idx_t segsumm_get_seglet( loc_t loc );
void segsumm_inverse( segsumm_t *result, const segsumm_t *s );
bool_t segsumm_intersects( const segsumm_t *s1, const segsumm_t *s2 );
bool_t segsumm_intersection_if_nonempty( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 );
void segsumm_intersection( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 );
void segsumm_union( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 );
bool_t segsumm_equal( const segsumm_t *s1, const segsumm_t *s2 );
void segsumm_init_full_seg( segsumm_t *s );
void segsumm_init_empty_seg( segsumm_t *s );
void segsumm_copy( segsumm_t *target, const segsumm_t *source );
bool_t segsumm_is_empty( const segsumm_t *s );

}  // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_SEGSUMM_H


