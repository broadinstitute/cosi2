#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <boost/swap.hpp>
#include <cosi/utils.h>
#include <cosi/segsumm.h>

namespace cosi {

#if 0
static void segsumm_print( const segsumm_t *s ) {
  printf( "\n" );
  for ( seglet_idx_t i = 0; i < NUM_SEGLETS; i++ ) 
		 printf( "%d", !segsumm_is_seglet_definitely_empty( s, i ) );
  printf( "\n" );
  for ( seglet_idx_t i = 0; i < NUM_SEGLETS; i++ ) 
		 printf( "%d", !segsumm_is_seglet_definitely_full( s, i ) );
  printf( "\n" );
}
#endif

ostream& operator<<( ostream& strm, const segsumm_t& s ) {
  for ( seglet_idx_t i = 0; i < NUM_SEGLETS; i++ ) 
		 strm << !segsumm_is_seglet_definitely_empty( &s, i );
	return strm;
}

static const int SEGSUMM_PRECOMP_DIST = 10;


segsumm_t segsumm_recomb[ NUM_SEGLETS ][2];
segsumm_t segsumm_gc[ NUM_SEGLETS][ SEGSUMM_PRECOMP_DIST ][2];

segsumm_t empty_segsumm, full_segsumm;

static void segsumm_set_bit( segsumm_word_t *w, seglet_idx_t i ) { *w |= ( ((segsumm_word_t)1) << i ); }

static void segsumm_set_bits( segsumm_word_t *w, seglet_idx_t i, seglet_idx_t j ) {
  while( i <= j )  segsumm_set_bit( w, i++ );
}

static void segsumm_clear_bit( segsumm_word_t *w, seglet_idx_t i ) { *w &= ~( ((segsumm_word_t)1) << i ); }
static void segsumm_clear_bits( segsumm_word_t *w, seglet_idx_t i, seglet_idx_t j ) {
  while( i <= j )  segsumm_clear_bit( w, i++ );
}

void segsumm_init_module(void) {
  static bool_t initialized;
  if ( initialized ) return;
  initialized = True;

  /* Precompute the segsumms for [] and [0,1] */
  segsumm_clear_bits( &( empty_segsumm.summ ), 0, NUM_SEGLETS-1 );
  segsumm_set_bits( &( empty_segsumm.inv_summ ), 0, NUM_SEGLETS-1 );
  
  segsumm_set_bits( &( full_segsumm.summ ), 0, NUM_SEGLETS-1 );
  segsumm_clear_bits( &( full_segsumm.inv_summ ), 0, NUM_SEGLETS-1 );

  /* Precompute bitmasks for quickly computing <segsumm_split()> and <segsumm_split_on_seg()> */
  
  for ( seglet_idx_t i = 0; i < NUM_SEGLETS; i++ ) {

		/* Precompute bitmasks for computing the left and the right segsumm after a seglist is */
		/* split at a point falling into seglet i. */
		/* The right result is definitely empty at seglets strictly before i; */
		 
		segsumm_init_empty_seg( &( segsumm_recomb[i][DIR_L] ) );
		segsumm_set_bits( &( segsumm_recomb[i][DIR_L].summ ), 0, i );
		segsumm_recomb[i][DIR_L].inv_summ = 0;
		segsumm_set_bits( &( segsumm_recomb[i][DIR_L].inv_summ ), i, NUM_SEGLETS-1 );

		segsumm_inverse( &( segsumm_recomb[i][DIR_R] ), &( segsumm_recomb[i][DIR_L] ) );

		/*******************/

		for ( int dist = 0; dist < SEGSUMM_PRECOMP_DIST && ( ( i+dist ) < NUM_SEGLETS ); dist++ ) {
			segsumm_init_empty_seg( &( segsumm_gc[i][dist][GC_INNER] ) );
			segsumm_set_bits( &( segsumm_gc[i][dist][GC_INNER].summ ), i, i+dist );
			segsumm_inverse( &( segsumm_gc[i][dist][GC_OUTER] ), &( segsumm_gc[i][dist][GC_INNER] ) );
		}
  }
}

bool_t segsumm_is_empty( const segsumm_t *s ) { return s->summ == 0; } 

void segsumm_copy( segsumm_t *target, const segsumm_t *source ) { *target = *source; }

void segsumm_init_empty_seg( segsumm_t *s ) { segsumm_copy( s, &empty_segsumm ); }
void segsumm_init_full_seg( segsumm_t *s ) { segsumm_copy( s, &full_segsumm ); }

bool_t segsumm_equal( const segsumm_t *s1, const segsumm_t *s2 ) { return s1->summ == s2->summ; }

void segsumm_union( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 ) {
  result->summ = s1->summ | s2->summ;
  /* so, if a bit is completely full in one or completely full in the other, it is definitely completely full in the union. */
  result->inv_summ = s1->inv_summ & s2->inv_summ;
  
}

void segsumm_intersection( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 ) {
  result->summ = s1->summ & s2->summ;
  result->inv_summ = s1->inv_summ | s2->inv_summ;
}

bool_t segsumm_intersection_if_nonempty( segsumm_t *result, const segsumm_t *s1, const segsumm_t *s2 ) {
  segsumm_word_t inters = s1->summ & s2->summ;
  if ( inters ) {
		result->summ = inters;
		result->inv_summ = s1->inv_summ | s2->inv_summ;
		return True;
  } else return False;
}

bool_t segsumm_intersects( const segsumm_t *s1, const segsumm_t *s2 ) { return ( s1->summ & s2->summ ) != 0; }
void segsumm_inverse( segsumm_t *result, const segsumm_t *s ) {
  if( s != result ) { result->summ = s->inv_summ; result->inv_summ = s->summ; }
  else boost::swap( result->summ, result->inv_summ );
}

seglet_idx_t segsumm_get_seglet( loc_t loc ) {
  return loc >= MAX_LOC ? NUM_SEGLETS-1 : ( loc <= MIN_LOC ? 0 : ((seglet_idx_t)(get_loc( loc ) * ((double)NUM_SEGLETS))) );
}


bool_t segsumm_is_seglet_definitely_empty( const segsumm_t *s, seglet_idx_t i ) {
  segsumm_word_t mask = ( ((segsumm_word_t)1ULL) << i );
  segsumm_word_t the_and = ( s->summ & mask );
  return the_and == 0;
}

bool_t segsumm_is_seglet_definitely_full( const segsumm_t *s, seglet_idx_t i ) {
  segsumm_word_t mask = ( ((segsumm_word_t)1ULL) << i );
  segsumm_word_t the_and = ( s->inv_summ & mask );
  return the_and == 0;
}


void segsumm_split( const segsumm_t *s, loc_t loc, segsumm_t *s_left, segsumm_t *s_right ) {
  extern segsumm_t segsumm_recomb[ NUM_SEGLETS ][2];
  seglet_idx_t i = segsumm_get_seglet( loc );
  segsumm_t s_segsumm = *s;
  segsumm_intersection( s_left, &( s_segsumm ), &(segsumm_recomb[i][DIR_L]) );
  segsumm_intersection( s_right, &( s_segsumm ), &(segsumm_recomb[i][DIR_R]) );
  
}

bool_t segsumm_split_on_seg( const segsumm_t *s, loc_t loc1, loc_t loc2,
														 segsumm_t *s_inside, segsumm_t *s_outside ) {

  segsumm_t s_segsumm;
  segsumm_copy( &s_segsumm, s );
  seglet_idx_t i = segsumm_get_seglet( loc1 );
  seglet_idx_t j = segsumm_get_seglet( loc2 );

  seglet_idx_t dist = j - i;
  if ( dist < SEGSUMM_PRECOMP_DIST ) {
		if ( !segsumm_intersection_if_nonempty( s_inside,
																						&( s_segsumm ), &(segsumm_gc[i][ dist ][GC_INNER]) ) )
			 return False;
	 
		segsumm_intersection( s_outside,
													&( s_segsumm ), &(segsumm_gc[i][ dist ][GC_OUTER]) );
  } else {
		segsumm_t inner, outer;
		segsumm_init_empty_seg( &inner );
		segsumm_set_bits( &( inner.summ ), i, j );
		segsumm_inverse( &outer, &inner );
	 
		if ( !segsumm_intersection_if_nonempty( s_inside,
																						&( s_segsumm ), &inner ) )
			 return False;
		segsumm_intersection( s_outside,
													&( s_segsumm ), &outer );
	 
  }
  return True;
}

}  // namespace cosi
