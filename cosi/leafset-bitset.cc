/*
 * File: leafset.c
 *
 * Code for manipulating sets of leaves of the ARG (corresponding to present-day <chroms>).
 */

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <boost/scoped_array.hpp>
#include <cosi/utils.h>
#include <cosi/leafset-bitset.h>

namespace cosi {

using boost::scoped_array;

static Mempool mempool_leafsets;

/*static void leafset_clear( leafset_p  );*/

/* Static FuncP: leafset_add */
/* Add one leaf to the given leafset. */
void leafset_add( leafset_p leafset, leaf_id_t leaf );

/* Static var: max_leaf_id */
/* The maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
static leaf_id_t max_leaf_id = NULL_LEAF_ID;

/* Static var: leafset_nwords */
/* Number of words in the bitset representation of leafsets. */
size_t leafset_nwords;

/* Static var: nbits_per_word */
/* Number of bits in one <leafset_word_t>. */
static size_t nbits_per_word;

/* Static var: leafset_bits_size */
/* Size of <leafset_t::bits>, in bytes */
static size_t leafset_bits_size;

/* Static var: leafset_size */
/* Leafset size in bytes */
static size_t leafset_size;

/* Static var: leaf_word_idx */
/* For each possible leaf id, the index of the word in <leafset_t::bits> containing the bit representing this leaf. */
static scoped_array<int> leaf_word_idx;

/* Var: leaf_set_mask */
/* The mask used for setting a particular leaf, once the corresponding word is known. */
scoped_array<leafset_word_t> leaf_set_mask;

/* Var: leaf_clear_mask */
/* The mask used for clearing a particular leaf, once the corresponding word is known. */
scoped_array<leafset_word_t> leaf_clear_mask;

/* Var: singleton_leafsets */
/* All the possible singleton leafsets */
static vector<leafset_p> singleton_leafsets;

/* Func: set_max_leaf_id */
/* Sets the maximum leaf id that can be used in a leafset, */
/* as determined by the sample size of our simulated samples. */
void leafset_set_max_leaf_id( leaf_id_t max_leaf_id_ ) {
  max_leaf_id = max_leaf_id_;


  /* precompute the sizes of various structures that we'll repeatedly allocate */ 
  nbits_per_word = sizeof( leafset_word_t ) * BITS_PER_BYTE;
  leafset_nwords = ((int)( (max_leaf_id+1) / nbits_per_word )) + 1;
  leafset_bits_size = (leafset_nwords) * sizeof( leafset_word_t );
  leafset_size = sizeof( leafset_t ) + leafset_bits_size;

  mempool_leafsets.mempool_init( leafset_size, 65536, 4096 );
  {
	 leaf_id_t leaf;

	 /* precompute bitmasks that will help up set and clear particular bits of the leafset */
	 leaf_word_idx.reset( new int[ max_leaf_id ] );
	 leaf_set_mask.reset( new leafset_word_t[ max_leaf_id ] );
	 leaf_clear_mask.reset( new leafset_word_t[ max_leaf_id ] );
	 for ( leaf = 0; leaf < max_leaf_id; leaf++ ) {
		int bit_id = leaf % nbits_per_word;
		leafset_word_t mask = 1;
		
		leaf_word_idx[ leaf ] = (int)(leaf / nbits_per_word);
		leaf_set_mask[ leaf ] = ( mask << bit_id );
		leaf_clear_mask[ leaf ] = ~( leaf_set_mask[ leaf ] );
	 }
  }

  for ( leaf_id_t leaf = 0; leaf < max_leaf_id; leaf++ )
	 singleton_leafsets.push_back( make_singleton_leafset( leaf ) );
}  /* leafset_set_max_leaf_id() */

leaf_id_t leafset_get_max_leaf_id(void) { return max_leaf_id; }

void* leafset_struct::operator new (size_t size) {
  return mempool_leafsets.mempool_alloc_item();
}

void leafset_struct::operator delete (void *p) {
  mempool_leafsets.mempool_free_item( p );
}


leafset_struct::leafset_struct() {
  memset( bits, 0, leafset_bits_size );
}

leafset_p leafset_alloc(void) {
  return leafset_p( new leafset_struct );
}

leafset_p make_empty_leafset(void) {
  return leafset_alloc();
}

leafset_p make_singleton_leafset( leaf_id_t leaf ) {
  leafset_p leafset = make_empty_leafset();
  leafset_add( leafset, leaf );
  return leafset;
}

/* FuncP: make_range_leafset */
/* Returns the leafset containing the given range of leaves. */
leafset_p make_range_leafset( leaf_id_t fromLeaf, leaf_id_t toLeaf ) {
  leafset_p leafset = make_empty_leafset();
  for ( leaf_id_t leaf = fromLeaf; leaf < toLeaf; leaf++ )
	 leafset_add( leafset, leaf );
  return leafset;
}


bool_t leafset_contains( leafset_p leafset, leaf_id_t leaf ) {
  return ( leafset->bits[ leaf_word_idx[ leaf ] ] & leaf_set_mask[ leaf ] ) > 0;
}

bool_t leafset_is_subset( leafset_p L1, leafset_p L2 ) {
  int i = leafset_nwords;
  const leafset_word_t *L1p = L1->bits, *L2p = L2->bits;
  while( i-- ) {
	 if ( ( *L2p | *L1p ) != *L2p ) return False;
	 L1p++; L2p++;
  }
  return True;
}

leafset_p leafset_union( leafset_p s1, leafset_p s2) {
  assert( s1 );
  assert( s2 );
  leafset_p s_ret = leafset_alloc();

  leafset_word_t *s = s_ret->bits;
  const leafset_word_t *s1p = s1->bits, *s2p = s2->bits;
  int i = leafset_nwords;
  /*bool_t same_as_s1 = True, same_as_s2 = True;*/
  /* if ( s1 == s2 ) return (leafset_t *)s1;*/
  while( i-- ) {
	 *s = *s1p | *s2p;
	 /*if ( *s1 & *s2 ) { printf( "chego??\n" ); exit(2); }*/
#if 0
	 same_as_s1 = same_as_s1 && ( *s == *s1 );
	 same_as_s2 = same_as_s2 && ( *s == *s2 );
#endif	 
	 s++; s1p++; s2p++;
  }
#if 0
  if ( same_as_s1 ) { /*leafset_delete( s_ret );*/ return (leafset_t *) s1; }
  if ( same_as_s2 ) { /*leafset_delete( s_ret );*/ return (leafset_t *) s2; }
#endif

  /*printf( "%s | %s = %s\n", leafset_str( s1 ), leafset_str( s2 ), leafset_str( s_ret ) );*/
  return s_ret;
}

leafset_p leafset_intersection( leafset_p s1, leafset_p s2) {
  assert( s1 );
  assert( s2 );
  leafset_p s_ret = leafset_alloc();

  leafset_word_t *s = s_ret->bits;
  const leafset_word_t *s1p = s1->bits, *s2p = s2->bits;
  int i = leafset_nwords;
  /*bool_t same_as_s1 = True, same_as_s2 = True;*/
  /* if ( s1 == s2 ) return (leafset_p )s1;*/
  while( i-- ) {
	 *s = *s1p & *s2p;
	 /*if ( *s1 & *s2 ) { printf( "chego??\n" ); exit(2); }*/
#if 0
	 same_as_s1 = same_as_s1 && ( *s == *s1 );
	 same_as_s2 = same_as_s2 && ( *s == *s2 );
#endif	 
	 s++; s1p++; s2p++;
  }
#if 0
  if ( same_as_s1 ) { /*leafset_delete( s_ret );*/ return (leafset_p ) s1; }
  if ( same_as_s2 ) { /*leafset_delete( s_ret );*/ return (leafset_p ) s2; }
#endif

  /*printf( "%s | %s = %s\n", leafset_str( s1 ), leafset_str( s2 ), leafset_str( s_ret ) );*/
  return s_ret;
}



bool_t leafset_intersects( leafset_p s1, leafset_p s2) {
  assert( s1 );
  assert( s2 );

  const leafset_word_t *s1p = s1->bits, *s2p = s2->bits;
  int i = leafset_nwords;
  while( i-- )
	 if ( *s1p++ & *s2p++ ) return True;

  return False;
}


leafset_p leafset_complement( leafset_p sa ) {

  assert( sa );
  leafset_p s_ret = leafset_alloc();

  leafset_word_t *s = s_ret->bits;
  const leafset_word_t *sp = sa->bits;
  int i = leafset_nwords;

  while( i-- ) {
	 *s = ~(*sp);
	 s++; sp++;
  }

  return s_ret;
}

leafset_p leafset_difference( leafset_p s1, leafset_p s2) {
  assert( s1 );
  assert( s2 );
  leafset_p s2c = leafset_complement( s2 );
  leafset_p result = leafset_intersection( s1, s2c );
  //leafset_delete( s2c );
  return result;
}


void leafset_add( leafset_p leafset, leaf_id_t leaf ) {
  assert( 0 <= leaf );
  assert( leaf < max_leaf_id );
  leafset->bits[ leaf_word_idx[ leaf ] ] |= leaf_set_mask[ leaf ];
}

void leafset_print( leafset_p leafset ) {
  int i;
  for ( i = 0; i < max_leaf_id; i++ )
	 if ( leafset_contains( leafset, i ) )
		printf( "%d ", i );
  printf( "\n" );
}

int leafset_cardinality( leafset_p leafset )  {
  leaf_id_t i;
  int n = 0;
  for ( i = 0; i < max_leaf_id; i++ )
	 if ( leafset_contains( leafset, i ) ) n++;
  return n;
}
bool_t leafset_is_empty( leafset_p leafset ) {
  const leafset_word_t *sp = leafset->bits;
  int i = leafset_nwords;
  while( i-- )
	 if ( *sp++ ) return False;
  return True;
}

leafset_p leafset_from_str( const char *s ) {
  leafset_p leafset = make_empty_leafset();
  char *str = cosi_strdup( s );
  char *saveptr;
  for ( char *token = cosi_strtok_r( str, ",", &saveptr ); token; token = cosi_strtok_r( NULL, ",", &saveptr ) ) {
	 leaf_id_t leaf;
	 int num_scanned = sscanf( token, "%d", &leaf );
	 chkCond( num_scanned == 1, "leafset_from_str: bad set representation - %s", s );
	 leafset_add( leafset, leaf );
  }
  free( str );
  /*printf( "leafset from str: %s is ", s ); leafset_print( leafset );*/
  return leafset;
}

void leafset_write_bin( FILE *f, leafset_p leafset ) {
  cosi_fwrite_val( leafset_nwords, f );
  cosi_fwrite( leafset->bits, sizeof( leafset_word_t ), leafset_nwords, f );
}
leafset_p leafset_read_bin( FILE *f ) {
  leafset_p leafset = make_empty_leafset();
  size_t nwords;
  cosi_fread_val( nwords, f );
  chkCond( nwords == leafset_nwords, "leafset_read_bin: leafset size mismatch -- read %u expected %u", nwords, leafset_nwords );
  cosi_fread( leafset->bits, sizeof( leafset_word_t ), leafset_nwords, f );
  return leafset;
}

std::string leafset_str( leafset_p leafset ) {
  std::ostringstream buf;

  bool_t first = True;
  leaf_id_t range_start = NULL_LEAF_ID;
  leaf_id_t leaf = 0;
  for ( ; leaf < max_leaf_id; leaf++ ) {
	 if ( leafset_contains( leafset, leaf ) ) {
		if ( range_start == NULL_LEAF_ID ) {
		  range_start = leaf;
		  if ( !first ) buf << ",";
		  
		  first = False; 
		  buf << leaf;
		}
	 } else if ( range_start != NULL_LEAF_ID ) {
		if ( leaf-1 > range_start )
			 buf << "-" << ( leaf-1 );
		range_start = NULL_LEAF_ID;
	 }
  }
  if ( range_start != NULL_LEAF_ID ) {
	 if ( leaf-1 > range_start )
			buf << "-" << ( leaf-1 );
	 range_start = NULL_LEAF_ID;
  }  
  
  return buf.str();
}

bool_t leafset_is_singleton( leafset_p leafset, leaf_id_t leaf ) { return leafset_equal( leafset, singleton_leafsets[ leaf ] ); }

bool_t leafset_equal( leafset_p leafset1, leafset_p leafset2) {
  extern size_t leafset_bits_size;
  return leafset1 == leafset2 ||
		!memcmp( leafset1->bits, leafset2->bits, leafset_bits_size );
}


}  // namespace cosi
