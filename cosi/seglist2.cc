
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <cosi_rand/random.h>
#include <cosi_rand/mtwist.h>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/leafset.h>
#include <cosi/mempool.h>
#include <cosi/seglist2.h>

namespace cosi {

namespace seglist2 {

namespace seglist_detail {

using util::chkCond;
using util::Date;
using util::cosi_strdup;
using util::cosi_strtok_r;

int n_seglist_empty, n_left_or_right, n_containing, n_in_gap, n_nonempty, n_segsumm, n_oneseg;

#ifdef COSI_DEV_SEGLIST_DEBUG
static Seglist *all_seglists;
#endif

clock_t time_gc;

/*
 * Var: mt_state
 *
 * The random state used to generate the randomly chosen <levels>
 * of the <Segs>.  Kept separate from the random generator of the simulator,
 * because the random choices of the Seg levels should only affect the
 * performance but never the results of seglist operations.  If we change
 * the random seed or the generator used for seg levels, this should have no effect
 * on otherwise deterministic simulations that use the seglists.
 */
static mt_state seglist_random_state;

/* Func: seglist_random_bit */
/* Generate one random bit from <seglist_random_state>. */
static inline bool_t seglist_random_bit(void) {
  static unsigned long val;
  static int bits_left;

  if ( !bits_left ) {
		val = mts_lrand( &seglist_random_state );
		bits_left = 31;
  } else {
		val >>= 1;
		bits_left--;
  }
  return ( val & 1 ) != 0;
}

/* Static func: seglist_new_level */
/* Choose a level for a newly allocated <Seg>. */

#ifdef COSI_STATS
int level_hist[1024];
#endif

/*
 * Static func: seglist_new_level
 *
 * Choose a level for a new <Seg> being added to a <Seglist>.  The level is chosen from
 * a geometric distribution, currently with p=.25 .  We use a source of randomness that is
 * separate from the one used to guide the simulations.
 */
static level_t seglist_new_level(void) {
  level_t level = seglist_min_level;
  while ( level < seglist_max_level-1 &&
					/* Experimentally, seglist performance seems best when p=.25 */
					/* that is, when about one-quarter of segs in level i are also in level i+1. */
					seglist_random_bit() && seglist_random_bit() )
		 level++;
  IF_COSI_STATS( level_hist[ level ]++ );
  return level;
}

/* Var: mempool_seglist */
/* <Mempool> from which <Seglist> objects are allocated. */
static Mempool mempool_seglist;

/* Var: mempool_seg */
/* <Mempool>s from which <Seg> objects are allocated; a separate Mempool */
/* for each possible Seg level, since the memory size of Segs differ by level. */
static boost::scoped_array<Mempool> mempool_seg;

static Mempool mempool_finger;

int n_free_seg, n_alloc_seg, n_free_adj, n_free_doneseg;

void* seg::operator new (size_t /*size*/, level_t level) {
  // size_t item_size = sizeof( Seg ) + ( level + seglist_min_level ) * sizeof( Segptr );
  // return malloc( item_size );
  return mempool_seg[ level ].mempool_alloc_item( );
}
void seg::operator delete (void *p) {
  //free( p );
  Seg *seg = (Seg *)p;
  mempool_seg[ seg->level ].mempool_free_item( seg );
}


/* Static func: seglist_alloc_seg */
/* Allocate a new <Seg> with the given level. */
static Seg *seglist_alloc_seg( level_t level ) {
  Seg *seg = new (level) Seg;
  seg->level = level;
	seg->leafset = LEAFSET_NULL;
	seg->lastCoalGen = NULL_GEN;
	
#ifdef COSI_DEV_SEGLIST_DEBUG  
  seg->id = get_obj_id();
#endif
	
#ifdef COSI_TRACK_LAST_COAL
	seg->lastCoalGen = genid(0);
#endif	
	
  n_alloc_seg++;
  return seg;
}

/*
 * Static func: seglist_alloc_header
 *
 * Allocate a <Seglist> header.  The header is a special element, kept in each Seglist
 * (including empty ones), that for each level points to the leftmost <Seg> at that level.
 * The <Seg::level> field of the header stores the level of the Seglist: one higher than
 * the highest level of any Seg in the Seglist.
 */
static Seg *seglist_alloc_header( void ) {
  Seg *seg = new (seglist_max_level) Seg;
  seg->level = seglist_max_level;
  seg->beg = seg->end = MIN_LOC;
  seg->leafset = LEAFSET_NULL;
#ifdef COSI_DEV_SEGLIST_DEBUG  
  seg->id = get_obj_id();
#endif
  n_alloc_seg++;
  return seg;
}


/* Static global var: NIL */
/* A special <Seg> that ends every <Seglist>.  Its (beg,end) is at +infinity. */
/* It is shared by all seglists and is never deallocated. */
Seg *NIL;

/*
 * Static func: seglist_free_seg
 * Free one Seg of a Seglist.
 *
 * Params:
 *
 *    s - the seg we're freeing
 */
static void seglist_free_seg_helper( Seg *seg ) {
  IF_COSI_STATS( n_free_seg++ );
  assert( seg != NIL );
  delete seg;
}
static inline void seglist_free_seg( Seg *&s) { seglist_free_seg_helper(s); s = NULL; }

static void seglist_free_seglist_helper( Seglist *seglist );


void seglist_free_seglist_helper( Seglist *seglist ) {
  Seg *seg = seglist->header;
	seg->level = seglist_max_level;  // make sure the header is freed from the right mempool
  while( seg != NIL ) {
		Seg *nxt = seg->fwd[ seglist_min_level ].seg;
		seglist_free_seg( seg );
		seg = nxt;
  }
#ifdef COSI_DEV_SEGLIST_DEBUG
  if ( seglist->prev ) {
		seglist->prev->next = seglist->next;
		if ( seglist->next ) 
			 seglist->next->prev = seglist->prev;
  } else {
		all_seglists = seglist->next;
		if ( seglist->next ) seglist->next->prev = NULL;
  }
#endif
  mempool_seglist.mempool_free_item( seglist );
}
static inline void seglist_free_seglist(Seglist *& s) { seglist_free_seglist_helper(s); s = NULL; }

void seglist_delete( Seglist **seglistp ) {
  if ( *seglistp ) seglist_free_seglist( *seglistp );
}


/* Static func: finger_rightmost_seg */
/* Gives the rightmost seg at a given level, as an l-value. */
static inline Seg *& finger_rightmost_seg( Finger *finger, level_t level ) { return (finger)[ (level) ].seg; }

/* Static func: finger_travel_len */
/* Gives the length traveled along a given level, as an l-value. */
static inline len_t& finger_travel_len( Finger *finger, level_t level ) { return (finger)[ (level) ].skip_len; }

//static size_t seglist_seg_size( level_t level ) { return sizeof( Seg ) + ( level + seglist_min_level ) * sizeof( Segptr ); }

static size_t finger_size;

/*
 * Func: seglist_init_module
 * Initializes the seglist module; must be called once, before any operations in this module.
 */
void seglist_init_module(void) {
  segsumm_init_module();
  if ( !NIL ) {
	 
		mts_seed32new( &seglist_random_state, getenv( "COSI_SEGLIST_SEED" ) ? atol( getenv( "COSI_SEGLIST_SEED" ) ) : 373737 );
		/*mts_goodseed( &seglist_random_state );*/

		/*
		 * Allocate Mempools for Seglists, seglist headers, and for Segs of each level.
		 */

		mempool_seglist.mempool_init( sizeof( Seglist ), /* block_size_in_items= */ 4096,
																	/* initial_freelist_size_in_items= */ 1024 );

		mempool_seg.reset( new Mempool[ seglist_max_level+seglist_min_level+1 ] );
		int pool_size = 512000;
		int min_pool_size = 1024;
		int max_level_pool_size = 131072;
		//PRINT5( seglist_min_level, seglist_max_level, pool_size, min_pool_size, max_level_pool_size );
		for ( level_t level = seglist_min_level; level <= seglist_max_level; level++ ) {
			size_t item_size = sizeof( Seg ) + ( level + seglist_min_level ) * sizeof( Segptr );
			//PRINT2( level, item_size );
			mempool_seg[ level ].mempool_init( item_size,
																				 /* block_size_in_items= */ level < seglist_max_level ? pool_size : max_level_pool_size,
																				 /* initial_freelist_size_in_items= */ 1024 );
			if ( pool_size > min_pool_size )
				 pool_size >>= 1;
		}

		/*
		 * Set up the NIL element that terminates all seglists.
		 */
	 
		NIL = seglist_alloc_seg( seglist_max_level );
		NIL->beg = make_loc( MAX_PLOC + plen_t( 1.0 ), MAX_GLOC + glen_t( 1.0 ) );
		NIL->end = make_loc( MAX_PLOC + plen_t( 2.0 ), MAX_GLOC + glen_t( 2.0 ) );
		NIL->leafset = LEAFSET_NULL;
		NIL->fwd[ seglist_min_level ].seg = NULL;
		NIL->fwd[ seglist_min_level ].skip_len = ZERO_LEN;
		for ( level_t level = seglist_min_level+1; level <= seglist_max_level; level++ )
			 NIL->fwd[ level ] = NIL->fwd[ seglist_min_level ];
	 
	 
		finger_size = ( seglist_max_level + seglist_min_level ) * sizeof( Segptr );

		mempool_finger.mempool_init( finger_size, 128, 64 );
	 
  }
}  /* seglist_init_module() */

struct seglist_initializer {
	 seglist_initializer() { seglist_init_module(); }
	 ~seglist_initializer() {
		 if ( NIL ) {
			 delete NIL;
			 NIL = NULL;
		 }
	 }
};
static seglist_initializer seglist_initializer_dummy;

/* Func: seglist_make_empty */
/* Create an empty seglist. */
Seglist *seglist_make_empty( ) {
  Seglist *seglist = (Seglist *)mempool_seglist.mempool_alloc_item();
  
#ifdef COSI_DEV_SEGLIST_DEBUG
  seglist->id = get_obj_id();
  
  if ( all_seglists ) all_seglists->prev = seglist;
  seglist->next = all_seglists;
  all_seglists = seglist;
	all_seglists->prev = NULL;
#endif  
  
  seglist->header = seglist_alloc_header();
  seglist->header->level = seglist_min_level;
  seglist->header->fwd[ seglist_min_level ].seg = NIL;
  seglist->header->fwd[ seglist_min_level ].skip_len = ZERO_LEN;
  seglist->header->leafset = LEAFSET_NULL;
  seglist->setEnd( MIN_LOC );
  
  segsumm_init_empty_seg( &( seglist->segsumm ) );
  
  return seglist;
}


/*
 * Func: seglist_compact_level
 *
 * Reduce the level of a seglist to the smallest necessary.
 * The header's forward pointer at the highest level needs to point to NIL;
 * any levels above that, are extraneous.
 *
 * Seglists in which more than one forward pointer from the header points to NIL
 * can arise because for some operations on two seglists, we raise the
 * level of one seglist to make the two seglists have the same level
 * (see <seglist_increase_seglist_level_to()>).
 */
static void seglist_compact_level( Seglist *seglist ) {
  Seg *h = seglist->header;
  while( h->level >= seglist_min_level+1 && h->fwd[ h->level-1 ].seg == NIL )
		 h->level--;
}

/*
 * Func: seglist_init_finger
 *
 * Initialize the cursor/iterator for traversing a Seglist to the beginning of the seglist.
 *
 * Input params:
 *
 *   seglist - the seglist we want to traverse
 *   init_lens - whether to zero out the finger_travel_len value for each level.
 *
 * Output params:
 *
 *   finger - the finger we're initializing 
 */
void seglist_init_finger( const Seglist *seglist, Finger *finger, bool_t /*init_lens*/ ) {
  for ( level_t level = seglist_min_level; level <= seglist->header->level; level++ ) {
		finger_rightmost_seg( finger, level ) = seglist->header;
		/*if ( init_lens )*/ finger_travel_len( finger, level ) = ZERO_LEN;
  }
}

#if 0
/*
 * Func: seglist_advance_finger_to_next
 *
 * Move the lowest level of the finger to the next Seg in the Seglist (assumes we're not yet at the last seg).
 * Note that after this, some higher level pointers in the finger may no longer point to the rightmost
 * seg we have traveled to.
 */
static void seglist_advance_finger_to_next( Seglist * /*seglist*/, Finger *finger ) {
  assert( finger_rightmost_seg( finger, seglist_min_level ) != NIL );
  finger_rightmost_seg( finger, seglist_min_level ) = seglist_next_seg( finger_rightmost_seg( finger, seglist_min_level ) );
}

/*
 * Func: seglist_advance_finger_to_loc
 *
 * Set each level of a Finger to the rightmost Seg of a Seglist that is still fully to the left of loc.
 */
static void seglist_advance_finger_to_loc( Seglist *seglist, Finger *finger, loc_t loc ) {
  Seg *seg = finger_rightmost_seg( finger, seglist->header->level );
  for ( level_t level = seglist->header->level; level >= seglist_min_level; level-- ) {
		assert( seg->end <= loc );
		Seg *nxt;
		while( ( nxt = seg->fwd[ level ].seg )->end <= loc )
			 seg = nxt;
		finger_rightmost_seg( finger, level ) = seg;
  }
}
#endif

/*
 * Func: seglist_advance_by_len
 *
 * Advance a Finger by the given amount of length along the seglist, counting only
 * the length traveled along the seglist's segments and not the gaps between segs.
 * Used by <mutate_put_muts_on_seglist()>.  Making this routine fast is the whole reason
 * for keeping <Segptr::skip_len> for each forward pointer.
 *
 * Params:
 *
 *   seglist - the seglist over which we're moving
 *   finger - current position in the seglist; more specifically, for each level, the rightmost seg that we have
 *      seen at that level
 *   cur_loc - the current location along the seglist.  updated as we move along the list.
 *   cur_len - the total length we have traveled so far along the seglist (counting only the segs not the gaps)
 *   cur_leafset - set to the leafset of the seg on which cur_loc is located
 *   advance_by - on input, the length by which we need to advance along the list.  on output,
 *     if advancing by this much would put us beyond the end of the list, the distance left to travel along the segs of
 *     the list from loc to its end is subtracted from *advance_by.
 *
 * Returns:
 *   if we're still within the seglist, returns True, and updates loc to the new location; also updates finger.
 *   if advancing this much will move past the end of the list, returns False, and reduces *advance_by by the
 *   length traveled along the list.
 *
 * Invariants and state:
 *
 *    finger_rightmost_seg(level) points to the rightmost seg at that level that we have reached;
 *    finger_travel_len(level) gives the total length of the segments preceding finger_rightmost_seg(level).
 */
bool_t seglist_advance_by_len( const Seglist *seglist, Finger *finger,
															 loc_t *cur_loc, len_t *cur_len,  leafset_p *cur_leafset, len_t *advance_by ) {

  len_t target_len = *cur_len + *advance_by;

  len_t left_in_list = ( seglist_tot_len( seglist ) - *cur_len );
  
  if ( get_phys_len( *advance_by ) > get_phys_len( left_in_list ) ) {
		*advance_by -= left_in_list;
		return False;
  }

  const Seg *seg = finger_rightmost_seg( finger, seglist->header->level );
  len_t len_through_seg = finger_travel_len( finger, seglist->header->level );
  for ( level_t level = seglist->header->level; seg != NIL && level >= seglist_min_level; level-- ) {

		/* Advance as far as possible along this level without moving by more than *advance_by length */
		bool_t keep_moving;
		do {

			keep_moving = False;
			const Seg *nxt = seg->fwd[ level ].seg;
			if ( nxt != NIL ) {
				len_t advance_here = seg->fwd[ level ].skip_len + ( nxt->end - nxt->beg );
				len_t new_travel_len = len_through_seg + advance_here;
				if ( get_phys_len( new_travel_len ) < get_phys_len( target_len )  ) {
					seg = nxt;
					len_through_seg = new_travel_len;
					keep_moving = True;
				} 
			}
		} while( keep_moving );
		finger_rightmost_seg( finger, level ) = const_cast<Seg *>(seg);
		finger_travel_len( finger, level ) = len_through_seg;
  }

  Seg *our_seg = seglist_next_seg( finger_rightmost_seg( finger, seglist_min_level ) );

  if ( our_seg == NIL ) {
		*advance_by -= left_in_list;
		return False;
  } else {
		*cur_len = target_len;
		*cur_loc = loc_t( get_loc( our_seg->beg ) + get_phys_len( target_len - finger_travel_len( finger, seglist_min_level ) ) );
		*cur_leafset = our_seg->leafset;
		return True;
  }
}  // seglist_advance_by_len

inline glen_t seg_glen( const Seg *seg ) { return get_gloc( seg->end ) - get_gloc( seg->beg ); }

//
// Func: seglist_find_glen
//
// Find a point on the given seglist corresponding to traveling a given glen along
// the seglist's segs, counting only the seg and not the gaps between them.
//
// Params:
//
//    seglist - the seglist
//    glen - how far to travel along the seglist's segs; must be between 0 and
//       the total glen of the seglist.
//
// See also: <seglist_advance_by_len()>
//
gloc_t seglist_find_glen( const Seglist *seglist, glen_t glen ) {

	assert( !seglist_is_empty( seglist ) );
	assert( glen_t(0.0) <= glen && glen <= seglist_tot_glen( seglist ) );

	const Seg *seg = seglist_header( seglist );
	glen_t glen_through_seg( 0.0 );
	
  for ( level_t level = seglist->header->level; seg != NIL && level >= seglist_min_level; level-- ) {

		/* Advance as far as possible along this level without moving by more than 'glen' */
		bool_t keep_moving;
		do {

			keep_moving = False;
			const Seg *nxt = seg->fwd[ level ].seg;
			if ( nxt != NIL ) {
				glen_t advance_here = get_glen( seg->fwd[ level ].skip_len ) + seg_glen( nxt );
				glen_t new_travel_len = glen_through_seg + advance_here;
				if ( new_travel_len < glen ) {
					seg = nxt;
					glen_through_seg = new_travel_len;
					keep_moving = True;
				}
			}
		} while( keep_moving );
  }

  const Seg *our_seg = seglist_next_seg_const( seg );
	chkCond( our_seg != NIL, "can't happen" );
	return get_gloc( our_seg->beg ) + ( glen - glen_through_seg );
}  // seglist_find_glen

#ifdef __GNUC__
#if ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) )
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#endif
	
/*
 * Func: seglist_increase_seglist_level_to
 *
 * Increase the level of the given seglist to a given level.  Used to make the levels of two
 * seglists equal before doing an operation on them.
 *
 * Params:
 *
 *    seglist - the seglist whose level we're updating
 *    target_level - the level to which we want to raise the seglist's level
 *    finger - if non-NULL, a Finger pointing to Segs of this Seglist;
 *      we ensure that the rightmost seg at any levels we add, is the seglist's header.
 */
static void seglist_increase_seglist_level_to( Seglist *seglist, Finger *finger, level_t target_level ) {
	Seg *header = seglist->header;
	level_t *header_level = &( header->level );
	
	while ( *header_level < target_level ) {
	
		(*header_level)++;
		header->fwd[ *header_level ] = header->fwd[ (*header_level)-1 ];
		if( finger ) {
			finger_rightmost_seg( finger, *header_level ) = header;
			finger_travel_len( finger, *header_level ) = ZERO_LEN;
		}
	}
}

#if defined(__GNUC__) && ( ( __GNUC__ > 4 ) || ( ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ > 5 ) ) )
#pragma GCC diagnostic pop
#endif

/*
 * Func: seglist_split
 *
 * Split a seglist in two at a given location.  The original seglist is destroyed.
 * This routine is used to compute the two seglists resulting
 * from a recombination event, as well as a subroutine in several other seglist operations.
 *
 * Input params:
 *
 *   seglistp - the seglist we're splitting.  It is destroyed by this operation.
 *   loc - the location at which we're splitting the seglist
 *   split_seg - affects what happens when loc falls in the middle of a Seg of this Seglist.  If split_seg is True,
 *     we break up this Seg at loc, adding its left part to the left output seglist and its right part to the right output seglist.
 *     if split_seg is False, we include the segment entirely in the right output seglist, in effect pretending that loc fell
 *     at the beginning of the seg rather than in its middle.
 *
 * Output params:
 *
 *   s_left, s_right - seglists representing the portions of the input seglist to the left and to the right of loc, respectively.
 *
 * Returns: a Finger pointing to the rightmost seg of *s_right at every level.
 *  The same static Finger is always returned, overwritten by each call to seglist_split().
 */
Finger *seglist_split( Seglist **seglistp, loc_t loc, Seglist **s_left, Seglist **s_right, bool_t split_seg, bool *did_split ) {

	
  assert( seglistp );
  Seglist *seglist = *seglistp;
  *seglistp = NULL;
  
  seglist_chk( seglist );

  segsumm_t segsumm_orig = seglist->segsumm;

  Seg *seg = seglist_header( seglist );

  loc_t orig_end = seglist_end( seglist );

  static Finger *finger;
  if ( !finger ) finger = seglist_finger_alloc();

  /*
	 * For each level, find the rightmost seg at that level, that is still fully to the left of loc (except that loc may equal its end).
	 * Also, record the total length of segments over which we travel at each level: this will be needed to set <Seg::skip_len> of the
	 * rightmost Seg of *s_left in each level (i.e. the number of Segs that the pointer going to NIL skips over).
	 */
  for ( level_t level = seglist->header->level; level >= seglist_min_level; level-- ) {
		len_t travel_len = ZERO_LEN;
		
		Seg *nxt;
		while ( ( nxt = seg->fwd[ level ].seg )->end <= loc ) {
			travel_len += ( seg->fwd[ level ].skip_len + ( nxt->end - nxt->beg ) );
			seg = nxt;
		}
		finger_travel_len( finger, level ) = travel_len;
		finger_rightmost_seg( finger, level ) = seg;
  }

  Seg *next_seg = seglist_next_seg( finger_rightmost_seg( finger, seglist_min_level ) );

  /* Invariant: seg is now the rightmost seg at level seglist_min_level, that is still fully to the left of loc. */
  /* So, next_seg->end > loc */

  /* does loc break the successor of seg, or fall in the gap between seg and its successor? */
  if ( split_seg && next_seg->beg < loc ) {

		/*
		 * loc breaks a seg of seglist, next_seg.
		 * Rewrite Seglist to split this seg into two segs, with a boundary at loc --
		 * effectively reducing the analysis to the case where loc falls in a gap between two segs.
		 * Subsequent code can then assume that loc fell in a gap.
		 */

		/* Local var: left_part */
		/* The part of next_seg to the left of loc. */
		Seg *left_part = next_seg;
	 
		/* Local var: right_part */
		/* The part of next_seg to the right of loc. */
		Seg *right_part = seglist_alloc_seg( seglist_new_level() );

		level_t left_level = left_part->level, right_level = right_part->level;
	 
		seglist_increase_seglist_level_to( seglist, finger, right_level+1 );

		loc_t next_seg_end = next_seg->end;
		left_part->end = right_part->beg = loc;
		right_part->end = next_seg_end;
		right_part->leafset = left_part->leafset;
#ifdef COSI_TRACK_LAST_COAL
		right_part->lastCoalGen = left_part->lastCoalGen;
#endif		

		/* If left_part->level > right_part->level, the forward pointers from left_part */
		/* above right_level now skip over right_part.  Record this in their skip_len fields. */
		for ( level_t level = right_level+1; level <= left_level; level++ )
			 left_part->fwd[ level ].skip_len += ( right_part->end - loc );

		/* Update finger to indicate that we now travel through left_part before stopping to the left of loc. */
		len_t skip_len_tot = ZERO_LEN;
		level_t seg_min_level = std::min( left_level, right_level );
		for ( level_t level = seglist_min_level; level <= seg_min_level; level++ )
			 skip_len_tot += finger_travel_len( finger, level );
	 
		finger_travel_len( finger, left_level ) +=
			 ( finger_rightmost_seg( finger, left_level )->fwd[ left_level ].skip_len + ( loc - left_part->beg ) );
		finger_rightmost_seg( finger, left_level ) = left_part;

		for ( level_t level = seglist_min_level; level < left_level; level++ ) {
			finger_rightmost_seg( finger, level ) = left_part;
			finger_travel_len( finger, level ) = ZERO_LEN;
		}

		/* Set forward pointers of right_part to what were the forward pointers of the original next_seg that we broke up, */
		/* up to right_part->level. */
		memcpy( &( right_part->fwd[ seglist_min_level ] ), &( left_part->fwd[ seglist_min_level ] ),
						( seg_min_level - seglist_min_level + 1 ) * sizeof( Segptr ) );

		/* Forward pointers of left_part, up to min(left_level,right_level), point directly to right part. */
		for ( level_t level = seglist_min_level; level <= seg_min_level; level++ ) {
			left_part->fwd[ level ].seg = right_part;
			left_part->fwd[ level ].skip_len = ZERO_LEN;
		}

		for ( level_t level = seg_min_level+1; level <= right_level; level++ ) {

			/* If right_part sticks out above left_part (and the former next_seg), */
			/* adjust the pointers at the higher levels that used to skip over next_seg */
			/* to point to right_part. */
		
			Seg *rightmost_seg = finger_rightmost_seg( finger, level );
			Segptr *rightmost_seg_fwd = &( rightmost_seg->fwd[ level ] );
			right_part->fwd[ level ].seg = rightmost_seg_fwd->seg;
			rightmost_seg_fwd->seg = right_part;

			len_t prev_skip_len = rightmost_seg_fwd->skip_len;

			rightmost_seg_fwd->skip_len = ( skip_len_tot + ( loc - left_part->beg ) );
			right_part->fwd[ level ].skip_len = prev_skip_len - skip_len_tot - ( right_part->end - left_part->beg );

			skip_len_tot += finger_travel_len( finger, level );  
		}

		if ( did_split ) *did_split = True;

		seglist_chk( seglist );
  }  /* if ( split_seg && next_seg->beg < loc ) */
	else { if ( did_split ) *did_split = False; }

  segsumm_t segsumm_left, segsumm_right;
  segsumm_split( &segsumm_orig, loc, &segsumm_left, &segsumm_right );

  len_t skip_len_tot = ZERO_LEN;

  *s_left = seglist;
  *s_right = seglist_make_empty( );
  seglist_increase_seglist_level_to( *s_right, (Finger *)NULL, seglist->header->level );

  /*
	 * Make rightmost pointers of *s_left at each level point to NIL, and make the header
	 * pointers of *s_right point to the former tarets of the rightmost pointers of *s_left.
	 */
  Seg *right_header = (*s_right)->header;
  for ( level_t level = seglist_min_level; level <= seglist->header->level; level++ ) {
		Seg *rseg = finger_rightmost_seg( finger, level );
		right_header->fwd[ level ].seg = rseg->fwd[ level ].seg;
		rseg->fwd[ level ].seg = NIL;
		
		right_header->fwd[ level ].skip_len = rseg->fwd[ level ].skip_len - skip_len_tot;
		rseg->fwd[ level ].skip_len = skip_len_tot;

		skip_len_tot += finger_travel_len( finger, level );
  }
  
  seglist_end( *s_right ) = orig_end;
  seglist_end( *s_left ) = finger_rightmost_seg( finger, seglist_min_level )->end;

  seglist_compact_level( *s_left );
  seglist_compact_level( *s_right );

  if ( seglist_is_empty( *s_left ) ) seglist_end( *s_left ) = MIN_LOC;
  if ( seglist_is_empty( *s_right ) ) seglist_end( *s_right ) = MIN_LOC;

  seglist_chk( *s_left );
  seglist_chk( *s_right );

  (*s_left)->segsumm = segsumm_left;
  (*s_right)->segsumm = segsumm_right;

  return finger;

}  /* seglist_split() */

// Func: seglist_find_seg
// Return the seg within the seglist that contains the given location,
// or NULL if no such seg.
const Seg *seglist_find_seg( const Seglist *seglist, loc_t loc ) {
  const Seg *seg = seglist_header( seglist );
  /*
	 * For each level, find the rightmost seg at that level, that is still fully to the left of loc (except that loc may equal its end).
	 */
  for ( level_t level = seglist->header->level; level >= seglist_min_level; level-- ) {
		const Seg *nxt;
		while ( ( nxt = seg->fwd[ level ].seg )->end <= loc )
			 seg = nxt;
  }
  const Seg *next_seg = seglist_next_seg_const( seg );
	
  /* Invariant: seg is now the rightmost seg at level seglist_min_level, that is still fully to the left of loc. */
  /* So, next_seg->end > loc */
	assert( next_seg && next_seg->end > loc );
	return next_seg->beg < loc ? next_seg : NULL;
}

Finger *seglist_finger_alloc(void) { return (Finger *)mempool_finger.mempool_alloc_item(); /*COSI_MALLOC( Finger_elt, seglist_max_level + seglist_min_level );*/ }
void seglist_finger_free( Finger *finger ) { mempool_finger.mempool_free_item( finger ); /*free( finger );*/ }

/***************************
 * Section: Seglist builder
 **************************/

/*
 * Struct: seglist_builder
 *
 * An object used when building up a Seglist.   Keeps the fledgling Seglist, as well
 * as a <Finger> pointing to the rightmost <Seg> at each level.  Has methods for appending
 * one <Seg> or an entire <Seglist> to the end of the Seglist being built.
 */
typedef struct seglist_builder {
	 /* Field: seglist */
	 /* The Seglist we're building up. */
	 Seglist *seglist;
  
	 /* Field: finger */
	 /* For each level, the rightmost <Seg> in <seglist> at that level */
	 /* (which may be seglist->header if there are no segs at that level). */
	 /* If <finger_at_end> is False, the pointers at some level may be "behind" and may */
	 /* not point to the rightmost Seg at that level. */
	 Finger *finger;

	 /* Field: finger_at_end */
	 /* If True, <finger> correctly points to the rightmost Seg at each level. */
	 /* If False, <finger> may point to earlier segs and needs to be updated */
	 /* (by a call to <seglist_builder_ensure_finger_at_end()>) before append operations */
	 /* on this Seglist_builder will work correctly. */
	 bool_t finger_at_end;
} Seglist_builder;

ostream& operator<<( ostream& s, const Seglist_builder& b );

ostream& operator<<( ostream& s, const Seglist_builder& b ) {
	s << "{B:" << *b.seglist << "}";
	return s;
}

void seglist_chk_builder_helper( const Seglist_builder *builder, const char *fname, int line );

/*
 * Func: seglist_builder_init
 *
 * Initialize a Seglist_builder.  Must be called before any segs or seglists are appended to the builder.
 */
static void seglist_builder_init( Seglist_builder *builder ) {
  builder->seglist = seglist_make_empty( );
  builder->finger = seglist_finger_alloc();
  finger_rightmost_seg( builder->finger, seglist_min_level ) = builder->seglist->header;
  builder->finger_at_end = True;
}

/*
 * Func: seglist_builder_ensure_finger_at_end
 *
 * Ensure that builder->finger points to the rightmost seg of builder->seglist at each level.
 */
static void seglist_builder_ensure_finger_at_end( Seglist_builder *builder ) {
  if ( !builder->finger_at_end ) {
		Finger *finger = builder->finger;
		level_t the_level = builder->seglist->header->level;
		Seg *seg = finger_rightmost_seg( finger, the_level );
		assert( seg != NIL );
	 
		for( level_t level = the_level; level >= seglist_min_level; level-- ) {
			Seg *nxt;
			while( ( nxt = seg->fwd[ level ].seg ) != NIL )
				 seg = nxt;
			finger_rightmost_seg( finger, level ) = seg;
		}
		builder->finger_at_end = True;
  }
}

/*
 * Func: seglist_builder_free_and_get_result
 *
 * Free any auxiliary storage allocated for the builder (in particular, builder->finger), and return
 * the seglist the builder has been constructing.
 */
static Seglist *seglist_builder_free_and_get_result( Seglist_builder *builder ) {
  seglist_finger_free( builder->finger );
  return builder->seglist;
}


/*
 * Static func: seglist_builder_add_existing_seg
 *
 * Append a <Seg> to the end of the given list.  The seg must be to the right of the end of the list.
 * This routine takes an existing Seg object and incorporates it into the list, rather than allocating a new Seg.
 */
static void seglist_builder_add_existing_seg( Seglist_builder *builder, Seg *new_seg ) {
  seglist_chk_builder( builder );

  assert( new_seg != NIL );
  assert( new_seg->beg >= builder->seglist->getEnd() );

  seglist_builder_ensure_finger_at_end( builder );

  bool_t free_new_seg = False;

  /* See if the new seg can be merged with the last seg of builder->seglist */
  if ( !seglist_is_empty( builder->seglist ) && new_seg->beg == builder->seglist->getEnd() &&
			  !leafset_may_differ( new_seg->leafset,
														 finger_rightmost_seg( builder->finger, seglist_min_level )->leafset )
				  ) {
		Seg *last_seg = finger_rightmost_seg( builder->finger, seglist_min_level );
		builder->seglist->setEnd( last_seg->end = new_seg->end );
		free_new_seg = True;
  } else {
		/* Ok, cannot merge the new seg with the last seg of builder->seglist, so add it at the end of the seglist */

		/* if new_seg has higher level than the list, increase the list level */
		seglist_increase_seglist_level_to( builder->seglist, builder->finger, new_seg->level+1 );
 
		/* for rightmost segs up to new_seg->level, set them to point to new_seg instead of NIL */
		/* and move the finger at that level to new_seg */
		level_t level = seglist_min_level;
		while ( level <= new_seg->level ) {
			new_seg->fwd[ level ].seg = NIL;
			new_seg->fwd[ level ].skip_len = ZERO_LEN;
	 
			finger_rightmost_seg( builder->finger, level )->fwd[ level ].seg = new_seg;
			finger_rightmost_seg( builder->finger, level ) = new_seg;
			level++;
		}
  }

  {
		/* for rightmost segs above new_seg->level, they still point to NIL, but add new_seg's length */
		/* to their skip_len. */
		len_t new_seg_len = ( new_seg->end - new_seg->beg );
		level_t header_level = builder->seglist->header->level;
		level_t level = finger_rightmost_seg( builder->finger, seglist_min_level )->level+1;
		while( level <= header_level ) {
			finger_rightmost_seg( builder->finger, level )->fwd[ level ].skip_len += new_seg_len;
			level++;
		}
  }

  seglist_end( builder->seglist ) = new_seg->end;
  if ( free_new_seg ) { seglist_free_seg( new_seg ); IF_COSI_STATS( n_free_adj++ ); }
  seglist_chk_builder( builder );
}  /* seglist_builder_add_existing_seg() */

/*
 * Func: seglist_builder_add_seg
 *
 * Add a specified seg to the end of the builder's seglist.
 */
static void seglist_builder_add_seg( Seglist_builder *builder, loc_t beg, loc_t end, leafset_p leafset, genid lastCoalGen ) {
  seglist_chk_builder( builder );
  if ( end > beg ) {
		Seg *new_seg = seglist_alloc_seg( seglist_new_level() );
	 
		new_seg->beg = beg;
		new_seg->end = end;
		seglist_seg_set_leafset( new_seg, leafset );
		new_seg->lastCoalGen = lastCoalGen;
	 
		seglist_builder_add_existing_seg( builder, new_seg );
  }
}

/*
 * Func: seglist_builder_add_seglist
 *
 * Append a given seglist to builder->seglist.  If seglist_end_finger is given, it
 * must point to the rightmost seg of seglist at each level.
 */
static void seglist_builder_add_seglist( Seglist_builder *builder, Seglist *seglist,
																				 Finger *seglist_end_finger ) {

  if ( !seglist_is_empty( seglist ) ) {

		seglist_builder_ensure_finger_at_end( builder );
	 
		seglist_chk_builder( builder );
		seglist_chk( seglist );
		Seglist *s1 = builder->seglist, *s2 = seglist;

		/* make the two lists the same level */
		seglist_increase_seglist_level_to( s1, builder->finger, s2->header->level );
		seglist_increase_seglist_level_to( s2, /* finger= */ seglist_end_finger, s1->header->level );
		seglist_chk( s1 );
		seglist_chk( s2 );

		level_t the_level = builder->seglist->header->level;
		Finger *finger = builder->finger;
		for ( level_t level = seglist_min_level; level <= the_level; level++ ) {
			finger_rightmost_seg( finger, level )->fwd[ level ].seg = s2->header->fwd[ level ].seg;

			finger_rightmost_seg( finger, level )->fwd[ level ].skip_len += s2->header->fwd[ level ].skip_len;
		}
		if ( seglist_end_finger ) {
			for ( level_t level = seglist_min_level; level <= the_level && finger_rightmost_seg( seglist_end_finger, level ) != seglist->header;
						level++ )
				 finger_rightmost_seg( finger, level ) = finger_rightmost_seg( seglist_end_finger, level );
		} else
			 builder->finger_at_end = False;
		seglist_end( builder->seglist ) = seglist_end( seglist );

		seglist_chk_builder( builder );
  
		seglist->header->level = seglist_min_level;
		seglist->header->fwd[ seglist_min_level ].seg = NIL;
  }
  seglist_free_seglist( seglist );
  seglist_chk_builder( builder );	 
}

/************************************************************************/

/*
 * Func: seglist_remove_head_seg
 *
 * Remove the first seg from a non-empty seglist, and return it.
 */
static Seg *seglist_remove_head_seg( Seglist *seglist ) {
  seglist_chk( seglist );
  Seg *head = seglist_first_seg( seglist );
  assert( head != NIL );
  level_t level = seglist_min_level;
  while ( level <= head->level ) {
		seglist->header->fwd[ level ] = head->fwd[ level ];
		level++;
  }
  {
		len_t head_len = head->end - head->beg;
		while( level <= seglist->header->level ) {
			seglist->header->fwd[ level ].skip_len -= head_len;
			level++;
		}
  }
  seglist_compact_level( seglist );
  if ( seglist_is_empty( seglist ) ) seglist->setEnd( MIN_LOC );
  seglist_chk( seglist );
  return head;
}

/*
 * Func: seglist_move_beg
 *
 * Move the beginning of the first seg of a non-empty seglist either forward or backwards.
 * If moving forward, must not move past the end of the first segment.
 */
void seglist_move_beg( Seglist *seglist, loc_t new_beg ) {
  seglist_chk( seglist );

  Seg *seg = seglist_first_seg( seglist );
  assert( seg != NIL && new_beg < seg->end );
  len_t advance = new_beg - seg->beg;  // may be negative!
  if ( !is_zero( advance ) ) {
		seg->beg = new_beg;
		// For any forward pointers from the seglist header that skipped over the first
		// segment, adjust the amount of length that they skip.
		for ( level_t level = seg->level+1; level <= seglist->header->level; level++ )
			 seglist->header->fwd[ level ].skip_len -= advance;
  }

  seglist_chk( seglist );
}

/*
 * Func: seglist_move_end
 *
 * Move the end of the last seg of a non-empty seglist either forward or backwards.
 * If moving backward, must not move past the beginning of the last segment.
 *
 * Params:
 *
 *   seglist - the seglist
 *   finger - a <Finger> pointing to the rightmost seg of 'seglist' at every level,
 *      such as returned by <seglist_split()>.
 *   new_end - the new endpoint
 */
void seglist_move_end( Seglist *seglist, Finger *finger, loc_t new_end ) {
  seglist_chk( seglist );

  Seg *last_seg = finger_rightmost_seg( finger, seglist_min_level );
  assert( last_seg != NIL && new_end > last_seg->beg );
  len_t advance = new_end - last_seg->end;  // may be negative!
  if ( !is_zero( advance ) ) {
		seglist->setEnd( last_seg->end = new_end );
		
		// For any forward pointers that skipped over the last
		// segment, adjust the amount of length that they skip.
		for ( level_t level = last_seg->level+1; level <= seglist->header->level; level++ )
			 finger_rightmost_seg( finger, level )->fwd[ level ].skip_len += advance;
  }

  seglist_chk( seglist );
}  // seglist_move_end

/*
 * Func: seglist_union_or_inters
 *
 * Computes the union or the intersection of two seglists.  The input seglists are destroyed.
 *
 * Input params:
 *
 *    seglist1p, seglist2p - the input seglists
 */
Seglist *seglist_union( Seglist **seglist1p, Seglist **seglistp,
												genid gen,
												seglist_seg_union_callback_t seg_union_callback ) {
	

  /*printf( "union_or_inters: %s got %p %p %s %s\n",
		is_union ? "union" : "inters", seglist1p, seglistp, seglist_str( *seglist1p ), seglist_str( *seglistp ) );*/

  assert( seglist1p && seglistp );
  Seglist *seglist1 = *seglist1p, *seglist = *seglistp;
  *seglist1p = *seglistp = NULL;
  
  Seglist_builder builder;

  seglist_chk( seglist1 );
  seglist_chk( seglist );

  //loc_t orig_end = !is_union ? MIN_LOC : std::max( seglist1->end, seglist->end );
  
  seglist_builder_init( &builder );

  segsumm_union( &( builder.seglist->segsumm ), &( seglist1->segsumm ), &( seglist->segsumm ) );

  seglist_chk_builder( &builder );

  Seglist **seglists[2] = { &seglist1, &seglist };

  if ( !segsumm_is_empty( &( builder.seglist->segsumm ) ) ) {

		while( !seglist_is_empty( *seglists[0] ) && !seglist_is_empty( *seglists[1] ) ) {
			seglist_chk( *seglists[0] );
			seglist_chk( *seglists[1] );

			Seg *first_seg[2] = { (*seglists[0])->header->fwd[ seglist_min_level ].seg, (*seglists[1])->header->fwd[ seglist_min_level ].seg };
			/* Local var: which */
			/* The index (0 or 1) of the input seglist that starts earlier (leftmost) on the chromosome. */
			int which = first_seg[0]->beg < first_seg[1]->beg ? 0 :
				 ( first_seg[0]->beg > first_seg[1]->beg ? 1 : ( first_seg[0]->end > first_seg[1]->end ? 0 : 1 ) );
			int other = 1 - which;

			/* add to the output any segs in seglists[which] that lie completely to the left of seglists[other]. */
			Seglist *bef, *aft;

			Finger *bef_end_finger = seglist_split( seglists[ which ],
																							seglist_beg( *seglists[ other ] ), &bef, &aft, /* split_seg= */ False );


			seglist_builder_add_seglist( &builder, bef, bef_end_finger );

			*seglists[ which ] = aft;

			if ( !seglist_is_empty( aft ) ) {

				if ( seglist_beg( *seglists[ which ] ) > seglist_beg( *seglists[ other ] ) ) continue;

				{  /* is_union */
					bool_t in_seg = True;
					Seg *cur_seg = seglist_first_seg( *seglists[ which ] );
					loc_t end = cur_seg->beg;
		
					/* process any segs in other list lying fully within this seg */
					Seg *other_seg;
					while ( !seglist_is_empty( *seglists[ other ] ) &&
									( other_seg = seglist_first_seg( *seglists[ other ] ) )->end <= cur_seg->end ) {
						assert( other_seg->beg >= cur_seg->beg );
						seglist_chk( *seglists[which] );
						seglist_chk( *seglists[other] );

						seglist_chk_builder( &builder );
						seglist_builder_add_seg( &builder, end, other_seg->beg, cur_seg->leafset, cur_seg->lastCoalGen );
						seglist_chk_builder( &builder );
						other_seg = seglist_remove_head_seg( *seglists[ other ] );
						bool_t dropSeg = False;
						/* if ( other_seg->leafset && cur_seg->leafset && !leafset_is_subset( other_seg->leafset, cur_seg->leafset ) ) */ {
							{
								if ( seg_union_callback ) {
									seg_union_callback( other_seg->beg, other_seg->end, other_seg->lastCoalGen, gen, other_seg->leafset );
									seg_union_callback( other_seg->beg, other_seg->end, cur_seg->lastCoalGen, gen, cur_seg->leafset );
								}
								 seglist_seg_set_leafset( other_seg, leafset_union( other_seg->leafset, cur_seg->leafset ) );
								 if ( leafset_is_full( other_seg->leafset ) ) {
									 dropSeg = True;
#ifdef COSI_DEV_PRINT_COALESCED_LEN									 
									 static len_t totCoalLen(0.0);
									 static len_t lastPrintedLen(0.0);
									 static const len_t printInterval(0.2);
									 totCoalLen += ( other_seg->end - other_seg->beg );
									 if ( totCoalLen - lastPrintedLen > printInterval ) {
										 PRINT2( totCoalLen, gen );
										 lastPrintedLen = totCoalLen;
									 }
#endif									 
								 }
#ifdef COSI_TRACK_LAST_COAL
								 other_seg->lastCoalGen = gen;
#endif								 
								 
							}

							seglist_chk_builder( &builder );
							end = other_seg->end;
							if ( dropSeg ) { seglist_free_seg( other_seg ); }
							else {
								seglist_builder_add_existing_seg( &builder, other_seg );
								seglist_chk_builder( &builder );
							}
						} /*else
								seglist_free_seg( other_seg );*/

						/* so, actually, if not dealing with leafsets, can also do this step in one step -- just cut the relevant part
							 of seglists[other].
						*/
					}  /* while other_seg lies fully within cur_seg */

					/* If there is an overlapping seg in other list, switch to that and continue; else, done with this segment. */
					if ( ( in_seg = ( !seglist_is_empty( *seglists[ other ] ) && other_seg->beg < cur_seg->end ) ) ) {
						seglist_chk_builder( &builder );
						seglist_builder_add_seg( &builder, end, other_seg->beg, cur_seg->leafset, cur_seg->lastCoalGen );
						seglist_chk_builder( &builder );

						seglist_move_beg( *seglists[ which ], other_seg->beg );
		  
						which = 1-which; other = 1-other;
						end = cur_seg->end;
						cur_seg = other_seg;
					} else {
						seglist_chk_builder( &builder );
						seglist_builder_add_seg( &builder, end, cur_seg->end, cur_seg->leafset, cur_seg->lastCoalGen );  /* reuse if can */
						seglist_chk_builder( &builder );
#ifndef NDEBUG						
						Seg *a_seg =
#endif							 
							 seglist_remove_head_seg( *seglists[ which ] );
						assert( a_seg == cur_seg );

						seglist_free_seg( cur_seg );
						IF_COSI_STATS( n_free_doneseg++ );
					}  /* no overlapping seg */
				} /* is_union */
			}  /* if ( !seglist_is_empty( aft ) ) */
		}  /* while both input seglists are non-empty */

		seglist_chk_builder( &builder );
  }
  
  {

		/* add leftover segs from one of the lists, if any */
		for ( int which = 0; which < 2; which++ ) {
			if ( !seglist_is_empty( *seglists[which] ) ) {
				seglist_chk_builder( &builder );
				seglist_builder_add_seglist( &builder, *seglists[which], /* end_finger= */ (Finger *)NULL );
				seglist_chk_builder( &builder );
			} else
				 seglist_free_seglist( *seglists[ which ] );
		}
  } 

  seglist_chk_builder( &builder );

  return seglist_builder_free_and_get_result( &builder );
}  // seglist_union

/*
 * Func: seglist_dup
 *
 * Create a separate, independent copy of a given seglist.
 */
Seglist *seglist_dup( const Seglist *seglist ) {
  Seglist_builder builder;
  seglist_builder_init( &builder );

  for ( const Seg *seg = seglist_first_seg_const( seglist ); seg != NIL; seg = seg->fwd[ seglist_min_level ].seg )
		 seglist_builder_add_seg( &builder, seg->beg, seg->end, seg->leafset, seg->lastCoalGen );

  segsumm_copy( &( (builder.seglist)->segsumm ), &( seglist->segsumm ) );
  
  return seglist_builder_free_and_get_result( &builder );
}

void seglist_canonicalize( Seglist **seglist ) {
  Seglist *result = seglist_dup( *seglist );
  seglist_free_seglist( *seglist );
  *seglist = result;
}

/*
 * Func: seglist_make_one_seg_list
 *
 * Construct a Seglist consisting of one specified seg.
 */
static Seglist *seglist_make_one_seg_list( loc_t beg, loc_t end, leafset_p leafset, genid lastCoalGen ) {
  assert( end > beg );
  Seglist_builder builder;
  seglist_builder_init( &builder );
  seglist_builder_add_seg( &builder, beg, end, leafset, lastCoalGen );
  return seglist_builder_free_and_get_result( &builder );
}

/* Func: seglist_make_full */
/* Create a seglist consisting one one [0,1] seg covering the entire chromosome. */
Seglist *seglist_make_full( leafset_p leafset, genid lastCoalGen ) {
  Seglist *seglist = seglist_make_one_seg_list( MIN_LOC, MAX_LOC, leafset, lastCoalGen );
  segsumm_init_full_seg( &( seglist->segsumm ) );
  return seglist;
}

/*
 * Params:
 *
 *    for_sentinels - if False, if one or both resulting lists are empty, seglist is left untouched;
 *       if True, seglist is always consumed and split into inside and outside.
 */
bool_t seglist_intersect_for_gc( Seglist **seglist, loc_t loc1, loc_t loc2, Seglist **inside, Seglist **outside,
																 bool_t for_sentinels ) {
  assert( MIN_LOC <= loc1 && loc1 < loc2 && loc2 <= MAX_LOC );
  *inside = *outside = NULL;
  if ( seglist_is_empty( *seglist ) ) {
		if ( for_sentinels ) {
			*inside = *seglist;  *seglist = NULL;
			*outside = seglist_make_empty(  );
		}
		n_seglist_empty++;
		return False;
  } else if ( loc2 <= seglist_beg( *seglist ) || seglist_end( *seglist ) <= loc1 ) {
		if ( for_sentinels ) {
			*inside = seglist_make_empty( );
			*outside = *seglist;
			*seglist = NULL;
		}
		n_left_or_right++;
		return False;
  } else if ( loc1 <= seglist_beg( *seglist )  &&  loc2 >= seglist_end( *seglist ) ) {
		if ( for_sentinels ) {
			*inside = *seglist;
			*outside = seglist_make_empty( );
			*seglist = NULL;
		}
		n_containing++;
		return False;
  } else {

		segsumm_t s_inside, s_outside;
#if 0	 
		if ( !segsumm_split_on_seg( &( (*seglist)->segsumm ), loc1, loc2, &s_inside, &s_outside ) ) {
			if ( for_sentinels ) { *inside = seglist_make_empty( ); *outside = *seglist; *seglist = NULL; }
			return False;
		}
#endif
		segsumm_init_empty_seg( &s_inside );
		segsumm_init_empty_seg( &s_outside );
		segsumm_split_on_seg( &( (*seglist)->segsumm ), loc1, loc2, &s_inside, &s_outside );

		if ( segsumm_is_empty( &s_inside ) ) {
			if ( for_sentinels ) { *outside = *seglist; *seglist = NULL; *inside = seglist_make_empty(  ); }
			return False;
		}

		{
			Seg *seg = ( *seglist )->header, *nxt = NULL;
			static Seg *rightmost_segs[ seglist_max_level+1 ];
			for ( level_t level = (*seglist)->header->level; level >= seglist_min_level; level-- ) {
				while ( seg != NIL && ( nxt = seg->fwd[ level ].seg )->end <= loc1 )
					 seg = nxt;
				rightmost_segs[ level ] = seg;
			}
			if ( seg == NIL  ||  loc2 <= ( nxt = seglist_next_seg( seg ) )->beg ) {
				if ( for_sentinels ) { *outside = *seglist; *seglist = NULL; *inside = seglist_make_empty(  ); }
				return False;
			} else if ( loc2 < nxt->end && loc1 <= nxt->beg ) {
				n_oneseg++;
				*inside = seglist_make_one_seg_list( std::max( nxt->beg, loc1 ), loc2,
																						 nxt->leafset,
																						 nxt->lastCoalGen );

				len_t inner_len = loc2 - nxt->beg;
				for ( level_t level = nxt->level+1; level <= (*seglist)->header->level; level++ ) {
					assert( rightmost_segs[ level ] != NIL );
					assert( rightmost_segs[ level ]->level > nxt->level );
					rightmost_segs[ level ]->fwd[ level ].skip_len -= inner_len;
				}
		  
				nxt->beg = loc2;
				assert( nxt->beg < nxt->end );

				/* now must subtract ( loc2 - nxt->beg ) from any skip edges that skip over nxt, including at least the top skip. */
				/* only if !for_sentinels, of course, but that's the most time-consuming step. */
		  
				*outside = *seglist; *seglist = NULL;
				segsumm_copy( &( (*inside)->segsumm ), &s_inside );
				segsumm_copy( &( (*outside)->segsumm ), &s_outside );

				seglist_chk( (*outside) );
		  
				return True;
			}
		}

		/*clock_t t_start = clock();*/

		/* [   loc1   loc2    ] */
		/*   a      b      c    */
		Seglist *seglist_ab, *seglist_a, *seglist_c;
	 
		seglist_split( seglist, loc2, &seglist_ab, &seglist_c, /* split_seg= */ True );
		Finger *a_end_finger = seglist_split( &seglist_ab, loc1, &seglist_a, inside, /* split_seg= */ True );

		Seglist_builder builder;
		seglist_builder_init( &builder );
		seglist_builder_add_seglist( &builder, seglist_a, a_end_finger );
		seglist_builder_add_seglist( &builder, seglist_c, (Finger *)NULL );
	 
		*outside = seglist_builder_free_and_get_result( &builder );

		/*time_gc += ( clock() - t_start );*/

		if ( seglist_is_empty( *inside ) ) {
			if ( !for_sentinels ) { *seglist = *outside; seglist_delete( inside ); }
			n_in_gap++;
			if ( segsumm_is_empty( &s_inside ) ) n_segsumm++;
			return False;
		}
		assert( !seglist_is_empty( *outside ) );
#if 0	 
		if ( seglist_is_empty( *outside ) ) {
			if ( !for_sentinels ) { *seglist = *inside; seglist_delete( outside ); }
			return False;
		}
#endif

		segsumm_copy( &( (*inside)->segsumm ), &s_inside );
		segsumm_copy( &( (*outside)->segsumm ), &s_outside );

		n_nonempty++;
	 
		seglist_chk( (*outside) );
	 
		return True;
  }
}  /* seglist_intersect_for_gc() */

/*
 * Func: seglist_tot_len
 *
 * Return the sum of lengths of all segs in a seglist.
 */
len_t seglist_tot_len( const Seglist *seglist ) {
  return seglist->header->fwd[ seglist->header->level ].skip_len;
}

/*
 * Func: seglist_tot_glen
 *
 * Return the sum of lengths of all segs in a seglist.
 */
glen_t seglist_tot_glen( const Seglist *seglist ) {
  return get_glen( seglist->header->fwd[ seglist->header->level ].skip_len );
}

/***********************************************************************************************/

/* Section: Debugging and testing code */

/*
 * Func: seglist_equal
 *
 * Returns True if the two Seglists are identical.
 */
bool_t seglist_equal( const Seglist *seglist1, const Seglist *seglist ) {
  seglist_chk( seglist1 );
  seglist_chk( seglist );

  if ( !segsumm_equal( &( seglist1->segsumm ), &( seglist->segsumm ) ) ) return False;
  const Seg *seg1 = seglist_first_seg_const( seglist1 );
  const Seg *seg2 = seglist_first_seg_const( seglist );

  while( seg1 != NIL && seg2 != NIL ) {
		if ( !( seg1->beg == seg2->beg && seg1->end == seg2->end && ( !leafset_may_differ( seg1->leafset, seg2->leafset ) ) ) )
			 return False;
		seg1 = seg1->fwd[ seglist_min_level ].seg;
		seg2 = seg2->fwd[ seglist_min_level ].seg;
  }
  return seg1 == NIL && seg2 == NIL;
}

/*
 * Func: seglist_equal_dbg
 *
 * Returns True if the two Seglists are identical; print any differences between them.
 */
void seglist_equal_dbg( const Seglist *seglist1, const Seglist *seglist ) {
  seglist_chk( seglist1 );
  seglist_chk( seglist );
  
  if ( !segsumm_equal( &( seglist1->segsumm ), &( seglist->segsumm ) ) )
		 std::cerr << "\nsegsumms differ: " << seglist1->segsumm << " " << seglist->segsumm;
  
  const Seg *seg1 = seglist_first_seg_const( seglist1 );
  const Seg *seg2 = seglist_first_seg_const( seglist );
  
  int seg_idx = 0;
  while( seg1 != NIL && seg2 != NIL ) {
		if ( !( seg1->beg == seg2->beg ) ) { PRINT3( seg_idx, seg1->beg, seg2->beg ); }
		if ( !( seg1->end == seg2->end ) ) { PRINT3( seg_idx, seg1->end, seg2->end ); }
	 
		seg1 = seg1->fwd[ seglist_min_level ].seg;
		seg2 = seg2->fwd[ seglist_min_level ].seg;
		seg_idx++;
  }
  
  if( seg1 != NIL ) printf( "list 1 not done\n" );
  if( seg2 != NIL ) printf( "list 2 not done\n" );
}

#if 0	
static unsigned long seglist_ptr( const void *p ) { return ((unsigned long)(p ? p : (const Seg *)1)); }

/*
 * Func: seglist_write_dot_helper
 *
 * Write out a representation of a seglist as a dot file for Graphviz rendering.  If a non-NULL finger is given, includes a representation
 * of the finger and the segs it points to.
 */
void seglist_write_dot_helper( const Seglist *seglist, const Finger *finger, const char *fname, const char *fn, int line ) {
  printf( "writing dot to %s at %s:%d\n", fname, fn, line );
  /*seglist_chk( seglist );*/
  FILE *out = fopen( fname, "wt" );
  if ( out ) {
		fprintf( out, "digraph seglist {\n" );
		fprintf( out, "  rankdir=LR;\n" ); 

		Seg *seg = seglist->header;
		while ( True ) {
			/*		printf( "seg: level=%d\n", seg->level );*/
		
			fprintf( out, "%lu [ shape=box,label=<<TABLE>", seglist_ptr( seg ) );
			for ( level_t level = seg->level; level >= seglist_min_level; level-- )
				 fprintf( out, "<TR><TD PORT=\"%d\"></TD></TR>", level );
			fprintf( out, "<TR><TD>" );
			if ( seg == seglist->header ) fprintf( out, "HEAD" );
			else if ( seg == NIL ) fprintf( out, "NIL" );
			else fprintf( out, "" LOC_FMT "," LOC_FMT "", seg->beg, seg->end );
			fprintf( out, "</TD></TR>");
			if ( seg == seglist->header ) fprintf( out, "<TR><TD>" LOC_FMT "</TD></TR>", seglist->getEnd() );
			fprintf( out, "</TABLE>> ]\n" );

			if ( seg == NIL ) break;
			for ( level_t level = seg->level; level >= seglist_min_level; level-- )
				 fprintf( out, "%lu:%d -> %lu:%d [ label=\"" LOC_FMT "\" ]\n", seglist_ptr( seg ), level, seglist_ptr( seg->fwd[ level ].seg ), level,
									seg->fwd[ level ].skip_len );
		
			seg = seg->fwd[ seglist_min_level ].seg;
		}

		if ( finger ) {

			fprintf( out, "%lu [ shape=box,label=<<TABLE>", seglist_ptr( finger ) );
			for ( level_t level = seglist->header->level; level >= seglist_min_level; level-- )
				 fprintf( out, "<TR><TD PORT=\"%d\"></TD></TR>", level );
			fprintf( out, "<TR><TD>FINGER</TD></TR></TABLE>>]\n" );

			for ( level_t level = seglist->header->level; level >= seglist_min_level; level-- ) {
				unsigned long from_node = seglist_ptr( finger );
				unsigned long to_node = seglist_ptr( finger[ level ].seg );
				len_t skip_len = finger[ level ].skip_len;

				fprintf( out, "%lu", from_node );
				fprintf( out, ":" );
				fprintf( out, "%d", level );
				fprintf( out, " -> " );
				fprintf( out, "%lu", to_node );
				fprintf( out, "[ label=\"");
				fprintf( out, "" LOC_FMT "", skip_len );
				fprintf( out, "\" ]\n" );
#if 0
				fprintf( out, "%lu:%d -> %lu [ label=\"" LOC_FMT "\" ]\n", from_node, level, to_node, 
								 skip_len );
#endif		  
			}
		}
	 
		fprintf( out, "}\n" );
		fclose( out );
  }
}
#endif	// #if 0

ostream& operator<<( ostream& buf, const Seg *seg ) {
	buf << "[Seg: " << seg->beg << "-" << seg->end << "|" << seg->leafset;
	IF_COSI_TRACK_LAST_COAL( buf << "|" << seg->lastCoalGen );
	buf << "]";
	return buf;
}

ostream& operator<<( ostream& buf, const Seglist *seglist ) {
  buf << "[";
  const Seg *seg = seglist_first_seg_const( seglist );
  bool_t is_first = True;
  while( seg != NIL ) {
		if ( !is_first ) buf << " ";
		buf << seg;
		seg = seg->fwd[ seglist_min_level ].seg;
		is_first = False;
  }
	buf << "]";
  return buf;
}



void seglist_segsumm_chk_helper( const Seglist *my_s, const char *fname, int line ) {

  if ( !my_s ) printf( "what? %s:%d\n", fname, line );
  assert( my_s );
  if ( my_s ) {

		/*	 printf( "checking: " ); seglist_print( s );*/
		BOOST_FOREACH( const Seg& seg, *my_s ) {
			seglet_idx_t i_beg = segsumm_get_seglet( seg.getBeg() );
			seglet_idx_t i_end = segsumm_get_seglet( seg.getEnd() );
			for ( int i = i_beg; i <= i_end; i++ ) {
				segsumm_t my_s_copy = my_s->segsumm;
				if ( segsumm_is_seglet_definitely_empty( &my_s_copy, i ) ) {
					printf( "\nchego? i=%d i_beg=%d i_end=%d beg=" LOC_FMT " end=" LOC_FMT " %s:%d\n", i, i_beg, i_end, get_loc( seg.getBeg() ), get_loc( seg.getEnd() ),
									fname, line );
					std::cerr << "seglist: " << my_s;
					exit(1);
				}
			}
		}
  }
}

/*
 * Func: seglist_chk_helper
 *
 * Checks the sanity of a seglist, causes an assert failure if an inconsistency is found.
 */
#ifndef NDEBUG
void seglist_chk_helper( const Seglist *seglist, const char *fname, int line ) {
  asrt( seglist );
  asrt( seglist->header );
  asrt( seglist->header->level <= seglist_max_level );
	
  /* Local var: seg_end */
  /* Keeps the rightmost end of a Seg we have seen -- used to check that <Seglist::end> field is set correctly. */
  loc_t seg_end = MIN_LOC;
  for ( Seg *seg = seglist->header; seg != NIL; seg = seg->fwd[ seglist_min_level ].seg ) {
		asrt( seg );
		if ( seg != seglist->header ) {
			asrt( seglist_min_level <= seg->level );
			asrt( seg->level < seglist->header->level );
			asrt( MIN_LOC <= seg->beg );
			asrt( seg->end <= MAX_LOC );
			if ( !( seg->beg < seg->end ) ) {
				PRINT3( seg->beg, seg->end, seg->end - seg->beg );
			}
			asrt( seg->beg < seg->end );
			asrt( seg->end <= seglist_next_seg( seg )->beg );
			asrt( seg->leafset != LEAFSET_NULL );
			seg_end = seg->end;
		}
	 
		for ( level_t level = seglist_min_level; level <= seg->level; level++ ) {
			/* check that the forward pointer at this level does point to the next seg at this level or higher. */
			/* check that skip_len is correct */
			len_t skip_len = ZERO_LEN;
			bool_t points_to_next = False;
			for ( Seg *seg2 = seg->fwd[ seglist_min_level ].seg; seg2 != NIL &&
							 !( points_to_next = ( seg2 == seg->fwd[ level ].seg ) ); seg2 = seg2->fwd[ seglist_min_level ].seg ) {
				asrt( seg2->level < level );
				skip_len += ( seg2->end - seg2->beg );
			}
			asrt( ( seg->fwd[ level ].seg == NIL ) || points_to_next );
			asrt( equal_eps( skip_len, seg->fwd[ level ].skip_len ) );
		}
	 
  }
  
  if ( !( seglist->getEnd() == seg_end ) ) {
		PRINT2( seglist->getEnd(), seg_end );
  }
  
  asrt( seglist->getEnd() == seg_end );
  
  /*
	  invariants:

		- segsumm is a correct overapproximation, both the summ and inv_summ parts.
		  
		- header level is greater than level of any node in the list
		- no cycles
		- each list ends at the NIL elt
		- segs go in order and do not overlap
		- there are no zero-length segs

		- skip_len correctly represents the sum of skipped items

		- each level is a sublist of the level below
	  
		- in the header, everything above level points to NIL,
		at the level it points to (not) NIL (unless the list is empty?)
		below -- depending on the policy.

		- no empty segs, though segs may be completely adjacent (but only if their leafsets aren't equal)
	*/
  
}

/*
 * Func: seglist_chk_builder_helper
 *
 * Doa sanity-check on a <Seglist_builder>.
 */
void seglist_chk_builder_helper( const Seglist_builder *builder, const char *fname, int line ) {
  seglist_chk_helper( builder->seglist, fname, line );
  /* check that at each level, the finger points to the rightmost seg at that level */
  if ( builder->finger_at_end ) {
		for ( level_t level = seglist_min_level; level <= builder->seglist->header->level; level++ ) {
			const Seg *seg = builder->seglist->header;
			const Seg *nxt;
			while( ( nxt = seg->fwd[ level ].seg ) != NIL ) seg = nxt;
			asrt( seg == finger_rightmost_seg( builder->finger, level ) );
		}
  }
}
#endif // #ifndef NDEBUG


#if 0

static Seglist *seglist_from_str( const char *str );

	 Seglist *seglist_from_str( const char *str ) {
  bool_t starts_with_L = ( strlen( str ) >= 2 && str[0] == 'L' && str[1] == ':' );
  if ( starts_with_L ) str += 2;
  Seglist_builder builder;
  seglist_builder_init( &builder );
  char *s = cosi_strdup( str );

  char *saveptr;
  char *token1 = cosi_strtok_r( s, " ", &saveptr );
  while( token1 ) {
		{
			double beg_f, end_f;
			loc_t beg, end;
			char leafset_str[1024];
			int num_scanned = sscanf( token1, LOC_FMT "-" LOC_FMT ":%s", &beg_f, &end_f, leafset_str );
			beg = loc_t( beg_f );
			end = loc_t( end_f );
			chkCond( num_scanned == 3, "seglist_from_str: error scanning range with leafset - %s", token1 );
			seglist_builder_add_seg( &builder, beg, end, leafset_from_str( leafset_str ), NULL_GEN );
		} 
		token1 = cosi_strtok_r( NULL, " ", &saveptr );
  }
  seglist_chk_builder( &builder );
  free( s );
  return seglist_builder_free_and_get_result( &builder );
}


#define seglist_chk_eq(seglist,str) seglist_chk_eq_helper(seglist,str,#seglist,__FILE__,__LINE__)
void seglist_chk_eq_helper( const Seglist *seglist, const char *str, const char * /*expr*/, const char * /*fname*/, int /*line*/ ) {
  Seglist *a_seglist = seglist_from_str( str );

  if ( !seglist_equal( seglist, a_seglist ) ) {
		//printf( "error at %s:%d: %s is %s should be %s", fname, line, expr, seglist_str( seglist ), str );
		assert(0);
  }
  
  seglist_free_seglist( a_seglist );
}

void seglist_test0(void ) {
  set_rng_seed(220);

  Seglist *s1 = seglist_from_str( "0-.2 .2-.4 .4-.6 .6-.8 .8-.9 .9-1" );

  seglist_write_dot( s1, /* finger= */ NULL, "seglist.dot" );

  Seglist *left, *right;
  seglist_split( s1, .5, &left, &right, /* split_seg= */ True );
  seglist_chk( s1 );
  seglist_chk( left );
  seglist_chk( right );

  seglist_chk_eq( left, "0-.2 .2-.4 .4-.5" );
  seglist_chk_eq( right, ".5-.6 .6-.8 .8-.9 .9-1" );

  seglist_write_dot( left, NULL, "left.dot" );
  seglist_write_dot( right, NULL, "right.dot" );

  Seglist *back = seglist_union_or_inters( left, right, /* is_union= */ True, /* has_leafsets= */ False );
  seglist_write_dot( back, NULL, "back.dot" );

  seglist_chk_eq( back, "0-.2 .2-.4 .4-.5 .5-.6 .6-.8 .8-.9 .9-1" ); 

  Seglist_builder builder2;
  seglist_builder_init( &builder2, /* has_leafsets= */ False );
  seglist_builder_add_seg( &builder2, .3, .65, LEAFSET_NULL, NULL_GEN );

  seglist_write_dot( builder2.seglist, builder2.finger, "three.dot" );

  printf( "CCCCCCCCCCCCCCCComputing union\n" );
  Seglist *with_three = seglist_union_or_inters( back, builder2.seglist, /* is_union= */ True, /* has_leafsets= */ False );
  seglist_write_dot( with_three, NULL, "withthree.dot" );

  seglist_chk_eq( seglist_union_or_inters( seglist_from_str( ".1-.6" ), seglist_from_str( ".5-.6" ), /* is_union= */ False,
																					 /* has_leafsets= */ False ), "" );

  printf( "\n\nDone.\n" );
}
#endif

#ifdef COSI_DEV_SEGLIST_DEBUG
void seglist_chk_all( void ) {
  for ( Seglist *seglist = all_seglists; seglist; seglist = seglist->next ) {
		assert( !seglist->next || seglist->next->prev == seglist );
		assert( !seglist->prev || seglist->prev->next == seglist );
		for ( Seg *seg = seglist->header; seg != NIL; seg = seglist_next_seg( seg ) ) {
			assert( seg );
			seg->seglist = NULL;
		}
  }
  for ( Seglist *seglist = all_seglists; seglist; seglist = seglist->next ) {
		for ( Seg *seg = seglist->header; seg != NIL; seg = seglist_next_seg( seg ) ) {
			assert( seg );
			assert( !seg->seglist );
			seg->seglist = seglist;
		}
  }
}
#endif

#ifdef COSI_DEV_NODE_DEBUG
void seglist_set_node( Seglist *seglist, struct node *node ) { /*seglist->node = node;*/ }
struct node *seglist_get_node( const Seglist *seglist ) { return NULL; /*seglist->node;*/ }
#endif


ostream& operator<<( ostream& s, const Seg& seg ) {
	if ( &seg == NIL )
		 s << "(NIL)";
	else
		 s << "(" << seg.beg << "-" << seg.end << ")";
	return s;
}

ostream& operator<<( ostream& s, const Seglist& seglist ) {
	s.precision( 16 );
	s << "[";
	for ( const Seg *seg = seglist_first_seg_const( &seglist ); seg!= NIL; seg = seglist_next_seg_const( seg ) )
		 s << *seg;
	s << "]";
	return s;
}

#if 0
void seglist_test_1(void ) {
  printf("seed: %ld\n", -1 * seed_rng());
  
  /*set_rng_seed(34220);*/

  leafset_set_max_leaf_id( 360 );

  {
		FILE *f = fopen( "tests/segtestintersgc.dat", "rb" );
		assert(f);
		int nchecked = 0;
		while( !feof( f ) ) {
			if ( ( nchecked % 100 ) == 0 ) printf( "nchecked=%d\n", nchecked );
			long file_pos = ftell( f );
			Seglist *sl_input = seglist_read_bin( f );
			if ( !sl_input ) {
				assert( feof( f ) );
				break;
			}

			loc_t loc1, loc2;
			cosi_fread_val( loc1, f );
			cosi_fread_val( loc2, f );

			bool_t result;
			cosi_fread_val( result, f );

			Seglist *inside_chk = NULL, *outside_chk = NULL;
			if ( result ) {
				inside_chk = seglist_read_bin( f );
				outside_chk = seglist_read_bin( f );
		  
				/*s1->has_leafsets = s2->has_leafsets = False;*/
		  
				assert( !sl_input->has_leafsets || ( ( inside_chk->has_leafsets || seglist_is_empty( inside_chk ) ) &&
																						 ( outside_chk->has_leafsets || seglist_is_empty( outside_chk ) ) ) );
		  
		  
				inside_chk->has_leafsets = outside_chk->has_leafsets = sl_input->has_leafsets;
				seglist_canonicalize( &inside_chk );
				seglist_canonicalize( &outside_chk );
		  
				/*printf( "%s union %s gives %s\n", seglist_str( s1 ), seglist_str( s2 ), seglist_str( result_chk ) );*/
				{
			 
					Seglist *inside, *outside;
					seglist_intersect_for_gc( &sl_input, loc1, loc2, &inside, &outside, /* for_sentinels= */ True );
					seglist_canonicalize( &inside );
					seglist_canonicalize( &outside );
			 
					if( !seglist_equal( inside, inside_chk ) || !seglist_equal( outside, outside_chk ) ) {
				
						fseek( f, file_pos, SEEK_SET );
				
						sl_input = seglist_read_bin( f );
				
						// printf( "\n\nsl_input=%s\n", seglist_str( sl_input ) );
						// printf( "\n\n    inside=%s\ninside_chk=%s\n\n", seglist_str( inside ), seglist_str( inside_chk ) ); 
						// printf( "\n\n    outside=%s\noutside_chk=%s\n\n", seglist_str( outside ), seglist_str( outside_chk ) ); 
						chkCond( False, "incorrect result!" );
					}
					nchecked++;
					seglist_free_seglist( inside );
					seglist_free_seglist( outside );
				} /*else {
						seglist_free_seglist( s1 );
						seglist_free_seglist( s2 );
						}*/
		  
				seglist_free_seglist( inside_chk );
				seglist_free_seglist( outside_chk );
		  
			} 
		  
		}
		fclose( f );
		printf( "\n\nNCHECKED=%d\n\n", nchecked );
  }

  const char *UNION = "union";
  const char *INTERS = "inters";
  
  const char *tests[]= {
		"", UNION, "", "",
		".2-.5", UNION, ".6-.7", ".2-.5 .6-.7",
		"0-.6", UNION, ".4-1", "0-1",
		"0-1:1", UNION, "0-1:2", "0-1:1,2",

		".2-.5:1,2", UNION, ".6-.7:3,4", ".2-.5:1,2 .6-.7:3,4",
		"0-.6:1", UNION, ".4-1:2", "0-.4:1 .4-.6:1,2 .6-1:2",
		".1-.4:1 .4-.5:2", UNION, "0-1:3", "0-.1:3 .1-.4:1,3 .4-.5:2,3 .5-1:3",

		"0-1", INTERS, "0-1", "0-1",
		"0-.6", INTERS, ".4-1", ".4-.6",
		NULL
  };

  const char **s = tests;
  int num_checked = 0;
  while( *s ) {
		printf( "s0=%s s1=%s s2=%s s3=%s\n", s[0], s[1], s[2], s[3] );
	 
		const char *op = s[1];
		Seglist *s1 = seglist_from_str( *s );
		Seglist *s2 = seglist_from_str( s[2] );
		Seglist *result_chk = seglist_from_str( s[3] );

		Seglist *result = seglist_union_or_inters( &s1, &s2, /* is_union= */ op == UNION, /* inters= */ NULL );
	 
		if( !seglist_equal( result, result_chk ) ) {
			//printf( "s1=%s s2=%s op=%s result=%s expect=%s\n", seglist_str( s1 ), seglist_str( s2 ), op, seglist_str( result ), seglist_str( result_chk ) );
			assert(0);
		}
		seglist_free_seglist( result );
		seglist_free_seglist( result_chk );
		num_checked++;
		s += 4;
  }

  printf( "\n\nDone: checked %d tests.\n", num_checked );
}
#endif

int seglist_run_tests() {
	std::cout.precision(18);
	std::cerr.precision(18);
	leafset_set_max_leaf_id(9);

	Seglist_builder builder;
	seglist_builder_init( &builder );
	seglist_builder_add_seg( &builder, loc_t( ploc_t( 0.1 ), gloc_t( 0.1 ) ), loc_t( ploc_t( .2 ), gloc_t( .2 ) ), make_singleton_leafset( 1 ), NULL_GEN );
	seglist_builder_add_seg( &builder, loc_t( ploc_t( 0.3 ), gloc_t( 0.3 ) ), loc_t( ploc_t( .4 ), gloc_t( .4 ) ), make_singleton_leafset( 1 ), NULL_GEN );
	
	Seglist *s1 = seglist_builder_free_and_get_result( &builder );
	chkCond( equal_eps( seglist_find_glen( s1, glen_t(.05) ), gloc_t(.15) ), "err" );
	chkCond( equal_eps( seglist_find_glen( s1, glen_t(.15) ), gloc_t(.35), 1e-16 ), "err" );
	seglist_delete( &s1 );
	
	return EXIT_SUCCESS;
}

}  // namespace seglist_detail
}  // namespace seglist
}  // namespace cosi
