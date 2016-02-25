/*
 * Header: seglist.h
 *
 * Defines Seglists -- lists of segments of the <simulated region>.
 * The main class defined here is <Seglist>.
 * Example of a seglist is: ([.2-.3],[.5-.8]).  With each segment
 * of a seglist, we keep the leafset -- the set of leaves of the ARG
 * (present-day samples) inheriting that segment.
 *
 * Seglists are used in representing <Node::segs> and
 * <Node::sentinel_for>.
 *
 * Implementation notes:
 *
 * Seglists are represented using skip lists.  Some references on skiplists include:
 *
 * "Skip Lists: A Probabilistic Alternative to Balanced Trees" by William Pugh
 *  ftp://ftp.cs.umd.edu/pub/skipLists/skiplists.pdf
 *
 * "A Skip List Cookbook" by William Pugh
 * http://tinyurl.com/6jhmram
 */

#ifndef __INCLUDE_COSI_SEGLIST_H
#define __INCLUDE_COSI_SEGLIST_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <iterator>
#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <cosi/general/utils.h>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>
#include <cosi/segsumm.h>
#include <cosi/seglistfwd.h>

namespace cosi {
namespace seglist2 {
namespace seglist_detail {

#if 0
template <class Impl>
class SeglistPool {
public:
	 SeglistPool();
	 ~SeglistPool();

	 typedef typename Impl::Seglist Seglist;
	 typedef typename Impl::loc_t loc_t;

	 Seglist makeSeglist();
	 Seglist joinSeglists( Seglist s1, Seglist s2 );

};

/////

template <class Impl>
 inline SeglistPool<Impl>::SeglistPool() { }

 template <class Impl>
 inline SeglistPool<Impl>::~SeglistPool() { }

template <class Impl> inline
typename SeglistPool<Impl>::Seglist SeglistPool<Impl>::makeSeglist() { return static_cast<Impl*>(this)->makeSeglist(); }

template <class Impl> inline
typename SeglistPool<Impl>::Seglist SeglistPool<Impl>::joinSeglists( Seglist s1, Seglist s2 ) { return static_cast<Impl*>(this)->joinSeglists( s1, s2 ); }

//////

class SeglistPoolImpl: public SeglistPool<SeglistPoolImpl> {
	 typedef int Seglist;
	 typedef double loc_t;
	 
	 Seglist makeSeglist();
	 Seglist joinSeglists( Seglist s1, Seglist s2 );
};

inline
SeglistPoolImpl::Seglist SeglistPoolImpl::makeSeglist() { return Seglist(0); }

inline
SeglistPoolImpl::Seglist SeglistPoolImpl::joinSeglists( Seglist s1, Seglist s2 ) { return s1 + s2; }

/////
#endif  // comment out experimental code

struct seglist;
typedef struct seglist Seglist;

struct seg;
typedef struct seg Seg;

struct segptr;
typedef struct segptr Segptr;


#ifdef COSI_TRACK_LAST_COAL
#define IF_COSI_TRACK_LAST_COAL(...) __VA_ARGS__
#else
#define IF_COSI_TRACK_LAST_COAL(...)
#endif

/* Type: level_t */
/* Represents the skiplist level of a <Seg>: the height of its tower.  */
typedef int level_t;

const level_t seglist_min_level = 1;
const level_t  seglist_max_level = 15;

const level_t NULL_LEVEL = -1;

/*
 * Struct: Segptr
 *
 * Forward pointer from one <Seg> to the seg Seg at that or higher level.
 * A Seg of level k has k such forward pointers.
 */
struct segptr {
	 /* Field: seg */
	 /* Pointer to the next seg in this seglist with the same or higher level as this seg, */
	 /* or to <NIL> if there are no more such segs in this seglist. */
	 Seg *seg;
  
	 /* Field: skip_len */
	 /* Total length of segs skipped over by the <Segptr::seg> pointer above. */
	 /* Note that this counts only the length of the segs, not the gaps between them. */
	 /* Does not include the lengths of this seg or the target seg of the link. */
	 /* So, if the seg link above points to the immediate successor in the seglist, */
	 /* skip_len is zero. */
	 len_t skip_len;
};

/*
 * Struct: Seg
 *
 * One segment of a <Seglist>.  
 */
struct seg: private boost::noncopyable {
#ifdef COSI_DEV_SEGLIST_DEBUG
	 /* Field: seglist */
	 /* The seglist to which this seg belongs. */
	 Seglist *seglist;

	 /* Field: id */
	 /* Unique id of this seg -- its order in the creation of all objs. */
	 obj_id_t id;
#endif  
  
	 /* Fields: Endpoints */
	 /*    beg - the start of the seg on the chromosome */
	 /*    end - the end of the seg on the chromosome */
	 loc_t beg, end;

	 /* Field: leafset */
	 /* The set of leaves that will inherit any mutation placed on this seg. */
	 leafset_p leafset;

	 loc_t getBeg() const { return beg; }
	 loc_t getEnd() const { return end; }
	 leafset_p getLeafset() const { return leafset; }
  
	 /* Field: level */
	 /* The height of this seg's tower, in the skiplist representation of its <Seglist>. */
	 /* If this Seg is a header of its Seglist, then the height of the tower is seglist_max_level */
	 /* and this field stores the level of the seglist (>= max level of any Seg in it, plus one). */
	 level_t level;

#ifdef COSI_TRACK_LAST_COAL
	 // Field: lastCoalGen;
	 // The last generation at which this chunk coalesced.
	 genid lastCoalGen;
#endif	 
  
	 /* Field: fwd */
	 /* This seg's tower of forward pointers to the elements later in this <Seglist>. */
	 /* The size of this array is <level>, except if it's the header of the list; in that case, */
	 /* the size is <seglist_max_level>. */
	 Segptr fwd[];

	 void* operator new (size_t size, level_t);
	 void operator delete (void *p);

};  // struct Seg

ostream& operator<<( ostream&, const Seg& );

inline Seg *seglist_next_seg( Seg *a_seg ) { return ((a_seg)->fwd[ seglist_min_level ].seg); }
inline const Seg *seglist_next_seg_const( const Seg *a_seg ) { return ((a_seg)->fwd[ seglist_min_level ].seg); }

/*
 * Define mutator routines for setting fields of structures.
 * Ordinarily, these translate to simple assignment statements.
 * But for debugging purposes, we can use these to track the provenance of values,
 * and catch the moment when a particular field of a particular object gets modified.
 */

// inline void seglist_seg_set_beg(Seg *a_seg, loc_t new_beg) { (a_seg)->beg = (new_beg); } 
// inline void seglist_seg_set_end(Seg *a_seg, loc_t new_end) { (a_seg)->end = (new_end); }
inline void seglist_seg_set_leafset(Seg *a_seg, leafset_p new_leafset) { (a_seg)->leafset = (new_leafset); }
//inline void seglist_seg_set_level(Seg *a_seg, level_t new_level) { (a_seg)->level = (new_level); }
//inline void seglist_seg_dec_level(Seg *a_seg) { (((a_seg)->level)--); }
// inline void seglist_seg_set_fwd(Seg *a_seg, level_t level, Seg *new_fwd) { (a_seg)->fwd[ (level) ] = new_fwd; }
// inline void seglist_seg_set_fwd_seg(Seg *a_seg, level_t level, new_fwd_seg) (a_seg)->fwd[ (level) ].seg = new_fwd_seg
// #define seglist_seg_set_fwd_skip_len(a_seg, level, new_fwd_skip_len) (a_seg)->fwd[ (level) ].skip_len = new_fwd_skip_len

//#define seglist_seg_cpy_fwd( trg_seg, trg_level, src_seg, src_level, n_levels ) memcpy( &( (trg_seg)->fwd[ (trg_level) ] ), &( (src_seg)->fwd[ (src_level) ] ), ( n_levels ) )

class Seglist_iterator;
class Seglist_const_iterator;

/*
 * Struct: Seglist
 *
 * The list of segments representing portions of the chromosome inherited by at least some present-day chromosomes.
 * Semantically, a seglist looks like this: [ [.2,.3], [.45,.66], [.7,.93] ].  That is, a seglist is a list of segments
 * (segs for short), each of which is defined by a pair of <locations>.  The segs in a seglist a sorted left to right
 * along the chromosome and may touch, but never overlap.
 *
 * Represented as a single-linked skip list, augmented as described below.
 * Some references on skiplists include:
 *
 * "Skip Lists: A Probabilistic Alternative to Balanced Trees" by William Pugh
 * ftp://ftp.cs.umd.edu/pub/skipLists/skiplists.pdf
 *
 * "A Skip List Cookbook" by William Pugh
 * http://tinyurl.com/6jhmram
 *
 * In brief, a skip list is an ordinary linked list, in which each element points to the next; but also, roughly
 * half the elements also store a pointer to the element-after-next; of these, roughly half also store a pointer
 * to the fourth element down, and so on.  By following these forward pointers, one can get to a given point
 * on the list in logarithmic time.  Skip lists can also be cut and merged in logarithmic time, which is important
 * for our application.
 *
 * The skip lists used here are modified somewhat compared to the ones described in the literature.
 * For one, we have lists of segments rather than lists of points.  
 * Additionally, our skip lists are augmented as follows: each forward pointer stores the sum of lengths of the segs it skips over.
 *
 * There are actually two types of seglists in use: with leafsets, and without.  The former are used for <Node::segs>, the latter for
 * <Node::sentinel_for>.  In the former, for each <Seg> we keep, besides its beginning and end, its <leafset>: the set of leaves
 * (present-day sample nodes) that inherit any mutation placed on that Seg.  We also keep, for each forward link the total length
 * of segments it skips over (for speeding up <mutate_put_muts_on_seglist()>).  For seglists stored in <Node::sentinel_for>, we
 * keep neither of those things, although we still allocate space for them for simplicity.
 *
 */
struct seglist {
	 /* Field: header */
	 /* The header element of the list.  Its level one higher than the max level of all segs in the list, and its (beg,end) are -infinity.  */
	 /* It is part of every list, even the empty one; and its own (beg,end) segment is _not_ part of the list.  */
	 Seg *header;

	 /* Field: _end */
	 /* The end of the rightmost Seg of this list; stored here so we can check when a recombination loc is */
	 /* completely to the right of the list. */
	 loc_t _end;

	 /* Field: segsumm */
	 /* A coarse over-approximation of this seglist, indicating which <seglets> (divisions of the [0,1] interval) it may touch. */
	 /* Used for fast intersection checks.  */
	 /* While segsumm should always be a correct over-approximation of the seglist when viewed from outside the seglist module, */
	 /* this invariant is sometimes broken for efficiency's sake during internal seglist operations. */
	 segsumm_t segsumm;

	 ///////////
	 
	 typedef Seglist_const_iterator const_iterator;
	 typedef Seglist_iterator iterator;

	 const_iterator begin() const;
	 const_iterator end() const;
	 iterator begin();
	 iterator end();

	 loc_t getEnd() const { return _end; }
	 void setEnd( loc_t end_ ) { this->_end = end_; }
};  // class Seglist


#if 0 

template <typename LOC_T>
class SeglistPool: private boost::noncopyable {
public:
	 SeglistPool();

	 typedef ... SeglistP;

	 SeglistP makeSeglistFull();
	 SeglistP makeSeglistEmpty();

	 SeglistP join( SeglistP seglist1, SeglistP seglist2 );
	 pair<SeglistP,SeglistP> split( SeglistP seglist );
	 
};

#endif

void seglist_init_module(void);
void seglist_chk( const Seglist *seglist );

#ifndef COSI_DEV_SEGLIST_DEBUG
#define seglist_chk_all()
#else
void seglist_chk_all( void );
#endif

inline Seg *seglist_header( const Seglist *s ) { return s->header; }

Seglist *seglist_make_empty();
Seglist *seglist_make_full( leafset_p leafset, genid lastCoalGen );

bool_t seglist_equal( const Seglist *seglist1, const Seglist *seglist );

/*
 * Type: Finger_elt
 *
 * Represents one level of a <Finger> (see below).
 */
typedef Segptr Finger_elt;

/*
 * Type: Finger
 *
 * An iterator that keeps our place in a <Seglist> as we iterate through it.
 * More specifically, for each level of a seglist it keeps a pointer to
 * the rightmost segment at that level that is still to the left of a given <loc>.
 * It may also keep, for each level, the total length we have traveled along that level
 * while searching for the rightmost Seg that's still to the left of a loc.
 * Because this is the same information kept in a <Segptr>, we reuse that physical type.
 */
typedef Finger_elt Finger;


Finger *seglist_split( Seglist **seglistp, loc_t loc, Seglist **s_left, Seglist **s_right, bool_t split_seg,
											 bool_t *did_split = NULL );

// FuncP: seglist_find_seg
// Return the seg within the seglist that contains the given location,
// or NULL if no such seg.
const Seg *seglist_find_seg( const Seglist *seglist, loc_t loc );


// Func: seglist_contains_loc
// Test whether a Seglist contains a given location in one of its segs.
inline bool_t seglist_contains_loc( const Seglist *seglist, loc_t loc ) { return seglist_find_seg( seglist, loc ); }

void seglist_delete( Seglist **seglistp );

typedef boost::function<void (loc_t beg, loc_t end, genid genMoreRecent, genid genLessRecent, leafset_p leafset)> seglist_seg_union_callback_t;

Seglist *seglist_union( Seglist **seglist1p, Seglist **seglistp, genid gen = NULL_GEN, seglist_seg_union_callback_t seg_union_callback = 0 );

Seglist *seglist_concat( Seglist **seglist1, Seglist **seglist );

Seglist *seglist_dup( const Seglist *seglist );
len_t seglist_tot_len( const Seglist *seglist );
glen_t seglist_tot_glen( const Seglist *seglist );

bool_t seglist_intersect_for_gc( Seglist **seglist, loc_t loc1, loc_t loc2, Seglist **inside, Seglist **outside,
																 bool_t for_sentinels );

void seglist_equal_dbg( const Seglist *seglist1, const Seglist *seglist );

extern Seg *NIL;

inline Seg *seglist_first_seg( const Seglist *seglist ) { return seglist->header->fwd[ seglist_min_level ].seg; }
inline const Seg *seglist_first_seg_const( const Seglist *seglist ) { return seglist->header->fwd[ seglist_min_level ].seg; }
inline const loc_t seglist_beg( const Seglist *seglist ) { return seglist_first_seg_const( seglist )->beg; }
inline bool_t seglist_is_empty( const Seglist *seglist ) { return seglist_first_seg_const( seglist ) == NIL; }

inline loc_t& seglist_end(Seglist *s) { return ((s)->_end); }

Finger *seglist_finger_alloc(void);
void seglist_init_finger( const Seglist *seglist, Finger *finger, bool_t init_lens );
void seglist_finger_free( Finger *finger );
bool_t seglist_advance_by_len( const Seglist *seglist, Finger *finger,
															 loc_t *cur_loc, len_t *cur_len,  leafset_p *cur_leafset, len_t *advance_by );

//
// Func: seglist_find_glen
//
// Find a point on the given seglist corresponding to traveling a given glen along
// the seglist's segs.
//
// Params:
//
//    seglist - the seglist
//    glen - how far to travel along the seglist's segs; must be between 0 and
//       the total glen of the seglist.
//    genMap - the genetic map of the region
//
gloc_t seglist_find_glen( const Seglist *seglist, glen_t glen ); 

/*
 * Func: seglist_move_beg
 *
 * Move the beginning of the first seg of a non-empty seglist either forward or backwards.
 * If moving forward, must not move past the end of the first segment.
 */
void seglist_move_beg( Seglist *seglist, loc_t new_beg );

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
void seglist_move_end( Seglist *seglist, Finger *finger, loc_t new_end );


void seglist_segsumm_chk_helper( const Seglist *my_s, const char *fname, int line );
#define seglist_segsumm_chk(seglist)

ostream& operator<<( ostream& s, const Seglist& seglist );

// Class: Seglist_const_iterator
// Iterates over <Segs> of a <Seglist>
class Seglist_const_iterator: public std::iterator<std::forward_iterator_tag, const Seg, /* Distance= */ void, const Seg*, const Seg&>  {
	 typedef Seglist_const_iterator this_type;
public:

	 explicit Seglist_const_iterator( const Seglist *seglist ):
		 _seg( seglist_first_seg_const( seglist ) ) { }
	 Seglist_const_iterator(): _seg( NULL ) {}
	 explicit Seglist_const_iterator( const Seg *seg_ ):
		 _seg( seg_ ) { }

   Seglist_const_iterator& operator++()
   { _increment();   return *this;   }
  
   Seglist_const_iterator operator++(int)
   {
      this_type result (*this);
      _increment();
      return result;
   }

   const Seg& operator*() const
   { return _dereference(); }

   const Seg* operator->() const
   { return &(_dereference()); }

   friend bool operator== (const this_type& i, const this_type& i2)
   { return i._equal(i2); }

   friend bool operator!= (const this_type& i, const this_type& i2)
   { return !(i == i2); }
	 

private:
	 const Seg *_seg;

	 const Seg& _dereference() const { assert( _seg && _seg != NIL ); return *_seg; }
	 void _increment() { assert( _seg && _seg != NIL ); _seg = seglist_next_seg_const( _seg ); }
	 bool _equal( const this_type& other ) const { return this->_seg == other._seg; }
};  // class Seglist_const_iterator

inline Seglist::const_iterator Seglist::begin() const { return Seglist_const_iterator( this ); }
inline Seglist::const_iterator Seglist::end() const { return Seglist_const_iterator( NIL ); }


// Class: Seglist_iterator
// Iterates over <Segs> of a <Seglist>
class Seglist_iterator: public std::iterator<std::forward_iterator_tag, Seg>  {
	 typedef Seglist_iterator this_type;
public:

	 explicit Seglist_iterator( Seglist *seglist ):
		 _seg( seglist_first_seg( seglist ) ) { }
	 Seglist_iterator(): _seg( NULL ) {}
	 explicit Seglist_iterator( Seg *seg_ ):
		 _seg( seg_ ) { }

   Seglist_iterator& operator++()
   { _increment();   return *this;   }
  
   Seglist_iterator operator++(int)
   {
      this_type result (*this);
      _increment();
      return result;
   }

   Seg& operator*() const
   { return _dereference(); }

   Seg* operator->() const
   { return &(_dereference()); }

   friend bool operator== (const this_type& i, const this_type& i2)
   { return i._equal(i2); }

   friend bool operator!= (const this_type& i, const this_type& i2)
   { return !(i == i2); }
	 

private:
	 Seg *_seg;

	 Seg& _dereference() const { assert( _seg && _seg != NIL ); return *_seg; }
	 void _increment() { assert( _seg && _seg != NIL ); _seg = seglist_next_seg( _seg ); }
	 bool _equal( const this_type& other ) const { return this->_seg == other._seg; }
};  // class Seglist_iterator

inline Seglist::iterator Seglist::begin() { return Seglist_iterator( this ); }
inline Seglist::iterator Seglist::end() { return Seglist_iterator( NIL ); }


#ifdef NDEBUG
#define seglist_cmp(s2,s,FN,LN)
#else
#define seglist_cmp(s2,s,FN,LN) seglist_cmp_helper(s2,s,#s2,#s,__FILE__,__LINE__,FN,LN)
#endif

ostream& operator<<( ostream& buf, const Seg *seg );
ostream& operator<<( ostream& buf, const Seglist *seglist );

void seglist_chk_helper( const Seglist *seglist, const char *fname, int line );

/* Func: seglist_chk */
/* Check the validity of the representation of the given seglist. */
#ifdef NDEBUG
#define seglist_chk(s)
#define seglist_chk_builder(b)
#define asrt(c)
#else
#define seglist_chk(s) seglist_chk_helper(s,__FILE__,__LINE__)
#define asrt(c) if(!(c)){ printf( "error at %s:%d - %s\n", fname, line, #c ); assert(0); }
#define seglist_chk_builder(b) seglist_chk_builder_helper(b, __FILE__, __LINE__ )
#endif

void seglist_write_dot_helper( const Seglist *seglist, const Finger *finger, const char *fname, const char *fn, int line );
#define seglist_write_dot(seglist,finger,fname) seglist_write_dot_helper(seglist,finger,fname,__FILE__,__LINE__)

void seglist_canonicalize( Seglist **seglist );

int seglist_run_tests();
}  // namespace seglist_detail

using seglist_detail::Seglist;
using seglist_detail::Seg;

using seglist_detail::seglist_init_module;

using seglist_detail::seglist_make_full;
using seglist_detail::seglist_split;
using seglist_detail::seglist_union;
using seglist_detail::seglist_delete;
using seglist_detail::seglist_chk;

using seglist_detail::Finger;
using seglist_detail::seglist_advance_by_len;
using seglist_detail::seglist_init_finger;
using seglist_detail::seglist_finger_alloc;

using seglist_detail::seglist_run_tests;

} // namespace seglist

}  // namespace cosi
  
#endif  /* __INCLUDE_COSI_SEGLIST_H */

