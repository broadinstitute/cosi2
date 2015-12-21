/*
 * Header: mempool.h
 *
 * Code for maintaining a memory pool of fixed-size blocks that are never deallocated.
 */

#ifndef __INCLUDE_COSI_MEMPOOL_H
#define __INCLUDE_COSI_MEMPOOL_H

#include <cstdlib>
#include <vector>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <cosi/general/utils.h>

namespace cosi {

/*
 * Struct: Mempool
 *
 * A memory pool for allocating many items of a particular size.
 */
class Mempool: boost::noncopyable {
public:
	 Mempool();
	 Mempool( size_t item_size_, unsigned block_size_in_items_,
						unsigned initial_freelist_size_in_items_ );
	 ~Mempool();

	 void mempool_init( size_t item_size_, unsigned block_size_in_items_,
											unsigned initial_freelist_size_in_items_ );
	 
	 void *mempool_alloc_item();
	 void mempool_free_item( void *obj );
	 
private:
	 
	 /* Field: item_size */
	 /* Size of one item. */
	 size_t item_size;
	 
	 /* Field: block_size */
	 /* Size of a block of memory that we allocate at once. */
	 size_t block_size;
	 
	 /* Field: blocks */
	 /* Pointers to the memory blocks that we have allocated so far. */
	 vector<void *> blocks;
	 
	 /* Field: cur_block */
	 /* Pointer to the current memory block, from which we're allocating items. */
	 char *cur_block;
	 
	 /* Field: cur_block_end */
	 /* The end of the current block. */
	 const char *cur_block_end;
	 
	 /* Field: next_item */
	 /* The next item to be returned */
	 char *next_item;
	 
	 /* Field: free_items */
	 /* A freelist, containing pointers to items that have been freed. */
	 vector<void *> free_items;

	 /* Field: num_live_items */
	 /* Current number of live items in the pool */
	 long num_live_items;
	 
	 /* Field: max_live_items */
	 /* Max number of live items in the pool */
	 long max_live_items;

	 // Field: num_reused
	 // Number of times we were able to return a freelist object instead of
	 // allocating a new block
	 long num_reused;


////////////////////

	 void add_block();
	 
};  // class Mempool


} // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_MEMPOOL_H


