#include <boost/foreach.hpp>
#include <cosi/general/mempool.h>

namespace cosi {

using util::chk;	

#define ForEach BOOST_FOREACH  
  
Mempool::Mempool():
	item_size( 0 ), block_size( 0 ),
	cur_block( NULL ), cur_block_end( NULL ),
	next_item( NULL ), num_live_items( 0 ),
	max_live_items( 0 ), num_reused( 0 ) {
}

Mempool::Mempool( size_t item_size_, unsigned block_size_in_items_, unsigned initial_freelist_size_in_items_ ) {
	mempool_init( item_size_, block_size_in_items_, initial_freelist_size_in_items_ );
}

Mempool::~Mempool() {
#ifdef COSI_MEMPOOL_STATS	
	PRINT6( item_size, num_live_items, max_live_items, num_reused, blocks.size(), free_items.size() );
#endif	
	ForEach( void *p, blocks )
	  free( p );
}

void Mempool::mempool_init( size_t item_size_, unsigned block_size_in_items_, unsigned initial_freelist_size_in_items_ ) {
	item_size = item_size_;
#ifndef COSI_MEMPOOL_DISABLE  
  block_size = block_size_in_items_ * item_size_;
	blocks.reserve( 4 );
	free_items.reserve( initial_freelist_size_in_items_ );
  add_block();
  num_live_items = 0;
  max_live_items = 0;
	num_reused = 0;
#ifdef COSI_MEMPOOL_STATS	
	PRINT4( "mempool_init", item_size, block_size, blocks.size() );
#endif	
#endif  
}

void Mempool::add_block() {
#ifdef COSI_MEMPOOL_STATS	
	PRINT7( "add_block", item_size, block_size, blocks.size(), num_live_items, max_live_items, num_reused );
#endif	
  cur_block = (char *)malloc( block_size );
  chk( cur_block, "error allocating mempool block" );
  cur_block_end = cur_block + block_size;
  next_item = cur_block;
  blocks.push_back( cur_block );
}


void *Mempool::mempool_alloc_item() {
#ifdef COSI_MEMPOOL_DISABLE
  return malloc( item_size );
#else
#ifdef COSI_MEMPOOL_STATS	
  num_live_items++;
	if ( num_live_items > max_live_items )
		 max_live_items = num_live_items;
#endif	
  if ( !free_items.empty() ) {
#ifdef COSI_MEMPOOL_STATS		
		num_reused++;
#endif		
		void *lastElem = free_items.back();
		free_items.pop_back();
		return lastElem;
  }
  if ( next_item == cur_block_end )
		 add_block();
  void *result = next_item;
  next_item += item_size;
  return result;
#endif  
}

void Mempool::mempool_free_item( void *obj ) {
#ifdef COSI_MEMPOOL_DISABLE
  free( obj );
#else
#ifdef COSI_MEMPOOL_STATS	
  --num_live_items;
#endif	
  free_items.push_back( obj );
#endif  
}

}  // namespace cosi

