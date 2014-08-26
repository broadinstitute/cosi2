/* $Id: pop.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */
//#define COSI_DEV_PRINT
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cosi/pop.h>
#include <cosi/node.h>
#include <cosi/hullmgr.h>
#include <cosi/seglist.h>
#include <cosi/utils.h>

namespace cosi {

using std::string;
using util::chkCond;

Pop::Pop(popid name_, int popsize_, const string& label_) :
	name( name_ ), popsize( popsize_ ), label( label_ ),
	isRestrictingCoalescence( false )
#ifdef COSI_SUPPORT_COALAPX
	, hullMgr( boost::make_shared<HullMgr>() )
#endif	
#ifdef COSI_DEV	
	, coalRate( 0.0 )
#endif
{
	chkCond( popsize >= 0, "creating pop with negative size" );
}

Pop::~Pop() { }

void Pop::pop_remove_node_by_idx ( int idx_in_pop) {
	Node *n = members[ idx_in_pop ];
	PRINT3( "pop_remove_node_by_idx", this->name, idx_in_pop );
	if ( popListener ) {
		popListener->nodeRemoved( this, idx_in_pop );
		if ( idx_in_pop < ((int)members.size())-1 ) {
			Node *last_node = members.back();
			popListener->nodeRemoved( this, last_node->get_idx_in_pop() );
			popListener->nodeAdded( this, last_node->getName(), idx_in_pop, seglist_beg( last_node->getSegs() ),
														 seglist_end( last_node->getSegs() ) );
		}
	}
  nodelist_remove_idx( &(members), idx_in_pop );
	Node::PopAccess::SetNodePop( n, NULL );
#ifdef COSI_SUPPORT_COALAPX
	if  ( isRestrictingCoalescence ) {
		hullMgr->removeHull( Node::PopAccess::GetNodeHullPtr( n ) );
		chkHullMgr();
	}
#endif	
}

void 
Pop::pop_remove_node ( Node *nodeptr) 
{
  assert( nodeptr );
	assert( nodeptr->getPop() == this );
	assert( members[ nodeptr->get_idx_in_pop() ] == nodeptr );
  pop_remove_node_by_idx( nodeptr->get_idx_in_pop() );
}

void 
Pop::pop_add_node ( Node *nodeptr) 
{
	Node::PopAccess::SetNodePop( nodeptr, this );
  nodeptr->set_idx_in_pop( nodelist_add(&(members), nodeptr) );
#ifdef COSI_SUPPORT_COALAPX	
	if ( restrictingCoalescence() ) {
		 hullMgr->addHull( seglist_beg( nodeptr->getSegs() ), seglist_end( nodeptr->getSegs() ),
											 Node::PopAccess::GetNodeHullPtr( nodeptr ) );

		 chkHullMgr();
	}
#endif
	if ( popListener ) popListener->nodeAdded( this, nodeptr->getName(), nodeptr->get_idx_in_pop(),
																						 seglist_beg( nodeptr->getSegs() ),
																						 seglist_end( nodeptr->getSegs() ) );
}

#ifdef COSI_SUPPORT_COALAPX
void Pop::setCoalMargin( len_t margin_ ) {
	hullMgr->setMargin( margin_ );
	this->isRestrictingCoalescence = ( margin_ < ( MAX_LOC - MIN_LOC ) );
}

bool Pop::restrictingCoalescence() const { return this->isRestrictingCoalescence; }

bool apxWithTrajOk = True;
factor_t apxMinFactor = 0;

bool Pop::useCoalApx() const {
	return isRestrictingCoalescence && ( popsize > members.size() * apxMinFactor ) &&
		 ( apxWithTrajOk || !coalArrivalProcess );
}

nchromPairs_t Pop::getNumCoalesceableChromPairs() const {
	return !useCoalApx() ? ( members.size() * ( members.size() - 1 ) / 2 ) :
		 hullMgr->getNumIntersections();
}

std::pair< Node *, Node * > Pop::chooseRandomIntersection( RandGenP randGen ) {
	if ( !useCoalApx() ) {
		nchroms_t node1idx = randGen->random_idx( members.size() );
		nchroms_t node2idx = randGen->random_idx( members.size() - 1 );
		if ( node2idx >= node1idx ) node2idx++;

		return std::make_pair( members[ node1idx ], members[ node2idx ] );

	} else {
		std::pair< const HullMgr::Hull *, const HullMgr::Hull * > p = 
			 hullMgr->chooseRandomIntersection( randGen );
		return std::make_pair( Node::PopAccess::GetNodeFromHullPtr( p.first ),
													 Node::PopAccess::GetNodeFromHullPtr( p.second ) );
													 
	}
}



void Pop::chkHullMgr() const {

#if !defined(NDEBUG) && defined(COSI_SUPPORT_COALAPX)

	if ( !isRestrictingCoalescence ) return;
	
	//typedef HullMgr::ninters_t ninters_t;

	nchromPairs_t npairs = 0;
	nchromPairs_t npairs_tot = 0;
	len_t maxCoalDist = hullMgr->getMargin();

	//PRINT2( maxCoalDist, members.size() );
	for ( size_t ii = 0; ii < members.size(); ii++ ) {
		const Node *node1 = members[ ii ];
		for ( size_t jj = ii+1; jj < members.size(); jj++ ) {
			const Node *node2 = members[ jj ];
			npairs_tot++;
			if (	!( seglist_beg( node1->getSegs() ) - seglist_end( node2->getSegs() ) > maxCoalDist ||
							 seglist_beg( node2->getSegs() ) - seglist_end( node1->getSegs() ) > maxCoalDist ) ) {
				npairs++;
			}
		}
	}
			
	//hullMgr.chkMap( loc2count );

	typedef util::order_statistics_tree<ploc_t> ost_t;
	typedef ost_t::iterator ost_iter_t;
	ost_t begs, ends;
			
	nchroms_t n_beg0 = 0;
	for ( size_t ii = 0; ii < members.size(); ii++ ) {
		const Node *n = members[ ii ];
		ploc_t beg( get_ploc( seglist_beg( n->getSegs() ) ) );
		if ( beg == ploc_t(0.0) )
			 n_beg0++;
		else
			 begs.insert( beg );

		ploc_t end_ext = get_ploc( seglist_end( n->getSegs() ) ) + maxCoalDist;
		if (  end_ext < ploc_t(1.0) )
			 ends.insert( end_ext );
	}
	nchromPairs_t n_pairs2 = n_beg0 * (n_beg0-1) / 2;
	for ( ost_iter_t bi = begs.begin(); bi != begs.end(); bi++ ) {
		ost_iter_t closestEnd = ends.upper_bound( *bi );
		nchroms_t n_end_before = closestEnd.position();
		nchroms_t n_beg_after = begs.size() - bi.position();

		n_pairs2 += ( members.size() - n_end_before - n_beg_after );
	}

//	PRINT4( npairs_tot, npairs, n_pairs2, hullMgr->getNumIntersections() );
	static size_t nchecked = 0;

	nchecked++;
	if (  n_pairs2 != npairs  || hullMgr->getNumIntersections() != n_pairs2 ) {
		PRINT4( npairs_tot, npairs, n_pairs2, hullMgr->getNumIntersections() ); 
//				PRINT2( n_pairs2, lastIntersCount );
		throw std::runtime_error( "inters count mismatch" );
	} else {
		if ( !( nchecked % 10000 ) ) {
			PRINT6( nchecked, npairs_tot, npairs, n_pairs2, hullMgr->getNumIntersections(), hullMgr->getMargin() );
		}
	}

			
	//coalRate = ((double)n_pairs2) / ( 2 * std::max(popsize,1) );
	//PRINT3( npairs, npairs_tot, coalRate );

#endif  // #ifndef NDEBUG
	
}  // Pop::chkHullMgr()

#endif // #ifdef COSI_SUPPORT_COALAPX


}
