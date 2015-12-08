/* $Id: recomb.c,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

/*
	File: recomb.c

	Maintains the genetic map (the recombination rate at each genomic position), and provides a method to choose
	a recombination point according to this map.
*/


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/next_prior.hpp>
#include <boost/foreach.hpp>
#include <cosi/utils.h>
#include <cosi/node.h>
#include <cosi/seglist.h>
#include <cosi/pop.h>
#include <cosi/recomb.h>
#include <cosi/genmap.h>
#include <cosi/demography.h>

namespace cosi {

Recomb::Recomb( DemographyP demography_, GenMapP genMap_ ):
	demography( demography_ ), genMap( genMap_ ), nrecombs( 0 ),
	ignoreRecombsInPop( NULL_POPID ) {
}

Recomb::~Recomb() {
}

// MethodP: getAllNodesRecombRate
// Returns the probability, per generation, of a recombination that splits one of the currently active nodes.
glen_t Recomb::getAllNodesRecombRate() const {
	return demography->getNodePool()->getAllNodesRecombRate();
}

void Recomb::recomb_execute (genid gen, int popindex, loc_t *location, Node**nodes_out) {
	Pop *popptr = demography->dg_get_pop_by_index( popindex );
  nchroms_t nodeindex( (int) (random_double() * ToInt(popptr->pop_get_num_nodes())) );
  Node *n = popptr->pop_get_node (nodeindex);
	*location = genMap->getLoc( gloc_t( random_double() ) );
	demography->dg_recombine( n, gen, *location, nodes_out );
	nrecombs++;
}

void Recomb::recomb_execute (genid gen, frac_t frac) {
	loc_t loc;
	Node *n = demography->getNodePool()->findRecomb( frac, &loc );
	if ( this->ignoreRecombsInPop == NULL_POPID || n->getPop()->pop_get_name() != this->ignoreRecombsInPop ) {
		demography->dg_recombine( n, gen, loc, (Node **)NULL );
		nrecombs++;
	}
}

// End class impl: Recomb

// Class impl: RecombRecorder

RecombRecorder::RecombRecorder() { }
RecombRecorder::~RecombRecorder() {}

// Virtual method: handle_recomb
// Called after each recombination.
void RecombRecorder::handle_recomb( Node *, Node *, loc_t loc, genid) {
	recombLocs.push_back( loc );
}

// End class impl: RecombRecorder

}  // namespace cosi

