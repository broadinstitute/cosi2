//
// Header: mutlist.h
//
// Represent mutation information for individual mutations (<Mut>), and
// for all mutations in a simulation (<Mutlist>).
//

/* $Id: mutlist.h,v 1.2 2011/05/03 18:50:54 sfs Exp $ */
#ifndef __INCLUDE_COSI_MUTLIST_H
#define __INCLUDE_COSI_MUTLIST_H

#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/operators.hpp>
#include <boost/container/flat_map.hpp>
#ifndef COSI_NO_CPU_TIMER
#include <boost/timer/timer.hpp>
#endif
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/leafset.h>
#include <cosi/mutate.h>
#include <cosi/cosirand.h>

namespace cosi {

using std::vector;
using boost::container::flat_map;
using util::chkCond;

class TreeStatsHook;
typedef boost::shared_ptr<TreeStatsHook> TreeStatsHookP;

// Struct: Mut
// Represents one mutation.
struct Mut: public boost::totally_ordered<Mut> {
	 
	 // Field: loc
	 // The physical location of the mutation on the chromosome.
	 loc_t loc;

	 // Field: leaves
	 // The leaves that inherit the (derived allele of the) mutation.
	 leafset_p leaves;

	 // Field: mutIdOrig
	 // The number of this mutation, numbered from left to right of simulated region.
	 int mutIdOrig;

	 // Field: mutId
	 // The id of this mutation, with any filtered out mutations removed.
	 int mutId;

	 // Field: gen
	 // The generation in which this mutation was born.
	 genid gen;

	 // Field: popName
	 // The population id (given in config file) of the population in which this mutation was born
	 popid popName;

	 Mut();
	 Mut( loc_t loc_, leafset_p leaves_, genid gen_, popid popName_ );

	 // Method: get_int_loc
	 // Return the location of this mut rounded to the nearest integer basepair.
	 int get_int_loc( len_bp_t length  ) const { return (int)( get_loc( loc ) * length ); }

	 bool_t isMonomorphicIn( leafset_p core ) const {
		 leafset_p inters = leafset_intersection( core, leaves );
		 return leafset_is_empty( inters ) || leafset_equal( inters, core );
	 }

};  // struct Mut

inline bool operator<( const Mut& m1, const Mut& m2 ) { return m1.loc < m2.loc; }
inline bool operator==( const Mut& m1, const Mut& m2 ) { return m1.loc == m2.loc; }

inline ostream& operator<<( ostream& s, const Mut& m ) {
  s << "[Mut: " << m.loc << "]";
  return s;
}

//
// Struct: Mutlist
//
// A list of all mutations generated during a simulation,
// and their locations.
// A singleton instance of this struct keeps a global list
// of all mutations we place.  The mutations are NOT sorted
// by position.
//
// This class is used as follows: first <addMut()> is called
// to add mutations one-by-one as they are generated during
// the simulation; then, after all muts have been added, <freeze()>
// is called to do various post-processing (sort the muts by order,
// create mapping from each leaf to the muts on it, etc).
//
class Mutlist {
public:

	 // Constructor: Mutlist
	 // Create an empty mutlist.
	 Mutlist(): frozen( False ) { }

	 // Method: addMut
	 // Add a mutation to the mutlist.  <freeze()> must not have been called yet.
	 void addMut( loc_t loc, leafset_p leaves, genid gen, popid popName );
	 
	 // Method: freeze
	 // Called after the simulation finishes and all muts have been placed.
	 void freeze( bool_t inf_sites, len_bp_t length );

	 // Method: size
	 // Return the number of muts in the mutlist. 
	 size_t size() const { return muts.size(); }

	 typedef vector<Mut>::const_iterator const_iterator;

	 // Method: getLeafMuts
	 // Return list of muts falling on a given leaf (present-day chromosome).
	 const vector< const_iterator >& getLeafMuts( leaf_id_t leaf ) const {
		 assert( frozen );
		 return leaf2muts[ leaf ];
	 }

	 // Method: getMuts
	 // Return the vector of muts.   Note that this is unordered until <freeze()> has been called.
	 const vector<Mut>& getMuts() const { return muts; }

	 // Method: print_haps_ms
	 //
	 // Write the output of one simulation in ms format to the specified stream.
	 // Note that the standard format also includes one header describing all the simulation replicas:
	 // the command line of the simulator, the random seeds, the total number of replicas,
	 // and the number of samples in each replica.  Those are written in <CoSiMain::cosi_main()>;
	 // print_haps_ms only writes a single simulation (replica).
	 //
	 // Params:
	 //
	 //   strm - the stream to which the ms-format simulation is written
	 //   sampleSizes - for each sampled population, the number of sampled haplotypes in it
	 //   treeStatsHook - the hook object used to gather requested tree statistics; may be NULL
	 //   outputMutGens - whether to write out the time of each mutation
	 //   recombLocs - if non-NULL, vector (unsorted) of recombination locations
	 //   outputPrecision - number of decimal places in the output
	 void print_haps_ms( ostream& strm, const vector< nchroms_t >& sampleSizes, TreeStatsHookP treeStatsHook,
											 bool_t outputMutGens, const vector< loc_t > *recombLocs,
											 bool_t outputMutGlocs,
											 GenMapP genMap,
											 int outputPrecision,
#ifndef COSI_NO_CPU_TIMER											 
											 boost::timer::cpu_timer
#else											 
											 void
#endif											 
											 *cpuTimer,
											 const genid *endGen,
											 const std::vector< leaf_id_t > *leafOrder = NULL
		 ) const;

	 // Method: loadFromMs
	 // Load a Mutlist from the output of the ms simulator
	 static MutlistP loadFromMs( istream& );

	 // Method: forMutPairs
	 //
	 // Iterate over all pairs of muts meeting specified conditions.
	 // Useful for gathering LD stats.
	 //
	 // Params:
	 //
	 //    maxSep - max distance between muts.  We will consider all mut pairs separated
	 //     by at most this distance.
	 //
	 //    binFn - the binary function to be called for all mut pairs within 'maxSep'.
	 //
	 template <class BinFn>
	 void forMutPairs( plen_t maxSep, BinFn binFn );
	 
private:
	 // Field: mutVec
	 // The vector of mutations added during the simulation.  Unordered until <freeze()> is called.
	 vector<Mut> muts;

	 // Field: frozen
	 // When True, no more muts can be added.
	 bool_t frozen;

	 // Field: haps
	 // For each leaf, the list of muts placed on it
	 vector< vector< const_iterator > > leaf2muts;

	 typedef vector<Mut>::iterator iterator;
	 
};  // class Mutlist

//
// Class: MutProcessor_AddToMutlist
//
// Implementation of <MutProcessor> which adds the mutation to a <Mutlist>
// This is the default implementation of MutProcessor used by cosi.
//
class MutProcessor_AddToMutlist: public MutProcessor {
public:
	 MutProcessor_AddToMutlist( MutlistP mutlist_ ):
		 mutlist( mutlist_ ) {}
	 virtual ~MutProcessor_AddToMutlist();

	 // Virtual method: processMut
	 // Adds the mutation to the mutlist.
	 virtual void processMut(loc_t loc, leafset_p leaves, genid gen, popid popName);

	 // Virtual method: postprocess
	 // Sorts mutations in order.
	 virtual void postprocess();

private:
	 // Field: mutlist
	 // The <Mutlist> to which we will be adding mutations sent to us.
	 MutlistP mutlist;
};
typedef boost::shared_ptr<MutProcessor_AddToMutlist> MutProcessor_AddToMutlistP;

class MutProcessor_AddToMutlist_WithAscertainment: public MutProcessor_AddToMutlist, private HasRandGen {
	 typedef MutProcessor_AddToMutlist PARENT;
public:
	 MutProcessor_AddToMutlist_WithAscertainment( MutlistP mutlist_, frac_t dropSingletonsFrac_,
																								RandGenP randGen_ ):
		 PARENT( mutlist_ ), HasRandGen( randGen_ ), dropSingletonsFrac( dropSingletonsFrac_ ) { }
	 virtual ~MutProcessor_AddToMutlist_WithAscertainment();

	 // Virtual method: processMut
	 // Adds the mutation to the mutlist.
	 virtual void processMut(loc_t loc, leafset_p leaves, genid gen, popid popName);

private:
	 frac_t dropSingletonsFrac;
	 
};  // class MutProcessor_AddToMutlist_WithAscertainment

//
// Implementations
//


// Method: forMutPairs
//
// Iterate over all pairs of muts meeting specified conditions.
// Useful for gathering LD stats.
//
// Params:
//
//    maxSep - max distance between muts.  We will consider all mut pairs separated
//     by at most this distance.
//
//    binFn - the binary function to be called for all mut pairs within 'maxSep'.
//
template <class BinFn>
void Mutlist::forMutPairs( plen_t maxSep, BinFn binFn ) {
	chkCond( frozen );

	for ( const_iterator m1 = muts.begin(); m1 != muts.end(); m1++ )
		 for ( const_iterator m2 = boost::next( m1 ); m2 != muts.end() && ( get_ploc( m2->loc ) - get_ploc( m1->loc ) < maxSep ); m2++ )
				binFn( *m1, *m2 );
	
}  // Mutlist::forMutPairs()


}  // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_MUTLIST_H

