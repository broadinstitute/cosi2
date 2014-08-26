//
// Header: condsnp.h
//
// Conditioning the simulations on the existence of a SNP at a given location with specified allele freqs.
//

#ifndef COSI_INCLUDE_CONDSNP_H
#define COSI_INCLUDE_CONDSNP_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <utility>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <cosi/defs.h>
#include <cosi/decls.h>
#include <cosi/hooks.h>
#include <cosi/seglist.h>

namespace cosi {

//
// Class: CondSnpDef
//
// Criteria for a SNP: its location, and allowed frequencies.
//
class CondSnpDef {
public:
	 CondSnpDef(): loc( NULL_LOC ) { }
	 CondSnpDef( std::string def_ ) { init( def_ ); }

	 // Method: init
	 // Parse the SNP condition, from a format such as ".5/1:35-35,4:10-15"
	 // where .5 is the SNP loc, and 1:35-35 means in pop 1 its sample freq must be
	 // 35 while in pop 4 it must be between 10 and 15 (inclusive).
	 void init( std::string def_ ) {
		 using std::vector;
		 using std::string;
		 using std::invalid_argument;
		 using std::make_pair;
		 using boost::algorithm::split;
		 using boost::algorithm::is_any_of;
		 using boost::lexical_cast;

		 try {
			 vector< string > locAndCondDefs;
			 split( locAndCondDefs, def_, is_any_of( "/" ) );
			 if ( locAndCondDefs.size() > 2 )
					throw invalid_argument( "invalid CondSnpDef: " + def_ + " - too many parts separated by /" );
			 this->loc = lexical_cast< ploc_t >( locAndCondDefs.front() );

			 pop2cond.clear();
			 if ( locAndCondDefs.size() > 1 ) {
				 vector< string > condDefs;
				 split( condDefs, locAndCondDefs[1], is_any_of( ",;" ) );
				 BOOST_FOREACH( string condDef, condDefs ) {
					 vector< string > condParts;
					 split( condParts, condDef, is_any_of( ":" ) );
					 if ( condParts.size() != 2 ) throw invalid_argument( "invalid part of conditional snp def: " + condDef );
					 pop2cond.insert( make_pair( lexical_cast< popid >( condParts.front() ), util::ValRange<nchroms_t>( condParts.back() ) ) );
				 }
			 }
		 } catch( std::exception& e ) {
			 std::cerr << "error: " << e.what() << "\n";
		 }
		 
	 }

	 loc_t getLoc() const { return loc; }

	 bool sampleFreqsMatch( const std::map< popid, nchroms_t >& pop2Counts ) const;

private:
	 // Field: loc
	 // The location of the SNP.
	 loc_t loc;


	 typedef std::map< popid, util::ValRange< nchroms_t > >	pop2cond_map_t;
	 
	 // Field: pop2cond
	 // Map from population to the frequency restriction in that population.
	 pop2cond_map_t pop2cond;
	 
	 friend std::ostream& operator<<( ostream& s, const CondSnpDef& condSnpDef ) {
		 s << condSnpDef.loc << "/" << condSnpDef.pop2cond;
		 return s; 
	 }

	 friend std::istream& operator>>( std::istream& s, CondSnpDef& csd ) {
		 std::string def;
		 s >> def;
		 csd.init( def );
		 return s;
	 }
	 
};  // class CondSnpDef

inline std::istream& operator>>( std::istream& s, boost::optional<CondSnpDef>& csd_opt ) {
	CondSnpDef csd;
	s >> csd;
	csd_opt = csd;
	return s;
}


//
// Class: CondSnpMgr
//
// Manage the conditioning of simulations on the presence of a SNP at a specified location with the
// specified features.
//
//
class CondSnpMgr: public Hook {
public:

	 CondSnpMgr( DemographyP demography_, CondSnpDef condSnpDef_);
	 virtual ~CondSnpMgr();
	 // Virtual method: handle_add_edge
	 //
	 // Called after adding an edge to the ARG.  The new edge may be the result of a coalescence,
	 // a recombination or a gene conversion.
	 //
	 // Params:
	 //
	 //   nodeId_moreRecent - <nodeid> of the more recent (lower gen) node of the edge
	 //   nodeId_lessRecent - <nodeid> of the less recent (higher gen) node of the edge.
	 //   genId_moreRecent - <genid> of the more recent (lower gen) node of the edge
	 //   genId_lessRecent - <genid> of the less recent (higher gen) node of the edge.
	 //   seglist - the segments inherited along the edge.  NOTE: this seglist may be destroyed
	 //       after the call, so make a copy if you need to save it.
	 //       
	 virtual void handle_add_edge( nodeid nodeId_moreRecent,
																 nodeid nodeId_lessRecent,
																 genid genId_moreRecent,
																 genid genId_lessRecent, const Seglist *seglist,
																 edge_kind_t edgeKind );

	 virtual void handle_sweep_end( leafset_p sel_leaves_ ) {
		 this->sel_leaves = sel_leaves_;
	 }

	 void printResults( MutlistP muts, GenMapP genMap ) const;
	 
private:
	 // Field: demography
	 // The demography object, for mapping leaves to pops.
	 DemographyP demography;
	 
	 // Field: condSnpDef
	 // If gathering stats for a snp at a given loc and given freqs, the loc and the freqs.
	 CondSnpDef condSnpDef;

	 // Field: condSnpTreeTimeTot
	 // Total tree time in the tree at the cond snp's loc.
	 gens_t condSnpTreeTimeTot;

	 // Field: condSnpTreeTimeMatching
	 // Tree time on edges on which a mutation would produce a snp matching <condSnpDef> criteria.
	 gens_t condSnpTreeTimeMatching;

	 std::vector< leafset_p > savedLeafsets;
	 std::vector< gens_t > savedEdgesCumLens;

	 leafset_p sel_leaves;

	 size_t nwrongEdges;
	 
	 
};  // class CondSnpMgr


}  // namespace cosi

#endif // #ifndef COSI_INCLUDE_CONDSNP_H
