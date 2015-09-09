#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <fstream>
#include <boost/format.hpp>
#include <boost/scoped_array.hpp>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>
#include <cosi/mutlist.h>
#include <cosi/utils.h>
#include <cosi/demography.h>
#include <cosi/output.h>

namespace cosi {

#define ForEach BOOST_FOREACH
  
// Method: print_haps
//
// Write the output of one simulation in cosi format.
//
// Params:
//
//    demography - the <Demography> object describing what populations are there, size of sample from each pop, etc
//    filebase - the common filename prefix of the output files.  For each population with <popid> of 'i', two files will be
//      written: a 'filebase'.pos-i file giving SNP positions and a 'filebase'.hap-i file giving the haplotypes.
//    length - the length of the simulated region, in bases
//    mutlist - the list of mutations
//    inf_sites - whether to use an infinite-sites model (if not, of mutations falling on the same integer coordinate only
//      the first is kept)
void print_haps(DemographyP demography, const string& filebase, len_bp_int_t length, MutlistP mutlist, bool_t inf_sites) {

  FILE *outf = NULL;
  freq_t freq;

  size_t nmuts = mutlist->size();
  size_t line_size = 64 + 2 * nmuts;
  char *line = (char *)malloc( line_size );
  char *blank_line = (char *)malloc( line_size );
	 
  char *p = blank_line;
  for ( size_t i = 0; i < nmuts; i++ ) {
		*p++ = '2';
		*p++ = ' ';
  }
  *p = 0;

  size_t the_buf_size = 16777216;
  char *the_buf = (char *)malloc( the_buf_size );
	leaf_id_t leaf = 0;
	const vector< popid >& popNames = demography->getPopNames();
	const vector< nchroms_t >& sampleSizes = demography->getSampleSizes();
  for (size_t ipop = 0; ipop < popNames.size(); ipop++) {
    if ( sampleSizes[ipop] > 0) {

			vector< nchroms_t > mutcount( nmuts );

			string filename( (boost::format( "%s.hap-%d" ) % filebase % popNames[ipop]).str() );
			outf = fopen(filename.c_str(), "wt");
			if (outf == NULL) {fprintf(stderr, "Could not open %s\n", filename.c_str());}
			if ( the_buf ) setvbuf( outf, the_buf, _IOFBF, the_buf_size );

#if 0
			static char buf[65536];
			setvbuf( outf, buf, _IOFBF, 65536 );
#endif		
			leaf_id_t popEndLeaf = leaf + sampleSizes[ipop];

			for ( ; leaf < popEndLeaf; leaf++) {
				fprintf(outf, "%d\t%d\t", leaf, ToInt( popNames[ipop] ) );
			
				memcpy( line, blank_line, line_size ); 
			
				const vector< Mutlist::const_iterator >& leafMuts = mutlist->getLeafMuts( leaf );
				ForEach( Mutlist::const_iterator m, leafMuts ) {
					int mutId = m->mutId;
					line[ 2 * mutId ] = '1';
					mutcount[ mutId ]++;
				}
			
				fputs( line, outf );
				fputs( "\n", outf);
			}  // write out haps for this pop

      fclose(outf);
			
			string pos_filename( (boost::format( "%s.pos-%d" ) % filebase % popNames[ipop]).str() );

      outf = fopen(pos_filename.c_str(), "wt");
			if ( the_buf ) setvbuf( outf, the_buf, _IOFBF, the_buf_size );
			
      if (outf == NULL) {fprintf(stderr, "Could not open %s\n", pos_filename.c_str());}
      fprintf(outf, "SNP     CHROM   CHROM_POS       ALLELE1 FREQ1   ALLELE2 FREQ2\n");
      BOOST_AUTO( it, mutlist->getMuts().begin() );
      for (size_t im = 0; im < nmuts; im++, it++) {
				freq = (freq_t) mutcount[im] / sampleSizes[ipop];
				if (inf_sites) {
					fprintf(outf, "%d\t1\t%.4f\t1\t%.4f\t2\t%.4f\n", (int)(it->mutIdOrig+1), double( length * get_loc( it->loc ) ), 
									double( freq ), double( 1 - freq ) );
				}
				else {
					fprintf(outf, "%d\t1\t%d\t1\t%.4f\t2\t%.4f\n", (int)(it->mutIdOrig+1), (int) (length * get_loc( it->loc ) ), 
									double( freq ), double( 1 - freq ) );
				}
			}
		}  // if this pop is nonempty
		fclose(outf);
  }  // for each pop

	using util::cosi_free;
  cosi_free( line );
  cosi_free( blank_line );
  cosi_free( the_buf );
}

void print_mut_contexts( DemographyP demography, const string& filebase, len_bp_int_t length, const mutcontext::mutContexts_t& mutContexts ) {

	using mutcontext::mutContexts_t;

	const vector< popid >& popNames = demography->getPopNames();
	const vector< nchroms_t >& sampleSizes = demography->getSampleSizes();
	leaf_id_t leaf = 0;
  for (size_t ipop = 0; ipop < popNames.size(); ipop++) {
    if ( sampleSizes[ipop] > 0) {
			string filename( (boost::format( "%s.mutContexts-%d.tsv" ) % filebase % popNames[ipop]).str() );
			std::ofstream out( filename.c_str() );

			out << "leafId" << "\t" << "pop";
			ForEach( mutContexts_t::value_type mc, mutContexts )
				 out << "\tbeg_" << static_cast<len_bp_int_t>( ToDouble( mc.first ) * length ) << "\t" << "end_" << static_cast<len_bp_int_t>( ToDouble( mc.first ) * length );
			out << "\n";

			leaf_id_t popEndLeaf = leaf + sampleSizes[ipop];
			for ( ; leaf < popEndLeaf; leaf++) {
				out << leaf << "\t" << ToInt( popNames[ipop] );
				ForEach( mutContexts_t::value_type mc, mutContexts )
					 out <<
					 "\t" << static_cast<int>( length * ToDouble( mc.second[ leaf ].getBeg() ) ) <<
					 "\t" << static_cast<int>( length * ToDouble( mc.second[ leaf ].getEnd() ) );
				out << "\n";
			}
		}
	}
}  // print_mut_contexts

//
// Class impl: ARGOutputHook
//

void ARGOutputHook::handle_add_edge( nodeid nodeId_moreRecent,
																		 nodeid nodeId_lessRecent,
																		 genid genId_moreRecent,
																		 genid genId_lessRecent, const Seglist *seglist,
																		 edge_kind_t edgeKind ) {
	std::cout << "E " << edgeKind2code( edgeKind ) << " " << nodeId_moreRecent << " " << nodeId_lessRecent << " "
						<< genId_moreRecent << " " << genId_lessRecent;

	BOOST_FOREACH( const seglist::Seg& s, *seglist )
		std::cout << " " << s.getBeg() << " " << s.getEnd(); 

	std::cout << "\n";
}

char ARGOutputHook::edgeKind2code( edge_kind_t edgeKind ) {
	assert( edgeKind == EDGE_RECOMB || edgeKind == EDGE_GC || edgeKind == EDGE_COAL );
	assert( int(EDGE_RECOMB) == 0  && int(EDGE_GC) == 1 && int(EDGE_COAL) == 2 ); 
	return "RGC"[ int(edgeKind) ];
}


ARGOutputHook::~ARGOutputHook() { }

#ifdef COSI_TREE_OUTPUT
namespace {
static double genFactor;

void write_tree( leafset_p t, nodeid nodeId ) {
	using std::cout;
	if ( leafset_is_singleton( t ) )
		 cout << ( leafset_get_singleton_leaf( t ) + 1 );
	else {
		cout << "(";
		write_tree( t->childA, t->nodeIdA );
		cout << ":" << ( t->gen - t->childA->gen ) * genFactor << ",";
		write_tree( t->childB, t->nodeIdB );
		cout << ":" <<  ( t->gen - t->childA->gen ) * genFactor
				 << ")" << nodeId;
	}
}
}  // anonymous namespace

void output_trees( len_bp_t region_len, nchroms_t N0_tot ) {
	typedef std::map< loc_t, leafset_p > m_t;
	extern m_t loc2tree;

	//std::cerr << "output_trees: " << loc2tree.size() << " segs.\n";

	genFactor = 1. / ( 2. * N0_tot );
	//std::cerr.precision(15);

	loc_t lastEnd( MIN_LOC );
	//bool isFirst = True;
	for ( m_t::const_iterator it = loc2tree.begin(); it != loc2tree.end(); ++it ) {
		len_bp_t seg_len_bp = region_len * ToDouble( it->first - lastEnd );
		// std::cerr << "seg: " << it->first << ": " << it->second << " seg_len_bp=" << seg_len_bp <<
		// 	 " lastEnd=" << lastEnd << " isFirst=" << isFirst << "\n";
		if ( seg_len_bp > .5 ) {
			std::cout << "[" << (long)( round( ToDouble( seg_len_bp ) ) ) << "]";
			write_tree( it->second, it->second->nodeId );
			std::cout << ";\n";
		}
		//util::chkCond( isFirst || it->first == lastEnd );
		lastEnd = it->first;
		//isFirst = False;
		//std::cout << ":" << get_ploc( it->first ) << " ;\n";
	}
	//util::chkCond( isFirst || lastEnd == MAX_LOC );
	loc2tree.clear();
}
#endif  // #ifdef COSI_TREE_OUTPUT

}  // namespace cosi
