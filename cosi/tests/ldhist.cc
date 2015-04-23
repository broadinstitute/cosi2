//
// Program: ldhist.cc
//
// Plot the decay of LD as a function of distance between SNPs.
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <TH1I.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TStyle.h>

#include <cosi/mutlist.h>

namespace cosi {

using boost::make_shared;
using boost::shared_ptr;

struct AddToRootHist {
	 TH2I *hist;
	 
	 AddToRootHist( TH2I *hist_ ): hist( hist_ ) {}
	 
	 void operator()( const Mut& m1, const Mut& m2 ) {
		 // PRINT4( m1.loc, m2.loc, ToDouble( get_plen( m2.loc - m1.loc ) ),
		 // 				 leafset_t::compute_r2( m1.leaves, m2.leaves ) );
		 hist->Fill( ToDouble( get_plen( m2.loc - m1.loc ) ),
								 leafset_t::compute_r2( m1.leaves, m2.leaves ) );
	 }
};

int ldhist_main( int argc, char *argv[] ) {

	std::ifstream s( "tst.ms" );

	MutlistP mutlist = Mutlist::loadFromMs( s );
	PRINT( leafset_get_max_leaf_id() );
	// vector< nchroms_t > sampleSizes( 1, leafset_get_max_leaf_id() );
	// mutlist->print_haps_ms( std::cout, sampleSizes );

	TH1::SetDefaultSumw2();
	gStyle->SetOptStat( "neMRou" );
	gStyle->SetOptStat(0);
	
	plen_t maxSep(0.05);

	TH2I *ldHist = new TH2I( "LDvsDist", "LD vs distance", 20, 0.0, ToDouble( maxSep ),
													 40, 0.0, 1.01 );
	ldHist->GetXaxis()->SetTitle( "plen" );
	ldHist->GetYaxis()->SetTitle( "r^2" );
	//ldHist->SetBit( TH1::kCanRebin );

	TH1I *ldOnlyHist = new TH1I( "LDonly", "LD only", 40, 0.0, 1.01 );
	TH1I *sepOnlyHist = new TH1I( "seponly", "sep only", 20, 0.0, ToDouble( maxSep ) );

	TProfile *ldProf = new TProfile( "ldProf", "LD profile", 20, 0.0, ToDouble( maxSep ), 0.0, 1.0 );

	
	const vector<Mut>& muts = mutlist->getMuts();
	for ( auto m1 = muts.begin(); m1 != muts.end(); m1++ ) {
		for ( auto m2 = boost::next( m1 ); m2 != muts.end() && ( get_ploc( m2->loc ) - get_ploc( m1->loc ) < maxSep ); m2++ ) {

			plen_t pairSep = get_plen( m2->loc - m1->loc );
			cosi_double pair_r2 = leafset_t::compute_r2( m1->leaves, m2->leaves );

			if ( !(boost::math::isnan)( pair_r2 ) ) {

				//PRINT4( m1->loc, m2->loc, pairSep, pair_r2 );
			
				ldHist->Fill( ToDouble( pairSep ), pair_r2 );
				ldProf->Fill( ToDouble( pairSep ), pair_r2 );
				ldOnlyHist->Fill( pair_r2 );
				sepOnlyHist->Fill( ToDouble( pairSep ) );
				
			}

		}
	}
	
	//	mutlist->forMutPairs( plen_t( 1.0 ), AddToRootHist( ldHist ) );

	shared_ptr<TCanvas> c1 = make_shared<TCanvas>( "c1", "canvas 1", 1200, 2000);
	c1->Divide( 1, 4 );
	c1->cd(1);
	ldHist->Draw( "TEXT"  );
	c1->cd(2);
	ldOnlyHist->Draw();
	c1->cd(3);
	sepOnlyHist->Draw();
	c1->cd(4);
	ldProf->Draw();

	c1->Update();
	c1->Print( "/home/unix/ilya/public_html/cc.svg" );
	
	return EXIT_SUCCESS;
}

}

int main( int argc, char *argv[] ) { return cosi::ldhist_main( argc, argv ); }
