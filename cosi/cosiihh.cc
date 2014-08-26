
#include <cosi/coalescent.h>

int 
main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
	using std::ofstream;
	using namespace cosi;

	
	
	ofstream ihhout( getenv( "IHHVALSFN" ) ? getenv( "IHHVALSFN" ) : "ihhvals.tsv" );
  ihhout << "rseed\tnder\tnanc\tDleft\tAleft\tDright\tAright\n";
  for ( int i = 0; i < 10; i++ ) {
	 cout << "\n\nRUN " << i << endl;
	 CoSi cosi;
	 
	 cosi.parseArgs( argc, argv );
	 cosi.setUpSim();

	 ParamFileReaderP params = cosi.getParams();
	 MutlistP muts = cosi.getMuts();
	 MutateP mutate = cosi.getMutate();
	 SweepP sweep = cosi.getSweep();
	 DemographyP demography = cosi.getDemography();
	 
	 sweep->sweep_set_fix_freq( True );
	 
	 cosi.runSim();

	 mutate->addSamplingLocMut();

	 leafset_p derLeaves = muts->getMut( sweep->getSelLoc() ).leaves,
			ancLeaves = leafset_difference( demography->pop2leaves[0], derLeaves );
	 ehh::ihh_t derR = ehh::computeIHH( muts, loc_t( 0.5 ), DIR_R,
																			derLeaves,
																			demography->pop2leaves[0] );
	 ehh::ihh_t derL = ehh::computeIHH( muts, loc_t( 0.5 ), DIR_L,
																			derLeaves,
																			demography->pop2leaves[0] );
	 ehh::ihh_t ancR = ehh::computeIHH( muts, loc_t( 0.5 ), DIR_R,
																			ancLeaves,
																			demography->pop2leaves[0] );
	 ehh::ihh_t ancL = ehh::computeIHH( muts, loc_t( 0.5 ), DIR_L,
																			ancLeaves,
																			demography->pop2leaves[0] );

	 /*
		 so, we can save, for each sim:
		    the freqs
		    the freqs of the four haps between this and the causal (DD, DA, AD, AA)
				and the left and right ihh values, at this SNP and the causal.
		*/
				 
	 PRINT6( leafset_cardinality( derLeaves ), leafset_cardinality( ancLeaves ), derL, ancL, derR, ancR );
	 ihhout <<
			params->getRandomSeed() << "\t" <<
			leafset_cardinality( derLeaves ) << "\t" <<
			leafset_cardinality( ancLeaves ) << "\t" <<
			derL << "\t" << ancL << "\t" << derR << "\t" << ancR << "\n";
  }
  return 0;
}
