//
// Program: cosimain
//
// The top-level executable for the cosi simulator.
// Parses command-line args, and invokes the cosi library to run the simulation.
//

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <ios>
#include <string>
#include <boost/program_options.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/optional/optional_io.hpp>
#ifndef COSI_NO_CPU_TIMER
#include <boost/timer/timer.hpp>
#endif
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>

#include <cosi/defs.h>
#include <cosi/mutlist.h>
#include <cosi/mutate.h>
#include <cosi/cositop.h>
#include <cosi/file.h>
#include <cosi/demography.h>
#include <cosi/genmap.h>
#include <cosi/recomb.h>
#include <cosi/condsnp.h>
#include <cosi/cosicfg.h>
#include <cosi/customstats.h>

namespace cosi {

using util::ToString;

void on_print(const std::string& str);

void on_print(const std::string& str)
{
  std::cout << str;
}

CoSiMain::CoSiMain():
	deltaTfactor( 1.0 ),
	msOutput( False ), outputMutGens( False ), outputRecombLocs( False ), segfp( NULL ), logfp( NULL ), len_length(-1),
	sim_only( False ), showProgress( 0 ), verbose( False ), treeSizeOnly( False ), nsims( 1 ),
	showNumRecombs( False ), outputTreeStats( False ), outputPrecision( 6 ), randSeed( 0 )
#ifdef COSI_SUPPORT_COALAPX	
	, maxCoalDist( 1.0 ), maxCoalDistCvxHull( False )
#endif
	, genMapShift( 0 ), sweepFracSample( False ), outputSimTimes( False ), outputEndGens( False ), stopAfterMinutes( 0 ),
							 outputARGedges( False ), freqsOnly( False ),
							 dropSingletonsFrac( 0 ), genmapRandomRegions( False ), outputPopInfo( False ),
							 outputGenMap( False ), customStats( False )
{
}

void CoSiMain::printCompileOptions() {
	using std::cerr;

	cerr << "cosi coalescent simulator, version 2.0\n\n";
	cerr << "Compile-time options:\n";

	cerr << "   COSI_SUPPORT_COALAPX (support for approximating the coalescent): " <<
		 IFELSE_COSI_SUPPORT_COALAPX( 1, 0 ) << "\n";
	
}  // CoSiMain::printCompileOptions

double poisPrec = 1e-5;
unsigned int poisMaxSteps = 100000;

int
CoSiMain::parse_args( int argc, char *argv[] ) {
	using std::cout;
  using std::vector;
	using std::string;
	using std::istringstream;

	extern bool apxWithTrajOk;
	extern factor_t apxMinFactor;
	extern double coalApxRejSamp_u;
	extern bool sweep3_no_oneSidedRecombs;
	extern bool showOneSimProgress;

	namespace po = boost::program_options;

	string prog( argv[0] );

	po::options_description main_options( "Specifying the model" );
	main_options.add_options()
		 ( "paramfile,p", po::value(&paramfile)->required(), "parameter file" )
		 ( "recombfile,R", po::value(&recombfileFN), "genetic map file (if specified, overrides the one in paramfile)" )
		 ( "genmapRandomRegions", po::bool_switch(&genmapRandomRegions),
			 "for each simulation use a randomly chosen subregion of the genetic map" )
		 ( "trajfile,J", po::value(&trajFN), "file from which to read sweep trajectory" )
		 ( "nsims,n", po::value(&nsims)->default_value(1), "number of simulations to output" )
		 ( "seed,r", po::value(&randSeed)->default_value(0), "random seed (0 to use current time)" ) 
#ifdef COSI_SUPPORT_COALAPX
		 ( "max-coal-dist,u", po::value(&maxCoalDist), "max dist betw segs for coalescence" )
		 ( "max-coal-dist-cvx-hull,U", po::bool_switch(&maxCoalDistCvxHull), "use convex hull for max coal dist" )
#endif
		 ( "genmapshift,G", po::value(&genMapShift), "shift all genmap locations by this delta" )
		 ( "sweep-frac-sample,E", po::bool_switch(&sweepFracSample), "sweep end freq specifies exact sample fraction")
		 ( "apx-with-traj-ok", po::bool_switch(&apxWithTrajOk), "enable coalapx even when pop size traj specified (experimental)" )
		 ( "apx-min-factor", po::value(&apxMinFactor), "enable coalapx only when pop size exceeds active sample size by this factor (experimental)" )
		 ( "apx-rej-samp", po::value(&coalApxRejSamp_u)->default_value(1.0), "use rejection sampling (experimental)" )
		 ( "sweep3-no-one-sided-recombs", po::bool_switch(&sweep3_no_oneSidedRecombs), "handle one-sided recombs (experimental)" )
		 ( "show-one-sim-progress", po::bool_switch(&showOneSimProgress), "show progress of coalescence in each sim" )
		 ( "pois-max-steps", po::value(&poisMaxSteps)->default_value(100000), "max # of steps when evaluating waiting times" ) 
		 ( "pois-prec", po::value(&poisPrec)->default_value(1e-5), "precision when evaluating waiting times" )
#ifdef COSI_FREQONLY		 
		 ( "freqs-only", po::bool_switch(&freqsOnly), "output frequencies only" )
#endif		 
#ifdef COSI_CONDSNP
		 ( "condsnp,c", po::value(&this->condSnpDef), "condition sims on a SNP at this loc with these freqs" )
#endif		 
		 ;

	po::options_description output_options( "Specifying the output format" );
	output_options.add_options()
		 ( "outfilebase,o", po::value(&outfilebase), "base name for output files in cosi format" )
		 ( "outms,m", po::bool_switch(&msOutput), "write output to stdout in ms format" )
		 ( "output-pop-info", po::bool_switch(&outputPopInfo), "output pop info in ms format output" )
		 ( "output-gen-map", po::bool_switch(&outputGenMap), "output genetic map in ms format output" )
		 ;

	po::options_description output_details_options( "Specifying output details" );
	output_details_options.add_options()
		 ( "output-precision,P", po::value( &outputPrecision ), "number of decimal places used for floats in the outputs" )
		 ( "trajoutfile,j", po::value(&trajOutFN), "file to which to output sweep trajectory" )
		 
		 ( "write-tree-stats,T", po::bool_switch(&outputTreeStats), "output tree stats" )
		 ( "write-mut-ages,M", po::bool_switch(&outputMutGens), "output mutation ages" )
		 ( "write-recomb-locs,L", po::bool_switch(&outputRecombLocs), "output recombination locations" )
		 ( "output-ARG,e", po::bool_switch(&outputARGedges), "output ARG edges.  Edges are written in ms output format (see -m option"
			 "), one per line, after the '//' line but before the 'segsites: ' line of each simulation. "
			 "Format is:\n\tE <edgeKind> <node_1_id> <node_2_id> <node_1_generation> <node_2_generation> <seg_1_beg> <seg_1_end> ... <seg_k_beg> <seg_k_end>.\n"
			 "Edge kinds are: R, recombination; G, gene conversion; C, coalescence. "
			 "seg_i_beg, seg_i_end give chromosomal segments inherited along the edge; locations are values in [0.0,1.0] representing locations "
			 "within the simulated region.")
		 ( "write-mut-contexts,C", po::value(&outputMutContextsFor)->value_name( "position" ),
			 "output mutation contexts for these locations" )
		 ( "drop-singletons", po::value(&dropSingletonsFrac)->value_name("fraction"),
			 "drop this fraction of singleton SNPs" )
		 ;

	po::options_description misc_options( "Misc options" );
	misc_options.add_options()
		 ( "help,h", "produce help message" )
		 ( "version,V", "print version info and compile-time options" )
		 ( "verbose,v", po::bool_switch(&verbose), "verbose output" )
		 ( "show-progress,g", po::value( &showProgress )->implicit_value( 10 ), "show progress" );

	po::options_description dev_options( "Developer options" );
	dev_options.add_options()
		 ( "logfile,l", po::value(&logfile), "log file" )
		 ( "segfile,s", po::value(&segfile), "seg file" )

		 
		 ( "deltaTfactor,d", po::value(&deltaTfactor)->default_value(factor_t(1.0)), "delta factor for sweep" )

		 ( "tree-size-only,t", po::bool_switch(&treeSizeOnly), "compute tree size only" )
		 ( "show-num-recombs,k", po::bool_switch(&showNumRecombs), "print number of recombs" )
		 ( "sim-and-quit,S", po::bool_switch(&sim_only), "do simulation and quit" )
		 ( "output-sim-times", po::bool_switch(&outputSimTimes), "for each sim output the time it took" )
		 ( "output-end-gens", po::bool_switch(&outputEndGens), "for each sim output the generation at which it ended" )
		 ( "stop-after-minutes", po::value(&stopAfterMinutes)->default_value(0.0), "stop simulation after this many minutes" )
		 ( "custom-stats", po::bool_switch(&customStats), "compute custom stats" )
		 ;

	po::options_description cosi_options;
	cosi_options.add( main_options ).add( output_options ).add( output_details_options ).add( misc_options ).add( dev_options );

	if ( argc == 1 ) {
    cout << cosi_options << "\n";
    return EXIT_FAILURE;
	}

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cosi_options), vm);

	if (vm.count("help")) {
    cerr << cosi_options << "\n";
    return EXIT_FAILURE;
	}
	if (vm.count("version")) {
		printCompileOptions();
    return EXIT_FAILURE;
	}
	
	po::notify(vm);

#ifdef COSI_FREQONLY
	freqsOnly = true;
#endif	


	if (!segfile.empty()) {
		if ((segfp = cosi_fopen(segfile, "w")) == NULL) {
			fprintf(stderr, 
							"ERROR %s: cannot open file %s for writing\n",
							prog.c_str(), 
							segfile.BOOST_FILESYSTEM_C_STR);
			exit (EXIT_FAILURE);
		}
	}
			 
	else segfp = NULL;
			 
	if (!logfile.empty()) {
		if ((logfp = cosi_fopen(logfile, "w")) == NULL) {
			fprintf(stderr, 
							"ERROR %s: cannot open file %s for writing\n",
							prog.c_str(), logfile.BOOST_FILESYSTEM_C_STR);
			exit(EXIT_FAILURE);
		}
	}
	else logfp = NULL;

	return EXIT_SUCCESS;
}  // CoSiMain::parse_args

int 
CoSiMain::cosi_main(int argc, char *argv[]) {
	using std::cout;
	using std::cerr;

	if ( parse_args( argc, argv ) == EXIT_FAILURE ) return EXIT_FAILURE;

	PRINT( nsims );
	RandGenP randGen;
	if ( this->randSeed ) randGen = boost::make_shared<RandGen>( this->randSeed );
#ifndef COSI_NO_CPU_TIMER
	boost::timer::cpu_timer overallTimer;
	double stopAfterNs = stopAfterMinutes * 1e9 * 60.0;
#endif	
	for ( int simNum = 0; simNum < nsims; simNum++ ) {
#ifndef COSI_NO_CPU_TIMER		
		if ( stopAfterNs > 0 && overallTimer.elapsed().wall > stopAfterNs ) {
			cout << "// cosi-early-exit\n";
			cerr << "cosi: exiting after " << overallTimer.elapsed().wall << " ns; completed " <<
				 simNum << " of " << nsims << " sims.\n";
			break;
		}
		
		boost::timer::cpu_timer cpuTimer;
#endif		
		if ( showProgress && !( simNum % showProgress ) ) { cerr << " sim " << simNum << " of " << nsims << endl; }
		CoSi cosi;

		cosi.set_segfp( segfp );
		cosi.set_logfp( logfp );
		cosi.set_verbose( verbose );
		cosi.set_trajFN( trajFN );
		cosi.set_trajOutFN( trajOutFN );
		cosi.set_outputTreeStats( outputTreeStats );
		cosi.set_outputMutGens( outputMutGens );
		cosi.set_outputRecombLocs( outputRecombLocs );
		cosi.set_deltaTfactor( deltaTfactor );
#ifdef COSI_SUPPORT_COALAPX		
		cosi.set_maxCoalDist( plen_t( maxCoalDist ) );
		cosi.set_maxCoalDistCvxHull( maxCoalDistCvxHull );
#endif
		cosi.set_sweepFracSample( this->sweepFracSample );

#ifdef COSI_CONDSNP		
		if ( condSnpDef ) cosi.set_condSnpDef( boost::make_shared<CondSnpDef>( *this->condSnpDef ) );
#endif
		cosi.set_recombfileFN( this->recombfileFN );
		cosi.set_outputARGedges( this->outputARGedges );
		cosi.set_genmapRandomRegions( this->genmapRandomRegions );

		cosi.setUpSim( paramfile, randGen );

		if ( simNum == 0 ) {
			randGen = cosi.getRandGen();
			if ( !msOutput ) cerr << "coalescent seed: " << randGen->getSeed() << "\n";
			if ( msOutput ) {
				cout.precision( outputPrecision );
				DemographyP dem = cosi.getDemography();
				cout << "ms " << dem->getTotSamples() << " " << nsims << "\n";
				if ( outputPopInfo ) {
					cout << "pops " << dem->getPopNames().size();
					for ( size_t popNum = 0; popNum < dem->getPopNames().size(); ++popNum )
						 cout << " " << dem->getPopNames()[ popNum ] << " " << dem->getSampleSizes()[ popNum ];
					cout << "\n";
				}
				cout << "cosi_rand " << randGen->getSeed() << "\n\n";
			}
			customstats::init( cosi.getDemography(), nsims );
		}

		using boost::make_shared;
		MutlistP muts = make_shared<Mutlist>();
		if ( dropSingletonsFrac < 1e-10 )
			 cosi.setMutProcessor( make_shared<MutProcessor_AddToMutlist>( muts ) );
		else
			 cosi.setMutProcessor( make_shared<MutProcessor_AddToMutlist_WithAscertainment>( muts,
																																											 dropSingletonsFrac, randGen ) );

		if ( msOutput ) { cout << "// seed=" << randGen->getSeed() << "\n"; }
	
		ParamFileReaderP params = cosi.getParams();
		cosi.getMutate()->setFreqsOnly( freqsOnly );
		genid endGen = cosi.runSim();

		if ( showNumRecombs ) { PRINT( cosi.getRecomb()->getNumRecombs() ); }

		if ( freqsOnly ) cosi.getMutate()->writeTreeSize();
		if ( msOutput || !outfilebase.empty() || cosi.getCondSnpMgr() || customStats ) {
			//PRINT( "freezing" );
			muts->freeze( params->getInfSites() || msOutput || cosi.getCondSnpMgr(),
										cosi.getGenMap()->recomb_get_length() );
			//PRINT( "frozen" );

			if ( cosi.getCondSnpMgr() ) cosi.getCondSnpMgr()->printResults( muts, cosi.getGenMap() );
			
			if (!msOutput && !outfilebase.empty()) {
			  std::ostringstream fbase;
			  fbase << outfilebase.c_str();
			  if ( nsims > 1 )
			    fbase << "_" << simNum;
			  print_haps( cosi.getDemography(), fbase.str(),
										 params->getLength(), muts,
										 params->getInfSites() );
#ifdef COSI_DEV_MUTCONTEXT				 
				 if ( !outputMutContextsFor.empty() )
				   print_mut_contexts( cosi.getDemography(), fbase.str(), params->getLength(),
																mutcontext::getSavedMutContexts() );
#endif				 
			}

			if ( customStats ) {
				customstats::record_sim( cosi.getDemography(), cosi.getGenMap(), params->getLength(),
																 muts, params->getInfSites() );
			}
			
			if ( msOutput ) 
				 muts->print_haps_ms( cout, cosi.getDemography()->getSampleSizes(), cosi.getTreeStatsHook(),
															cosi.get_outputMutGens(),
															outputRecombLocs ? &(cosi.getRecombRecorder()->getRecombLocs()) : NULL,
															outputGenMap,
															cosi.getGenMap(),
															outputPrecision,
#ifndef COSI_NO_CPU_TIMER
															outputSimTimes ? &cpuTimer : NULL,
#else															
															/*outputSimTimes ? &cpuTimer : */NULL,
#endif															
															outputEndGens ? &endGen : NULL );
		}  // output simulation results
		
	}  // for each simulation

	if ( segfp ) fclose( segfp );
	if ( logfp ) fclose( logfp );

	if ( customStats ) customstats::finish();

  return EXIT_SUCCESS;
}  // CoSiMain::cosi_main()

}  // namespace cosi
