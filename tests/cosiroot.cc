
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <cosi/general/math/cosirand.h>
#include <cosi/coalescent.h>
#include <cosi/output.h>
#include <cosi/file.h>
#include <cosi/mutate.h>

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>

static void
printusage (void)
{
	fprintf(stderr, "usage: ");
	fprintf(stderr, "coalescent -p param_file ");
	fprintf(stderr, "[-l logfile] [-s segfile]\n");
  
	exit (EXIT_FAILURE);
}

namespace cosi {

class MutProcessor_Histogrammer: public MutProcessor {
public:
	 MutProcessor_Histogrammer();
	 virtual ~MutProcessor_Histogrammer();

	 virtual void processMut(loc_t loc, leafset_p leaves, genid gen, popid popName);

	 void saveTo( const char *fname ) {
		 using boost::make_shared;
		 using boost::shared_ptr;
		 shared_ptr<TCanvas> c1 = make_shared<TCanvas>("c1","c1",800,1000);
		 gPad->SetLogy(1);
		 hist->Draw();
		 c1->Print( fname );
	 }

private:
	 TH2F *hist;
	 
};  // class MutProcessor_Histogrammer

MutProcessor_Histogrammer::MutProcessor_Histogrammer() {
	hist = new TH2F( "freqages", "freq+age histogram", 100, 0.0, 1.0, 100, 0.0, 10000 );
	hist->GetXaxis()->SetTitle("total freq");
	hist->GetYaxis()->SetTitle("age");
	hist->SetBit(TH1::kCanRebin);
}

MutProcessor_Histogrammer::~MutProcessor_Histogrammer() {
}

void MutProcessor_Histogrammer::processMut(loc_t loc, leafset_p leaves, genid gen, popid popName) {
	nchroms_t leafsetSize = leafset_size( leaves );
	if ( leafsetSize > 360 ) PRINT( leafsetSize );
	hist->Fill( double(leafsetSize) / 360.0, ToDouble( gen ) );
}


int 
do_main(int argc, char *argv[]) {
	
  using namespace std;
	using boost::make_shared;
	using boost::shared_ptr;
	using namespace cosi;

	shared_ptr<MutProcessor_Histogrammer> histogrammer = make_shared<MutProcessor_Histogrammer>();

	string paramfile, segfile, logfile, outfilebase;

	bool_t msOutput = False;
	FILE *segfp = NULL, *logfp = NULL;
	len_bp_int_t len_length(-1);
	bool_t sim_only = False;
	bool_t verbose = False;
	bool_t treeSizeOnly = False;
	int nsims = 1;
		 
	bool_t usage = False;
	string prog = argv[0];
	vector<string> args( argv, argv + argc );
	size_t i = 0;
		 
	while ( i+1 < args.size() && (args[++i][0] == '-') ) {
		char c = args[i][1];
		switch (c) {
				 
		case 'l':      /* set logfile */
			logfile = args[++i];
			if (DEBUG)
				 printf("logfile: %s\n", logfile.c_str());
			break;
				 
		case 'p':      /* set param file */
			paramfile = args[++i];
			if (DEBUG) 
				 printf ("paramfile: %s\n", paramfile.c_str());
			break;
				 
		case 's':      /* set segsites file */
			segfile = args[++i];
			if (DEBUG)
				 printf("segfile: %s\n", segfile.c_str());
			break;
				 
		case 'o':      /* set haplo/marker base file name */
			outfilebase = args[++i];
			break;

		case 'n':
			istringstream( args[++i] ) >> nsims;
			break;

		case 'm': msOutput = True; break;
		case 't': treeSizeOnly = True; break;

		case 'v':      /* set verbose flag. does nothing yet. */
			verbose = True;
			if (DEBUG)
				 printf("set to verbose\n");
			break;
				 
		case 'S':      /* do the simulation and quit -- for profiling/development purposes only */
			sim_only = True;
			break;
				 
		default:
			fprintf(stderr, 
							"WARNING coalescent: illegal option %c\n", 
							c);
			usage = True;
			break;				
		}
	}
		 
	if (usage)
		 printusage();
		 
	else if (paramfile.empty()) {
		fprintf(stderr, "coalescent: must use -p option\n");
		printusage();
	}
	else {
		if (!segfile.empty()) {
			if ((segfp = fopen(segfile.c_str(), "w")) == NULL) {
				fprintf(stderr, 
								"ERROR %s: cannot open file %s for writing\n",
								prog.c_str(), 
								segfile.c_str());
				exit (EXIT_FAILURE);
			}
		}
			 
		else segfp = NULL;
			 
		if (!logfile.empty()) {
			if ((logfp = fopen(logfile.c_str(), "w")) == NULL) {
				fprintf(stderr, 
								"ERROR %s: cannot open file %s for writing\n",
								prog.c_str(), logfile.c_str());
				exit(EXIT_FAILURE);
			}
		}
		else logfp = NULL;
	}  // usage ok, param file given

	RandGenP randGen;
	for ( int simNum = 0; simNum < nsims; simNum++ ) {
		PRINT( simNum );
		CoSi cosi;

		cosi.set_segfp( segfp );
		cosi.set_logfp( logfp );
		cosi.set_verbose( verbose );
	

		if ( simNum == 0 ) {
			cosi.setUpSim( paramfile );
			randGen = cosi.getRandGen();
			if ( !msOutput ) std::cerr << "coalescent seed: " << randGen->getSeed() << "\n";
		} else
			 cosi.setUpSim( paramfile, randGen );

		cosi.setMutProcessor( histogrammer );
	
		cosi.runSim();
	}

	histogrammer->saveTo( "/home/unix/ilya/public_html/t.jpg" );

	if ( segfp ) fclose( segfp );
	if ( logfp ) fclose( logfp );
	
  return EXIT_SUCCESS;
}

}  // namespace cosi

int main( int argc, char *argv[] ) { return cosi::do_main( argc, argv ); }
