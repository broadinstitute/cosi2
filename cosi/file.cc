/* $Id: file.c,v 1.4 2011/05/06 15:35:43 sfs Exp $ */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <sstream>
#include <fstream>
#include <utility>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/exception/all.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/triangle_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/core/demangle.hpp>
#include <cosi/general/utils.h>
#include <cosi/file.h>
#include <cosi/recomb.h>
#include <cosi/demography.h>
#include <cosi/simulator.h>
#include <cosi/output.h>
#include <cosi/geneconversion.h>
#include <cosi/historical.h>
#include <cosi/sweep.h>
#include <cosi/basemodel.h>

namespace cosi {

struct cosi_param_file_error: virtual cosi_io_error {};
typedef boost::error_info<struct errinfo_err_detail_,std::string> errinfo_err_detail;
typedef boost::error_info<struct errinfo_file_name2_, filename_t> errinfo_file_name2;
typedef boost::error_info<struct errinfo_param_name_,std::string> errinfo_param_name;
typedef boost::error_info<struct errinfo_param_details_,std::string> errinfo_param_details;

static const size_t BUF_MAX = 4096;

ParamFileReader::ParamFileReader( DemographyP demography_, RandGenP randGen_ ):
	HasRandGen( randGen_ ), demography( demography_ ) {
	init();
}

void ParamFileReader::init() {
//  seeded = False;
  popsize = 0;
  sampsize = 0;
  length = 0;
  mu = prob_per_bp_per_gen_t( 0.0 );
  recombrate = 0;
	geneConv2RecombRateRatio = factor_t(0.0);
	geneConversionMeanTractLength = 500 /* bp */;
	geneConversionMinTractLength = 4 /* bp */;
	geneConversionModel = GeneConversion::GCM_LEFT_AND_RIGHT_FROM_ORIGIN;
//  rseed = 0;

	histEvents = boost::make_shared<HistEvents>( demography );
	baseModel = boost::make_shared<BaseModel>();
	infSites = False;
//	printSeed = True;
	ignoreRecombsInPop = NULL_POPID;
	if ( !getRandGen() )
		 setRandGen( boost::make_shared<RandGen>() );
}


void ParamFileReader::file_read(boost::filesystem::path filename, FILE *segfp) 
{
	init();
	this->paramFileName = filename;

	try {
		FILE *infileptr = cosi_fopen (filename.c_str(), "r");
		if ( !infileptr ) {
			BOOST_THROW_EXCEPTION( cosi_param_file_error()
														 << boost::errinfo_errno(errno)
														 << error_msg( "Error opening param file" ) );
		}
		if (segfp != NULL) {
		  std::string s = filename.string();
		  fprintf(segfp, "paramfile: %s\n", s.c_str());
		}
		file_get_data (infileptr, segfp);
		fclose(infileptr);
		if ( getenv( "COSI_NEWSIM" ) ) {
			histEvents->constructBaseModel( baseModel );
			// std::cerr.precision(16);
			// std::cerr << *baseModel << "\n";
		}
	} catch( boost::exception& e ) {
		e << boost::errinfo_file_name( filename.string() );
		throw;
	}
}

/*********************************************************/

void ParamFileReader::sample_distribution_values( char *buf ) {
	bool sampled = false;
	while ( char *beg = strstr( buf, "N(" ) ) {
		double mean = NAN, std = NAN;
		if ( sscanf( beg, "N(%lf,%lf)", &mean, &std ) != 2 )
			 BOOST_THROW_EXCEPTION( cosi_param_file_error() <<
															error_msg( "invalid parameter distribution spec" ) );
		boost::random::normal_distribution<double> nd( mean, std );
		std::string s( buf );
		s.replace( beg-buf, strchr( beg, ')' ) -  beg + 1,
							 boost::lexical_cast<std::string>( nd( *getRandGen() ) ) );
		strcpy( buf, s.c_str() );
		sampled = true;
	}

	while ( char *beg = strstr( buf, "U(" ) ) {
		double dmin = NAN, dmax = NAN;
		if ( sscanf( beg, "U(%lf,%lf)", &dmin, &dmax ) != 2 ||
				 !( dmin < dmax ) )
			 BOOST_THROW_EXCEPTION( cosi_param_file_error() <<
															error_msg( "invalid parameter distribution spec" ) );
		boost::random::uniform_real_distribution<double> ud( dmin, dmax );
		std::string s( buf );
		s.replace( beg-buf, strchr( beg, ')' ) -  beg + 1,
							 boost::lexical_cast<std::string>( ud( *getRandGen() ) ) );
		strcpy( buf, s.c_str() );
		sampled = true;
	}

	while ( char *beg = strstr( buf, "T(" ) ) {
		double dmin = NAN, dmid = NAN, dmax = NAN;
		if ( sscanf( beg, "T(%lf,%lf,%lf)", &dmin, &dmid, &dmax ) != 3 ||
				 !( dmin < dmid && dmid < dmax) )
			 BOOST_THROW_EXCEPTION( cosi_param_file_error() <<
															error_msg( "invalid parameter distribution spec" ) );
		boost::random::triangle_distribution<double> td( dmin, dmid, dmax );
		std::string s( buf );
		s.replace( beg-buf, strchr( beg, ')' ) -  beg + 1,
							 boost::lexical_cast<std::string>( td( *getRandGen() ) ) );
		strcpy( buf, s.c_str() );
		sampled = true;
	}

	while ( char *beg = strstr( buf, "E(" ) ) {
		double lambda_arg = NAN;
		if ( sscanf( beg, "E(%lf)", &lambda_arg ) != 1 )
			 BOOST_THROW_EXCEPTION( cosi_param_file_error() <<
															error_msg( "invalid parameter distribution spec" ) );
		boost::random::exponential_distribution<double> ed( lambda_arg );
		std::string s( buf );
		s.replace( beg-buf, strchr( beg, ')' ) -  beg + 1,
							 boost::lexical_cast<std::string>( ed( *getRandGen() ) ) );
		strcpy( buf, s.c_str() );
		sampled = true;
	}
	if (sampled && getenv("COSI_SAVE_SAMPLED")) {
		std::ofstream f(getenv("COSI_SAVE_SAMPLED"), std::ofstream::app);
		f << buf;
	}
	
}

/*********************************************************/

int 
ParamFileReader::file_proc_buff(char *var, char* buffer, FILE* segfp) 
{

	popid     popname;
	int  intarg;

	sample_distribution_values( buffer );
		
	/*	printf("%s\n", var); */
	if (strcmp(var,"length") == 0) {
		length = atoi(buffer);
	}
	else if (strcmp(var,"recomb_file") == 0) {
		if (buffer[strlen(buffer)-1] == '\n')
			 buffer[strlen(buffer)-1] = '\0';
		chkCond( length > 0, "cosi paramfile error: must specify length before recombmap" );

		namespace fs = boost::filesystem;

		fs::path recombFilePath( buffer );
		if ( !this->recombfileFN.empty() ) recombFilePath = this->recombfileFN;
		fs::path paramFileDir( fs::canonical( this->paramFileName ).parent_path() );
		recombFilePath = fs::absolute( recombFilePath /*, paramFileDir */ );

		recombfileFN = recombFilePath;
		
//		genMap = boost::make_shared<GenMap>( recombFilePath, length, this->genMapShift );
//		assert( genMap.get() );
	}
	else if (strcmp(var, "mutation_rate") == 0) {
		mu = prob_per_bp_per_gen_t( atof(buffer) );
	}
	else if (strcmp(var, "infinite_sites") == 0) {
	  if (buffer[strlen(buffer)-1] == '\n')
			 buffer[strlen(buffer)-1] = '\0';
	  if (strcmp(buffer, "yes") == 0 || strcmp(buffer, "Yes") == 0 || strcmp(buffer, "YES") == 0)
			 infSites = True;
	}
	else if (strcmp(var, "gene_conversion_rate") == 0) {
		//geneConversionRate = atof(buffer);
		file_exit( "deprecated_parameter", "please specify gene_conversion_relative_rate giving ratio of gene conversion to recomb rate" );
	}
	else if (strcmp(var, "gene_conversion_relative_rate") == 0) {
		geneConv2RecombRateRatio = factor_t( atof(buffer) );
	}
	else if (strcmp(var, "gene_conversion_mean_tract_length") == 0) {
		geneConversionMeanTractLength = atoi(buffer);
	}
	else if (strcmp(var, "gene_conversion_min_tract_length") == 0) {
		geneConversionMinTractLength = atoi(buffer);
	}
	else if (strcmp(var, "gene_conversion_model") == 0) {
		geneConversionModel = GeneConversion::parseGCModel( boost::trim_copy( string( buffer ) ) );
	}
	else if (strcmp(var, "pop_size") == 0) {
		popname = popid( atoi(strtok (buffer, " ")) );
		intarg = atoi(strtok (NULL, " " ));
		if (FILE_DEBUG)
			 printf("popsize: %d\n", intarg);
		/* 
		 * Throw a fatal error if pop [popname] does not exist.
		 */		
		if (! demography->dg_set_pop_size_by_name (ZERO_GEN, popname, intarg))
			 file_exit("file_proc_buff", 
								 "parameter file - pop specified does not exist.");
		
		baseModel->popInfos[ popname ].setSizeFrom( ZERO_GEN, popsize_float_t( intarg ) );
	}
	else if (strcmp(var, "sample_size") == 0) {		
		popname = popid( atoi(strtok (buffer, " ")) );
		intarg = atoi(strtok (NULL, " " ));
		if (FILE_DEBUG)
			 printf("sampsize: %d\n", intarg);
		demography->dg_populate_by_name (popname, intarg);
		if (segfp != NULL) {fprintf(segfp, "A %d %d\n", ToInt( popname ), intarg);}
	}
	else if ( !strcmp( var, "pop_ignore_recombs" ) ) {
		this->ignoreRecombsInPop = popid( atoi(strtok (buffer, " ")) );
	}
	else if (strcmp(var, "pop_define") == 0) {
		popname = popid( atoi(strtok (buffer, " ")) );
		string pop_label( strtok (NULL, " ") );
		if ( pop_label.empty() ) pop_label = "some_pop";
		if ( pop_label.at( pop_label.size()-1 ) == '\n' )
			 pop_label.resize( pop_label.size()-1 );
		demography->dg_create_pop (popname, pop_label, /* genid=*/ ZERO_GEN);
	}
	else if (strcmp(var, "pop_event") == 0) {
		histEvents->addEvent( histEvents->parseEvent( buffer ) );
	}
	else if (strcmp(var, "random_seed") == 0 ) {
		if ( !getRandGen() ) {
			std::istringstream is( buffer );
			is.exceptions( std::ios_base::failbit | std::ios_base::badbit );

			unsigned long rseed = -1;
			try { is >> rseed; }
			catch( const std::ios_base::failure& e ) {
				BOOST_THROW_EXCEPTION( cosi_param_file_error()
															 << error_msg( "invalid random seed" )
															 << boost::errinfo_errno( errno ) );
			}
			setRandGen( boost::make_shared<RandGen>( rseed ) );
		} else
			 std::cerr << "ignoring random_seed from param file\n";
	} else if (strcmp(var,"sweep_traj_file") == 0) {
		if (buffer[strlen(buffer)-1] == '\n')
			 buffer[strlen(buffer)-1] = '\0';
		//sweep_set_traj_file(buffer);
	} else if (!strcmp(var, "pop_traj_file" )) {
		std::istringstream is( buffer );
		is.exceptions( std::ios_base::failbit | std::ios_base::badbit );
		filename_t sizeTrajFN;
		popid popname( NULL_POPID );
			
		try { is >> popname >> sizeTrajFN; }
		catch( const std::ios_base::failure& e ) {
			BOOST_THROW_EXCEPTION( cosi_param_file_error()
														 << error_msg( "invalid pop size traj line" )
														 << boost::errinfo_errno( errno ) );
		}
		try {
			math::Function< genid, popsize_float_t, math::Piecewise< math::Const<> > > sizeTraj;
			boost::filesystem::ifstream sizeTrajFile;
			try { sizeTrajFile.open( sizeTrajFN ); }
			catch( const std::ios_base::failure& e ) {
				BOOST_THROW_EXCEPTION( cosi_param_file_error()
															 << error_msg( "could not open pop size traj file" )
															 << boost::errinfo_errno( errno ) );
			}
			loadFrom( sizeTrajFile, sizeTraj );
			this->pop2sizeTraj.insert( std::make_pair( popname, sizeTraj ) );
		} catch( boost::exception& e ) {
			e << errinfo_file_name2( sizeTrajFN );
			throw;
		}
	} else {
		BOOST_THROW_EXCEPTION( cosi_param_file_error() << error_msg( "unknown param file line" ) );
	}

	return 1;
}

int 
ParamFileReader::file_get_data (FILE *fileptr, FILE *segfp) 
{
	char c,
		 buffer[BUF_MAX],
		 var[50];
	unsigned lineNum = 1;

	c = getc(fileptr);
	while (c != EOF) {
		switch (c) {
		case '#':
			c = getc(fileptr);
			while (c != '\n' && c != EOF)
				 c = getc(fileptr);
			lineNum++;
			break;
		case '\n':
			c = getc(fileptr);
			lineNum++;
			break;
		case ' ':
			c = getc(fileptr);
			break;
		default: 
			{
				ungetc(c, fileptr);
				int n_read = fscanf(fileptr, "%s", var);
				chkCond(n_read == 1, "error reading parameter file");
				file_killwhitespace(fileptr);
				char *fgets_result = fgets(buffer, BUF_MAX, fileptr);
				chkCond(fgets_result == buffer, "error reading parameter file");
				//strcat( buffer, " " );  // make sure no error if eof
			}
			try {
				try { file_proc_buff(var, buffer, segfp); }
				catch( const std::ios_base::failure& ) {
					BOOST_THROW_EXCEPTION( cosi_param_file_error()
																 << boost::errinfo_errno( errno ) );
				}
				catch( const std::exception& e ) {

					std::cerr << "genMap::readFrom - caught UNKNOWN exception of type " <<
						 typeid( e ).name() << " demangled " << 
						 ( boost::core::demangle( typeid( e ).name() ) ) << " ; exception is " << e.what() << "\n";
					
					BOOST_THROW_EXCEPTION( cosi_param_file_error()
																 << boost::errinfo_errno( errno ) );
				}
			} catch( boost::exception& e ) {
				e << boost::errinfo_at_line( lineNum )
					<< errinfo_param_name( var )
					<< errinfo_param_details( boost::trim_copy( std::string( buffer ) ) );
				throw;
			}
			c = getc(fileptr);
			lineNum++;
			break;
		}
	}

	// if (!seeded)
	// 	 std::cerr << "coalescent seed: " << ( rseed = seed_rng()) << endl;

	if (segfp != NULL) {fprintf(segfp, "params: length %d mu %.9f\n", length, (double)ToDouble( mu ) );}
	if (segfp != NULL) {fprintf(segfp, "L %d\n", length);}

	return 1;
}

int 
ParamFileReader::file_killwhitespace(FILE * fileptr) 
{
	char c;
	int i = 0;
	c = getc(fileptr);
	while (isspace((int) c)) {
		c = getc(fileptr);
		i++;
	}
	ungetc(c, fileptr);
	return i;
}

/*
 * FILE_EXIT
 */
void 
ParamFileReader::file_exit(const char* funct_name, const char* error_string) 
{
	fprintf(stderr, "file.c | %s: %s\n",
					funct_name, error_string);
	exit(EXIT_FAILURE);
}

void 
ParamFileReader::file_error_nonfatal(const char* funct_name, const char* error_string) 
{
	fprintf(stderr, "file.c | %s: %s\n",
					funct_name, error_string);
}

}  // namespace cosi
