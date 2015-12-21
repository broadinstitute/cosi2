#include <ios>
#include <boost/filesystem/fstream.hpp>
#include <cosi/general/utils.h>
#include <cosi/traj.h>

namespace cosi {

using util::chkCond;
	
FreqTraj::~FreqTraj() { }

//
// Class impl: TrajFromFile
//

TrajFromFile::TrajFromFile( filename_t fname, genid firstGen_, freq_t firstFreq_ )
{
	using std::ios;
	boost::filesystem::ifstream trajStrm( fname );
	trajStrm.exceptions( ios::failbit | ios::badbit | ios::eofbit );
	traj.insert( make_pair( firstGen_, firstFreq_ ) );
	while ( True ) {
		try {
			genid a_gen;
			freq_t a_freq;
			trajStrm >> a_gen >> a_freq;
			traj.insert( make_pair( a_gen, a_freq ) );
		} catch( std::ifstream::failure e ) {
			break;
		}
	}
	if ( traj.rbegin()->second > static_cast<freq_t>(0.0) )
		traj.insert( make_pair( traj.rbegin()->first + gens_t( 1.0 ), freq_t( 0.0 ) ) );
	curTrajPt = traj.begin();
}

TrajFromFile::~TrajFromFile() { }

// End class impl: TrajFromFile


////////////////

DeterministicSweepTraj::DeterministicSweepTraj( popid sweepPop_, genid sweepEndGen_,
																								gensInv_t selCoeff_, freq_t final_sel_freq_, nchroms_t final_pop_size_, double deltaTfactor_  ):
	sweepPop( sweepPop_ ), selCoeff( selCoeff_ ), final_sel_freq( final_sel_freq_ ), epsilon( 1 / (2. * final_pop_size_) ) {
	
	end_shift =  (final_sel_freq > 1-epsilon) ? ZERO_GENS :
		 ( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / selCoeff );
	deltaT = gens_t( .01 / selCoeff ) * deltaTfactor_;
	tend = sweepEndGen_ - ( 2 * log(epsilon) / selCoeff );
	t = sweepEndGen_;
	freq = final_sel_freq;
}
DeterministicSweepTraj::~DeterministicSweepTraj() { }




}  // namespace cosi


