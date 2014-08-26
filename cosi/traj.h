//
// Header: traj.h
//
// Representation of population size trajectories, and methods for
// generating trajectories under various conditions.
//

#ifndef __INCLUDE_COSI_TRAJ_H
#define __INCLUDE_COSI_TRAJ_H

#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <cosi/defs.h>
#include <cosi/utils.h>

namespace cosi {

using std::vector;
using std::pair;
using util::chkCond;

// Abstract class: FreqTraj
// Specifies the frequency trajectory of an allele in a set of pops.
// Keeps track of our current point in the trajectory, and allows stepping through
// the trajectory, once.
class FreqTraj {
public:
	 FreqTraj() { }
	 virtual ~FreqTraj();

	 // MethodP: getCurGen
	 // Returns the generation of our current point on the trajectory.
	 virtual genid getCurGen() const = 0;

	 // MethodP: getCurFreq
	 // Returns the frequency of the given pop at the current point on the trajectory.
	 virtual freq_t getCurFreq( popid popname ) const = 0;

	 // MethodP: next
	 // Moves to the next point on the trajectory.
	 virtual void next() = 0;

	 // MethodP: done
	 // Checks whether the trajectory is completed.
	 virtual bool_t done() const = 0;
};

//
// Class: TrajFromFile
//
// Single-population trajectory loaded from file.
//
class TrajFromFile: public FreqTraj {
public:
	 TrajFromFile( filename_t fname, genid firstGen_, freq_t firstFreq_ );
	 virtual ~TrajFromFile();

	 // MethodP: getCurGen
	 // Returns the generation of our current point on the trajectory.
	 virtual genid getCurGen() const { chkCond( !done(), "read past end of traj" ); return curTrajPt->first; }

	 // MethodP: getCurFreq
	 // Returns the frequency of the given pop at the current point on the trajectory.
	 virtual freq_t getCurFreq( popid /*popname*/ ) const { chkCond( !done(), "read past end of traj" ); return curTrajPt->second; }

	 // MethodP: next
	 // Moves to the next point on the trajectory.
	 virtual void next() {
		 chkCond( !done(), "reading past the end of trajectory" );
		 ++curTrajPt;
	 }

	 // MethodP: done
	 // Checks whether the trajectory is completed.
	 virtual bool_t done() const { return curTrajPt == traj.end(); }

private:
	 typedef std::map< genid, freq_t > traj_t;
	 typedef traj_t::const_iterator const_iterator;

	 // Field: traj
	 // The entire trajectory vector
	 traj_t traj;

	 // Field: curTrajPt
	 // The current traj point
	 const_iterator curTrajPt;
	 
};  // class TrajFromFile

//
// Class: DeterministicSweepTraj
// 
// Deterministic trajectory of a selective sweep in a single constant-size population,
// without migration to any other populations.  This corresponds to cosi's original
// sweep implementation.
//
struct DeterministicSweepTraj: public FreqTraj {
public:
	 DeterministicSweepTraj( popid sweepPop_, genid sweepEndGen_, gensInv_t selCoeff_, freq_t final_sel_freq_, nchroms_t final_pop_size, double deltaTfactor_ = 1.0 );
	 virtual ~DeterministicSweepTraj();

	 virtual genid getCurGen() const { return t; }
	 virtual freq_t getCurFreq( popid popname ) const {
		 return popname == sweepPop ? freq : 0.0;
	 }

	 virtual void next() {
		 t += deltaT;
		 freq = epsilon / (epsilon + (1 - epsilon) * exp(selCoeff*(t + end_shift - tend)));				
	 }
			
	 // Method: done
	 // Checks whether the trajectory is completed.
	 virtual bool_t done() const { return freq <= epsilon; }
	 
private:
	 // Field: sweepPop
	 // The sole population in which the sweep happens.
	 // The frequency of the derived allele in other populations always remains zero.
	 popid sweepPop;

	 // Field: selCoeff
	 // The selection coefficient.
	 gensInv_t selCoeff;

	 // Field: final_sel_freq
	 // The frequency of the derived allele at the end of the sweep.
	 freq_t final_sel_freq;

	 // Field: epsilon
	 // 1/2N.
	 freq_t epsilon;
	 
	 // Field: t
	 // The current generation of the trajectory.
	 genid t;
			
	 // Field: deltaT
	 // The step we take each time
	 gens_t deltaT;

	 // Field: tend
	 // Time when sweep ends?
	 genid tend;

	 // Field: end_shift
	 // The end shift
	 gens_t end_shift;

	 // Field: freq
	 // Current frequency of the allele
	 freq_t freq;

};  // class DeterministicSweepTraj




}  // namespace cosi

#endif // #ifndef __INCLUDE_COSI_TRAJ_H

