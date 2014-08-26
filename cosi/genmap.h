/* $Id: recomb.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_GENMAP_H
#define __INCLUDE_COSI_GENMAP_H

#include <cstdlib>
#include <functional>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/generalmath.h>

namespace cosi {

// Class: GenMap
//
// Keeps the correspondence between physical locations and genetic map locations. 
class GenMap: boost::noncopyable {
public:
	 GenMap( const boost::filesystem::path& fname, len_bp_t length_, ploc_bp_diff_t genMapShift_ = 0 );
	 GenMap( istream&, len_bp_t length_, ploc_bp_diff_t genMapShift_ = 0 );

	 // Method: recomb_get_length
	 // Returns the physical length of the simulated region, in basepairs.
	 len_bp_t recomb_get_length () const {
		 return rec_length;
	 }

	 // MethodP: get_one_chrom_recomb_rate
	 // Return the probability, per generation, of recombination _somewhere_ on one chromosome
	 // (i.e. the total absolute genetic length of the simulated region).
	 glen_cM_t getRegionRecombRateAbs(void) const {
		 return rec_recombrate;
	 }

	 // Func: recomb_get_gdPos_at
	 // Return the genetic map position corresponding to the given physical position.
	 gloc_t getGdPos( loc_t loc ) const {
		 bool is_nan = (boost::math::isnan)( ToDouble( loc.gdVal ) );
		 //PRINT5( loc, loc.gdVal, is_nan, this->loc2cumRate( get_ploc( loc ) ), this->loc2cumRate.getNumPts() );
		 return is_nan ? this->loc2cumRate( get_ploc( loc ) ) : loc.gdVal;
	 }

	 loc_t getLoc( gloc_t gdPos ) const {
		 return make_loc( get_ploc_from_gloc( gdPos ), gdPos );
	 }

	 loc_t add( const loc_t& loc, const plen_t& plen ) const {
		 ploc_t result_ploc = std::min( get_ploc( loc ) + plen, MAX_PLOC );
		 return make_loc( result_ploc, loc2cumRate( result_ploc ) );
	 }
	 loc_t sub( const loc_t& loc, const plen_t& plen ) const {
		 ploc_t result_ploc = std::max( get_ploc( loc ) - plen, MIN_PLOC );
		 return make_loc( result_ploc, loc2cumRate( result_ploc ) );
	 }
																																					 
	 
private:
	 
	 // Field: rec_recombrate
	 // The total per-chromosome recombination rate per generation.
	 glen_cM_t rec_recombrate;

	 // Field: rec_length
	 // Length in basepairs of the simulated region.
	 const len_bp_t rec_length;

	 math::InterpFn<ploc_t,gloc_t> loc2cumRate;
	 math::InterpFn<gloc_t,ploc_t> cumRate2loc;

	 void readFrom( istream& recombfp, ploc_bp_diff_t genMapShift_ );

	 ploc_t forceLocOntoBpBoundary( ploc_t loc ) const {
		 return loc;
		 //return  ploc_t( ( (double)( (int)( ToDouble( loc ) * ((double)recomb_get_length() ) ) ) ) / ((double)recomb_get_length()) );
	 }

	 ploc_t get_ploc_from_gloc( gloc_t gdPos ) const {
		 //return cumRate2loc( gdPos );
		 return forceLocOntoBpBoundary( cumRate2loc( gdPos ) );
	 }
	 
};  // class GenMap

}  // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_RECOMB_H
