/* $Id: recomb.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_GENMAP_H
#define __INCLUDE_COSI_GENMAP_H

#include <cstdlib>
#include <functional>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/join.hpp>
#include <cosi/defs.h>
#include <cosi/utils.h>
#include <cosi/decls.h>

namespace cosi {

template <typename TScale> struct scale_loc_type;
template <typename TScale> struct scale_len_type;
template <typename TScale> typename scale_loc_type<TScale>::type scale_min( const TScale& );
template <typename TScale> typename scale_loc_type<TScale>::type scale_max( const TScale& );
template <typename TScale> typename scale_len_type<TScale>::type scale_len( const TScale& s ) {
	return scale_max( s ) - scale_min( s );
}

template <typename TScaleIn, typename TScaleOut>
typename scale_loc_type<TScaleOut>::type cvt( TScaleIn& scaleIn, TScaleOut& scaleOut,
																							typename scale_loc_type<TScaleIn>::type x ) {
	return scale_min( scaleOut ) + ( ( x - scale_min( scaleIn ) ) / scale_len( scaleIn ) ) * scale_len( scaleOut );
}


// Class: GenMap
//
// Keeps the correspondence between physical locations and genetic map locations. 
class GenMap: boost::noncopyable {
	 // Types for representing physical distance (pd) and genetic distance (pd),
	 // as locations (loc) or lengths (len), and in the original untis (orig: bp for pd,
	 // cM for gd) or in normalized units (norm: fraction of the physical or genetic length
	 // of the whole simulated region).

	 typedef loc_bp_t  pd_orig_loc_t;
	 typedef len_bp_t  pd_orig_len_t;
	 typedef gloc_cM_t gd_orig_loc_t;
	 typedef glen_cM_t gd_orig_len_t;
	 typedef gloc_t gd_norm_loc_t;
	 typedef glen_t gd_norm_len_t;
	 typedef ploc_t pd_norm_loc_t;
	 typedef plen_t pd_norm_len_t;

public:
	 GenMap( const boost::filesystem::path& fname, len_bp_t length_ );
	 GenMap( istream&, len_bp_t length_ );

	 // Method: recomb_get_length
	 // Returns the physical length of the simulated region, in basepairs.
	 len_bp_t recomb_get_length () const {
		 return pd_range_len;
	 }

	 // MethodP: get_one_chrom_recomb_rate
	 // Return the probability, per generation, of recombination _somewhere_ on one chromosome
	 // (i.e. the total absolute genetic length of the simulated region).
	 glen_cM_t getRegionRecombRateAbs(void) const {
		 return gd_range_len;
	 }

	 // Func: recomb_get_gdPos_at
	 // Return the genetic map position corresponding to the given physical position.
	 gloc_t getGdPos( loc_t loc ) const {
		 bool is_nan = (boost::math::isnan)( ToDouble( loc.gdVal ) );
		 //PRINT5( loc, loc.gdVal, is_nan, this->loc2cumRate( get_ploc( loc ) ), this->loc2cumRate.getNumPts() );
		 return is_nan ? orig2norm_gd( pd2gd( norm2orig_pd( get_ploc( loc ) ) ) ) : loc.gdVal;
	 }

	 loc_t getLoc( gloc_t gdPos ) const {
		 return make_loc( get_ploc_from_gloc( gdPos ), gdPos );
	 }

	 loc_t add( const loc_t& loc, const plen_t& plen ) const {
		 ploc_t result_ploc = std::min( get_ploc( loc ) + plen, MAX_PLOC );
		 return make_loc( result_ploc, orig2norm_gd( pd2gd( norm2orig_pd( result_ploc ) ) ) );
	 }
	 loc_t sub( const loc_t& loc, const plen_t& plen ) const {
		 ploc_t result_ploc = std::max( get_ploc( loc ) - plen, MIN_PLOC );
		 return make_loc( result_ploc, orig2norm_gd( pd2gd( norm2orig_pd( result_ploc ) ) ) );
	 }

//	 void pickRandomRegion( RandGenP randGen_ );

	 void setStart( pd_orig_loc_t start );
	 
private:
	 
	 // Representation of the entire genetic map, as two parallel arrays of physical and genetic distance
	 // from the start of the map.
	 
	 std::vector< pd_orig_loc_t > pd_locs;
	 std::vector< gd_orig_loc_t > gd_locs;

	 // The sub-ranges of the genetic map, used for the currently simulated region.
	 std::vector< pd_orig_loc_t > pd_locs_range;
	 std::vector< gd_orig_loc_t > gd_locs_range;

	 // The length of the current simulated region, in cM
	 gd_orig_len_t gd_range_len;
	 // The length of the current simulated region, in bp
	 pd_orig_len_t pd_range_len;

	 gd_orig_loc_t pd2gd( pd_orig_loc_t pd_orig_loc ) const
			{ return util::interp( pd_locs_range, gd_locs_range, pd_orig_loc ); }
	 pd_orig_loc_t gd2pd( gd_orig_loc_t gd_orig_loc ) const
			{ return util::interp( gd_locs_range, pd_locs_range, gd_orig_loc ); }
	 
	 // math::InterpFn<ploc_t,gloc_t> loc2cumRate;
	 // math::InterpFn<gloc_t,ploc_t> cumRate2loc;

	 void readFrom( istream& recombfp );

	 ploc_t forceLocOntoBpBoundary( ploc_t loc ) const {
		 return loc;
		 //return  ploc_t( ( (double)( (int)( ToDouble( loc ) * ((double)recomb_get_length() ) ) ) ) / ((double)recomb_get_length()) );
	 }

	 gd_orig_loc_t norm2orig_gd( gd_norm_loc_t gdPos ) const
			{  return *boost::begin( gd_locs_range ) + ToDouble( gdPos ) * gd_range_len; }

	 pd_orig_loc_t norm2orig_pd( pd_norm_loc_t pdPos ) const
			{ return *boost::begin( pd_locs_range ) + ToDouble( pdPos ) * pd_range_len; }
	 gd_norm_loc_t orig2norm_gd( gd_orig_loc_t gdPos ) const
			{ return static_cast<gd_norm_loc_t>( ( gdPos - *boost::begin( gd_locs_range ) ) / gd_range_len ); }
	 pd_norm_loc_t orig2norm_pd( pd_orig_loc_t pdPos ) const
			{ return static_cast<pd_norm_loc_t>( ( pdPos - *boost::begin( pd_locs_range ) ) / pd_range_len ); }

	 ploc_t get_ploc_from_gloc( gloc_t gdPos ) const {
		 //return cumRate2loc( gdPos );
		 return forceLocOntoBpBoundary( orig2norm_pd( gd2pd( norm2orig_gd( gdPos )  ) ) );
	 }
	 
};  // class GenMap

}  // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_RECOMB_H
