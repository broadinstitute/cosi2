/* $Id: defs.h,v 1.3 2011/05/31 20:32:29 sfs Exp $ */

/**
 * Header: defs.h
 *
 * Defines types and constants used throughout cosi, as well as conditional defines.
 */

#ifndef __COSI_INCLUDE_DEFS_H
#define __COSI_INCLUDE_DEFS_H

#include <cmath>
#include <limits>
#include <iostream>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/operators.hpp>
#include <cosi/typedval.h>

namespace cosi {

/*
 * Section: Debugging and profiling defines
 */

#define MUTATE_DEBUG 0
#define DEMOG_DEBUG 0
#define POP_DEBUG 0
#define DEBUG 0
#define FILE_DEBUG 0

/* #define DO_LOGGING */

/* #define COSI_STATS */

#ifdef COSI_STATS
#define IF_COSI_STATS(...) __VA_ARGS__
#else
#define IF_COSI_STATS(...)
#endif

/* #define NODE_DEBUG */

/* Type: bool_t */
/* A Boolean value.  */
typedef bool bool_t;

const bool_t True = true;
const bool_t False = false;

// Logical type: nchroms_t
// A count of chromosome instances: for example, the size of a population or a sample,
// or a subset of either.
typedef int nchroms_t;

typedef nchroms_t popsize_t;


// // Type: frac_t
// // A fraction, in the interval [0.0,1.0]
// struct frac_t: public factor_t {
// 	 frac_t() { }
// 	 explicit frac_t( double val_ ): factor_t( val_ ) { }
// };

// inline frac_t operator+( const frac_t& f1, const frac_t& f2 ) { return frac_t( f1.val + f2.val ); }

// /**
//  * Type: prob_t
//  *
//  * A probability/frequency/rate value.
//  */
// struct prob_t: public frac_t {
// 	 prob_t() { }
// 	 explicit prob_t( double val_ ): frac_t( val_ ) { }
// };

// // Type: freq_t
// // Represents the frequency of an allele, in the [0.0,1.0] range.
// struct freq_t: public frac_t {
// 	 freq_t() { }
// 	 explicit freq_t( double val_ ): frac_t( val_ ) { }
// };

typedef double frac_t;
typedef double freq_t;
typedef double prob_t;

COSI_DEFINE_TYPEDVAL_REL(prob_per_gen_t);
COSI_DEFINE_TYPEDVAL_REL(prob_per_chrom_per_gen_t);

COSI_DEFINE_TYPEDVAL_MULT(nchroms_t, prob_per_chrom_per_gen_t, prob_per_gen_t);

/**
 * Type: len_bp_t
 *
 * Length of a <region>, measured in basepairs (possibly fractional).
 */
typedef cosi_double len_bp_t;

COSI_DEFINE_TYPEDVAL_REL(prob_per_bp_per_gen_t);
COSI_DEFINE_TYPEDVAL_REL(prob_per_len_per_gen_t);

COSI_DEFINE_TYPEDVAL_MULT(prob_per_bp_per_gen_t, len_bp_t, prob_per_len_per_gen_t);

// Type: idx_t
// An index into an array.
typedef int idx_t;

/**
 * Logical Type: gens_t
 *
 * Length of time, in generations.
 */
COSI_DEFINE_TYPEDVAL_ABSREL(genid, gens_t, cosi_double);

/**
 * Logical type: genid
 *
 * A particular generation.  Since the coalescent is a continuous approximation of the Wright-Fischer model, generations are
 * floating-point rather than integer.
 */

// Const: NULL_GEN
// A deliberately invalid value for a generation
const genid NULL_GEN(NAN);

const genid ZERO_GEN(0.0);

const gens_t ZERO_GENS(0.0);
const gens_t NULL_GENS(NAN);

inline bool_t is_null( const gens_t& gens ) { return (boost::math::isnan)( ToDouble( gens ) ); }
inline bool_t is_null( const genid& gen ) { return (boost::math::isnan)( ToDouble( gen ) ); }

/**
 * Type: nodeid
 * 
 * Identifier of one ARG <Node>.
 */
typedef int nodeid;

/**
 * Type: popid
 *
 * Identifier of one population.   This is the numeric id assigned to the population in the <config file>.
 */
COSI_DEFINE_TYPEDVAL_ID( popid );

typedef int pop_idx_t;

const pop_idx_t NULL_POP_IDX( -1 );

// Const: NULL_POPID

// An invalid sentinel value for population name.

const popid NULL_POPID(-1);

inline bool_t is_null( const popid& p ) { return p == NULL_POPID; }

/****************************/

//
// Section: Representation of locations
//
// Here we describe the types used for representing locations in the simulated region,
// and distances between these locations.  Each point in the region has a physical
// location (represented by ploc_t) and a genetic map location (represented by gloc_t);
// there is a bijection between physical locations and genetic locations.  (So, a segment
// of non-zero physical distance must have non-zero genetic distance.)
// Both physical and genetic locations are represented by a real value in the interval [0.0,1.0],
// representing the fraction of the region's
// total physical or genetic length at that location.  Physical distance between two points
// is represented by plen_t; genetic distance, by glen_t.
//
// The type used throughout most of cosi to represent locations is loc_t.  This type includes
// the physical location and the genetic location.

/////////////////////////////////////////////

// Struct: ploc_t
// A physical location in the simulated region.
// Represented as a fraction of the total physical length of the region.

// Logical type: plen_t
// Physical distance between two points, represented as a fraction
// of the total physical length of the region.

COSI_DEFINE_TYPEDVAL_ABSREL(ploc_t, plen_t, cosi_double);

/* Const: MIN_PLOC */
/* Leftmost possible physical location. */
const ploc_t MIN_PLOC(0.0);

/* Const: MAX_PLOC */
/* Rightmost possible physical location. */
const ploc_t MAX_PLOC(1.0);

// Const: NULL_PLOC
// A deliberately invalid value for a location.
const ploc_t NULL_PLOC(std::numeric_limits<cosi_double>::quiet_NaN());

const plen_t ZERO_PLEN(0.0);

///////////////////////////

// Logical type: gloc_cM_t
// Location on the genetic map, in centimorgans.  Throughout most of cosi, genetic map locations are represented not
// by this type but by gloc_t, which is the fraction at this location of the region's total genetic distance.
//typedef double gloc_cM_t;

// Logical type: glen_cM_t
// Genetic distance between two locations, in centimorgans.
COSI_DEFINE_TYPEDVAL_ABSREL(gloc_cM_t,glen_cM_t,cosi_double);

const gloc_cM_t ZERO_GLOC_CM(0.0);

COSI_DEFINE_TYPEDVAL_ABSREL(gloc_t,glen_t,cosi_double);

/* Leftmost possible genetic map location. */
const gloc_t MIN_GLOC(0.0);

/* Const: MAX_GLOC */
/* Rightmost possible genetic map location. */
const gloc_t MAX_GLOC(1.0);

// Const: NULL_GLOC
// A deliberately invalid value for a location.
const gloc_t NULL_GLOC(std::numeric_limits<cosi_double>::quiet_NaN());

inline glen_cM_t operator*( const glen_t& glen, const glen_cM_t& glenAbs ) { return glen_cM_t( glen.val * glenAbs.val ); }
inline glen_cM_t operator*( const glen_cM_t& glenAbs, const glen_t& glen ) { return glen * glenAbs; }

const glen_t ZERO_GLEN(0.0);
const glen_t NULL_GLEN(std::numeric_limits<cosi_double>::quiet_NaN());

inline glen_t operator*( const frac_t& frac, const glen_t& glen ) { return glen_t( frac * glen.val ); }


#ifdef COSI_LONG_DOUBLE
#define LOC_FMT "%Lf"
#else
#define LOC_FMT "%lf"
#endif

inline bool_t is_null( const ploc_t& ploc ) { return (boost::math::isnan)(ploc.val); }
inline bool_t is_null( const gloc_t& gloc ) { return (boost::math::isnan)(gloc.val); }

// Logical type: len_t
// Represents the length of a particular stretch of the simulated region.
// Stores the genetic length as well as the
// physical length.

struct len_t: public plen_t {
  glen_t gdVal;
  len_t() {}
  explicit len_t( cosi_double val_ ): plen_t( val_ ), gdVal( NULL_GLEN ) {}
	explicit len_t( const plen_t& plen_ ): plen_t( plen_ ), gdVal( NULL_GLEN ) {}
	explicit len_t( const plen_t& plen_, const glen_t& gdVal_ ): plen_t( plen_ ), gdVal( gdVal_ ) {}
	len_t& operator-=( const len_t& len ) { plen_t::operator-=( len ); gdVal -= len.gdVal; return *this; }
	len_t& operator+=( const len_t& len ) { plen_t::operator+=( len ); gdVal += len.gdVal; return *this; }

  const plen_t& as_plen() const { return *this; } 	 
};

struct loc_t: public ploc_t {
   gloc_t gdVal;

	 loc_t() { }
   explicit loc_t( cosi_double val_ ): ploc_t( val_ ), gdVal(NULL_GLOC) { }
   explicit loc_t( const ploc_t& val_ ): ploc_t( val_ ), gdVal(NULL_GLOC) { }
   loc_t( const ploc_t& val_, const gloc_t& gdVal_ ): ploc_t( val_ ), gdVal( gdVal_ ) {}
   const ploc_t& as_ploc() const { return *this; } 	 
};

inline len_t operator+( const len_t& len1, const len_t& len2 ) { return len_t( len1.as_plen() + len2.as_plen(),
len1.gdVal + len2.gdVal ); }
inline len_t operator-( const len_t& len1, const len_t& len2 ) { return len_t( len1.as_plen() - len2.as_plen(), len1.gdVal - len2.gdVal ); }

inline len_t operator-( const loc_t& loc1, const loc_t& loc2 ) { return len_t( loc1.as_ploc() - loc2.as_ploc(), loc1.gdVal - loc2.gdVal); }

inline std::ostream& operator<<( std::ostream& s, const len_t& len ) { s << len.val << ":" << len.gdVal ; return s; }
inline std::istream& operator>>( std::istream& s, len_t& len ) { s >> len.val; len.gdVal = NULL_GLEN; return s; }

const len_t ZERO_LEN(ZERO_PLEN, ZERO_GLEN);
const loc_t MIN_LOC( MIN_PLOC, MIN_GLOC ), MAX_LOC( MAX_PLOC, MAX_GLOC );

inline std::ostream& operator<<( std::ostream& s, const loc_t& loc ) { s << loc.val << ":" << loc.gdVal; return s; }
inline std::istream& operator>>( std::istream& s, loc_t& loc ) { s >> loc.val; loc.gdVal = NULL_GLOC; return s; }

inline loc_t make_loc( const ploc_t& ploc_, const gloc_t& gloc_ ) { return loc_t( ploc_, gloc_ ); }

const loc_t NULL_LOC(NULL_PLOC);

inline cosi_double get_loc( const loc_t& loc ) { return ToDouble( loc.val ); } 

inline cosi_double get_len( const len_t& len ) { return ToDouble( len.val ); }
inline cosi_double get_phys_len( const len_t& len ) { return ToDouble( len.val ); }

// Func: get_glen
// Return the genetic distance, from either a <glen_t> or a <len_t>.
inline glen_t get_glen( const glen_t& glen ) { return glen; }
inline glen_t get_glen( const len_t& len ) { return len.gdVal; }

// Func: get_plen
// Return the physical distance, from either a <plen_t> or a <len_t>.
inline plen_t get_plen( const plen_t& plen ) { return plen; }
inline plen_t get_plen( const len_t& len ) { return len; }

// Func: get_gloc
// Return the genetic distance, from either a <gloc_t> or a <loc_t>.
inline gloc_t get_gloc( const gloc_t& gloc ) { return gloc; }
inline gloc_t get_gloc( const loc_t& loc ) { return loc.gdVal; }

// Func: get_ploc
// Return the physical distance, from either a <ploc_t> or a <loc_t>.
inline ploc_t get_ploc( const ploc_t& ploc ) { return ploc; }
inline ploc_t get_ploc( const loc_t& loc ) { return loc; }


/////////////////////////////////////////////

/**
 * Type: loc_bp_t
 *
 * A point location on the chromosome, represented as a base-pair coordinate.
 */
typedef cosi_double loc_bp_t;

/**
 * Type: loc_bp_int_t
 *
 * A point location on the chromosome, represented as a base-pair integer coordinate.
 */
typedef int loc_bp_int_t;

/**
 * Type: len_t
 *
 * Length of a <region>, measured as a fraction of the total chromosome length; range [0-1].
 */
//typedef loc_t len_t;





/**
 * Type: len_bp_int_t
 *
 * Length of a <region>, measured in basepairs (integer).
 */
typedef int len_bp_int_t;


/* Enum: Direction */
/* Direction from a given point, e.g. to identify  */
/*   */
/*    DIR_L - leftward direction */
/*    DIR_R - rightward direction. */
enum dir_t { DIR_L, DIR_R };

/* Enum: Gene conversion part */
/* Identifies one of the two chromosome parts resulting from a gene conversion. */
/*     GC_INNER - the inner part resulting from gene conversion */
/*     GC_INNER - the outer part resulting from gene conversion */
enum gc_kind_t { GC_INNER, GC_OUTER };

}  // namespace cosi
  
#endif /* DEFS_H */

