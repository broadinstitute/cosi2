//
// Header: cosirand.h
//
// Defines a uniform interface for random number generation used by cosi;
// to change cosi to another random number generator, just change this interface.
//

#ifndef __INCLUDE_COSIRAND_H
#define __INCLUDE_COSIRAND_H

#include <cstdio>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/integer/integer_mask.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <cosi/general/utildefs.h>
#include <cosi/general/typedval.h>

namespace cosi {


//
// Class: RandGen
//
// A random number generator.
//
class RandGen {
	 
public:
	 typedef boost::mt19937 rand_engine_t;
	 
	 // Type: rseed_t
	 // A random seed.
	 typedef uint64_t rseed_t;

	 RandGen() {
		 boost::random_device rd;
		 setSeed( rd() );
	 }
	 //RandGen( unsigned long seed_ ) { setSeed( static_cast<uint32_t>( seed_ & 0xFFFFFFFF ) ); }
	 RandGen( rseed_t seed_ ) { setSeed( seed_ ); }

	 void setSeed( rseed_t seed_ ) {
		 std::cerr << "orig seed=" << seed_ << "\n";
		 // seed_ = seed_ & 0xFFFFFFFF;
		 // uint32_t 
		 randEngine.seed( static_cast<boost::uint32_t>( seed_ & boost::low_bits_mask_t<32>::sig_bits ) ); u01.reset();
		 rseed = seed_;
		 std::cerr << "got seed" << seed_ << "\n";
	 }

	 /* Func: random_bit */
	 /* Generate one random bit from <random_state_t>. */
	 bool_t random_bit_batched() {

		 if ( !bits_left_in_batch ) {
			 random_bits_batch = randEngine();
			 bits_left_in_batch = 31;
		 } else {
			 random_bits_batch >>= 1;
			 bits_left_in_batch--;
		 }
		 return ( random_bits_batch & 1 ) != 0;
	 }

	 // Method: seed_rng
	 // Initialize this random number generator using the system time or other random source.
	 //rseed_t seed_rng () { return rseed = mts_seed( &state ); }

	 // Method: set_rng_seed
	 // Initialize this random number generator using the specified seed.
	 //void set_rng_seed(rseed_t rseed_) { mts_seed32new( &state, rseed_ ); rseed = rseed_; }

	 // Method: getSeed
	 // Returns the seed last used to initialize this random number generator.
	 rseed_t getSeed() const { return rseed; }
	 
	 // Method: random_double
	 // Return a uniformly distributed random value in the [0,1) range.
	 frac_t random_double() {
		 return u01( *this );
	 }

	 //rand_engine_t& getEngine() { return randEngine; }
	 
	 // Method: random_idx
	 // Return a random unsigned integer uniformly distributed in [0,N) range.
	 size_t random_idx( size_t N ) { return static_cast<size_t>( random_double() * ((double)N) ); }

	 // Method: random_bit
	 // Return a random bit from flipping a fair coin.
	 bool_t random_bit() { return random_double() < .5; } 

	 // Method: expdev
	 // Return an exponentially distributed waiting time.
	 factor_t expdev (void);

	 // void printState( filename_t fname ) const {
	 //   FILE *f = cosi_fopen( fname, "wt" );
	 //   mts_savestate( f, const_cast<mt_state *>(&this->state) );
	 //   fclose( f );
	 // }
	 // void loadState( filename_t fname ) {
	 //   FILE *f = cosi_fopen( fname.c_str(), "rt" );
	 // 	 mts_loadstate( f, &state );
	 // 	 fclose( f );
	 // }

	 // Method: ranbinom
	 // Return a sample from the binomial distribution.
	 //
	 // Params:
	 //
	 //    n - the total number of coin flips
	 //    p - the probability of success in each flip
	 //
	 // Returns:
	 //
	 //    the number of succeses
	 //
	 int ranbinom(int n, prob_t p) {
		 boost::binomial_distribution<> d( n, p );
		 return d( *this );
	 }

	 double poisson_get_next (double rate);

	 // Group: Uniform Random Number Generator members
	 // Members required for RandGen object to satisfy the concept requirements of
	 // Uniform Random Number Generator defined as
	 // http://www.boost.org/doc/libs/1_52_0/doc/html/boost_random/reference.html#boost_random.reference.concepts

	 typedef rand_engine_t::result_type result_type;
	 result_type min() const { return randEngine.min(); }
	 result_type max() const { return randEngine.max(); }
	 result_type operator()() {
		 result_type r = randEngine();
		 //std::cerr << "[" << r << "]\n";
		 return r;
	 }

	 // End group: Uniform Random Number Generator members
	 
private:
	 // Private field: state
	 // State of the random number generator.
	 rand_engine_t randEngine;

	 boost::uniform_01<> u01;

	 // Private field: seed
	 // The seed used to initialize the generator
	 rseed_t rseed;

	 // Field: random_bits_batch
	 // The last batch of random bits fetched; see random_bit_batched()
	 rand_engine_t::result_type random_bits_batch;

	 // Field: bits_left_in_batch
	 // Number of random bits still left in <random_bits_batch>.
	 unsigned short bits_left_in_batch;
#if 0
	 //
	 // Group: Auxiliary methods for ranbinom.
	 //
	 int ranb1(int, double);
	 double ranbeta(double, double);
	 double rangam(double);
	 double randev1(double);
	 double randev0(double);
	 double ranexp(void);

	 // Group: Other aux functions
	 int poisson( double xm );
#endif	 
	 
};  // class RandGen

typedef boost::shared_ptr<RandGen> RandGenP;

//
// Class: HasRandGen
//
// A class that includes a pointer to a random number generator.
//
class HasRandGen {
public:
	 HasRandGen() {}
	 HasRandGen( RandGenP randGen_ ): randGen( randGen_ ) {}
	 void setRandGen( RandGenP randGen_ ) { randGen = randGen_; }
	 RandGenP getRandGen() const { return randGen; }

	 frac_t random_double() const { return randGen->random_double(); }
	 bool_t random_bit() const { return randGen->random_bit(); }
	 factor_t expdev() const { return randGen->expdev(); }
	 int ranbinom(int n, prob_t p) const { return randGen->ranbinom( n, p ); }

	 double poisson_get_next (double rate) { return randGen->poisson_get_next( rate ); }
	 
private:
	 RandGenP randGen;
	 
};  // class HasRandGen

//bool_t random_bit(void);  

}  // namespace cosi

#endif  // #ifndef __INCLUDE_COSIRAND_H
