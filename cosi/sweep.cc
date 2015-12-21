/* $Id: sweep.c,v 1.5 2011/06/07 15:29:25 sfs Exp $ */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <ios>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/filesystem/fstream.hpp>
#include <cosi/general/utils.h>
#include <cosi/decls.h>
#include <cosi/node.h>
#include <cosi/pop.h>
#include <cosi/demography.h>
#include <cosi/recomb.h>
#include <cosi/geneconversion.h>
#include <cosi/mutlist.h>
#include <cosi/mutate.h>
#include <cosi/migrate.h>
#include <cosi/sweep.h>
#include <cosi/historical.h>
#include <cosi/genmap.h>
#include <cosi/leafset.h>
#include <cosi/seglist.h>
#include <cosi/traj.h>
#include <cosi/hooks.h>
#include <cosi/coalrate.h>

namespace cosi {

#define ForEach BOOST_FOREACH  
  
using std::vector;
using std::map;
using std::set;
using std::ifstream;
using std::ofstream;
using util::map_get;
using util::STLContains;
using node::Node;
using node::NodeList;

//
// Class impl: Sweep
//


Sweep::Sweep( DemographyP demography_ ):
	demography( demography_ ),
	finalFreq( std::numeric_limits<freq_t>::quiet_NaN() ),
	selLoc( NULL_LOC ),
	fix_freq( False ), deltaTfactor( 1.0 ), verbose( False ),
	coalesce_rate_nonsel(0.0), migrate_rate_nonsel( 0.0 ), recombination_rate_nonsel( 0.0 ),
	geneconv_rate_nonsel( 0.0 ), poisson_rate_nonsel( 0.0 )  {
}

void Sweep::sweep_setGenMap( GenMapP genMap_ ) { genMap = genMap_; }
void Sweep::sweep_setRecomb( RecombP recomb_ ) { recomb = recomb_; }
void Sweep::sweep_setGeneConversion( GeneConversionP geneConversion_ ) { geneConversion = geneConversion_; }
void Sweep::sweep_setMutate( MutateP mutate_ ) { mutate = mutate_; }
void Sweep::sweep_setNodePool( node::NodePoolP nodePool_ ) { nodePool = nodePool_; }

filename_t theTrajFN;

// Method: set_trajFN
// Set the name of the file from which to load the frequency trajectory of the selected allele
// in the selected pop, for use during sweep simulation.
void Sweep::setTrajFN( filename_t trajFN_ ) {
  trajFN = trajFN_;
	theTrajFN = trajFN_;
}


// Method: sweep_load_traj
// Load the causal allele trajectory from <trajFN> into <selAlleleFreqTraj>, if
// trajFN was given.
void Sweep::sweep_load_traj() {
  if ( !trajFN.empty() ) {
    boost::filesystem::ifstream f( trajFN );
		f.exceptions( ofstream::failbit | ofstream::badbit );
		
		while ( f ) {
			genid gen;
			freq_t freq;
			try { f >> gen >> freq; }
			catch( ios::failure fe ) { break; }
			//PRINT2( gen, freq );
			selAlleleFreqTraj.addPt( gen, freq );
		}
  }
}  // Sweep::sweep_load_traj()


void Sweep::sweep_set_fix_freq( bool_t fix_freq_ ) {
  fix_freq = fix_freq_;
}

unsigned long n_switch_recomb = 0, n_tot_recomb = 0;
  
//
// Method: sweep_execute
//
// Execute the selected sweep on one selected allele in one population,
// with a deterministically chosen selected allele trajectory.
//
// Params:
//
//   popname - the population in which the sweep happens
//   selCoeff - the selection coefficient.  We assume the heterozygous effect is .5.
//   gen - the generation when the sweep ends (moving forward) or starts (moving pastward).
//   sel_pos - the location of the selected allele within the simulated region
//   final_sel_freq - the frequency to which the selected allele rises by the end of the sweep
//      (moving forward), i.e. at time 'gen'.
//
genid Sweep::sweep_execute(popid popname, gensInv_t sel_coeff, genid gen, loc_t sel_pos, 
													 freq_t final_sel_freq) {
  Pop *popptr;
  int inode, nsel, nunsel;
  int newn = -1, foundit;
  Node **sel_nodes, **unsel_nodes, *rec_nodes[3], *gc_nodes[3];
  Node *old_allele_node, *new_allele_node;
	gens_t deltaT;
	genid t;
	genid tend;
  freq_t sel_freq;
	freq_t epsilon;
	double x, /*alpha,*/ rate;
  prob_t prob_coal_unsel, prob_coal_sel, prob_recomb, all_prob;
  genid tsave, t_nonsweep;
	gens_t poisson_event_time;
	gens_t end_shift=ZERO_GENS;
  loc_t loc;
  prob_t prob_gc;
  loc_t loc1, loc2;
  int selsize, unselsize, dorandom=1;

  /*  sw_trajectory(popname, sel_coeff); */
	finalFreq = final_sel_freq;
	selLoc = sel_pos;

  popptr = demography->dg_get_pop_by_name(popname);
  nunsel = nsel = 0;
  selsize = unselsize = popptr->pop_get_num_nodes();
  sel_nodes = (Node **)malloc(selsize * sizeof(Node*));
  unsel_nodes = (Node **)malloc(unselsize * sizeof(Node*));
  assert(sel_nodes != NULL && unsel_nodes != NULL);

  epsilon = 1 / (2. * popptr->pop_get_size());
  if (final_sel_freq > 1-epsilon) {
    end_shift = ZERO_GENS;
    dorandom = 0;
  }
  else {
    end_shift = ( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / sel_coeff );
  }
	//PRINT( deltaTfactor );
  deltaT = gens_t( ( .01 / sel_coeff ) * ToDouble( deltaTfactor ) ) ;
  tend = gen - ( 2 * log(epsilon) / sel_coeff );


	//PRINT4( deltaT, sel_coeff, deltaTfactor, tend );
	
  //alpha = 2 * popptr->pop_get_size() * sel_coeff;

	leafset_p sel_leaves = make_empty_leafset();
  for (inode = 0; inode < (int)popptr->pop_get_num_nodes(); inode++) {
    if (!dorandom || ( fix_freq && inode < ((int)( final_sel_freq * ((double)popptr->pop_get_num_nodes())) ) ) ||
				( !fix_freq && random_double() < final_sel_freq) ) {
      sel_nodes[nsel] = popptr->pop_get_node(inode);

			ForEach( const seglist::Seg& seg, *sel_nodes[nsel]->getSegs() )
				 if ( seg.getBeg() <= sel_pos && sel_pos <= seg.getEnd() )
						sel_leaves = leafset_union( sel_leaves, seg.getLeafset() );
			
      nsel++;
    }
    else {
      unsel_nodes[nunsel] = popptr->pop_get_node(inode);
      nunsel++;
    }
  }
  sel_freq = 1;
  tsave = gen;
  
  rate = sw_get_poisson_rate_nonsel(popname); 
  poisson_event_time = gens_t( 1e50 );
  if (rate > 0) {
    poisson_event_time = gens_t( poisson_get_next(rate) );
  }
  t_nonsweep = gen + poisson_event_time;

	bool_t trajSpecified = False;
	if ( !trajFN.empty() ) {
		sweep_load_traj();
		trajSpecified = True;
	}

#ifdef COSI_DEV_SWEEP_REDUCE_DELTA		
	const freq_t maxFreqDiff = 0.01;
	const prob_t maxEventProb = 0.08;
#endif	

	boost::filesystem::ofstream trajOut;
	trajOut.exceptions( ofstream::failbit | ofstream::badbit );
	if ( !trajOutFN.empty() ) trajOut.open( trajOutFN, ifstream::out );

	static int sweepCount=0;
	sweepCount++;
	//PRINT( sweepCount );

#if 0
	prob_t probCoalSelTot = 0, probCoalUnselTot = 0;
	typedef math::Function< genid, popsize_float_t, math::SweepPopSizeTraj > sweepTraj_t;
	typedef boost::shared_ptr<sweepTraj_t> sweepTraj_p;
	sweepTraj_p sweepTraj;
	sweepTraj = boost::make_shared<sweepTraj_t>(
									popsize_float_t( final_pop_size ),
									selCoeff,
									epsilon,
									end_shift,
									tend
		);

//	genid tstart = 
	

	BOOST_AUTO( sweepPopSizeComplementTraj,
							(
					math::Function< genid, popsize_float_t, math::SweepPopSizeComplementTraj >(
						popsize_float_t( final_pop_size ),
						selCoeff,
						epsilon,
						end_shift,
						tend
						)
								));

	BOOST_AUTO( coalRate, math::coalRateFunction( *sweepTraj ) );
	BOOST_AUTO( coalRateComplement, math::coalRateFunction( sweepPopSizeComplementTraj ) );

	// PRINT( coalRate( genid( 370.0 ) ) ); 
	// PRINT( coalRateComplement( genid( 370.0 ) ) ); 

	BOOST_AUTO( coalRateIntegral, math::integralFromPoint( coalRate, t ) );
	BOOST_AUTO( coalRateComplementIntegral, math::integralFromPoint( coalRateComplement, t ) );
#endif

	unsigned long nsteps = 0;
  for (t = gen; sel_freq >= epsilon || nsel > 1; t += deltaT) {
		nsteps++;
    while (t > t_nonsweep) {
      /* Process event in a nonsweep population, and update time for next nonsweep event */
      sw_do_poisson_nonsel(popname, t_nonsweep); 
      poisson_event_time = gens_t( poisson_get_next(sw_get_poisson_rate_nonsel(popname)) );
      t_nonsweep += poisson_event_time;  
    }

    sel_freq = trajSpecified ? selAlleleFreqTraj( t ) :
			 epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*(t + end_shift - tend)));
		//if ( sweepCount >= 264 ) PRINT4( t, sel_freq, nsel, nunsel );
		if ( !trajOutFN.empty() ) trajOut << t << "\t" << sel_freq << "\n";

#ifdef COSI_DEV_SWEEP_REDUCE_DELTA		
		while ( True ) {
			freq_t next_sel_freq = trajSpecified ? selAlleleFreqTraj( t + deltaT ) :
				 epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*ToDouble(t + deltaT + end_shift - tend)));
			if ( fabs( next_sel_freq - sel_freq ) < maxFreqDiff ) break;
			deltaT /= factor_t(2);
			//PRINT2( "deltaT reduced for freq, now ", deltaT );
		}

		while ( True ) {
#endif		

			prob_coal_sel = ToDouble( deltaT ) * nsel * (nsel-1) / 4 / popptr->pop_get_size() / sel_freq;
			prob_coal_unsel = ToDouble( deltaT ) * nunsel * (nunsel-1) / 4 / popptr->pop_get_size() / (1 - sel_freq);
#if 0
			BOOST_AUTO( curCoalRate, 1.0 / 2 / popptr->pop_get_size() / sel_freq  );
			BOOST_AUTO( curCoalRateComplement, 1.0 / 2 / popptr->pop_get_size() / (1-sel_freq)  );
			assert( curCoalRate == math::eval( coalRate, t ) );
			assert( curCoalRateComplement == math::eval( coalRateComplement, t ) );

			probCoalSelTot += ToDouble( deltaT )  / 2 / popptr->pop_get_size() / sel_freq;
			probCoalUnselTot += ToDouble( deltaT )  / 2 / popptr->pop_get_size() / ( 1 - sel_freq );

			PRINT5( t, probCoalSelTot, coalRateIntegral(t), probCoalUnselTot, coalRateComplementIntegral(t) );
			// PRINT2( t, );
			// PRINT2( t, 1.0 / 2 / popptr->pop_get_size() / ( 1 - sel_freq ) );
#endif
			
//		PRINT9( t, sel_freq, nsel, nsel_pop, nunsel, nunsel_pop, popptr->pop_get_size(), prob_coal_sel, prob_coal_unsel );
			prob_recomb = ToDouble( deltaT ) * (nsel + nunsel) * ToDouble( genMap->getRegionRecombRateAbs() );
			prob_gc = ToDouble( deltaT ) * (nsel + nunsel) * ToDouble( geneConversion->getRegionGeneConvRate() );
			all_prob = 0; 
			all_prob += prob_coal_unsel;
			all_prob += prob_coal_sel;
			all_prob += prob_recomb;
			all_prob += prob_gc;
			
#ifdef COSI_DEV_SWEEP_REDUCE_DELTA
			if ( all_prob < maxEventProb ) break;
			deltaT /= factor_t(2);
			//PRINT2( "deltaT reduced for prob, now ", deltaT );
		}
#endif		
		

		// nchroms_t nsel_pop = nchroms_t( popptr->pop_get_size() * sel_freq );
		// nchroms_t nunsel_pop = nchroms_t( popptr->pop_get_size() - nsel_pop );
		//PRINT13( deltaT, t, sel_freq, nsel, nsel_pop, nunsel, nunsel_pop, popptr->pop_get_size(), prob_coal_sel, prob_coal_unsel, prob_recomb, prob_gc, all_prob );

// #define COSI_DEV_SWEEP_ALT_RANDOM
// #ifndef COSI_DEV_SWEEP_ALT_RANDOM		
//     x = random_double();
// #else
		
// #endif

		(void)x;
		//static boost::mt19937 generator;
		//static boost::lagged_fibonacci44497 generator;
		enum choices_t { CHOICE_COAL_UNSEL, CHOICE_COAL_SEL, CHOICE_RECOMB, CHOICE_GC, CHOICE_NONE };
		double probs[] = { prob_coal_unsel, prob_coal_sel, prob_recomb, prob_gc,
											 1.0 - all_prob };
		boost::random::discrete_distribution<> dist(probs);
		choices_t rollResult( static_cast<choices_t>( dist( *getRandGen() ) ) );

    if (rollResult != CHOICE_NONE) {
      if (rollResult == CHOICE_RECOMB) {
				//PRINT5( "bef_recomb", newn, nsel, nunsel, t );
				
				/* Recombination.  Pick sel/unsel for new chrom piece */
				int selpop_index = demography->dg_get_pop_index_by_name(popname);
				recomb->recomb_execute(t, selpop_index, &loc, rec_nodes);

				/* Which segment carries old allele? */
				if (loc < sel_pos) {
					old_allele_node = rec_nodes[2];
					new_allele_node = rec_nodes[1];
				}
				else {
					old_allele_node = rec_nodes[1];
					new_allele_node = rec_nodes[2];
				}
				n_tot_recomb++;

				foundit = 0;
				bool_t oldWasSel = False;
				for (inode = 0; inode < nsel; inode++) {
					if (sel_nodes[inode] == rec_nodes[0]) {
						sel_nodes[inode] = old_allele_node ? old_allele_node : sel_nodes[ nsel-- - 1 ];
						foundit = 1;
						oldWasSel = True;
						break;
					}
				}
				if (!foundit) {
					for (inode = 0; inode < nunsel; inode++) {
						if (unsel_nodes[inode] == rec_nodes[0]) {
							unsel_nodes[inode] = old_allele_node ? old_allele_node : unsel_nodes[ nunsel-- - 1 ];
							foundit = 1;
							break;
						}
					}
				}
				chkCond( foundit, "sweep: Could not find node\n" );

				/* Decide whether the chrom that we've recombined onto is selected or unselected */
				if ( new_allele_node ) {
					if (random_double() < sel_freq) {
						if (nsel == selsize) {
							selsize *= 2;
							sel_nodes = (Node **)realloc(sel_nodes, selsize * sizeof(Node*));
							assert(sel_nodes != NULL);
						}
						sel_nodes[nsel] = new_allele_node;
						nsel++;
						if ( !old_allele_node && !oldWasSel ) n_switch_recomb++; 
					}
					else {
						if (nunsel == unselsize) {
							unselsize *= 2;
							unsel_nodes = (Node **)realloc(unsel_nodes, unselsize * sizeof(Node*));
							assert(unsel_nodes != NULL);
						}
						unsel_nodes[nunsel] = new_allele_node;
						nunsel++;
						if ( !old_allele_node && oldWasSel ) n_switch_recomb++;
					}
				}  // if ( new_allele_node )
      }  // if (x < prob_recomb / all_prob)
      else if (rollResult == CHOICE_COAL_UNSEL) {
				/* Coalescence among unselected chroms */
				//PRINT5( "bef_unsel_coal", newn, nsel, nunsel, t );
				newn = sw_coalesce(nunsel, unsel_nodes, popptr, t);
				//PRINT5( "aft_unsel_coal", newn, nsel, nunsel, t );
				assert(newn == nunsel - 1 || newn == nunsel - 2);
				nunsel = newn;
      }
      else if (rollResult == CHOICE_COAL_SEL) {
				/* Coalescence among selected chroms */
				//PRINT5( "bef_sw_coal", newn, nsel, nunsel, t );
				newn = sw_coalesce(nsel, sel_nodes, popptr, t);
				//PRINT5( "aft_sw_coal", newn, nsel, nunsel, t );
				assert(newn == nsel - 1 || newn == nsel - 2);
				nsel = newn;
      }
      else {
				assert( rollResult == CHOICE_GC );
				/* Gene conversion. Pick sel/unsel for new chrom piece */
				int selpop_index = demography->dg_get_pop_index_by_name(popname);
				geneConversion->gc_execute(t, selpop_index, &loc1, &loc2, gc_nodes);
				if (*gc_nodes != NULL) {
					/* Which segment carries old allele? */
					if (loc1 <= sel_pos && loc2 >= sel_pos) {
						old_allele_node = gc_nodes[1];
						new_allele_node = gc_nodes[2];
					}
					else {
						old_allele_node = gc_nodes[2];
						new_allele_node = gc_nodes[1];
					}

					foundit = 0;
					for (inode = 0; inode < nsel; inode++) {
						if (sel_nodes[inode] == gc_nodes[0]) {
							sel_nodes[inode] = old_allele_node;
							foundit = 1;
							break;
						}
					}
					if (!foundit) {
						for (inode = 0; inode < nunsel; inode++) {
							if (unsel_nodes[inode]== gc_nodes[0]) {
								unsel_nodes[inode] = old_allele_node;
								foundit = 1;
								break;
							}
						}
					}
					if (foundit != 1) {fprintf(stderr, "Could not find %d\n", gc_nodes[0]->getName());}
					assert(foundit == 1);

					/* Decide whether the chrom that we've recombined onto is selected or unselected */
					if (random_double() < sel_freq) {
						if (nsel == selsize) {
							selsize *= 2;
							sel_nodes = (Node **)realloc(sel_nodes, selsize * sizeof(Node*));
							assert(sel_nodes != NULL);
						}
						sel_nodes[nsel] = new_allele_node;
						nsel++;
					}
					else {
						if (nunsel == unselsize) {
							unselsize *= 2;
							unsel_nodes = (Node **)realloc(unsel_nodes, unselsize * sizeof(Node*));
							assert(unsel_nodes != NULL);
						}
						unsel_nodes[nunsel] = new_allele_node;
						nunsel++; 
					}
				}
      } 

      if (((int)popptr->pop_get_num_nodes()) != nsel + nunsel) {
				printf("mismatch: pop(n): %d  here: %d\n", (int)popptr->pop_get_num_nodes(), nsel + nunsel);
      }
    }
    /*    if (nsel == 1) {break;} */
  }
	
	if ( !leafset_is_empty( sel_leaves ) )
		 mutate->mutate_print_leafset(sel_pos, sel_leaves, t, popname);

  util::cosi_free( sel_nodes );
  util::cosi_free( unsel_nodes );

  //std::cerr << "selective sweep ends at " << t << " generations after " << nsteps << " steps.\n";
  if ( verbose ) std::cout << "selective sweep ends at " << t << " generations\n";
	//PRINT2( n_switch_recomb, n_tot_recomb );
  return t;
}  // Sweep::sweep_execute()


genid Sweep::sweep_stochastic_execute(popid popname, gensInv_t sel_coeff,
																			double /*heterozygousEffect*/, genid gen, loc_t sel_pos, 
																			freq_t final_sel_freq) {
  Pop *popptr;
  int nsel, nunsel, newn, foundit;
  Node **sel_nodes, **unsel_nodes, *rec_nodes[3], *gc_nodes[3];
  Node *old_allele_node, *new_allele_node;
	genid t, tend;
	gens_t deltaT;
	freq_t sel_freq;
  freq_t epsilon;
	double x, /*alpha,*/ rate;
  prob_t prob_coal_unsel, prob_coal_sel, prob_recomb, all_prob;
	genid tsave, t_nonsweep;
	gens_t end_shift = ZERO_GENS;
	gens_t poisson_event_time;
  loc_t loc;
  prob_t prob_gc;
  loc_t loc1, loc2;
  int selsize, unselsize, dorandom=1;

  /*  sw_trajectory(popname, sel_coeff); */

  popptr = demography->dg_get_pop_by_name(popname);
  nunsel = nsel = 0;
  selsize = unselsize = popptr->pop_get_num_nodes();
  sel_nodes = (Node **)malloc(selsize * sizeof(Node*));
  unsel_nodes = (Node **)malloc(unselsize * sizeof(Node*));
  assert(sel_nodes != NULL && unsel_nodes != NULL);

  epsilon = 1 / (2. * popptr->pop_get_size());
  if (final_sel_freq > 1-epsilon) {
    end_shift = ZERO_GENS;
    dorandom = 0;
  }
  else {
    end_shift = gens_t( log( (1-final_sel_freq) / (final_sel_freq * epsilon * (1-epsilon)) ) / sel_coeff );
  }
  deltaT = gens_t( .01 / sel_coeff );
  tend = gen - gens_t( 2 * log(epsilon) / sel_coeff );
  //alpha = 2 * popptr->pop_get_size() * sel_coeff;

  for (int inode = 0; inode < ((int)popptr->pop_get_num_nodes()); inode++) {
    if (!dorandom || random_double() < final_sel_freq) {
      sel_nodes[nsel] = popptr->pop_get_node(inode);
      nsel++;
    }
    else {
      unsel_nodes[nunsel] = popptr->pop_get_node(inode);
      nunsel++;
    }
  }
  sel_freq = 1;
  tsave = gen;
  
  rate = sw_get_poisson_rate_nonsel(popname); 
  poisson_event_time = gens_t( 1e50 );
  if (rate > 0) {
    poisson_event_time = gens_t( poisson_get_next(rate) );
  }
  t_nonsweep = gen + poisson_event_time;

  FILE *trajFile = fopen( "traj.tsv", "wt" );
  for (t = gen; sel_freq >= epsilon || nsel > 1; t += deltaT) {
    while (t > t_nonsweep) {
      /* Process event in a nonsweep population, and update time for next nonsweep event */
      sw_do_poisson_nonsel(popname, t_nonsweep); 
      poisson_event_time = gens_t( poisson_get_next(sw_get_poisson_rate_nonsel(popname)) );
      t_nonsweep += poisson_event_time;  
    }

    sel_freq = epsilon / (epsilon + (1 - epsilon) * exp(sel_coeff*(t + end_shift - tend)));
		fprintf( trajFile, "%f\t%f\n", double( ToDouble( t ) ), sel_freq );
    prob_coal_sel = ToDouble( deltaT ) * nsel * (nsel-1) / 4 / popptr->pop_get_size() / sel_freq;
    prob_coal_unsel = ToDouble( deltaT ) * nunsel * (nunsel-1) / 4 / popptr->pop_get_size() / (1 - sel_freq);
    prob_recomb = ToDouble( deltaT ) * (nsel + nunsel) * ToDouble( genMap->getRegionRecombRateAbs() );
    prob_gc = ToDouble( deltaT ) * (nsel + nunsel) * ToDouble( geneConversion->getRegionGeneConvRate() );
    all_prob = 0; 
    all_prob += prob_coal_unsel;
    all_prob += prob_coal_sel;
    all_prob += prob_recomb;
    all_prob += prob_gc; 

    x = random_double();

    if (all_prob >= x) {
      if (x < prob_recomb / all_prob) {
				/* Recombination.  Pick sel/unsel for new chrom piece */
				int selpop_index = demography->dg_get_pop_index_by_name(popname);
				recomb->recomb_execute(t, selpop_index, &loc, rec_nodes);
				if (*rec_nodes != NULL) {
					/* Which segment carries old allele? */
					if (loc < sel_pos) {
						old_allele_node = rec_nodes[2];
						new_allele_node = rec_nodes[1];
					}
					else {
						old_allele_node = rec_nodes[1];
						new_allele_node = rec_nodes[2];
					}

					foundit = 0;
					for (int inode = 0; inode < (int)nsel; inode++) {
						if (sel_nodes[inode] == rec_nodes[0]) {
							sel_nodes[inode] = old_allele_node;
							foundit = 1;
							break;
						}
					}
					if (!foundit) {
						for (int inode = 0; inode < (int)nunsel; inode++) {
							if (unsel_nodes[inode] == rec_nodes[0]) {
								unsel_nodes[inode] = old_allele_node;
								foundit = 1;
								break;
							}
						}
					}
					if (foundit != 1) {fprintf(stderr, "Could not find %d\n", rec_nodes[0]->getName());}
					assert(foundit == 1);

					/* Decide whether the chrom that we've recombined onto is selected or unselected */
					if (random_double() < sel_freq) {
						if (nsel == selsize) {
							selsize *= 2;
							sel_nodes = (Node **)realloc(sel_nodes, selsize * sizeof(Node*));
							assert(sel_nodes != NULL);
						}
						sel_nodes[nsel] = new_allele_node;
						nsel++;
					}
					else {
						if (nunsel == unselsize) {
							unselsize *= 2;
							unsel_nodes = (Node **)realloc(unsel_nodes, unselsize * sizeof(Node*));
							assert(unsel_nodes != NULL);
						}
						unsel_nodes[nunsel] = new_allele_node;
						nunsel++; 
					}
				}
      }

      else if (x < (prob_recomb + prob_coal_unsel) / all_prob) {
				/* Coalescence among unselected chroms */
				newn = sw_coalesce(nunsel, unsel_nodes, popptr, t);
				assert(newn == nunsel - 1 || newn == nunsel - 2);
				nunsel = newn;
      }
      else if (x < (prob_recomb + prob_coal_unsel + prob_coal_sel) / all_prob) {
				/* Coalescence among selected chroms */
				newn = sw_coalesce(nsel, sel_nodes, popptr, t);
				assert(newn == nsel - 1 || newn == nsel - 2);
				nsel = newn;
      }
      else {
				/* Gene conversion. Pick sel/unsel for new chrom piece */
				int selpop_index = demography->dg_get_pop_index_by_name(popname);
				geneConversion->gc_execute(t, selpop_index, &loc1, &loc2, gc_nodes);
				if (*gc_nodes != NULL) {
					/* Which segment carries old allele? */
					if (loc1 <= sel_pos && loc2 >= sel_pos) {
						old_allele_node = gc_nodes[1];
						new_allele_node = gc_nodes[2];
					}
					else {
						old_allele_node = gc_nodes[2];
						new_allele_node = gc_nodes[1];
					}

					foundit = 0;
					for (int inode = 0; inode < nsel; inode++) {
						if (sel_nodes[inode] == gc_nodes[0]) {
							sel_nodes[inode] = old_allele_node;
							foundit = 1;
							break;
						}
					}
					if (!foundit) {
						for (int inode = 0; inode < nunsel; inode++) {
							if (unsel_nodes[inode]== gc_nodes[0]) {
								unsel_nodes[inode] = old_allele_node;
								foundit = 1;
								break;
							}
						}
					}
					if (foundit != 1) {fprintf(stderr, "Could not find %d\n", gc_nodes[0]->getName());}
					assert(foundit == 1);

					/* Decide whether the chrom that we've recombined onto is selected or unselected */
					if (random_double() < sel_freq) {
						if (nsel == selsize) {
							selsize *= 2;
							sel_nodes = (Node **)realloc(sel_nodes, selsize * sizeof(Node*));
							assert(sel_nodes != NULL);
						}
						sel_nodes[nsel] = new_allele_node;
						nsel++;
					}
					else {
						if (nunsel == unselsize) {
							unselsize *= 2;
							unsel_nodes = (Node **)realloc(unsel_nodes, unselsize * sizeof(Node*));
							assert(unsel_nodes != NULL);
						}
						unsel_nodes[nunsel] = new_allele_node;
						nunsel++; 
					}
				}
      } 

      if (((int)popptr->pop_get_num_nodes()) != nsel + nunsel) {
				printf("mismatch: pop(n): %d  here: %d\n", ((int)popptr->pop_get_num_nodes()), nsel + nunsel);
      }
    }
    /*    if (nsel == 1) {break;} */
  }
  fclose( trajFile );
  {
		leafset_p sel_leaves = LEAFSET_NULL;
		for (int inode = 0; inode < nsel; inode++)
			 ForEach( const seglist::Seg& seg, *sel_nodes[inode]->getSegs() )
					if ( seg.getBeg() <= sel_pos && sel_pos <= seg.getEnd() )
						 sel_leaves = leafset_union( sel_leaves, seg.getLeafset() );
		
		if ( !leafset_is_empty( sel_leaves ) )
			 mutate->mutate_print_leafset(sel_pos, sel_leaves, t, popname);
  }

  util::cosi_free( sel_nodes );
  util::cosi_free( unsel_nodes );

  std::cout << "selective sweep ends at " << t << " generations\n";
  return t;
}  // Sweep::sweep_stochastic_execute()

//
// Method: coalesce
//
// Pick and coalesce two nodes either within the selected or the unselected subset of the
// population in which sweep is happening.
//
// Params:
//
//    nnode - the number of nodes in the selected/unselected subset ('nodes' below)
//    nodes - the nodes in the subset
//    popptr - the Pop object for the population undergoing sweep
//    t - the generation at which the coalescence happens
//
int Sweep::sw_coalesce(int nnode, Node **nodes, Pop *popptr, genid t) {
  /* The difference between this routine and the one in demography.c is that here 
     I have to update both the main population node list and the local list that I'm using. */
  int node1index, node2index, inode;
  Node *node1, *node2, *newnode;

  // Choose the pair of nodes to coalesce.
  node1index = (int) (random_double() * nnode);
  node2index = (int) (random_double() * (nnode - 1));
  if (node2index >= node1index) node2index++;
  node1 = nodes[node1index];
  node2 = nodes[node2index];

  mutate->mutate_put_muts_on_seglist( node1->getSegs(), t - node1->getGen(), node1->getGen(), popptr->pop_get_name() );
  mutate->mutate_put_muts_on_seglist( node2->getSegs(), t - node2->getGen(), node2->getGen(), popptr->pop_get_name() );

  /* STEP 3 */
  popptr->pop_remove_node (node1);
  popptr->pop_remove_node (node2);
  
  if (node1index < node2index) {
    for (inode = node2index; inode < nnode-1; inode++) {nodes[inode] = nodes[inode+1];}
    nnode--;
    for (inode = node1index; inode < nnode-1; inode++) {nodes[inode] = nodes[inode+1];}
  }
  else {
    for (inode = node1index; inode < nnode-1; inode++) {nodes[inode] = nodes[inode+1];}
    nnode--;
    for (inode = node2index; inode < nnode-1; inode++) {nodes[inode] = nodes[inode+1];}
  }
  nnode--;

  newnode = nodePool->node_coalesce (&node1, &node2, /* gen= */ t );
  
  /* STEP 4 */
	if ( newnode ) {
		popptr->pop_add_node (newnode);
		 nodes[nnode++] = newnode;

		 demography->dg_log (Demography::COALESCE, t, node1, node2, newnode, popptr);
	}

  return nnode;
}  // Sweep::sw_coalesce()

/* coalesce_get_rate, modified to ignore coalescence in target pop */
double Sweep::sw_coalesce_get_rate_nonsel (popid sel_popname) const {
  int numpops = demography->dg_get_num_pops();
  int i;
  double rate = 0;
  int numnodes;
  int popsize;
  
  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
    popsize = demography->dg_get_pop_size_by_index (i);
    if (numnodes > 1)
			 rate += (double) (numnodes * (numnodes - 1)) 
					/ (4 * popsize);
  }
  return rate;
}

/* Similar to preceeding function */
double Sweep::sw_recomb_get_rate_nonsel (popid sel_popname) const {
  int numpops = demography->dg_get_num_pops();
  int i;
  double rate = 0, r = ToDouble( genMap->getRegionRecombRateAbs() );
  int numnodes;

  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
    rate += ( numnodes * r);
  }
  return rate;  
}

/* Similar to preceeding function */
double Sweep::sw_gc_get_rate_nonsel (popid sel_popname) const {
  int numpops = demography->dg_get_num_pops();
  int i;
  double rate = 0, r = ToDouble( geneConversion->getRegionGeneConvRate() );
  int numnodes;

  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
    rate += ( numnodes * r);
  }
  return rate;  
}

int Sweep::sw_coalesce_pick_popindex (popid sel_popname) {
  double  randcounter, sumprob, totprob;
  int     popindex = -1,
		 numpops = demography->dg_get_num_pops(),
		 numnodes, popsize, i;

  totprob = 0;
  for (i = 0; i < numpops; i++) {
		popid thispop = demography->dg_get_pop_name_by_index(i);
		if (thispop == sel_popname) {continue;}
		numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
		popsize = demography->dg_get_pop_size_by_index (i);
		totprob += (double) (numnodes * (numnodes - 1)) 
			 / (4 * popsize);
  }
  randcounter = random_double() * totprob;

  sumprob = 0;
  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    numnodes = demography->dg_get_num_nodes_in_pop_by_index (i);
    popsize = demography->dg_get_pop_size_by_index (i);
    sumprob += (double) (numnodes * (numnodes - 1)) 
			 / (4 * popsize);
    if (randcounter <= sumprob) {
			popindex = i;
			break;
    }
  }		
  assert(popindex >= 0);
  return popindex;
}

int Sweep::sw_recomb_pick_popindex(popid sel_popname) {
  int popindex = -1, i, sumnodes = 0;
  int totnodes;
  int numpops = demography->dg_get_num_pops();
  double randcounter;

  totnodes = 0;
  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    totnodes += demography->dg_get_num_nodes_in_pop_by_index (i);
  }
  randcounter = totnodes * random_double();
  assert(randcounter <= totnodes);
  
  /* weight pops by numnodes. */
  for (i = 0; i < numpops; i++) {
    popid thispop = demography->dg_get_pop_name_by_index(i);
    if (thispop == sel_popname) {continue;}
    sumnodes += demography->dg_get_num_nodes_in_pop_by_index(i); 
    if (randcounter <= sumnodes) {
      popindex = i;
      break;
    }
  }
  assert(popindex >= 0);
  return popindex;
}

// Method: sw_get_poisson_rate_nonsel
//
// Returns the probability of a Poisson event (coalescence, gene conversion or recombination;
// migration is not included) in a single generation, in a non-selected population.
//
// Params:
//
//    sel_popname - the selected population
//
double Sweep::sw_get_poisson_rate_nonsel(popid sel_popname) {
  coalesce_rate_nonsel = sw_coalesce_get_rate_nonsel(sel_popname);
  /*  migrate_rate = migrate_get_rate(); */
  geneconv_rate_nonsel = sw_gc_get_rate_nonsel(sel_popname); 
  recombination_rate_nonsel = sw_recomb_get_rate_nonsel(sel_popname);
  poisson_rate_nonsel = (double) (coalesce_rate_nonsel + migrate_rate_nonsel + recombination_rate_nonsel + geneconv_rate_nonsel);

  return poisson_rate_nonsel;
}

void Sweep::sw_do_poisson_nonsel(popid sel_popname, genid gen) {
  prob_t randdouble = random_double();
  loc_t dum, dum2;
  int recomb_pop;
  
  if (randdouble < recombination_rate_nonsel / poisson_rate_nonsel) {
    recomb_pop = sw_recomb_pick_popindex(sel_popname);
    recomb->recomb_execute(gen, recomb_pop, &dum, (Node **)NULL);
  }
  else if (randdouble < (recombination_rate_nonsel + migrate_rate_nonsel) / poisson_rate_nonsel) {
    /* migrate_execute(gen);*/
  }
  else if (randdouble < (recombination_rate_nonsel + migrate_rate_nonsel + coalesce_rate_nonsel) / poisson_rate_nonsel) { 
    int popindex = sw_coalesce_pick_popindex(sel_popname);
    demography->dg_coalesce_by_index (popindex, gen);
  }
  else {
    int popindex = sw_recomb_pick_popindex(sel_popname);
    geneConversion->gc_execute(gen, popindex, &dum, &dum2, (Node **)NULL);  
  }
  return;

}

int Sweep::sw_trajectory(popid popname, double s) {
  int iter, nsel, t;
  double freq_thresh=.9, freq;
  Pop *popptr;  

  fprintf(stderr, "sel coeff: %f\n", s);
  popptr = demography->dg_get_pop_by_name(popname);
  for (iter = 0; iter < 100000; iter++) {

    nsel = 1;
    for (t = 0; t < 100000; t++) {
      if (nsel == 0) {break;}
      freq = (double) nsel / popptr->pop_get_size() / 2;
      if (freq > freq_thresh) {
				fprintf(stderr, "success, t = %d, iter = %d\n", t, iter);
				exit(0);
				return 1;
      }
      /*      freq = nsel * (1. + s) / (nsel * (1. + s) + 2*popptr->pop_get_size() - nsel);	 */
      freq = nsel * (1. + s) / (nsel * s + 2*popptr->pop_get_size());	
      nsel = ranbinom(2*popptr->pop_get_size(), freq);
      /*      fprintf(stderr, "nsel returned: %d\n", nsel); */
    }
  }
  fprintf(stderr, "failed\n");
  exit(0);
  return 0;
}

//
// Section: Simulating frequency trajectory
//
// Simulating the frequency trajectory of the causal allele,
// forward in time.
//


//
// Func: sweep_sim_freq_traj
//
// Simulate the frequency trajectory of causal allele.
//
void Sweep::sweep_sim_freq_traj( genid /*start_gen*/ ) {
  
  
}

// End class impl: Sweep


/////////////////

//
// Class: Event_SweepNew
//
// A selective sweep.
//
// The sweep is implemented as follows.  We take a joint allele frequency trajectory
// for the causal allele, specifying the allele's frequency in each population.
// (Or, we take parameters of the sweep and from them construct the frequency trajectory --
// how, is another subject).  The trajectory ends at <gen>, this sweep's ending generation,
// and starts at a generation at which the total frequency of the causal allele drops to zero
// (viewing time as going backwards.)
//
// We try to reuse as much of the neutral simulation machinery as possible when
// simulating sweeps.  At the end of the sweep, when we first encounter it during backwards
// simulation, we split each population into two: one will contain nodes from the original population
// that carry the derived allele at the causal mutation SNP, and the other nodes that carry the
// ancestral allele.  More specifically, for each population we create one new "companion" population
// and move the derived-allele nodes into it.  The original population, with the original <popid>,
// then contains nodes with the ancestral allele.  Keeping the original popid lets us catch events
// executed on the original population and forward them to the ancestral & derived subpops as appropriate.
// The number of nodes moved from the original pop to the companion derived-allele pop is determined by
// the frequency of the derived (causal) allele in the original pop at the end of the sweep,
// as determined by the frequency trajectory.
//
// (Let's call each original pop which we split into two 'orig-pop', and the two resulting pops
// 'anc-pop' and 'der-pop', with the understanding that 'anc-pop' is the same <Pop> object as 'orig-pop' --
// we just moved some of its <Nodes> into the companion 'der-pop'.)
//
// Then, we let the simulation run as usual, with the following adjustments:
//   - at each point of the frequency trajectory, for each orig-pop we adjust the relative sizes of the
//     anc-pop and der-pop to reflect the derived allele frequency in the orig-pop at that generation
//   - after each recombination or gene conversion, we take the resulting node that did NOT receive
//     the causal mutation location, and decide whether it has the derived or the ancestral allele
//     at the causal location by flipping a biased coin based on the frequency of the derived allele
//     (aka pop size of the der-pop); then, if the node is in der-pop and it was assigned the ancestral
//     allele, we move it to the anc-pop, and vice versa.
//   - when migration rate is set between two orig-pops, we set the same migration rate between their corresponding
//     der-pops.  since the anc-pops are the orig-pops, their migration rate has already been set.
//
class Event_SweepNew: virtual public HistEvents::Event {
public:
	 Event_SweepNew( HistEvents *histEvents_, const string& label_, genid gen_, popid sweepPop_, gensInv_t selCoeff_, loc_t selPos_,
									 freq_t final_sel_freq_ ):
		 Event( histEvents_, label_, gen_ ), sweepPop( sweepPop_ ), selCoeff( selCoeff_ ), selPos( selPos_ ), final_sel_freq( final_sel_freq_ )
			{ }
	 Event_SweepNew( HistEvents *histEvents_, istream& is ): Event( histEvents_, is )  {
		 is >> sweepPop >> gen >> selCoeff >> selPos >> final_sel_freq;
	 }
	 static const char *typeStr() { return "sweep_new"; }
	 virtual eventKind_t getEventKind() const { return E_SWEEP; }	 
				
	 virtual ~Event_SweepNew();

	 virtual genid execute();


private:
	 // Field: sweepPop
	 // The sole population in which the sweep happens.
	 // The frequency of the derived allele in other populations always remains zero.
	 popid sweepPop;

	 // Field: selCoeff
	 // The selection coefficient.
	 gensInv_t selCoeff;

	 // Field: selPos
	 // The location of the causal mutation.
	 loc_t selPos;
	 
	 // Field: final_sel_freq
	 // The frequency of the derived allele at the end of the sweep.
	 freq_t final_sel_freq;

	 // Field: pop2companion
	 // For each pop which we split when we first encountered the sweep (going backwards),
	 // the corresponding pop we created to hold the nodes carrying the derived allele at <selPos>;
	 // and the reverse map.
	 map<popid, popid> pop2companion;

	 // Field: origPops
	 // The original pops we have split when we first encountered the end of the sweep during
	 // backwards simulation.
	 set<popid> origPops;
	 
	 // Field: freqTraj
	 // The frequency trajectory of the causal allele
	 FreqTrajP freqTraj;

	 // Field: sel_leaves
	 // Leaves which inherit the causal allele.
	 leafset_p sel_leaves;

	 // Class: SweepHook
	 // Callbacks invoked when certain events happen during a sweep.  We're using the standard
	 // neutral simulation machinery for most of the sweep simulation, but a few things are
	 // different -- this class implements these differences.
	 class SweepHook: public Hook {
	 public:
			SweepHook( Event_SweepNew *evt_): evt( evt_ ), inHook( False ) {}
			virtual ~SweepHook();
			virtual void handle_recomb( Node *, Node *, loc_t, genid);
			virtual void handle_gc( Node *, Node *, loc_t, loc_t, genid);
			virtual void handle_coal( Node * );
			virtual void handle_set_migrate_rate( popid from, popid to, prob_per_chrom_per_gen_t rate );
			
	 private:
			// Private field: evt
			// The sweep event to which these hooks relate.
			Event_SweepNew *evt;

			// Private field: inHook
			// Whether we're already executing a handle_* method of this hook.
			// Used to prevent recursively invoking a handle_* method.
			bool_t inHook;

			// Private methodP: determineAlleleAtSelPos_
			// Given a node (resulting from a recombination or gene conversion) for which
			// <Node::segs> does not contain <Event_SweepNew::selPos>, choose the allele at selPos
			// and if needed move node to the appropriate partial pop (anc-pop or der-pop).
			void determineAlleleAtSelPos_( Node *node );
	 };  // class SweepHook

	 typedef boost::shared_ptr<SweepHook> SweepHookP;

	 // Field: hook
	 // Hook that makes adjustments to neutral coalescent machinery needed to simulate sweeps.
	 SweepHookP sweepHook;

	 // Private methodP: init_
	 // Called when the backward simulation first encounters the sweep
	 // (at the generation corresponding to the sweep's end).
	 // Splits each population existing at that time into anc-pop and der-pop
	 // containing nodes carrying the ancestral/derived allele at the causal site,
	 // and establishes hooks (callbacks) to be called during the simulation
	 // to handle events according to the sweep.
	 void init_();

};  // class Event_SweepNew

//
// Implementation of class: Event_SweepNew
//
const char *Event_SweepNew_typeStr();
const char *Event_SweepNew_typeStr() { return "sweep_new"; }
HistEvents::EventP make_shared_Event_SweepNew( HistEvents *histEvents, istream& is );
HistEvents::EventP make_shared_Event_SweepNew( HistEvents *histEvents, istream& is ) {
  HistEvents::EventP ep( new Event_SweepNew( histEvents, is ) );
  return ep;
}

Event_SweepNew::~Event_SweepNew() {}

// Private method: init_
// Called when the backward simulation first encounters the sweep
// (at the generation corresponding to the sweep's end).
// Splits each population existing at that time into anc-pop and der-pop
// containing nodes carrying the ancestral/derived allele at the causal site,
// and establishes hooks (callbacks) to be called during the simulation
// to handle events according to the sweep.
void Event_SweepNew::init_() {
	if ( theTrajFN.empty() ) {
		const nchroms_t popSize = getDemography()->dg_get_pop_by_name( sweepPop )->pop_get_size();
		freqTraj = boost::make_shared<DeterministicSweepTraj>( sweepPop, gen, selCoeff, final_sel_freq,
																													 popSize );
	} else 
		freqTraj = boost::make_shared<TrajFromFile>( theTrajFN, gen, final_sel_freq );
	
	vector<PopP> origPopsCopy( getDemography()->getPops() );  // save a copy of the current list of pops;
	// iterate over this copy which is unaffected when we add new pops.
	sel_leaves = make_empty_leafset();
	
	ForEach( PopP pop, origPopsCopy ) {
		//
		// Split pop into two pops, anc-pop and der-pop.
		//
		
		nchroms_t final_pop_size = pop->pop_get_size();
		nchroms_t derPopSize( double(final_pop_size) * freqTraj->getCurFreq( pop->pop_get_name() ) );
		std::string derPopLabel( pop->get_label() + "_sel" );

		// For nodes carrying the derived allele at selPos, create a new pop.
		// For nodes carrying the ancestral allele, reuse the original Pop.
		Pop *derPop =
			 getDemography()->dg_add_pop( boost::make_shared<Pop>( getDemography()->find_unused_popname(),
																																		derPopSize,
																																		derPopLabel ),
																					 gen );

		
		getDemography()->dg_move_nodes_by_name( pop->pop_get_name(), derPop->pop_get_name(), freqTraj->getCurFreq( pop->pop_get_name() ),
																									 gen
#ifdef COSI_DEV
																									 , /* exactFraction= */ True 																									 
#endif																									 
			);
		
		getDemography()->dg_set_pop_size_by_name( gen, pop->pop_get_name(), final_pop_size - derPop->pop_get_size() );

		ForEach (Node *it, derPop->getMembers() ) {
			ForEach( const seglist::Seg& seg, *it->getSegs() ) 
				 if ( seg.getBeg() <= selPos && selPos <= seg.getEnd() )
						sel_leaves = leafset_union( sel_leaves, seg.getLeafset() );
		}
		
		// Record the correspondence between the anc-pop and the der-pop.
		pop2companion.insert( make_pair( pop->pop_get_name(), derPop->pop_get_name() ) );
		pop2companion.insert( make_pair( derPop->pop_get_name(), pop->pop_get_name() ) );
		origPops.insert( pop->pop_get_name() );
	}  // for each pop

	// For any active migrations between origPops (now ancPops), add the same migration between their corresponding derPops.
	vector< const Migrate::MigrateRate * > existingMigrs;
	for ( const Migrate::MigrateRate *migration = getMigrate()->getMigrations(); migration; migration = migration->next )
		 existingMigrs.push_back( migration );

	// Set up a hook, so that we get called when certain events happen during the sweep
	// (specifically recombinations, gene conversions, and setting of migration rate.)
	Event_SweepNew *thisPtr = this;
	sweepHook = boost::make_shared<SweepHook>( thisPtr );
	getDemography()->getHooks()->addHook( sweepHook );
	ForEach( const Migrate::MigrateRate *migration, existingMigrs )
		 sweepHook->handle_set_migrate_rate( migration->frompop->pop_get_name(),
																				 migration->topop->pop_get_name(), migration->rate );

}  // Event_SweepNew::init()

unsigned nrecomb_sel = 0, nrecomb_unsel = 0, nnew_sel = 0, nnew_unsel = 0, ncoal_sel = 0, ncoal_unsel = 0;

//
// Method: execute
//
// Execute a selective sweep.  This method is invoked the first time the sweep is encountered during pastward simulation,
// and then each time the frequency of the selected allele in the population changes according to the frequency trajectory.
//
// When first called: splits each pop into der-pop and anc-pop; adds hooks to be called on recomb and gc events.
//
// When called for frequency changes:
// Updates the sizes of the der-pop and anc-pop making up each orig-pop
// (see <Event_SweepNew> for definitions of terms used here), based on the frequency of the selected allele in the orig-pop
// given by <freqTraj>.
//
genid Event_SweepNew::execute() {
	// If this is the first time this sweep is encountered,
	// initialize the processing of this sweep by splitting each orig-pop currently existing in the demographic model into
	// anc-pop and der-pop, etc.
	if ( pop2companion.empty() ) init_();

	bool_t someNotCoal = False;
	ForEach( popid origPop, origPops ) {
		Pop *ancPop = getDemography()->dg_get_pop_by_name( origPop );
		Pop *derPop = getDemography()->dg_get_pop_by_name( map_get( pop2companion, origPop ) );

		someNotCoal = someNotCoal || ( derPop->pop_get_num_nodes() > 1 );

		nchroms_t totSize = derPop->pop_get_size() + ancPop->pop_get_size();
		nchroms_t derSize = nchroms_t( double( totSize ) * freqTraj->getCurFreq( ancPop->pop_get_name() ) );
		//nchroms_t ancSize = totSize - derSize;

#ifdef COSI_DEV		
		{
			static ofstream trajOut;
			if ( !trajOut.is_open() ) {
				trajOut.open( "trajnew.tsv", ofstream::out );
				trajOut << "gen\tsel_freq\tnsel\tnsel_in_pop\tnunsel\tnunsel_in_pop\tnrecomb_sel\tnrecomb_unsel\tnnew_sel\tnnew_unsel\tncoal_sel\tncoal_unsel\tprob_coal_sel\tprob_coal_unsel\n";
			}
			frac_t sel_freq = freqTraj->getCurFreq( ancPop->pop_get_name() );
			nchroms_t nsel = derPop->pop_get_num_nodes();
			nchroms_t nunsel = ancPop->pop_get_num_nodes();
			int nsel_in_pop = ((int)(2*totSize * sel_freq));
			int nunsel_in_pop = 2*totSize - nsel_in_pop;
			if ( ToInt( origPop ) == 1 ) {
				trajOut << gen << "\t" << sel_freq << "\t" << nsel << "\t" << nsel_in_pop << "\t" << nunsel << "\t" << nunsel_in_pop
								<< "\t" << nrecomb_sel << "\t" << nrecomb_unsel << "\t" << nnew_sel << "\t" << nnew_unsel << "\t" << ncoal_sel
								<< "\t" << ncoal_unsel << "\t" << computeApproxCoalProb( derPop ) << "\t" <<
					 computeApproxCoalProb( ancPop ) << "\n";
			}
		}
#endif		


		// If derSize has dropped below the actual sample size, force sample size to be <= derSize by performing coalescences.
		// if ( derSize > 0 ) {
		// 	while ( derPop->pop_get_num_nodes() > derSize ) {
		// 		getDemography()->dg_coalesce_by_pop( derPop, gen );
		// 		gen += gens_t( 1e-12 );
		// 	}
		// }
			
		getDemography()->dg_set_pop_size_by_name( gen, derPop->pop_get_name(),
																										 std::max( derSize, derPop->pop_get_num_nodes() ) );
		
		getDemography()->dg_set_pop_size_by_name( gen, ancPop->pop_get_name(), totSize - derPop->pop_get_size() );
	}

	freqTraj->next();
	genid oldgen = gen;
	if ( someNotCoal || !freqTraj->done() ) {
		gen = freqTraj->getCurGen();
		addEvent( shared_from_this() );
	} else {

		//
		// The sweep has completed.  Merge each split population back into a single pop.
		// Remove the hooks that got called after each recomb or gc event.
		// Put the selected mutation on the leaves.  (We could have done it at the start -- in the pastwards sense --
		// of the sweep, but there we did not know the generation at which the mutation originated, and
		// we want to be able to output the generation of each mutation.)
		//
		
		std::pair<popid,popid> p;
		ForEach( popid origPop, origPops ) {
			Pop *ancPop = getDemography()->dg_get_pop_by_name( origPop );
			Pop *derPop = getDemography()->dg_get_pop_by_name( map_get( pop2companion, origPop ) );
			
			getDemography()->dg_move_nodes_by_name( derPop->pop_get_name(), ancPop->pop_get_name(), 1.00, gen, /* exactFraction= */ True );
		}  // for each split pop pair
		
		getDemography()->getHooks()->removeHook( sweepHook );

		getDemography()->getMutate()->mutate_print_leafset(selPos, sel_leaves, gen, sweepPop);		
	}  // sweep finished
	//
	// later: if an event such as a migration event gets executed, propagate this correctly.
	// but, initially migration is not an issue.
	//

	// if migration rates are unequal, wouldn't that change pop size?  but how would that affect effective pop size?
	return oldgen;
}  // Event_SweepNew::execute()

//
// Implementation of class: Event_SweepNew::SweepHook
//

Event_SweepNew::SweepHook::~SweepHook() {
}

// Private method: determineAlleleAtSelPos_
// Given a node (resulting from a recombination or gene conversion) for which
// <Node::segs> does not contain <Event_SweepNew::selPos>, choose the allele at selPos
// and if needed move node to the appropriate partial pop (anc-pop or der-pop).
void Event_SweepNew::SweepHook::determineAlleleAtSelPos_( Node *node ) {
	if ( node ) {
		Pop *curPop = node->getPop();
		Pop *companionPop = evt->getDemography()->dg_get_pop_by_name( map_get( evt->pop2companion, curPop->pop_get_name() ) );
		Pop *derPop, *ancPop;
		if ( STLContains( evt->origPops, curPop->pop_get_name() ) ) { ancPop = curPop; derPop = companionPop; nrecomb_unsel++; }
		else { ancPop = companionPop; derPop = curPop; nrecomb_sel++; }

		freq_t testFreq = evt->freqTraj->getCurFreq( ancPop->pop_get_name() );
		if ( testFreq > 0 ) {
			prob_t pval = evt->random_double();
			Pop * shouldBeInPop = ( pval < testFreq ) ? derPop : ancPop;
			( pval < testFreq ? nnew_sel : nnew_unsel )++;
			if ( curPop != shouldBeInPop ) {
				curPop->pop_remove_node( node );
				shouldBeInPop->pop_add_node( node );
			}
		}
	}
}
void Event_SweepNew::SweepHook::handle_recomb( Node *node1, Node *node2, loc_t loc, genid ) {
	determineAlleleAtSelPos_( evt->selPos < loc ? node2 : node1 );
}
	 void Event_SweepNew::SweepHook::handle_gc( Node *node1, Node *node2, loc_t loc1, loc_t loc2, genid) {
	determineAlleleAtSelPos_( loc1 <= evt->selPos && evt->selPos <= loc2 ? node2 : node1 );
}


void Event_SweepNew::SweepHook::handle_coal( Node *n ) {
	( STLContains( evt->origPops, n->getPop()->pop_get_name() ) ? ncoal_unsel : ncoal_sel )++;
}

void Event_SweepNew::SweepHook::handle_set_migrate_rate( popid from, popid to, prob_per_chrom_per_gen_t rate ) {
	if ( !inHook ) {
		inHook = True;
		evt->getMigrate()->migrate_set_rate( map_get( evt->pop2companion, from ),
																				 map_get( evt->pop2companion, to ), rate );
	}
}




}  // namespace cosi

