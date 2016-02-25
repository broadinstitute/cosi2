#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <cosi/coal_data_cosi.h>
//#define NPOPS 5 //2 is a dummy population 
//#define NPOPPAIRS 6 //(1, 3), (1, 4), (1, 5), (3, 4), (3, 5), (4, 5)

void calc_Fst_main( int numReps, int NPOPS,
										int excludePopIdx,
										void (*get_coal_data_from_cosi)( coal_data *, int /* irep */, int /* ipop */ ) ) ;

void calc_Fst_main( int numReps, int NPOPS,
										int excludePopIdx,
										void (*get_coal_data_from_cosi)( coal_data *, int /* irep */, int /* ipop */ ) ) {
    // char outfilename[264];
    // char basefile[264], filebase[264];
    // char repstr[3];
    // int numReps;


	int NPOPS_adj = NPOPS;
	if ( excludePopIdx >= 0 ) --NPOPS_adj;
	int NPOPPAIRS = NPOPS_adj * ( NPOPS_adj - 1 ) / 2;
	
	coal_data data;
    int nsnp = 0;
    int irep, ipop, jpop, isnp, isamp, ipopPair;
    int ni, nj;
    int nai[2], naj[2];
    double fst_total;
    int *nall1[NPOPS+1], *nall2[NPOPS+1];
    FILE *outf=NULL;
    double msp, msg, p[2], nc, num, pmean, nic, njc, denom;
  	double fst_vals[NPOPS+1][NPOPS+1], fst_numvals[NPOPS+1][NPOPS+1];
  	double fst_ave[NPOPPAIRS];
  	double sumsquaredif[NPOPPAIRS]; //for calculating variance

		for ( ipop = 0; ipop <= NPOPS; ++ipop ) {
			nall1[ ipop ] = NULL;
			nall2[ ipop ] = NULL;
			for ( jpop = 0; jpop <= NPOPS; ++jpop ) {
				fst_vals[ ipop ][ jpop ] = 0;
				fst_numvals[ ipop ][ jpop ] = 0;
			}
		}
		for ( int pp = 0; pp < NPOPPAIRS; ++pp ) {
			fst_ave[ pp ] = 0;
			sumsquaredif[ pp ] = 0;
		}

    // if (argc != 4) {
    //     fprintf(stderr, "Usage: calc_Fst_cosi <filebase> <# replicates> <outfilename>\n"); //same args as calc_pop_stats_cosi
    //     exit(0);
    //   }

    // strcpy(filebase, argv[1]);
    // numReps = atoi(argv[2]);
    // strcpy(outfilename, argv[3]);
    double fst_reps[NPOPPAIRS][numReps]; //(1, 3), (1, 4), (1, 5), (3, 4), (3, 5), (4, 5)

    outf = stdout; //fopen(outfilename, "w");
    assert(outf != NULL);

	/*******************
    LOOP OVER REPLICATES
    ********************/

	for (irep = 0; irep < numReps; irep++){

		/****************
    	GET GENOTYPE DATA
    	******************/

		for (ipop = 1; ipop <= NPOPS; ipop++){
	    	if (ipop-1 == excludePopIdx){continue;} //dummy population

		    // strcpy(basefile, filebase);
	      //   sprintf(repstr, "%d", irep);
	      //   strcat(basefile, repstr);

				get_coal_data_from_cosi(&data, irep, ipop);    

	    	if (data.nsample == 0) {continue;}
				nall1[ipop] = (int *)calloc(data.nsnp, sizeof(int));
				nall2[ipop] = (int *)calloc(data.nsnp, sizeof(int));
			for (isnp = 0; isnp < data.nsnp; isnp++) {
			    for (isamp = 0; isamp < data.nsample; isamp++) {
					if (data.genotype[isamp][isnp] == 1) {nall1[ipop][isnp]++;}
					else if (data.genotype[isamp][isnp] == 2) {nall2[ipop][isnp]++;}
			    } //end sample loop
			} //end snp loop

		//debugging
		//fprintf(stderr, "Replicate: %d\tPopulation: %d\tNSnps: %d\n", irep, ipop, data.nsnp);
		nsnp = data.nsnp;
		//free_coal_data(&data);
		} //end population loop

		
		/*****************************
	    CALC PAIRWISE FST AT EACH SITE
		*****************************/

	  	for (isnp = 0; isnp < nsnp; isnp++) {
	    	for (ipop = 1; ipop <= NPOPS; ipop++) {
	    		if (ipop-1 == excludePopIdx){continue;} //dummy pop
	      		if (nall1[ipop] == NULL) {continue;}
	      		nai[0] = nall1[ipop][isnp];
	      		nai[1] = nall2[ipop][isnp];
	      		ni = nai[0] + nai[1];
	      		p[0] = (double) nai[0] / ni;
	      		for (jpop = ipop+1; jpop <= NPOPS; jpop++) {
	      			if (jpop-1 == excludePopIdx){continue;} //dummy pop
					if (nall1[jpop] == NULL) {continue;}
					naj[0] = nall1[jpop][isnp];
					naj[1] = nall2[jpop][isnp];
					// na_both[0] = nai[0] + naj[0];
					// na_both[1] = nai[1] + naj[1];
					nj = naj[0] + naj[1];
					p[1] = (double) naj[0] / nj;
					if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) {continue;}
	        
					if (ni > 0 && nj > 0) {
		  				// Weir-Hill estimator
		  				pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
					  	nic = ni - (double) ni * ni / (ni + nj);
					  	njc = nj - (double) nj * nj / (ni + nj);
					  	nc = nic + njc;
					  	msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
					  	msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
					  	num = msp - msg;
					  	denom = msp + (nc - 1) * msg;
//					  	fwh = -99.;
					  	if (denom != 0) {
					    	fst_vals[ipop][jpop] += (num/denom);
					    	fst_numvals[ipop][jpop]++;
//					    	fwh = num / denom; 
					  	} // end if denom != 0
					}
	      		} // jpop
	    	} // ipop
	    } // end snp loop

	    //record Fst values for this replicate (flattening 2d array -> 1d; Steve's convention vs mine)
	    ipopPair = -1;
		for (ipop = 1; ipop <= NPOPS; ipop++) {
			if (nall1[ipop] == NULL) {continue;}
			if (ipop-1 == excludePopIdx){continue;} //dummy pop
			for (jpop = ipop+1; jpop <= NPOPS; jpop++) {
				if (jpop-1 == excludePopIdx){continue;} //dummy pop
				if (nall1[jpop] == NULL) {continue;}
				ipopPair++;
				fst_reps[ipopPair][irep] = fst_vals[ipop][jpop] / fst_numvals[ipop][jpop];
			} //jpop
       }//ipop

  	} //end replicates

  	//get average Fst for each pop pair across all replicates
  	for (ipopPair = 0; ipopPair < NPOPPAIRS; ipopPair++){
  		fst_total = 0;
  		for (irep = 0; irep < numReps; irep++){
			fst_total += fst_reps[ipopPair][irep];
  		}
  		fst_ave[ipopPair] = fst_total / numReps;
  	}

  	//get variance for each pop pair
  	for (ipopPair = 0; ipopPair < NPOPPAIRS; ipopPair++){
	  	for (irep = 0; irep < numReps; irep++){
  			sumsquaredif[ipopPair] += ((fst_ave[ipopPair] - fst_reps[ipopPair][irep])*(fst_ave[ipopPair] - fst_reps[ipopPair][irep]));
  		}
  	}

    //write to file and close file
		ipopPair = 0;
		for ( ipop = 1; ipop <= NPOPS; ++ipop ) {
			if ( ipop-1 == excludePopIdx ) continue;
			for ( jpop = ipop+1; jpop <= NPOPS; ++jpop ) {
				if ( jpop-1 == excludePopIdx ) continue;
				fprintf(outf, "%d, %d:\t%.15f\t%.15f\n", ipop, jpop, fst_ave[ipopPair], sumsquaredif[ipopPair]/numReps);
				++ipopPair;
			}
		}
//fclose(outf);
//return 0;
}
