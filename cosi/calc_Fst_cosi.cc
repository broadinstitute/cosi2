#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <cosi/coal_data_cosi.h>
//#define NPOPS 5 //2 is a dummy population 
#define NPOPPAIRS 6 //(1, 3), (1, 4), (1, 5), (3, 4), (3, 5), (4, 5)

void calc_Fst_main( int numReps, int NPOPS,
										void (*get_coal_data_from_cosi)( coal_data *, int /* irep */, int /* ipop */ ) ) {
    // char outfilename[264];
    // char basefile[264], filebase[264];
    // char repstr[3];
    // int numReps;
	coal_data data;
    int nsnp;
    int irep, ipop, jpop, isnp, isamp, ipopPair;
    int ni, nj;
    int nai[2], naj[2], na_both[2];
    double fst_total;
    int *nall1[NPOPS]={NULL}, *nall2[NPOPS]={NULL};
    FILE *outf=NULL;
    double msp, msg, p[2], nc, num, pmean, nic, njc, fwh, denom;
  	double fst_vals[NPOPS][NPOPS]={{0}}, fst_numvals[NPOPS][NPOPS] = {{0}};
  	double fst_ave[NPOPPAIRS];
  	double sumsquaredif[NPOPPAIRS]={0}; //for calculating variance


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
	    	if (ipop == 2){continue;} //dummy population

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
	    		if (ipop == 2){continue;} //dummy pop
	      		if (nall1[ipop] == NULL) {continue;}
	      		nai[0] = nall1[ipop][isnp];
	      		nai[1] = nall2[ipop][isnp];
	      		ni = nai[0] + nai[1];
	      		p[0] = (double) nai[0] / ni;
	      		for (jpop = ipop+1; jpop <= NPOPS; jpop++) {
	      			if (ipop == 2){continue;} //dummy pop
					if (nall1[jpop] == NULL) {continue;}
					naj[0] = nall1[jpop][isnp];
					naj[1] = nall2[jpop][isnp];
					na_both[0] = nai[0] + naj[0];
					na_both[1] = nai[1] + naj[1];
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
					  	fwh = -99.;
					  	if (denom != 0) {
					    	fst_vals[ipop][jpop] += (num/denom);
					    	fst_numvals[ipop][jpop]++;
					    	fwh = num / denom; 
					  	} // end if denom != 0
					}
	      		} // jpop
	    	} // ipop
	    } // end snp loop

	    //record Fst values for this replicate (flattening 2d array -> 1d; Steve's convention vs mine)
	    ipopPair = -1;
		for (ipop = 1; ipop <= NPOPS; ipop++) {
			if (nall1[ipop] == NULL) {continue;}
			if (ipop == 2){continue;} //dummy pop
			for (jpop = ipop+1; jpop <= NPOPS; jpop++) {
				if (ipop == 2){continue;} //dummy pop
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
  for (ipopPair = 0; ipopPair < NPOPPAIRS; ipopPair++){
  	if (ipopPair == 0){fprintf(outf, "1, 3:\t");}
  	if (ipopPair == 1){fprintf(outf, "1, 4:\t");}
  	if (ipopPair == 2){fprintf(outf, "1, 5:\t");}
  	if (ipopPair == 3){fprintf(outf, "3, 4:\t");}
  	if (ipopPair == 4){fprintf(outf, "3, 5:\t");}
  	if (ipopPair == 5){fprintf(outf, "4, 5:\t");}
  	fprintf(outf, "%.15f\t%.15f\n", fst_ave[ipopPair], sumsquaredif[ipopPair]/numReps);
  }

//fclose(outf);
//return 0;
}
