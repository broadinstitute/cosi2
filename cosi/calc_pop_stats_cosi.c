// Updated 03/07/15
// Now filters MAF > .2 for P(D' = 1)
// but uses all min_minor >= 3 for r2
// Updated 03/30/15 
// Now takes seqlen as a parameter (was prev. hard-coded; --> wonky pi values!)

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_cosi.h"
#define NPOPS 5

int main(int argc, char **argv){
    //const int maxdist = 70000; //only calculate LD for 70kb
    const int min_minor = 3; //don't calculate LD at singletons or doubletons
    const double min_freq = .2; 
    int seqlen; //must be consistent with cosi param file
    char outfilename[264];
    char filebase[100];
    char basefile[999];
    char repstr[3];
    coal_data data;
    //char filename[264];
    FILE *outf=NULL;
    int isnp, jsnp, isamp, ns, nderiv, nmarker;
    int ipop, irep, jbin;
    double pi, pi_sum=0.;
    double dprime, ddenom, r2denom, r2, genDist;
    int hap[2][2], dist, ai, aj, *nminor, *nmajor;
    int nanc, npoly, numReps;
    const int nhist = 6;
    const int nDistHist = 14;
    const int nGenDistHist = 17;
    int *mafHist, *ancHist, *physDistHist, *start_distBins, *end_distBins, *compLDHist, *genDistHist, ibin;
    double *start_fbin, *end_fbin, *start_genDistBins, *end_genDistBins, *r2sums, freq;
    double **reps_sfs, **reps_p_anc, **reps_r2, **reps_dprime, **reps_pi;
    double rep_sfs, rep_p_anc, rep_r2, rep_dprime, rep_pi;
    double *sfs_total, *sfs_ave, *sfs_var_total, *sfs_var;
    double *p_anc_total, *p_anc_ave, *p_anc_var_total, *p_anc_var;
    double *r2_total, *r2_ave, *r2_var_total, *r2_var;
    double *dprime_total, *dprime_ave, *dprime_var_total, *dprime_var;
    double pi_total, pi_ave, pi_var_total, pi_var;
    double p_a, p_b, q_a, q_b; //allele frequencies
    double p_00, p_11, p_01, p_10; //haplotype frequencies
    double d;
    double freqi, freqj;
    
    if (argc != 5) {
        fprintf(stderr, "Usage: ./calc_pop_stats_cosi <filebase> <# replicates> <outfilename> <seqlen>\n");
        exit(0);
      }

    strcpy(filebase, argv[1]);
    //fprintf(stderr, "filebase (arg): \t");
    //fprintf(stderr, filebase);
    numReps = atoi(argv[2]);
    strcpy(outfilename, argv[3]);
    seqlen = atoi(argv[4]);

    /****************************************
    MEMORY ALLOCATION AND ARRAY INITIALIZATION
    *****************************************/
    
    //Freq spectrum and P(ancestral|freq bin)
    mafHist = malloc(nhist * sizeof(int));
    ancHist = malloc(nhist * sizeof(int));
    start_fbin = malloc(nhist * sizeof(double));
    end_fbin = malloc(nhist * sizeof(double));
    start_fbin[0] = end_fbin[0] = -1;   // singletons
    start_fbin[1] = 1e-9;
    start_fbin[2] = end_fbin[1] = .1;
    start_fbin[3] = end_fbin[2] = .2;
    start_fbin[4] = end_fbin[3] = .3;
    start_fbin[5] = end_fbin[4] = .4;
    end_fbin[5] = .5;

    //r2 vs. phys distance
    physDistHist = malloc(nDistHist * sizeof(int));
    r2sums = malloc(nDistHist * sizeof(double));
    start_distBins = malloc(nDistHist * sizeof(int));
    end_distBins = malloc(nDistHist * sizeof(int));
    start_distBins[0] = 0;
    start_distBins[1] = end_distBins[0] = 5000;
    start_distBins[2] = end_distBins[1] = 10000;
    start_distBins[3] = end_distBins[2] = 15000;
    start_distBins[4] = end_distBins[3] = 20000;
    start_distBins[5] = end_distBins[4] = 25000;
    start_distBins[6] = end_distBins[5] = 30000;
    start_distBins[7] = end_distBins[6] = 35000;
    start_distBins[8] = end_distBins[7] = 40000;
    start_distBins[9] = end_distBins[8] = 45000;
    start_distBins[10] = end_distBins[9] = 50000;
    start_distBins[11] = end_distBins[10] = 55000;
    start_distBins[12] = end_distBins[11] = 60000;
    start_distBins[13] = end_distBins[12] = 65000;
    end_distBins[13] = 70000;

    //D'=1 vs gen distance
    compLDHist = malloc(nGenDistHist * sizeof(int));
    genDistHist = malloc(nGenDistHist * sizeof(int));
    start_genDistBins = malloc(nGenDistHist * sizeof(double));
    end_genDistBins = malloc(nGenDistHist * sizeof(double));
    start_genDistBins[0] = 0;
    start_genDistBins[1] = end_genDistBins[0] = .001;
    start_genDistBins[2] = end_genDistBins[1] = .002;
    start_genDistBins[3] = end_genDistBins[2] = .003;
    start_genDistBins[4] = end_genDistBins[3] = .004;
    start_genDistBins[5] = end_genDistBins[4] = .005;
    start_genDistBins[6] = end_genDistBins[5] = .006;
    start_genDistBins[7] = end_genDistBins[6] = .007;
    start_genDistBins[8] = end_genDistBins[7] = .008;
    start_genDistBins[9] = end_genDistBins[8] = .009;
    start_genDistBins[10] = end_genDistBins[9] = .01;
    start_genDistBins[11] = end_genDistBins[10] = .011;
    start_genDistBins[12] = end_genDistBins[11] = .012;
    start_genDistBins[13] = end_genDistBins[12] = .013;
    start_genDistBins[14] = end_genDistBins[13] = .014;
    start_genDistBins[15] = end_genDistBins[14] = .015;
    start_genDistBins[16] = end_genDistBins[15] = .016;
    end_genDistBins[16] = .017;
    
    /************
     OUTFILE PREP
     ************/
    outf = fopen(outfilename, "w");
    assert(outf != NULL);

    /*********************************
     LOOP ALL ANALYSIS OVER POPULATIONS
     *********************************/
    for (ipop = 1; ipop <= NPOPS; ipop++){
        if (ipop == 2){continue;} //ANI; dummy population

        //reset arrays to sum up summary statistics across replicates and get average, sd
        reps_sfs = malloc(nhist * sizeof(double*));
        reps_p_anc = malloc(nhist * sizeof(double*));
        reps_r2 = malloc(nDistHist * sizeof(double*));
        reps_dprime = malloc(nGenDistHist * sizeof(double*));
        reps_pi = malloc(NPOPS * sizeof(double*));
        
        for (ibin = 0; ibin < nhist; ibin++){
            reps_sfs[ibin] = calloc(numReps, sizeof(double));
            reps_p_anc[ibin] = calloc(numReps, sizeof(double));
        }
        
        for (ibin = 0; ibin < nDistHist; ibin++){
            reps_r2[ibin] = calloc(numReps, sizeof(double));
        }
        
        for (ibin = 0; ibin < nGenDistHist; ibin++){
            reps_dprime[ibin] = calloc(numReps, sizeof(double));
        }
        
        for (ibin = 0; ibin < (NPOPS); ibin++){
            reps_pi[ibin] = calloc(numReps, sizeof(double));
        }
        
        /*********************************
        LOOP OVER REPLICATES FOR A POP
        *********************************/
                         
        for (irep = 0; irep < numReps; irep++){
            
            /**********************************
             DATA ANALYSIS - ALLELE FREQUENCIES
             **********************************/
            //initialize sums to zero
            
            pi_sum = 0;
            
            //reset arrays for replicate
            for (ibin = 0; ibin < nhist; ibin++){
                mafHist[ibin] = 0;
                ancHist[ibin] = 0;
            }
            
            for (ibin = 0; ibin < nDistHist; ibin++){
                physDistHist[ibin] = 0;
                r2sums[ibin] = 0;
            }
            
            for (ibin = 0; ibin < nGenDistHist; ibin++){
                compLDHist[ibin] = 0;
                genDistHist[ibin] = 0;
            }
            
            //fprintf(stderr, basefile);
            //fprintf(stderr,"\n");
            //fprintf(basefile, "");
            strcpy(basefile, argv[1]);
            //fprintf(stderr, basefile);
                       // fprintf(stderr,"\n");
            sprintf(repstr, "%d", irep);
            strcat(basefile, repstr);
            //fprintf(stderr, basefile);
              //          fprintf(stderr,"\n");
            
            //fprintf(stderr, "getting coal data with basefile:\t");
            //fprintf(stderr, basefile);
            //fprintf(stderr, "\n");
            get_coal_data(&data, basefile, ipop);            
            //fprintf(stderr, "pop: %d  nsnp: %d\n  nsample: %d\n", ipop, data.nsnp, data.nsample);
            if (data.nsample == 0) {continue;}

            
            nminor = calloc(data.nsnp, sizeof(int));
            nmajor = calloc(data.nsnp, sizeof(int));
            nmarker = npoly = 0;
            
            //loop over SNPs
            for (isnp = 0; isnp < data.nsnp; isnp++) {
                //count anc, der
                if (data.nallele[isnp] != 2) {fprintf(stderr, "bloop bloop"); continue;}
                nderiv = 0;
                ns = 0;
                for (isamp = 0; isamp < data.nsample; isamp++) {
                    if (data.genotype[isamp][isnp] != 0) {
                        ns++;
                        if (data.genotype[isamp][isnp] != data.anc_base[isnp]) {nderiv++;}
                    }
                }
                //fprintf(stderr, "%d ", nderiv);
                //fprintf(stderr, "data.nsample is %d\n", data.nsample);
                //fprintf(stderr, "ns is %d\n", ns);
                assert(ns == data.nsample);
                nanc = ns - nderiv;
                //fprintf(stderr, "%d\t", nanc);
                //update pi_sum
                nminor[isnp] = (nderiv > nanc) ? nanc : nderiv;
                nmajor[isnp] = (nderiv > nanc) ? nderiv : nanc;
                if (nminor[isnp] != 0 && nmajor[isnp] != 0) {npoly++;}
                pi = 2. * nderiv * nanc / ns / (ns-1);
                pi_sum += pi;
                //fprintf(stderr, "%f\t", pi);
                nmarker++;
                
                //if polymorphic, record minor allele frequency, counting ancestrals
                if (nminor[isnp] > 0) {
                    freq = (double) nminor[isnp] / (nminor[isnp] + nmajor[isnp]);//(double) nderiv / ns;
                    if (nderiv == 1 || nanc == 1) { //singleton
                        mafHist[0]++;
                    }
                    else {
                        for (ibin = 1; ibin < nhist; ibin++) {
                            if (freq >= start_fbin[ibin] && freq < end_fbin[ibin]) {
                                mafHist[ibin]++;
                                break;
                            }
                        }
                    }
                    if (nanc < nderiv){ //if the minor allele is ancestral, update ancHist
                        
                        if (nanc == 1){ancHist[0]++;} //ancestral singleton
                        else
                        {
                            freq = (double) nminor[isnp] / (nminor[isnp] + nmajor[isnp]);
                            
                            for (ibin = 1; ibin < nhist; ibin++) {
                                if (freq >= start_fbin[ibin] && freq < end_fbin[ibin]) {
                                    ancHist[ibin]++;
                                    break;
                                }
                            }
                        }
                    }
                }// end if polymorphic
            } // end loop over snps
            //fprintf(stderr, "finished snp loop\n");

            //calculate summary statistics for this replicate
            for (ibin = 0; ibin < nhist; ibin++){
                //fprintf(stderr, "%d\t", mafHist[ibin]);
                rep_sfs = ((double) mafHist[ibin] / npoly);
                if (!isnan(rep_sfs)){
                    reps_sfs[ibin][irep] = rep_sfs;
                }

                rep_p_anc = ((double) ancHist[ibin] / (double) mafHist[ibin]);
                if (!isnan(rep_p_anc)){
                    reps_p_anc[ibin][irep] = rep_p_anc;
                }        
            }
            rep_pi = (pi_sum / (double)seqlen);
            //fprintf(stderr, "rep_pi: %f\n", rep_pi);
            if (!isnan(rep_pi)){
            reps_pi[ipop-1][irep] = rep_pi;}

            /**************************************
             DATA ANALYSIS - LINKAGE DISEQUILIBRIUM
             ***************************************/
            
            //outer loop over SNPs
            for (isnp = 0; isnp < data.nsnp; isnp++) {
                if (nminor[isnp] < min_minor){continue;}
                //inner loop over SNPs
                for (jsnp = isnp+1; jsnp < data.nsnp; jsnp++) {
                    if (nminor[jsnp] < min_minor) {continue;}
                    
                    dist = data.pos[jsnp] - data.pos[isnp];
                    //if (dist > maxdist) {break;}
                    genDist = getGenDist(&data, data.pos[isnp], data.pos[jsnp]);
                    //fprintf(stderr, "%d\t%d\t%f\n", data.pos[isnp], data.pos[jsnp], genDist);
                    
                    //loop over samples at these two SNPs and count haplotypes
                    hap[0][0] = hap[0][1] = hap[1][0] = hap[1][1] = 0;
                    for (isamp = 0; isamp < data.nsample; isamp++) {
                        ai = data.genotype[isamp][isnp];
                        aj = data.genotype[isamp][jsnp];
                        hap[ai-1][aj-1]++;
                    }
                    
                    p_00 = ((double)hap[0][0]/data.nsample);
                    p_01 = ((double)hap[0][1]/data.nsample);
                    p_10 = ((double)hap[1][0]/data.nsample);
                    p_11 = ((double)hap[1][1]/data.nsample);
                    
                    p_a = ((double)(hap[0][0] + hap[0][1])/data.nsample);//locus a, indexed by isnp
                    q_a = ((double)(hap[1][0] + hap[1][1])/data.nsample);//locus a, indexed by isnp
                    p_b = ((double)(hap[0][0] + hap[1][0])/data.nsample);//locus b, indexed by jsnp
                    q_b = ((double)(hap[0][1] + hap[1][1])/data.nsample);//locus b, indexed by isnp
                    
                    //only continue if both sites polymorphic in this pop
                    if (((p_a != 0) && (q_a !=0)) && ((p_b !=0) && (q_b !=0))){
                        d = (p_00 * p_11) - (p_01 * p_10);
                        r2denom = (p_a * q_a * p_b * q_b);
                        r2 = d * d / r2denom;
                        //fprintf(stderr, "%f\t", r2);
                        //d pos
                        if (d > 0.){
                            ddenom = ((q_a * p_b) > (p_a * q_b)) ? (p_a * q_b) : (q_a * p_b);
                        }
                        
                        //d neg
                        else{
                            ddenom = ((-p_a * p_b) > (-q_a*q_b)) ? (-p_a * p_b) : (-q_a*q_b);
                        }
                        
                        dprime = d / ddenom;
                        //fprintf(stderr, "%f\n", dprime);

                        freqi = (double)nminor[isnp]/data.nsample;
                        freqj = (double)nminor[jsnp]/data.nsample;
                        if (freqi > min_freq && freqj > min_freq){ 
                            for (ibin = 0; ibin < nGenDistHist; ibin++) {
                                if (genDist >= start_genDistBins[ibin] && genDist < end_genDistBins[ibin]) {
                                    genDistHist[ibin]++;
                                    if (dprime == 1) {compLDHist[ibin]++;}
                                    break;
                                }
                            }
                        }//end if thisfreq

                        //update values for r2 vs phys dist
                        for (ibin = 0; ibin < nDistHist; ibin++) {
                            if (dist >= start_distBins[ibin] && dist < end_distBins[ibin]) {
                                r2sums[ibin] += r2;
                                physDistHist[ibin]++;
                                break;
                            }
                        }
                        
                    } //end if polymorphic
                } // jsnp loop
            }  // end isnp loop
            
            //calculate summary statistics for this replicate
            for (ibin = 0; ibin < nDistHist; ibin++){
                rep_r2 = (r2sums[ibin] / (double) physDistHist[ibin]);
                if (!isnan(rep_r2)){
                reps_r2[ibin][irep] = rep_r2;}
            }
            
            for (ibin = 0; ibin < nGenDistHist; ibin++){
                rep_dprime = ((double) compLDHist[ibin] / (double) genDistHist[ibin]);
                if (!isnan(rep_dprime)){
                reps_dprime[ibin][irep] = rep_dprime;}

            }
            
            
            free(nminor);
            free(nmajor);
            free_coal_data(&data);
        } //end replicate
    
    
        /************************************
         CALC SUMMARY STATISTICS MEAN, VARIANCE
         ************************************/
        
        //PI
        pi_total = 0;
        pi_ave = 0;
        pi_var_total = 0;
        pi_var = 0;
        for (irep = 0; irep < numReps; irep++){
            pi_total += reps_pi[ipop-1][irep];
        }
        pi_ave = pi_total / (double) numReps;
        for (irep = 0; irep < numReps; irep++){
            pi_var_total += ((pi_total - reps_pi[ipop-1][irep]) * (pi_total - reps_pi[ipop-1][irep]));
        }
        pi_var = pi_var_total / (double) numReps;
        
        
        //SFS
        sfs_total = calloc(nhist, sizeof(double));
        sfs_ave = calloc(nhist, sizeof(double));
        sfs_var_total = calloc(nhist, sizeof(double));
        sfs_var = calloc(nhist, sizeof(double));
        for (ibin = 0; ibin < nhist; ibin++) {
            for (jbin = 0; jbin < numReps; jbin++){
                sfs_total[ibin] += reps_sfs[ibin][jbin];
            }
        }
        for (ibin = 0; ibin < nhist; ibin++){
            sfs_ave[ibin] = sfs_total[ibin]/(double)numReps;
        }
        for (ibin = 0; ibin < nhist; ibin++){
            for (irep = 0; irep < numReps; irep++){
                sfs_var_total[ibin] += ((sfs_ave[ibin] - reps_sfs[ibin][irep]) * (sfs_ave[ibin] - reps_sfs[ibin][irep]));
            }
        }
        for (ibin = 0; ibin < nhist; ibin++){
            sfs_var[ibin] = sfs_var_total[ibin] / (double) numReps;
        }

        //P(ANCESTRAL|FREQ)
        p_anc_total = calloc(nhist, sizeof(double));
        p_anc_ave = calloc(nhist, sizeof(double));
        p_anc_var_total = calloc(nhist, sizeof(double));
        p_anc_var = calloc(nhist, sizeof(double));
        for (ibin = 0; ibin < nhist; ibin++) {
            for (irep = 0; irep < numReps; irep++){
                p_anc_total[ibin] += reps_p_anc[ibin][irep];
            }
        }
        for (ibin = 0; ibin < nhist; ibin++){
            p_anc_ave[ibin] = p_anc_total[ibin]/(double)numReps;
        }
        for (ibin = 0; ibin < nhist; ibin++){
            for (irep = 0; irep < numReps; irep++){
                p_anc_var_total[ibin] += ((p_anc_ave[ibin] - reps_p_anc[ibin][irep]) * (p_anc_ave[ibin] - reps_p_anc[ibin][irep]));
            }
        }
        for (ibin = 0; ibin < nhist; ibin++){
            p_anc_var[ibin] = p_anc_var_total[ibin] / (double) numReps;
        }
        
        //R2 VS PHYS DIST
        r2_total = calloc(nDistHist, sizeof(double));
        r2_ave = calloc(nDistHist, sizeof(double));
        r2_var_total = calloc(nDistHist, sizeof(double));
        r2_var = calloc(nDistHist, sizeof(double));
        for (ibin = 0; ibin < nDistHist; ibin++) {
            for (jbin = 0; jbin < numReps; jbin++){
                r2_total[ibin] += reps_r2[ibin][jbin];
            }
        }
        for (ibin = 0; ibin < nDistHist; ibin++){
            r2_ave[ibin] = r2_total[ibin]/(double)numReps;
        }
        for (ibin = 0; ibin < nDistHist; ibin++){
            for (irep=0; irep < numReps; irep++){
                r2_var_total[ibin] += ((r2_ave[ibin] - reps_r2[ibin][irep]) * (r2_ave[ibin] - reps_r2[ibin][irep]));
            }
        }
        for (ibin = 0; ibin < nDistHist; ibin++){
            r2_var[ibin] = r2_var_total[ibin] / (double) numReps;
        }
        
        //DPRIME VS GEN DIST
        dprime_total = calloc(nGenDistHist, sizeof(double));
        dprime_ave = calloc(nGenDistHist, sizeof(double));
        dprime_var_total = calloc(nGenDistHist, sizeof(double));
        dprime_var = calloc(nGenDistHist, sizeof(double));
        for (ibin = 0; ibin < nGenDistHist; ibin++) {
            for (jbin = 0; jbin < numReps; jbin++){
                dprime_total[ibin] += reps_dprime[ibin][jbin];

            }
        }
        for (ibin = 0; ibin < nGenDistHist; ibin++){
            dprime_ave[ibin] = dprime_total[ibin]/(double)numReps;

        }
        for (ibin = 0; ibin < nGenDistHist; ibin++){
            for (irep = 0; irep < numReps; irep++){
                dprime_var_total[ibin] += ((dprime_ave[ibin] - reps_dprime[ibin][irep]) * (dprime_ave[ibin] - reps_dprime[ibin][irep]));
            }
        }
        for (ibin = 0; ibin < nGenDistHist; ibin++){
            dprime_var[ibin] = dprime_var_total[ibin] / (double) numReps;
        }
        
        
        for (ibin = 0; ibin < nGenDistHist; ibin++) {
            //fprintf(stderr, "%f\t", dprime_total[ibin]);
            //fprintf(stderr, "%f\t", dprime_ave[ibin]);
            //fprintf(stderr, "%f\t", dprime_var[ibin]);
            //fprintf(stderr, "%f\t\n", dprime_var_total[ibin]);
        }
        
        
        /************************************
         PRINT SUMMARY STATISTICS TO OUTFILE
         ************************************/
        fprintf(outf, "%d\n", ipop);
        fprintf(outf, "%.10e\n", pi_ave);           //pi
        for (ibin = 0; ibin < nhist; ibin++) {                          //SFS
            fprintf(outf, "%.10f\t", sfs_ave[ibin]);
        }
        fprintf(outf, "\n");
                    
        for (ibin = 0; ibin < nhist; ibin++) {                          //p(anc)
            fprintf(outf, "%.10f\t", p_anc_ave[ibin]);
        }
        fprintf(outf, "\n");

                    for (ibin = 0; ibin < nDistHist; ibin++) {                      //ave r2
            fprintf(outf, "%.10f\t", r2_ave[ibin]);
        }
        fprintf(outf, "\n");
        for (ibin = 0; ibin < nGenDistHist; ibin++) {                   //dprime
            fprintf(outf, "%.10f\t", dprime_ave[ibin]);
        }
        fprintf(outf, "\n");
                    
        fprintf(outf, "%.10f\n", pi_var);
        for (ibin = 0; ibin < nhist; ibin++) {                          //SFS
            fprintf(outf, "%.10f\t", sfs_var[ibin]);
        }
        fprintf(outf, "\n");
        for (ibin = 0; ibin < nhist; ibin++) {                          //p(anc)
            fprintf(outf, "%.10f\t", p_anc_var[ibin]);
        }
        fprintf(outf, "\n");
        
        for (ibin = 0; ibin < nDistHist; ibin++) {                          //r2
            fprintf(outf, "%.10f\t", r2_var[ibin]);
        }
        fprintf(outf, "\n");
        for(ibin = 0; ibin < nGenDistHist; ibin++){
            fprintf(outf, "%.10f\t", dprime_var[ibin]);
        }
        fprintf(outf, "\n");
                    
    }//end loop over pops
    

  return 0;
}
