#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "coal_data_cosi.h"

int main(int argc, char **argv){
	char filebase[264], basefile[1056];
	char repstring[8];
	char pos1filename[1056], pos3filename[1056], pos4filename[1056], pos5filename[1056];
	char hap1filename[1056], hap3filename[1056], hap4filename[1056], hap5filename[1056];
	FILE *posfile1=NULL, *posfile3=NULL, *posfile4=NULL, *posfile5=NULL;
	FILE *hapfile1=NULL, *hapfile3=NULL, *hapfile4=NULL, *hapfile5=NULL;
	coal_data data1, data3, data4, data5; //2 is dummy pop for legacy reasons
	int numReps, irep, isamp, isnp, singletonRate, nderiv=0;
	int *takeSNPs; //0=yes 1=no.
	int randint;
	int d1, d3, d4, d5;
	int nSNPstotake = 0;
	int **haps1, **haps3, **haps4, **haps5;
	int jsnp;

	if (argc != 4) {
		fprintf(stderr, "Usage: ./alleleDropout_cosi <filebase e.g. '/sims/test/rep_'> <# replicates> <singletonrate>\n");
		exit(0);
	}

	strcpy(filebase, argv[1]);
	numReps = atoi(argv[2]);
	singletonRate = atoi(argv[3]);
	srand (time(NULL));

	for (irep = 0; irep < numReps; irep++){
		strcpy(basefile, filebase); 
		sprintf(repstring, "%d", irep);
		strcat(basefile, repstring);
		get_coal_data(&data1, basefile, 1); 
		get_coal_data(&data3, basefile, 3); 
		get_coal_data(&data4, basefile, 4); 
		get_coal_data(&data5, basefile, 5); 

		assert(data1.nsnp == data3.nsnp);
		assert(data1.nsnp == data4.nsnp);
		assert(data1.nsnp == data5.nsnp);
		assert(data1.nsample == data3.nsample);
		assert(data1.nsample == data4.nsample);
		assert(data1.nsample == data5.nsample);

		strcpy(pos1filename, basefile);
		strcat(pos1filename, "_ad.pos-1"); //should maybe add allele dropout rate to the name of the file?
		strcpy(pos3filename, basefile);
		strcat(pos3filename, "_ad.pos-3");
		strcpy(pos4filename, basefile);
		strcat(pos4filename, "_ad.pos-4");
		strcpy(pos5filename, basefile);
		strcat(pos5filename, "_ad.pos-5");

		posfile1 = fopen(pos1filename, "w");
		assert(posfile1 != NULL);
		posfile3 = fopen(pos3filename, "w");
		assert(posfile3 != NULL);
		posfile4 = fopen(pos4filename, "w");
		assert(posfile4 != NULL);
		posfile5 = fopen(pos5filename, "w");
		assert(posfile5 != NULL);

		fprintf(posfile1, "SNP\tCHROM\tCHROM_POS\tALLELE1\tFREQ1\tALLELE2\tFREQ2\n");
		fprintf(posfile3, "SNP\tCHROM\tCHROM_POS\tALLELE1\tFREQ1\tALLELE2\tFREQ2\n");
		fprintf(posfile4, "SNP\tCHROM\tCHROM_POS\tALLELE1\tFREQ1\tALLELE2\tFREQ2\n");
		fprintf(posfile5, "SNP\tCHROM\tCHROM_POS\tALLELE1\tFREQ1\tALLELE2\tFREQ2\n");

		takeSNPs = calloc(data1.nsnp, sizeof(int));
		nSNPstotake = 0;
		for (isnp = 0; isnp < data1.nsnp; isnp++) {
			nderiv = 0;
			d1 = 0;
			d3 = 0;
			d4 = 0;
			d5 = 0;
			for (isamp = 0; isamp < data1.nsample; isamp++) {
				if (data1.genotype[isamp][isnp] != data1.anc_base[isnp]) {d1++;}
				if (data3.genotype[isamp][isnp] != data3.anc_base[isnp]) {d3++;}
				if (data4.genotype[isamp][isnp] != data4.anc_base[isnp]) {d4++;}
				if (data5.genotype[isamp][isnp] != data5.anc_base[isnp]) {d5++;}
			}			
			nderiv = d1 + d3 + d4 + d5;
			randint = rand() % 100; // between 0 and 99 

			if ( ((nderiv == 1) || nderiv == (4*data1.nsample - 1))  && (randint >= singletonRate)){ // singleton dropout
				takeSNPs[isnp] = 1; 
			} 

			else { //write to file
				nSNPstotake++;
				fprintf(posfile1, data1.snp_id[isnp]);
				fprintf(posfile1, "\t%d\t%d\t%d\t%f\t%d\t%f\n", data1.chrom[isnp], data1.pos[isnp],data1.snp_base[0][isnp],(double)d1/data1.nsample,data1.snp_base[1][isnp], (double)(data1.nsample - d1)/data1.nsample);
				fprintf(posfile3, data3.snp_id[isnp]);
				fprintf(posfile3, "\t%d\t%d\t%d\t%f\t%d\t%f\n", data3.chrom[isnp], data3.pos[isnp],data3.snp_base[0][isnp],(double)d3/data3.nsample,data3.snp_base[1][isnp], (double)(data3.nsample - d3)/data3.nsample);
				fprintf(posfile4, data4.snp_id[isnp]);
				fprintf(posfile4, "\t%d\t%d\t%d\t%f\t%d\t%f\n", data4.chrom[isnp], data4.pos[isnp],data4.snp_base[0][isnp],(double)d4/data4.nsample,data4.snp_base[1][isnp], (double)(data4.nsample - d4)/data4.nsample);
				fprintf(posfile5, data5.snp_id[isnp]);
				fprintf(posfile5, "\t%d\t%d\t%d\t%f\t%d\t%f\n", data5.chrom[isnp], data5.pos[isnp],data5.snp_base[0][isnp],(double)d5/data5.nsample,data5.snp_base[1][isnp], (double)(data5.nsample - d5)/data5.nsample);
			}
		} //end isnp loop

		//fprintf(stderr, "There are %d snps to take, out of an original total of %d.\n", nSNPstotake, data1.nsnp);
		fclose(posfile1);
		fclose(posfile3);
		fclose(posfile4);
		fclose(posfile5);

		//rewrite hap files 
		haps1 = malloc(data1.nsample * sizeof(int*));
		haps3 = malloc(data3.nsample * sizeof(int*));
		haps4 = malloc(data4.nsample * sizeof(int*));
		haps5 = malloc(data5.nsample * sizeof(int*));

		for (isamp = 0; isamp < data1.nsample; isamp++){
			haps1[isamp] = calloc(nSNPstotake, sizeof(int));
			haps3[isamp] = calloc(nSNPstotake, sizeof(int));
			haps4[isamp] = calloc(nSNPstotake, sizeof(int));
			haps5[isamp] = calloc(nSNPstotake, sizeof(int));
		}

		jsnp = 0;
		for (isnp = 0; isnp < data1.nsnp; isnp++) {
			if (takeSNPs[isnp] == 0){
				for (isamp = 0; isamp < data1.nsample; isamp++){
					assert(jsnp < nSNPstotake);
					haps1[isamp][jsnp] = data1.genotype[isamp][isnp];
					haps3[isamp][jsnp] = data3.genotype[isamp][isnp];
					haps4[isamp][jsnp] = data4.genotype[isamp][isnp];
					haps5[isamp][jsnp] = data5.genotype[isamp][isnp];				
				}
				jsnp++;
			} //end takeSNP condition
		} //end isnp loop

		strcpy(hap1filename, basefile);
		strcat(hap1filename, "_ad.hap-1"); //should maybe add allele dropout rate to the name of the file?
		strcpy(hap3filename, basefile);
		strcat(hap3filename, "_ad.hap-3");
		strcpy(hap4filename, basefile);
		strcat(hap4filename, "_ad.hap-4");
		strcpy(hap5filename, basefile);
		strcat(hap5filename, "_ad.hap-5");

		hapfile1 = fopen(hap1filename, "w");
		assert(hapfile1 != NULL);
		hapfile3 = fopen(hap3filename, "w");
		assert(hapfile3 != NULL);
		hapfile4 = fopen(hap4filename, "w");
		assert(hapfile4 != NULL);
		hapfile5 = fopen(hap5filename, "w");
		assert(hapfile5 != NULL);


		for (isamp = 0; isamp < data1.nsample; isamp++){
			fprintf(hapfile1, "%d\t%d\t", isamp, data1.chrom[0]);
			fprintf(hapfile3, "%d\t%d\t", isamp, data3.chrom[0]);
			fprintf(hapfile4, "%d\t%d\t", isamp, data4.chrom[0]);
			fprintf(hapfile5, "%d\t%d\t", isamp, data5.chrom[0]);						
			for (isnp = 0; isnp < nSNPstotake; isnp++){
				fprintf(hapfile1, "%d ", haps1[isamp][isnp]);
				fprintf(hapfile3, "%d ", haps3[isamp][isnp]);				
				fprintf(hapfile4, "%d ", haps4[isamp][isnp]);				
				fprintf(hapfile5, "%d ", haps5[isamp][isnp]);
			}
			fprintf(hapfile1, "\n");
			fprintf(hapfile3, "\n");
			fprintf(hapfile4, "\n");
			fprintf(hapfile5, "\n");

		}


		fclose(hapfile1);		
		fclose(hapfile3);
		fclose(hapfile4);
		fclose(hapfile5);

		free_coal_data(&data1);
		free_coal_data(&data3);
		free_coal_data(&data4);
		free_coal_data(&data5);

		free(takeSNPs);

		for (isamp = 0; isamp < data1.nsample; isamp++){
			free(haps1[isamp]);
			free(haps3[isamp]);
			free(haps4[isamp]);
			free(haps5[isamp]);
		}
		free(haps1);
		free(haps3);
		free(haps4);
		free(haps5);

	} // end replicate loop
	return 0;
}