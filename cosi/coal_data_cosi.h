#ifndef INCLUDE_COAL_DATA_COSI_H
#define INCLUDE_COAL_DATA_COSI_H

//
// * Struct coal_data
//
// Simulation output for one replica, for one population.
typedef struct coal_data {
	 // Number of samples from the population
    int nsample;
    char **samp_id;
	 // Number of SNPs in the replica
    int nsnp;
    char **snp_id;
	 // The position of each SNP
    int *pos;
	  double *gdPos;
    int *chrom;
    int **genotype;   // haploid
    int *snp_base[4];
    int *anc_base;
    //  int *snp_type;
    int *nallele;
    double *genPos; //for determining gen dist
    int *physPos; //parallel w above
    int nRecom; //# of lines in recombination file, as nsample or nsnp
    
} coal_data;

#define INTERGENIC 0
#define NONSYNON 1
#define SYNON 2
#define INTRONIC 3
#define OTHER 4

void get_coal_data(coal_data* data, char filebase[], int pop);
void free_coal_data(coal_data* data);
double getGenDist(coal_data* data, int pos_i, int pos_j);

#endif // #ifndef INCLUDE_COAL_DATA_COSI_H
