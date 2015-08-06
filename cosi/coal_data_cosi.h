typedef struct coal_data {
    int nsample;
    char **samp_id;
    int nsnp;
    char **snp_id;
    int *pos;
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