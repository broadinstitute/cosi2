#include <cosi/coal_data_cosi.h>

namespace cosi {
namespace customstats {

void get_coal_data(coal_data* data, char filebase[], int pop) {
  const int line_size = 200000;
  char filename[999], recomfilename[999];
  FILE *inf=NULL;
  char *newLine, *token, *running, popstr[4];
  int isamp, isnp, itoken, iRecom;

  newLine = malloc((line_size+1) * sizeof(char));
  assert(newLine != NULL);
    
  sprintf(popstr, "%d", pop);

  data->nsample = 0;
  data->nsnp = 0;
  data->samp_id = NULL;
  data->pos = NULL;
  data->snp_id = NULL;
  data->chrom = NULL;
  data->genotype = NULL;
  data->anc_base = NULL;
  data->nallele = NULL;
  data->genPos = NULL;
  data->physPos = NULL;
  data->nRecom = 0; 

  // Count number of samples in file
  strcpy(filename, filebase);
  strcat(filename, ".hap-");
  strcat(filename, popstr);
  inf = fopen(filename, "r");  
  if (inf == NULL) {
      fprintf(stderr, "Missing hap file;\n");
      fprintf(stderr, filename);
      fprintf(stderr, "\n");
      return;
  }
  while (fgets(newLine, line_size, inf) != NULL) {
    assert(strlen(newLine) < line_size - 2);
    data->nsample++;
  }
  fclose(inf);

  // Count number of SNPs in file
  strcpy(filename, filebase);
  strcat(filename, ".pos-");
  strcat(filename, popstr);
  inf = fopen(filename, "r");
  if (inf == NULL) {
      fprintf(stderr, "Missing pos file;\n");
      fprintf(stderr, filename);
      fprintf(stderr, "\n");
      return;
  }
  fgets(newLine, line_size, inf); // header
  while (fgets(newLine, line_size, inf) != NULL) {
    data->nsnp++;
  }
  fclose(inf);

  //count number of lines in recombination file 
  strcpy(recomfilename, filebase);
  strcat(recomfilename, ".recom");
  inf = fopen(recomfilename, "r");
  if (inf == NULL) {
      fprintf(stderr, "Missing recom file;\n");
      fprintf(stderr, recomfilename);
      fprintf(stderr, "\n");
      return;
  }
    
  while (fgets(newLine, line_size, inf) != NULL) {
      assert(strlen(newLine) < line_size - 2);
      data->nRecom++;
  }
  fclose(inf);
    
    
  //allocate memory; initialize
  data->samp_id = malloc(data->nsample * sizeof(char*));
  data->genotype = malloc(data->nsample * sizeof(int*));
  for (isamp = 0; isamp < data->nsample; isamp++) {
    data->samp_id[isamp] = malloc(64 * sizeof(char));
    assert(data->samp_id[isamp] != NULL);
    data->genotype[isamp] = malloc(data->nsnp * sizeof(int));
  }
  data->pos = malloc(data->nsnp * sizeof(int));
  data->chrom = malloc(data->nsnp * sizeof(int));
  data->snp_id = malloc(data->nsnp * sizeof(char*));
  data->anc_base = malloc(data->nsnp * sizeof(int));
  //  data->snp_type = malloc(data->nsnp * sizeof(int));
  data->snp_base[0] = malloc(data->nsnp * sizeof(int));
  data->snp_base[1] = malloc(data->nsnp * sizeof(int));
  data->snp_base[2] = calloc(data->nsnp, sizeof(int));
  data->snp_base[3] = calloc(data->nsnp, sizeof(int));
  data->nallele = malloc(data->nsnp * sizeof(int));
  data->genPos = malloc(data->nRecom * sizeof(double));
  data->physPos = malloc(data->nRecom * sizeof(int));
  assert(data->genPos != NULL);
  assert(data->physPos != NULL);
    
  for (isnp = 0; isnp < data->nsnp; isnp++) {
    data->snp_id[isnp] = malloc(64 * sizeof(char));
    assert(data->snp_id[isnp] != NULL);
  }

    
  // Reopen snp file
  inf = fopen(filename, "r");
  if (inf == NULL) {
      fprintf(stderr, "Missing pos file;\n");
      fprintf(stderr, filename);
      fprintf(stderr, "\n");
      return;
  }
  assert(inf != NULL);
  fgets(newLine, line_size, inf); // header
  isnp = 0;
  while (fgets(newLine, line_size, inf) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {
	       strcpy(data->snp_id[isnp], token);
	       data->nallele[isnp] = 2;
      }	
      else if (itoken == 1) {
	       data->chrom[isnp] = atoi(token);
      }
      else if (itoken == 2) {
	       data->pos[isnp] = atoi(token);
      }
      else if (itoken == 3) {
	       data->snp_base[0][isnp] = atoi(token);
      }
      else if (itoken == 5) {
      	data->snp_base[1][isnp] = atoi(token);
      	data->anc_base[isnp] = atoi(token);
      	isnp++;
      	break;
      }
    } //end for running=newLine
  } // end while fgets(newLine)
  fclose(inf);

  //switch to hap file
  strcpy(filename, filebase);
  strcat(filename, ".hap-");
  strcat(filename, popstr);
  inf = fopen(filename, "r");
  if (inf == NULL) {
      fprintf(stderr, "Missing hap file;\n");
      fprintf(stderr, filename);
      fprintf(stderr, "\n");
      return;
  }

  isamp = 0;
  while (fgets(newLine, line_size, inf) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t ")) != NULL; itoken++) {
      if (itoken >= 2) {
	     if (token[0] == '\n') {break;}
      	isnp = itoken - 2;
      	assert(isnp < data->nsnp);
      	data->genotype[isamp][isnp] = strtol(token, NULL, 10);
      }
    }
    isamp++;
  }
  fclose(inf);
    
    //now let's do recomb file
    inf = fopen(recomfilename, "r"); // NB assumes same name file in same folder as program !!
    if (inf == NULL) {
        fprintf(stderr, "Missing recom file;\n");
        fprintf(stderr, recomfilename);
        fprintf(stderr, "\n");
        return;
    }

    assert(inf != NULL);
    iRecom = 0;
    while (fgets(newLine, line_size, inf) != NULL) {
        for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
            if (itoken == 0) {
                data->physPos[iRecom] = atoi(token);
                }
            else if (itoken == 1) {
                data->genPos[iRecom] = 100*atof(token); //Morgans -> centimorgans
                iRecom++;
                break;
                }
        } // end for running = newLine
    } // end while fgets(newLine)
    fclose(inf);
  free(newLine);
}

void free_coal_data(coal_data* data) {
  int isamp, isnp;
  if (data == NULL) {return;}
  if (data->samp_id == NULL) {return;}
  for (isamp = 0; isamp < data->nsample; isamp++) {
    free(data->genotype[isamp]);
    free(data->samp_id[isamp]);
  }
  for (isnp = 0; isnp < data->nsnp; isnp++) {
    free(data->snp_id[isnp]);
  }
  free(data->pos);
  free(data->chrom);
  free(data->snp_id);
  free(data->anc_base);
  //  free(data->snp_type);
  free(data->snp_base[0]);
  free(data->snp_base[1]);
  free(data->snp_base[2]);
  free(data->snp_base[3]);
  free(data->nallele);
  free(data->samp_id);
  free(data->genotype);
    free(data->genPos);
    free(data->physPos);
  data->nsnp = 0;
  data->nsample = 0;
}

double getGenDist(coal_data* data, int pos_i, int pos_j){
    int pointer_i = (data->nRecom - 1);
    int pointer_j = (data->nRecom - 1); //last bin
    int dist, dist_i, dist_j, interval, ibin, jbin;
    double totaldist;
    
    //find correct bin for each pointer
    for (ibin = 0; ibin < data->nRecom; ibin++){
        if (data->physPos[ibin] > pos_i){
            pointer_i = ( ibin - 1 );
            break;
        }
    }

    for (jbin = 0; jbin < data->nRecom; jbin++){
        if (data->physPos[jbin] > pos_j){
            pointer_j = ( jbin - 1 );
            break;
        }
    }

    if (pointer_i == pointer_j){
        dist = pos_j - pos_i;
        return (dist*data->genPos[pointer_i]);
    }

    else{
        dist_i = (data->physPos[pointer_i+1] - pos_i);
        dist_j = (pos_j - data->physPos[pointer_j]);
        totaldist = ((dist_i*data->genPos[pointer_i]) + (dist_j*data->genPos[pointer_j])); //bookends

        if ((dist_i - dist_j) == 1){return totaldist;}
        else{
          pointer_i++;
            while(pointer_i < pointer_j){
            interval = data->physPos[pointer_i + 1] - data->physPos[pointer_i];
            totaldist += (interval * data->genPos[pointer_i]);
                pointer_i++;
            }
        }

    return totaldist; 
    }
    return 0.;
}

} // namespace customstats
} // namespace cosi
