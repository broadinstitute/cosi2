#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <cosi_rand/random.h>

int main(int argc, char **argv) {
  FILE *inf=NULL;
  int itoken, ichr, isite, *long_enough;
  const int line_size = 256;
  char filename[256], *new_line=NULL, *running=NULL, *token=NULL, chrstr[8], 
    command[128], seqstr[10];
  long min_pos, max_pos, *chrlen, seqtot=0, seqlen, sum=0;
  double x;
  
  if (argc != 4) {
    fprintf(stderr, "Usage: recomap_hapmap2 <1kG map directory> <sequence length (bp)> <output file>\n");
    return 0;
  }
  seqlen = strtol(argv[2], NULL, 10);
  chrlen = calloc(24, sizeof(long));
  long_enough = calloc(24, sizeof(int));
  fprintf(stdout, "random seed for choosing chromosome: %lu\n", seed_rng());
  new_line = malloc((line_size + 1) * sizeof(char));
  for (ichr = 1; ichr <= 22; ichr++) {
    strcpy(filename, argv[1]);
    strcat(filename, "/genetic_map_chr");
    sprintf(chrstr, "%d", ichr);
    strcat(filename, chrstr);
    strcat(filename, "_b36.txt");
    inf = fopen(filename, "r");
    if (inf == NULL) {
      fprintf(stderr, "Could not open map file %s\n", filename);
      return 0;
    }
    isite = max_pos = 0;
    fgets(new_line, line_size, inf);
    while (fgets(new_line, line_size, inf) != NULL) {
      for (itoken = 0, running = new_line; (token = strsep(&running, " \t")) != NULL; itoken++) {
	if (itoken == 0) {
	  if (isite == 0) {min_pos = strtol(token, NULL, 10);}
	  else {max_pos = strtol(token, NULL, 10);}
	  isite++;
	  break;
	}
      }
    }    
    if (max_pos - min_pos + 1 >= seqlen) {
      chrlen[ichr] = max_pos - min_pos + 1;
      seqtot += chrlen[ichr];
    }
    fclose(inf);
  }
  if (seqtot == 0) {
    fprintf(stderr, "Requested sequence is longer than any single chromosome\n");
    return 0 ;
  }
  x = random_double();
  for (ichr = 1; ichr <= 22; ichr++) {
    sum += chrlen[ichr];
    if ((double) sum / seqtot > x) {
      fprintf(stdout, "using chrom %d\n", ichr);
      break;
    }
    fclose(inf);
  }
  sprintf(seqstr, "%ld", seqlen);
  strcpy(command, "../../bin/get_recomap ");
  strcat(command, filename);
  strcat(command, " ");
  strcat(command, seqstr);
  strcat(command, " ");
  strcat(command, argv[3]);
  system(command);

  return 0;
}
