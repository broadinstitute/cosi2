/* Generate recombination rate map of a region, based on various models.  
   Output: file listing physical and genetic distances.*/
/* model 0 = regionally uniform, flat or drawn from histogram, 
    = model 0 + gamma-distributed hotspots */
/* Sequence numbering starts at 1. */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <cosi_rand/random.h>
#include <cosi_rand/gamma.h>

#define MAXBIN 200

int main(int argc, char **argv) {
  /* Program would require significant changes to handle larger spot sizes:
     1) include possibility of starting in the middle of a hotspot
     2) handle overlapping spots (including multiple overlaps)
  */
  char outname[256], distrname[256];
  FILE *parfile = NULL, *outfile = NULL, *distrfile = NULL;
  int nScan;
  char newLine[256], key[64], value[256];  
  int model = 0, size, nbin, ibin, spot_size = 1, local_size = 50000000, start_local;
  double bkgdFract = 0.1, baseRate = 1.0;
  double regRate = 0., x;
  double start[MAXBIN], end[MAXBIN], prob[MAXBIN];
  double meanInten, rspot, totDist, b, rhot, thisd;
  double shape = 1., local_shape = 0.3, intensity_shape = 0.3;
  double meanSpace = 9000., filler;
  unsigned long seed=0;

  if (argc != 3 && argc != 4) {
    printf("Usage: recosim <parameter file name> <size (bp)> [randomSeed]\n");
    printf("Parameter file entries (any order, all optional)\n");
    printf("  outfile <output file name>    [default=\"model.out\"]\n");
    printf("  model <0,1>    [default=0]\n");
    printf("     model 0: uniform recombination, constant or drawn from distribution.\n");
    printf("     model 1: model 0 + gamma-distributed hotspots.\n");
    printf("  baserate <mean recomb (cM/Mb)>    [default=1.0]\n");
    printf("  distribution <recomb distr. file name>    [default=none (const value)]\n");
    printf("     file format: bin_start bin_end cumulative_fraction\n");
    printf("  space <mean hotspot spacing (bp)>    [default=9000]\n");
    printf("  distance_shape <gamma function shape param>    [default=1.0]\n");
    printf("  intensity_shape <gamma function shape param>    [default=0.3]\n");
    printf("  local_shape <gamma function shape param, local variation>    [default=0.3]\n");
    printf("  local_size <size of region of local variation (bp) (e.g. 100000)>    [default=50000000]\n");
    printf("  bkgd <fraction in flat bkgd>    [default=0.1]\n");
    printf("  random_seed <unsignedinteger seed> (0=>picked by program based on time and PID) [default=0]\n");
    exit(EXIT_FAILURE);
  }

  if ( argc == 4 )
	 sscanf( argv[3], "%lu", &seed );
  
  strcpy(outname, "model.out");
  parfile = fopen(argv[1],"r");
  sscanf(argv[2], "%d", &size);
  if (parfile == NULL) {fprintf(stderr, "Error (recosim.c): Could not open parameter file %s.\n", argv[1]); exit(EXIT_FAILURE);}
  while (fgets(newLine, 256, parfile) != NULL) {
    if (strlen(newLine) == 1) {continue;}
    if (newLine[0] == '!' || newLine[0] == '#') {continue;}
    nScan = sscanf(newLine, "%s%s", key, value);
    if (nScan != 2) {fprintf(stderr, "Error (recosim.c): Illegal parameter file entry:\n%s\n", newLine); exit(EXIT_FAILURE);}
    if (!strcmp(key, "model")) {
      sscanf(value, "%d", &model);
    }
    else if (!strcmp(key, "outfile")) {
      strcpy(outname, value);
    }
    else if (!strcmp(key, "random_seed")) {
      sscanf(value, "%lu", &seed);

    }
    else if (!strcmp(key, "distribution")) {
      strcpy(distrname, value);
      distrfile = fopen(distrname, "r");
      if (distrfile == NULL) {fprintf(stderr, "Error (recosim.c): Could not open distribution file %s\n", distrname); exit(EXIT_FAILURE);}
    }
    else if (!strcmp(key, "baserate")) {
      sscanf(value, "%lf", &baseRate);
    }
    else if (!strcmp(key, "space")) {
      sscanf(value, "%lf", &meanSpace);
    }
    else if (!strcmp(key, "shape") || !strcmp(key, "distance_shape")) {
      sscanf(value, "%lf", &shape);
    }
    else if (!strcmp(key, "intensity_shape")) {
      sscanf(value, "%lf", &intensity_shape);
    }
    else if (!strcmp(key, "local_shape")) {
      sscanf(value, "%lf", &local_shape);
    }
    else if (!strcmp(key, "local_size")) {
      sscanf(value, "%d", &local_size);
    }
    else if (!strcmp(key, "bkgd")) {
      sscanf(value, "%lf", &bkgdFract);
    }
    else {
      fprintf(stderr, "Error (recosim.c): Unknown parameter name: %s\n", key);
      exit(EXIT_FAILURE);
    }
  }

  if (local_size > size) {local_size = size;}

  if (model > 0 || distrfile != NULL) {
    if (seed > 0) {
      set_rng_seed(seed);
    }
    else {
      seed = seed_rng();
    }
    printf("Recombination model seed: %ld\n", seed);
  }
  
  if (distrfile != NULL) {
    nbin = 0;
    while (fscanf(distrfile, "%lf%lf%lf", &start[nbin], &end[nbin], &prob[nbin]) == 3 && 
	   nbin < MAXBIN) {
      nbin++;
    }
    if (nbin == MAXBIN) {
      fprintf(stderr, "Warning (recosim.c): Too many recombination bins. Increase MAXBIN.\n");
    }
    /* Pick from this distribution */
    x = random_double();
    regRate = -1;
    for (ibin = 0; ibin < nbin; ibin++) {
      if (x <= prob[ibin]) {
	regRate = start[ibin] + random_double() * (end[ibin] - start[ibin]);
	break;
      }
    }
  }
  else {
    regRate = baseRate;
  }

  outfile = fopen(outname, "w");
  if (outfile == NULL) {fprintf(stderr, "Error (recosim.c): Unable to open output file %s\n", outname); exit(EXIT_FAILURE);}

  if (model == 0) {
    fprintf(outfile, "1\t%.5e\n", regRate*1.e-8);
  }
  else if (model == 1) {
    start_local = 0;
    rhot = (1 - bkgdFract) * regRate;
    b = (meanSpace-1) / shape;
    if (local_shape != 0) {
      rhot *= rndgamma(local_shape) / local_shape;
      /*      if (rhot > 10) {rhot = 10.;} */ /* trim extreme values */
    }
    meanInten = meanSpace * rhot;
    filler = 60 * meanSpace;
    totDist = -filler;
    thisd = b * rndgamma(shape);
    while (totDist + thisd < 1) {
      totDist += thisd;
      thisd = b * rndgamma(shape);
    }

    if ((int) (totDist + thisd + 0.5) > 1) {
      /* Starting region in cold spot */
      fprintf(outfile, "1\t%.5e\n", (bkgdFract * regRate)*1.e-8);
    }
    while (totDist + thisd < size) {
      if (totDist + thisd - start_local >= local_size) {
	rhot = (1 - bkgdFract) * regRate;
	if (local_shape != 0) {
	  rhot *= rndgamma(local_shape) / local_shape;
	  /*      if (rhot > 10) {rhot = 10.;} */ /* trim extreme values */
	}
	meanInten = meanSpace * rhot;
	start_local += local_size;
	continue;
      }
      totDist += thisd;
      rspot = meanInten;
      if (intensity_shape != 0) {
	rspot *= rndgamma(intensity_shape) / intensity_shape;
      }
      fprintf(outfile, "%d\t%.5e\n", (int) (totDist + 0.5), 
	      (rspot  / spot_size + bkgdFract * regRate)*1.e-8);  
      if (totDist + 0.5 + spot_size < size) {
	fprintf(outfile, "%d\t%.5e\n", (int) (totDist + 0.5 + spot_size), 
		(bkgdFract * regRate)*1.e-8);  /* hotspot is spot_size bp wide */
      }
      totDist += spot_size;
      thisd = b * rndgamma(shape);  /* distance to start of next hotspot */
    }
  }
  else {
    fprintf(stderr, "Error (recosim.c): Unknown model: %d\n", model);
  }
  return 0;
}
