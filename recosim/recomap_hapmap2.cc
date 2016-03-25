#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/random_device.hpp>
#include <boost/integer/integer_mask.hpp>
#include <boost/cstdint.hpp>

int main(int argc, char **argv) {
  int itoken, ichr, isite, *long_enough;
  const int line_size = 256;
  char filename[256], *new_line=NULL, *running=NULL, *token=NULL, chrstr[8], 
		 command[256], seqstr[10];
  long min_pos, max_pos, *chrlen, seqtot=0, seqlen, sum=0;
  double x;

	boost::mt19937 rand_engine;
	boost::uniform_01<> u01;
	boost::random_device rand_dev;

	boost::uint32_t seed = rand_dev() & boost::low_bits_mask_t<32>::sig_bits;
	rand_engine.seed( static_cast<boost::uint32_t>( seed ) );
  
  if (argc != 4) {
    fprintf(stderr, "Usage: recomap_hapmap2 <1kG map directory> <sequence length (bp)> <output file>\n");
    return EXIT_FAILURE;
  }
  seqlen = strtol(argv[2], NULL, 10);
  chrlen = (long *)calloc(24, sizeof(long));
	std::cerr << "seqlen=" << seqlen << " chrlen=" << chrlen << "\n";
  long_enough = (int *)calloc(24, sizeof(int));
  fprintf(stdout, "random seed for choosing chromosome: %lu\n", (unsigned long)seed);
  new_line = (char *)malloc((line_size + 1) * sizeof(char));
  for (ichr = 1; ichr <= 22; ichr++) {
    strcpy(filename, argv[1]);
    strcat(filename, "/genetic_map_chr");
    sprintf(chrstr, "%d", ichr);
    strcat(filename, chrstr);
    strcat(filename, "_b37.txt");
		std::cerr << "filename=" << filename << "\n";
    FILE *inf = fopen(filename, "r");
    if (inf == NULL) {
      fprintf(stderr, "Could not open map file %s\n", filename);
      return EXIT_FAILURE;
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
    return EXIT_FAILURE;
  }
  x = u01(rand_engine);
  for (ichr = 1; ichr <= 22; ichr++) {
    sum += chrlen[ichr];
    if ((double) sum / seqtot > x) {
      fprintf(stdout, "using chrom %d\n", ichr);
      break;
    }
  }
  sprintf(seqstr, "%ld", seqlen);
  strcpy(command, "get_recomap ");
  strcat(command, filename);
  strcat(command, " ");
  strcat(command, seqstr);
  strcat(command, " ");
  strcat(command, argv[3]);
	std::cerr << "now run: " << command << "\n";
  int exitCode = system(command);

	free(chrlen);
	free(new_line);

  return exitCode;
}
