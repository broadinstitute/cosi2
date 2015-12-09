#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>
#include <cosi/coal_data_cosi.h>
#include <cosi/mutlist.h>
#include <cosi/utils.h>
#include <cosi/genmap.h>
#include <cosi/demography.h>

int customstats_main(int numReps, int seqlen, int NPOPS, int excludePopIdx, void (*get_coal_data_from_cosi)( coal_data *, int, int )  );
void calc_Fst_main( int numReps, int NPOPS, int excludePopIdx,
										void (*get_coal_data_from_cosi)( coal_data *, int /* irep */, int /* ipop */ ) );

namespace cosi {

namespace customstats {

namespace {
coal_data *coal_datas;
size_t nsims;
//size_t simNum;
std::vector<popid> popNames;
coal_data *data;
int seqlen;
pop_idx_t excludePopIdx = NULL_POP_IDX;
}  // namespace anon

void init( DemographyP demography, size_t nsims_, int seqlen_, popid excludePop_ ) ;

void init( DemographyP demography, size_t nsims_, int seqlen_, popid excludePop_ ) {
	popNames = demography->getPopNames();
	data = coal_datas = new coal_data[ ( nsims = nsims_ ) * popNames.size() ];
	seqlen = seqlen_;
	if ( !is_null( excludePop_ ) ) excludePopIdx = demography->dg_get_pop_index_by_name( excludePop_ );
	std::cerr << "nsims is now " << nsims << "\n";
}


void my_get_coal_data_from_cosi( coal_data *d, int irep, int ipop );

void my_get_coal_data_from_cosi( coal_data *d, int irep, int ipop ) {
	//std::cerr << "getting coal data: irep=" << irep << " ipop=" << ipop << "popNames.size=" << popNames.size() << "\n";
	*d = coal_datas[ irep * popNames.size() + ipop - 1 ];
}

void finish();

void finish() {
	std::cerr << "computing stats...\n";
	customstats_main( nsims, seqlen, popNames.size(), excludePopIdx, my_get_coal_data_from_cosi );
	std::cerr << "computing Fst...\n";
	calc_Fst_main( nsims, popNames.size(), excludePopIdx, my_get_coal_data_from_cosi );
	std::cerr << "done computing stats...\n";
}

void record_sim(DemographyP demography, GenMapP genMap, len_bp_int_t length, MutlistP mutlist,
								bool_t inf_sites ); 

void record_sim(DemographyP demography, GenMapP genMap, len_bp_int_t length, MutlistP mutlist,
								bool_t /*inf_sites*/ ) {
	using std::fill;

//  freq_t freq;

  size_t nmuts = mutlist->size();
  // char *p = blank_line;
  // for ( size_t i = 0; i < nmuts; i++ ) {
	// 	*p++ = '2';
	// 	*p++ = ' ';
  // }
  // *p = 0;

	leaf_id_t leaf = 0;
	const vector< popid >& popNames = demography->getPopNames();
	const vector< nchroms_t >& sampleSizes = demography->getSampleSizes();
  for (size_t ipop = 0; ipop < popNames.size(); ipop++, ++data) {
		data->nsample = sampleSizes[ ipop ];
		data->samp_id = data->snp_id = NULL;
		fill( data->snp_base, data->snp_base+4, (int *)NULL );
		
    if ( sampleSizes[ipop] > 0) {
			data->nsnp = nmuts;
			data->genotype = (int **)malloc(data->nsample * sizeof(int*));
			for (int isamp = 0; isamp < data->nsample; isamp++) {
				data->genotype[isamp] = (int *)malloc(data->nsnp * sizeof(int));
				fill( data->genotype[isamp], data->genotype[isamp] + data->nsnp, 2 );
			}
			data->anc_base = (int *)malloc(data->nsnp * sizeof(int));
			fill( data->anc_base, data->anc_base + data->nsnp, 2 );
			data->pos = (int *)malloc(data->nsnp * sizeof(int));
			data->gdPos = (double *)malloc(data->nsnp * sizeof(double));
			data->nallele = (int *)malloc(data->nsnp * sizeof(int));
			fill( data->nallele, data->nallele + data->nsnp, 2 );
			data->genPos = NULL;//(double *)malloc(data->nRecom * sizeof(double));
			data->physPos = NULL; //(int *)malloc(data->nRecom * sizeof(int));

			//vector< nchroms_t > mutcount( nmuts );

			//string filename( (boost::format( "%s.hap-%d" ) % filebase % popNames[ipop]).str() );

			leaf_id_t popEndLeaf = leaf + sampleSizes[ipop];

			for ( int isamp = 0; leaf < popEndLeaf; leaf++, ++isamp) {
				
				//fprintf(outf, "%d\t%d\t", leaf, ToInt( popNames[ipop] ) );
				//memcpy( line, blank_line, line_size ); 
			
				const vector< Mutlist::const_iterator >& leafMuts = mutlist->getLeafMuts( leaf );
				BOOST_FOREACH( Mutlist::const_iterator m, leafMuts ) {
					int mutId = m->mutId;
					data->genotype[isamp][mutId] = 1;
					//mutcount[ mutId ]++;
				}
				// fputs( line, outf );
				// fputs( "\n", outf);
			}  // write out haps for this pop

			//string pos_filename( (boost::format( "%s.pos-%d" ) % filebase % popNames[ipop]).str() );

      //if (outf == NULL) {fprintf(stderr, "Could not open %s\n", pos_filename.c_str());}
      //fprintf(outf, "SNP     CHROM   CHROM_POS       ALLELE1 FREQ1   ALLELE2 FREQ2\n");
      BOOST_AUTO( it, mutlist->getMuts().begin() );
      for (size_t im = 0; im < nmuts; im++, it++) {
				int mutId = it->mutId;
				data->pos[mutId] = (int) (length * get_loc( it->loc ) );
				data->gdPos[mutId] = ToDouble( genMap->getGdPos( loc_t( util::getFrac( data->pos[mutId], length ) ) ) )
					 * 100.0 * ToDouble( genMap->getRegionRecombRateAbs() );
				
				// freq = (freq_t) mutcount[im] / sampleSizes[ipop];
				// if (inf_sites) {
				// 	fprintf(outf, "%d\t1\t%.4f\t1\t%.4f\t2\t%.4f\n", (int)(it->mutIdOrig+1), double( length * get_loc( it->loc ) ), 
				// 					double( freq ), double( 1 - freq ) );
				// }
				// else {
				// 	fprintf(outf, "%d\t1\t%d\t1\t%.4f\t2\t%.4f\n", (int)(it->mutIdOrig+1), (int) (length * get_loc( it->loc ) ), 
				// 					double( freq ), double( 1 - freq ) );
				// }
			}
		}  // if this pop is nonempty
  }  // for each pop
}

#if 0
void get_coal_data(coal_data* data, , int pop) {
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
#endif // if 0


} // namespace customstats
} // namespace cosi
