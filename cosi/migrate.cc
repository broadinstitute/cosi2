/* $Id: migrate.c,v 1.4 2011/05/26 21:36:22 sfs Exp $ */

/* MIGRATE.C
 * determines when a migration should next occur.
 *
 * c.f. migrate.h
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cosi/migrate.h>
#include <cosi/demography.h>
#include <cosi/pop.h>

namespace cosi {

Migrate::Migrate( DemographyP demography_ ):
  demography( demography_ ), migrations(NULL), mig_lastrate(0.0) {
}

// Method: migrate_set_rate
// Set the rate of migration between two specified pops, overriding any existing rate.
void Migrate::migrate_set_rate (popid from, popid to, prob_per_chrom_per_gen_t rate) {
	migrate_delete( from, to );
	migrate_add( from, to, rate );
	hooks->fire_set_migrate_rate( from, to, rate );
}

void 
Migrate::migrate_add (popid from, popid to, prob_per_chrom_per_gen_t rate) 
{
  /* old migration, if any, will have been replaced.*/
  if (is_zero( rate )) {return;}

  MigrateRate *newmigrate = (MigrateRate *)malloc(sizeof (MigrateRate));
  assert(newmigrate != NULL);
  
  newmigrate->frompop = demography->dg_get_pop_by_name(from);
  newmigrate->topop = demography->dg_get_pop_by_name(to);
  newmigrate->rate = rate;
  newmigrate->next = migrations;
  
  migrations = newmigrate;
}

void 
Migrate::migrate_delete (popid from, popid to) {
  for ( MigrateRate **cur = &migrations; *cur; cur = &( (*cur)->next ) ) {
		if ( (*cur)->frompop->pop_get_name() == from && (*cur)->topop->pop_get_name() == to ) {
			MigrateRate *to_del = *cur;
			*cur = (*cur)->next;
			free( to_del );
			return;
		}
  }
}


prob_per_gen_t
Migrate::migrate_get_all_nodes_rate () const
{
  nchroms_t numnodes;
  MigrateRate *tempmigrate = migrations;
  prob_per_gen_t rate(0.0);
  
  if (migrations == NULL)
		 mig_lastrate = prob_per_gen_t( 0.0 );
  else {
    while (tempmigrate != NULL) {
      numnodes = tempmigrate->frompop->pop_get_num_nodes();
      rate += numnodes * tempmigrate->rate;
      tempmigrate = tempmigrate->next;
    }
    mig_lastrate = rate;
  }
  return mig_lastrate;
}


void 
Migrate::migrate_execute (genid gen) 
{
  int numnodes;
  MigrateRate *tempmigrate = migrations;
  prob_per_gen_t rate( 0.0 );
  prob_per_gen_t randcounter( factor_t( random_double() ) * mig_lastrate );
  
  if (migrations == NULL)
		 fprintf(stderr, "ERROR in migrate.\n");
  
  else {
    
    while (tempmigrate != NULL && rate < randcounter) {
      numnodes = tempmigrate->frompop->pop_get_num_nodes();
      rate += numnodes * tempmigrate->rate;
      if (rate < randcounter)
				 tempmigrate = tempmigrate->next;
      else
				 demography->dg_migrate_one_chrom (tempmigrate->frompop,
																					 tempmigrate->topop, gen);
      
    }
  }
}

}  // namespace cosi

