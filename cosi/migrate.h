/* $Id: migrate.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_MIGRATE_H
#define __INCLUDE_COSI_MIGRATE_H

#include <boost/shared_ptr.hpp>
#include <cosi/decls.h>
#include <cosi/hooks.h>
#include <cosi/cosirand.h>

namespace cosi {

//
// Class: Migrate
//
// Keeps track of the current migration rates between populations.
// Calculates the total propability of migration in the current generation, based on these rates.
// Executes a migration.
class Migrate: public Hookable, public HasRandGen {
public:
	 Migrate( DemographyP demography_ );

	 // MethodP: migrate_set_rate
	 // Set the rate of migration between two specified pops, overriding any existing rate.
	 void migrate_set_rate (popid from, popid to, prob_per_chrom_per_gen_t rate);

	 void migrate_delete_all_for_pop( popid );
	 
	 struct MigrateRate {
			Pop* frompop;
			Pop* topop;
			prob_per_chrom_per_gen_t rate;
			MigrateRate* next;
	 };

	 const MigrateRate *getMigrations() const { return migrations; }
	 
	 // MethodP: migrate_get_all_nodes_rate
	 // Return the probability of one of the active nodes (chroms) migrating to another pop.
	 prob_per_gen_t migrate_get_all_nodes_rate (void) const;

	 // MethodP: migrate_execute
	 // Pick one node to migrate based on the current migration rates, and migrate it.
	 void migrate_execute (genid gen);

private:
	 DemographyP demography;

	 MigrateRate *migrations;
	 mutable prob_per_gen_t mig_lastrate;

	 void migrate_add (popid from, popid to, prob_per_chrom_per_gen_t rate);
	 void migrate_delete (popid from, popid to);

};  // class Migrate

typedef boost::shared_ptr<Migrate> MigrateP;
  
} // namespace cosi

#endif
// #ifndef __INCLUDE_COSI_MIGRATE_H

