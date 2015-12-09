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


void 
Migrate::migrate_delete_all_for_pop (popid pop) {
	PRINT2( "deleting_all_migrations_involving", pop );
	//std::cerr << "deleting_all_migrations_involving " << pop << "\n";
  for ( MigrateRate **cur = &migrations; *cur; ) {
		//std::cerr << " looking migration from " << (*cur)->frompop->pop_get_name() << " to " << (*cur)->topop->pop_get_name() << "\n";
		if ( (*cur)->frompop->pop_get_name() == pop || (*cur)->topop->pop_get_name() == pop ) {
			MigrateRate *to_del = *cur;
			*cur = (*cur)->next;
			//std::cerr << " deleting migration from " << to_del->frompop->pop_get_name() << " to " << to_del->topop->pop_get_name() << "\n";
			free( to_del );
		} else
			 cur = &( (*cur)->next );
  }
}


gensInv_t
Migrate::migrate_get_all_nodes_rate ( genid gen ) const
{
  popsize_float_t num_nodes;
  MigrateRate *tempmigrate = migrations;
  gensInv_t rate(0.0);

	if ( baseModel ) {
		for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
				 pi != baseModel->popInfos.end(); ++pi ) {
			num_nodes = static_cast<popsize_float_t>( demography->dg_get_pop_by_name( pi->first )->pop_get_num_nodes() );
			if ( num_nodes > static_cast<popsize_float_t>(0.) ) {
				BOOST_AUTO( const& popInfo, pi->second );
				for( BOOST_AUTO( mi, popInfo.migrRateTo.begin() );
						 mi != popInfo.migrRateTo.end(); ++mi ) {
					rate += num_nodes * mi->second( gen );
				}  // for each migration out of pop
			}  // if pop sample is nonempty
		}  // for each pop
	} // if basemodel
	else {

		if (migrations == NULL)
			 mig_lastrate = gensInv_t( 0.0 );
		else {
			while (tempmigrate != NULL) {
				// if ( tempmigrate->frompop->isInactive() || tempmigrate->topop->isInactive() )
				// 	 std::cerr << " looking migration from " << tempmigrate->frompop->pop_get_name() << " to " << tempmigrate->topop->pop_get_name() << "\n";
				util::chkCond( !tempmigrate->frompop->isInactive(), "migration on inactive node!" );
				util::chkCond( !tempmigrate->topop->isInactive(), "migration on inactive node!" );
				num_nodes = static_cast<popsize_float_t>( tempmigrate->frompop->pop_get_num_nodes() );
				rate += num_nodes * tempmigrate->rate;
				tempmigrate = tempmigrate->next;
			}
		}
	}  // if not basemodel
	mig_lastrate = rate;
  
  return mig_lastrate;
}


void 
Migrate::migrate_execute (genid gen) 
{
  popsize_float_t numnodes;
  MigrateRate *tempmigrate = migrations;
  gensInv_t rate( 0.0 );
  gensInv_t randcounter( factor_t( random_double() ) * mig_lastrate );
  
	{

		if ( baseModel ) {
			for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
					 pi != baseModel->popInfos.end(); ++pi ) {
				Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
				numnodes = static_cast<popsize_float_t>( srcPop->pop_get_num_nodes() );
				if ( numnodes > static_cast<popsize_float_t>(0.) ) {
					BOOST_AUTO( const& popInfo, pi->second );
					for( BOOST_AUTO( mi, popInfo.migrRateTo.begin() );
							 mi != popInfo.migrRateTo.end(); ++mi ) {
						rate += numnodes * mi->second( gen );
						if ( rate >= randcounter ) {
							Pop *dstPop = demography->dg_get_pop_by_name( mi->first );
							demography->dg_migrate_one_chrom (srcPop,
																								dstPop, gen);
							return;
						}
					}  // for each migration out of pop
				}  // if pop sample is nonempty
			}  // for each pop
		}  // if basemodel
		else {
			while (tempmigrate != NULL && rate < randcounter) {
				numnodes = static_cast< popsize_float_t >( tempmigrate->frompop->pop_get_num_nodes() );
				rate += numnodes * tempmigrate->rate;
				if (rate < randcounter)
					 tempmigrate = tempmigrate->next;
				else
					 demography->dg_migrate_one_chrom (tempmigrate->frompop,
																						 tempmigrate->topop, gen);
      
			}  // for each migration
		}  // if not basemodel
  }  // if migrations!=NULL
}

#if 0
class MigrationProcess;

void createMigrationProcesses( DemographyP demography, BaseModelP baseModel ) {
	std::vector< MigrationProcess > migrationProcesses;
	for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
			 pi != baseModel->popInfos.end(); ++pi ) {
		Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
		BOOST_AUTO( const& popInfo, pi->second );
		for( BOOST_AUTO( mi, popInfo.migrRateTo.begin() );
				 mi != popInfo.migrRateTo.end(); ++mi ) {
			Pop *dstPop = mi->first;
			migrationProcesses.push_back( makePoissonProcess( /* rateFn= */ mi->second,
																												MigrateProcess( srcPop, dstPop, demography ) ) )
		}
	}
}  // createMigrationProcesses

class MigrationProcess: public EventRunner< genid, RandGen
	 Pop *srcPop, *dstPop;
	 DemographyP demography;

	 virtual double getRateFactor() const { return static_cast<double>( pop->pop_get_num_nodes() ); }
	 virtual void executeNextEvent( genid gen ) { demography->dg_migrate_one_chrom( srcPop, dstPop, gen ); }
		 
	 }
};  // class MigrateProcess
#endif

}  // namespace cosi

#include <cosi/arrproc2.h>

namespace cosi {

class MigrationProcess: public arrival2::ArrivalProcessDef< genid, RandGen, popsize_float_t > {

	 DemographyP demography;
	 Pop *srcPop;
	 Pop *dstPop;

public:

	 MigrationProcess( DemographyP demography_, Pop *srcPop_, Pop *dstPop_ ):
		 demography( demography_ ), srcPop( srcPop_ ), dstPop( dstPop_ ) {
	 }

	 virtual popsize_float_t getRateFactor() const { return static_cast<popsize_float_t>( srcPop->pop_get_num_nodes() ); }
	 virtual void executeEvent( genid gen, RandGen& ) { demography->dg_migrate_one_chrom( srcPop, dstPop, gen ); }
};  // class MigrateProcess


boost::shared_ptr< Migrate::migr_processes_type >
Migrate::createMigrationProcesses() {
	using namespace arrival2;
	using math::Any;
	using math::Const;
	using math::Piecewise;

	boost::shared_ptr< migr_processes_type > migrProcs =
		 boost::make_shared<migr_processes_type>();
	for( BOOST_AUTO( pi, baseModel->popInfos.begin() );
			 pi != baseModel->popInfos.end(); ++pi ) {
		Pop *srcPop = demography->dg_get_pop_by_name( pi->first );
		BOOST_AUTO( const& popInfo, pi->second );
		for( BOOST_AUTO( mi, popInfo.migrRateTo.begin() );
				 mi != popInfo.migrRateTo.end(); ++mi ) {
			Pop *dstPop = demography->dg_get_pop_by_name( mi->first );

			boost::shared_ptr< MigrationProcess > mp = boost::make_shared<MigrationProcess>( demography, srcPop, dstPop );
			add( *migrProcs,
					 ArrivalProcess< genid, Stoch< RandGen, Poisson< Piecewise< Const<> >, popsize_float_t > > >
					 ( mi->second, genid(0.),  mp ) );
																												
		}
	}
	return migrProcs;
}  // createMigrationProcesses

}  // namespace cosi
