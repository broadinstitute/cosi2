//
// Header: decls.h
//
// Forward declarations of classes and structs used in cosi.
// Used to reduce explicit dependencies between header files.
//

#ifndef __INCLUDE_COSI_DECLS_H
#define __INCLUDE_COSI_DECLS_H

#include <boost/shared_ptr.hpp>

namespace cosi {

struct Sweep;  // defined in sweep.h
typedef boost::shared_ptr<Sweep> SweepP;

class Demography;  // defined in demography.h
typedef boost::shared_ptr<Demography> DemographyP;

class Migrate;  // defined in migrate.h
typedef boost::shared_ptr<Migrate> MigrateP;

struct Mutate;  // defined in mutate.h
typedef boost::shared_ptr<Mutate> MutateP;
class MutProcessor;  // defined in mutate.h
typedef boost::shared_ptr<MutProcessor> MutProcessorP;

class GeneConversion;  // defined in geneconversion.h
typedef boost::shared_ptr<GeneConversion> GeneConversionP;

class Recomb;  // defined in recomb.h
typedef boost::shared_ptr<Recomb> RecombP;

class GenMap;  // defined in genmap.h
typedef boost::shared_ptr<GenMap> GenMapP;

// namespace seglist {
// struct seglist;  // defined in seglist.h
// typedef struct seglist Seglist;

// struct seg;  // defined in seglist.h
// typedef struct seg Seg;

// struct segptr;  // defined in seglist.h
// typedef struct segptr Segptr;
// }

struct Pop;  // defined in pop.h
typedef boost::shared_ptr<Pop> PopP;

class HistEvents;  // defined in historical.h
typedef boost::shared_ptr<HistEvents> HistEventsP;

class Mutlist;  // defined in mutlist.h
typedef boost::shared_ptr<Mutlist> MutlistP;

class Mut;  // defined in mutlist.h

class ParamFileReader;  // defined in file.h
typedef boost::shared_ptr<ParamFileReader> ParamFileReaderP;

namespace coal {
class Coalesce;  // defined in coalesce.h
}
typedef boost::shared_ptr<coal::Coalesce> CoalesceP;

class Simulator;  // defined in simulator.h
typedef boost::shared_ptr<Simulator> SimulatorP;

class Hook;  // defined in hooks.h
typedef boost::shared_ptr<Hook> HookP;

class Hooks;  // defined in hooks.h
typedef boost::shared_ptr<Hooks> HooksP;

class Event_SweepNew;  // defined in sweep.cc
typedef boost::shared_ptr<Event_SweepNew> Event_SweepNewP;

class FreqTraj;  // defined in traj.h
typedef boost::shared_ptr<FreqTraj> FreqTrajP;

class DeterministicSweepTraj;  // defined in traj.h
typedef boost::shared_ptr<DeterministicSweepTraj> DeterministicSweepTrajP;

class RandGen; // defined in cosirand.h
typedef boost::shared_ptr<RandGen> RandGenP;

class RecombRecorder; // defined in recomb.h
typedef boost::shared_ptr<RecombRecorder> RecombRecorderP;

class RecomblessNeighborhoodsTracker; // defined in recomb.h
typedef boost::shared_ptr<RecomblessNeighborhoodsTracker> RecomblessNeighborhoodsTrackerP;

class HullMgr;
typedef boost::shared_ptr<HullMgr> HullMgrP;  // defined in hullmgr.h

class CondSnpDef;
class CondSnpMgr;
typedef boost::shared_ptr<CondSnpDef> CondSnpDefP;  // defined in condsnp.h
typedef boost::shared_ptr<CondSnpMgr> CondSnpMgrP;  // defined in condsnp.h

class BaseModel;  // defined in basemodel.h
typedef boost::shared_ptr<BaseModel> BaseModelP;

class MSweep; // defined in msweep.h
typedef boost::shared_ptr<MSweep> MSweepP;

}  // namespace cosi


#endif // #ifndef __INCLUDE_COSI_DECLS_H
