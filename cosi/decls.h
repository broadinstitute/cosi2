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

#define COSI_DECL(cls) \
	class cls;                                  \
	typedef boost::shared_ptr< cls > cls ## P ; \
	typedef boost::shared_ptr< const cls > cls ## CP ; \
	
COSI_DECL(Sweep); // defined in sweep.h
COSI_DECL(Demography);  // defined in demography.h
COSI_DECL(Migrate);  // defined in migrate.h
COSI_DECL(Mutate);  // defined in mutate.h
COSI_DECL(MutProcessor);  // defined in mutate.h
COSI_DECL(GeneConversion);  // defined in geneconversion.h
COSI_DECL(Recomb);  // defined in recomb.h
COSI_DECL(GenMap);  // genmap.h
COSI_DECL(Pop);  // defined in pop.h
COSI_DECL(HistEvents); // defined in historical.h
COSI_DECL(Mutlist); // defined in mutlist.h
COSI_DECL(Mut);  // defined in mutlist.h
COSI_DECL(ParamFileReader); // defined in file.h

namespace coal {
class Coalesce;  // defined in coalesce.h
}
typedef boost::shared_ptr<coal::Coalesce> CoalesceP;

COSI_DECL(Simulator);  // defined in simulator.h
COSI_DECL(Hook);  // defined in hooks.h
COSI_DECL(Hooks);  // defined in hooks.h
COSI_DECL(Event_SweepNew);  // defined in sweep.cc
COSI_DECL(FreqTraj);  // defined in traj.h
COSI_DECL(DeterministicSweepTraj);  // defined in traj.h
COSI_DECL(RandGen); // defined in cosirand.h
COSI_DECL(RecombRecorder); // defined in recomb.h
COSI_DECL(RecomblessNeighborhoodsTracker); // defined in recomb.h
COSI_DECL(HullMgr);
COSI_DECL(CondSnpDef); // defined in condsnp.h
COSI_DECL(CondSnpMgr); // defined in condsnp.h

COSI_DECL(BaseModel);  // defined in basemodel.h
COSI_DECL(SweepInfo);
COSI_DECL(MSweep); // defined in msweep.h

}  // namespace cosi


#endif // #ifndef __INCLUDE_COSI_DECLS_H
