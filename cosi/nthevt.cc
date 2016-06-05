#include <vector>
#include <iostream>
#include <string>
#include <boost/cstdint.hpp>
#include <cosi/general/utils.h>
#include <cosi/nthevt.h>

namespace cosi {
namespace nthevt {
namespace {
struct NthEvtRec {
	evt_kind_t kind;
	util::StatKeeper<gens_t> dt;

	NthEvtRec(): kind(EVT_NONE) { clear(); }
	void clear() { dt.clear(); }
};

static std::vector<NthEvtRec> evtStats;

static inline NthEvtRec& erec(evt_kind_t k){ return evtStats[(size_t)k]; }

static bool nthEvtActive = ( getenv( "COSI_NTHEVT" ) != NULL );

}  // anon namespace

bool skipRest() { return nthEvtActive; }

void record_evt( evt_kind_t k, genid gen ) {
	if ( !nthEvtActive ) return;
	if ( evtStats.empty() ) {
		evtStats.resize( 4 );
		erec(EVT_RECOMB).kind = EVT_RECOMB;
		erec(EVT_COAL).kind  = EVT_COAL;
		erec(EVT_MIGR).kind = EVT_MIGR;
	}
	erec(k).dt.add(gen-genid(0.));
}

namespace {
struct PrintStats {
	PrintStats() { }
	~PrintStats() {
		if ( !nthEvtActive ) return;
		std::cerr << "eventKind\tn\tfrac\tmean_dt\tstd_dt\n";
		evt_kind_t evtKinds[] = { EVT_NONE, EVT_RECOMB, EVT_COAL, EVT_MIGR };
		const char *evtKindNames[] = {"NONE", "RECOMB", "COAL", "MIGR" };
		double totEvts = 0.;
		for ( int i =  1; i < 4; ++i )
			totEvts += erec(evtKinds[i]).dt.getNumVals();

		for ( int i =  1; i < 4; ++i ) {
			NthEvtRec& rec=erec(evtKinds[i]);
			std::cerr.precision(3);
			std::cerr << evtKindNames[i] << "\t" << rec.dt.getNumVals() << "\t" 
								<< (rec.dt.getNumVals()/totEvts) << "\t" << rec.dt.getMean() << "\t" << rec.dt.getStd() << "\n";
		}
	}
} p;

}
}  // namespace nthevt
}  // namespace cosi
