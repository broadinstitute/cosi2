/* $Id: recomb.h,v 1.3 2011/05/04 15:46:19 sfs Exp $ */

#ifndef __INCLUDE_COSI_RECOMB_H
#define __INCLUDE_COSI_RECOMB_H

#include <cstdlib>
#include <vector>
#include <utility>
#include <boost/noncopyable.hpp>
#include <boost/multi_array.hpp>
#include <cosi/decls.h>
#include <cosi/cosirand.h>
#include <cosi/hooks.h>
#include <cosi/nodefwd.h>
#include <cosi/arrproc2.h>

namespace cosi {

using std::vector;
using std::pair;
using node::Node;

//
// Class: Recomb
//
// Executes recombination events.
//
class Recomb: boost::noncopyable, public HasRandGen {
public:
	 Recomb( DemographyP demography_,
					 GenMapP genMap_ );
	 ~Recomb();
	 
	 // MethodP: getAllNodesRecombRate
	 // Returns the probability, per generation, of a recombination that splits one of the currently active nodes.
	 glen_t getAllNodesRecombRate() const;

	 void recomb_execute (genid gen, int popindex, loc_t *location, Node**nodes_out);
	 void recomb_execute (genid gen, frac_t frac);

	 unsigned long getNumRecombs() const { return nrecombs; }

	 void setIgnoreRecombsInPop( popid ignoreRecombsInPop_ ) {
		 this->ignoreRecombsInPop = ignoreRecombsInPop_;
	 }

	 typedef
	 arrival2::ArrivalProcess< genid, arrival2::Stoch< RandGen, arrival2::Poisson< math::Const<>, double > > >
	 recomb_processes_type;

	 boost::shared_ptr< recomb_processes_type > createRecombProcesses();

private:
	 DemographyP demography;
	 GenMapP genMap;
	 unsigned long nrecombs;

	 // Field: ignoreRecombsInPop
	 // If not NULL_POPID, recomb events in this pop will be ignored.
	 popid ignoreRecombsInPop;
};  // class Recomb

// Class: RecombRecorder
// Records locations of recomb events, if the list of these locations is requested to be output.
class RecombRecorder: public Hook {
public:
	 RecombRecorder();
	 virtual ~RecombRecorder();

	 // Virtual methodP: handle_recomb
	 // Called after each recombination.
	 virtual void handle_recomb( Node *, Node *, loc_t, genid);

	 const vector<loc_t>& getRecombLocs() const { return recombLocs; }

private:
	 vector< loc_t > recombLocs;
};  // class RecombRecorder

}  // namespace cosi
  
#endif
// #ifndef __INCLUDE_COSI_RECOMB_H
