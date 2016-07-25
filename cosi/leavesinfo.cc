#include <cosi/leavesinfo.h>

namespace cosi {

LeavesInfoP computeStdLeavesInfo( std::vector<nchroms_t> const& sampleSizes, 
																	std::vector<popid> const& popNames ) {
	LeavesInfoP leavesInfo = boost::make_shared<LeavesInfo>();
	leavesInfo->sampleSizes = sampleSizes;
	leavesInfo->popNames = popNames;
	leaf_id_t leafId = 0;
	for ( size_t i = 0; i < sampleSizes.size(); ++i ) {
		for ( int k = 0; k < sampleSizes[i]; ++k ) {
			leavesInfo->leafOrder.push_back( leafId++ );
		}
	}
	return leavesInfo;
}

std::ostream& operator<<( std::ostream& s, LeavesInfo const& leavesInfo ) {
	s << "[LeavesInfo:";
	for ( size_t i = 0; i < leavesInfo.sampleSizes.size(); ++i ) 
		s << "  pop idx " << i << " name " << leavesInfo.popNames[i] << " sampleSize " << 
			leavesInfo.sampleSizes[i] << " ; ";

	s << " leafOrder: nleaves=" << leavesInfo.leafOrder.size() << " leaves: " ;
	for ( size_t j=0; j < leavesInfo.leafOrder.size(); ++j )
		s << " " << leavesInfo.leafOrder[j];

	s << "]";
	return s;
}

}  //  namespace cosi

																		
