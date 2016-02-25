#include <cosi/generalmath.h>

namespace cosi {
namespace math {

// cosi_double integrate( const vector<cosi_double>& f, const vector<cosi_double>& x, cosi_double from, cosi_double to ) {
//   assert( 0 <= from && from <= f.size()-1 );
//   assert( 0 <= to && to <= f.size()-1 );

//   cosi_double firstF = interpolate( f, from );
//   cosi_double firstX = interpolate( x, from );
//   cosi_double lastF = interpolate( f, to );
//   cosi_double lastX = interpolate( x, to );

//   if ( ((int)from) == ((int)to) )
// 		 // The integral doesn't cross any of the f samples,
// 		 // handle this case separately
// 		 return 0.5*(firstF+lastF)*(lastX-firstX);

//   int ceilFrom = (int)ceil( from );
//   int floorTo = (int)floor( to );

//   cosi_double I = 0.5*(f[ceilFrom]+firstF)*(x[ceilFrom]-firstX);
//   for ( int i = ceilFrom; i < floorTo; i++ ) {
// 		cosi_double incr = 0.5*(f[i]+f[i+1])*(x[i+1]-x[i]);
// 		I += incr;
//   }

//   cosi_double finalAdd = 0.5*(f[floorTo]+lastF)*(lastX-x[floorTo]);
//   I += finalAdd;
//   return I;
// }  // integrate()


}  // namespace math
}  // namespace cosi
