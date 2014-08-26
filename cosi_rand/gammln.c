/* $Id: gammln.c,v 1.2 2011/05/03 18:50:21 sfs Exp $ */
/* ln of gamma function.  From Numerical Recipes. */

#include <math.h>
#include <cosi_rand/gammln.h>

double gammln(double xx) {

  double x, tmp, ser;
  static double cof[6] = {76.18009173, -86.50532033, 24.01409822, 
       -1.231739516, 0.120858003e-2, -0.536382e-5};
  int j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j=0; j <= 5; j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp + log(2.50662827465 * ser);
}












