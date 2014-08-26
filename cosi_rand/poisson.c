/* $Id: poisson.c,v 1.2 2011/05/03 18:50:21 sfs Exp $ */

/* Return a random deviate drawn from Poisson distribution of mean xm. */
/* From Numerical Recipes. */

#include <math.h>
#include <cosi_rand/random.h>
#include <cosi_rand/poisson.h>
#include <cosi_rand/gammln.h>

#define PI 3.141592654

int poisson( double xm ) {

  /*double gammln(double argum);*/

  static double sq, alxm, g, oldm=(-1.0);
  double em, t, y;

  /**************************************************************************/

  if (xm < 12.0) {
    if (xm != oldm) {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1.;
    t = 1.0;
    do {
      em += 1.0;
      t *= random_double();
    } while (t > g);
  } 
  else {
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm * alxm - gammln(xm+1.0);
    }
    do {
      do {
	y = tan( PI * random_double() );
	em = sq*y + xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9 * (1.0 + y*y) * exp(em * alxm - gammln(em + 1.0) - g);
    } while (random_double() > t);
  }
  return (int) (em + 0.5);
}


double
poisson_get_next (double rate) {
  double ed;
  if (rate == 0) return -1;
  ed = expdev();
  return (ed / rate);
}

