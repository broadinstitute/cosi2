#include <cosi_rand/gammln.h>
#include <cosi_rand/random.h>
#include <cosi_rand/mtwist.h>
#include <cosi/cosirand.h>

namespace cosi {

factor_t RandGen::expdev (void) {
	double dum = 0;
	
	while (dum == 0.0)
		 dum = (double) 1 - random_double();
	return factor_t( -log (dum) );
}

int RandGen::ranbinom(int n, double p) 
{
/** 
		Knuth Vol 2, p 131
*/
	const int BTHRESH = 50 ;
	int a, b ;
	double x ; 
	if (p>=1) return n ;
	if (p<=0) return 0 ;
	if (n<=0) return 0 ;

	if (n<=BTHRESH) {
    return ranb1(n,p) ;  /** small case */
	}

	a  = 1 + n/2 ;  
	b  = n + 1 - a ;
	x = ranbeta((double) a, (double) b) ;
	if (x>=p) return ranbinom(a-1, p/x) ;
	return (a + ranbinom(b-1, (p-x)/(1.0-x)) ) ;
}

double 
RandGen::ranexp( void)
{
  /**
		 exponential mean 1
  */
  double          x, t;
  t = random_double();
  x = -log(1.0 - t);
  return x;
}


double 
RandGen::randev0(double a)
{
  /**
		 algorithm G6: Gamma for a < 1
  */
  double          r1, r2, x, w;
  double t = 1.0 - a;
  double p = t / (t + a * exp(-t));
  double s = 1.0 / a;
  for (;;) {
    r1 = random_double();
    if (r1 <= p) {
      x = t * pow(r1 / p, s);
      w = x;
    } else {
      x = t + log((1.0 - p) / (1.0 - r1));
      w = t * log(x / t);
    }
    r2 = random_double();
    while (r2 == 0) {r2 = random_double();}
    if ((1.0 - r2) <= w) {
      if ((1.0 / r2 - 1.0) <= w)
				 continue;
      if (-log(r2) <= w)
				 continue;
    }
    return x;

  }

}

double 
RandGen::randev1(double a)
{
  /**
		 Random gamma deviate:  a>=1
		 GBEST algorithm  (D.J. BEST: Appl. Stat. 29 p 181 1978
  */
  double          x, d, e, c, g, f, r1, r2;

  e = a - 1.0;
  c = 3.0 * a - 0.75;


  for (;;) {
    r1 = random_double();
    g = r1 - (r1 * r1);
    if (g <= 0.0)
			 continue;
    f = (r1 - 0.5) * sqrt(c / g);
    x = e + f;
    if (x <= 0.0)
			 continue;
    r2 = random_double();
    while (r2 == 0) {r2 = random_double();}
    d = 64.0 * r2 * r2 * g * g * g;
    if ((d >= 1.0 - 2.0 * f * f / x) && (log(d) >= 2.0 * (e * log(x / e) - f)))
			 continue;
    return (x);
  }

}


double RandGen::ranbeta(double a, double b) 
{
	double xa, xb ;

	xa = rangam(a) ;
	xb = rangam(b) ;
	return xa/(xa+xb) ;
}

double
RandGen::rangam(double a)
{
  /**
		 generate gamma deviate mean a
  */
  if (a < 1.0) {
    return( randev0(a));
  }
  if (a == 1.0) {
    return( ranexp());
  }
  return( randev1(a));
}

int RandGen::ranb1 (int n, double p) 
/** 
		binomial dis. 
		Naive routine
*/
{ 
	int cnt = 0, i ;

	for (i=0 ; i< n ; i++)  {
		if (random_double() <= p) ++ cnt ;
	}

	return cnt ;

}


int RandGen::poisson( double xm ) {
	const double  PI = 3.141592654;

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
      g = xm * alxm - ::gammln(xm+1.0);
    }
    do {
      do {
				y = tan( PI * random_double() );
				em = sq*y + xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9 * (1.0 + y*y) * exp(em * alxm - ::gammln(em + 1.0) - g);
    } while (random_double() > t);
  }
  return (int) (em + 0.5);
}


double
RandGen::poisson_get_next (double rate) {
	assert( rate > 0.0 );  
  double ed;
  if (rate == 0) return -1;
  ed = ToDouble( expdev() );
  return (ed / rate);
}

bool_t random_bit(void) {
  return ( mt_lrand() & 0x01 );
}


}  // namespace cosi
