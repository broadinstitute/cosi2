/* $Id: gamma.c,v 1.3 2011/05/06 23:13:09 sfs Exp $ */

#include <math.h>
#include <cosi_rand/random.h>
#include <cosi_rand/gamma.h>

double rndgamma1(double);
double rndgamma2(double);

/* gamma function routine #1 from Yang's PAML package  */
/* Assumes random seed already initialized */
/* s = shape param = k = mean**2/Var  */
/* multiply output by mean/k to get correct mean and Var */
double rndgamma (double s)
{
        double rndgamma1(double s1);
        double rndgamma2(double s2);
        double  r=0.0, x;

        if (s <= 0.0)
                return 0;
        else if (s < 1.0)
                r = rndgamma1 (s);
        else if (s > 1.0)
                r = rndgamma2 (s);
        else {
	        x = random_double();
	        while (x == 0) {x = random_double();} 
	        r = -log(x);	
	}
        return (r);
}

/* gamma function routine #2 from Yang's PAML package */
double rndgamma1 (double s)
{
        double                  r, x=0.0, small=1e-37, w;
        static double   a, p, uf, ss=10.0, d;

        if (s!=ss)
                {
                a  = 1.0-s;
                p  = a/(a+s*exp(-a));
                uf = p*pow(small/a,s);
                d  = a*log(a);
                ss = s;
                }
        for (;;)
                {
                r = random_double();
                if (r > p)
                        x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
                else if (r>uf)
                        x = a*pow(r/p,1/s), w=x;
                else
                        return (0.0);
                r = random_double();
                if (1.0-r <= w && r > 0.0)
                        if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                                continue;
                break;
                }
        return (x);
}

/* gamma function routine #3 from Yang's PAML package */
double rndgamma2 (double s)
{
        double                  r ,d, f, g, x;
        static double   b, h, ss=0;

        if (s!=ss)
                {
                b  = s-1.0;
                h  = sqrt(3.0*s-0.75);
                ss = s;
                }
        for (;;)
                {
                r = random_double();
                g = r-r*r;
                f = (r-0.5)*h/sqrt(g);
                x = b+f;
                if (x <= 0.0)
                        continue;
		r = random_double();
		while (r == 0) {r = random_double();}
                d = 64*r*r*g*g*g;
                if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
                        break;
                }
        return (x);
}
