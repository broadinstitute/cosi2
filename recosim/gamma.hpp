#ifndef INCLUDE_COSI_RECOSIM_GAMMA_HPP
#define INCLUDE_COSI_RECOSIM_GAMMA_HPP

#include <cmath>
#include <boost/random/uniform_01.hpp>

template <typename URNG> double rndgamma(double s, URNG& urng);
template <typename URNG> double rndgamma1(double s, URNG& urng);
template <typename URNG> double rndgamma2(double s, URNG& urng);

/* gamma function routine #1 from Yang's PAML package  */
/* Assumes random seed already initialized */
/* s = shape param = k = mean**2/Var  */
/* multiply output by mean/k to get correct mean and Var */
template <typename URNG>
double rndgamma (double s, URNG& urng)
{
	boost::uniform_01<> u01;
        double  r=0.0, x;

        if (s <= 0.0)
                return 0;
        else if (s < 1.0)
					 r = rndgamma1 (s,urng);
        else if (s > 1.0)
					 r = rndgamma2 (s,urng);
        else {
	        x = u01(urng);
	        while (x == 0) {x = u01(urng);} 
	        r = -log(x);	
	}
        return (r);
}

/* gamma function routine #2 from Yang's PAML package */
template <typename URNG>
double rndgamma1 (double s, URNG& urng)
{
        double                  r, x=0.0, small=1e-37, w;
        static double   a, p, uf, ss=10.0, d;
				boost::uniform_01<> u01;

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
									r = u01(urng);
                if (r > p)
                        x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
                else if (r>uf)
                        x = a*pow(r/p,1/s), w=x;
                else
                        return (0.0);
                r = u01(urng);
                if (1.0-r <= w && r > 0.0)
                        if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                                continue;
                break;
                }
        return (x);
}

/* gamma function routine #3 from Yang's PAML package */
template <typename URNG>
double rndgamma2 (double s, URNG& urng)
{
        double                  r ,d, f, g, x;
        static double   b, h, ss=0;
				boost::uniform_01<> u01;

        if (s!=ss)
                {
                b  = s-1.0;
                h  = sqrt(3.0*s-0.75);
                ss = s;
                }
        for (;;)
                {
									r = u01(urng);
                g = r-r*r;
                f = (r-0.5)*h/sqrt(g);
                x = b+f;
                if (x <= 0.0)
                        continue;
								r = u01(urng);
								while (r == 0) {r = u01(urng);}
                d = 64*r*r*g*g*g;
                if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
                        break;
                }
        return (x);
}

#endif // #ifndef INCLUDE_COSI_RECOSIM_GAMMA_HPP
