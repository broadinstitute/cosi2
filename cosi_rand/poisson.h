/* $Id: poisson.h,v 1.2 2011/05/03 18:50:21 sfs Exp $ */

#ifndef COSI_RAND_POISSON_H
#define COSI_RAND_POISSON_H

#ifdef __cplusplus
extern "C"
    {
#endif

int poisson(double xm);
double poisson_get_next (double rate);

#ifdef __cplusplus
    }
#endif
		
#endif  /* #ifndef COSI_RAND_POISSON_H */


