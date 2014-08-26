/* $Id: random.h,v 1.3 2011/05/06 23:13:09 sfs Exp $ */

#ifndef COSI_RAND_RANDOM_H
#define COSI_RAND_RANDOM_H

#ifdef __cplusplus
extern "C"
    {
#endif


double random_double (void);
double expdev (void);
unsigned long seed_rng (void);  /* pick "random" seed based on /dev/urandom or time */
void set_rng_seed(unsigned long newseed);

#ifdef __cplusplus
    }
#endif

#endif
