#ifndef COSI_RAND_EXPDEV_H
#define COSI_RAND_EXPDEV_H
/* $Id: expdev.h,v 1.2 2011/05/03 18:50:21 sfs Exp $ */

#ifdef __cplusplus
extern "C"
    {
#endif

/* calculates exponential deviates. used for interarrival
 * times for poisson processes.
 */

double expdev (void);

#ifdef __cplusplus
    }
#endif
		
#endif  /* #ifndef COSI_RAND_EXPDEV_H */
