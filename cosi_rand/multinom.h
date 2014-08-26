/* $Id: multinom.h,v 1.2 2011/05/03 18:50:21 sfs Exp $ */

#ifndef COSI_RAND_MULTINOM_H
#define COSI_RAND_MULTINOM_H

#ifdef __cplusplus
extern "C"
    {
#endif

void multinom(int nclass, int nitem, double prob[], int nbybin[]);

#ifdef __cplusplus
    }
#endif

			
#endif  /* #ifndef COSI_RAND_MULTINOM_H */

