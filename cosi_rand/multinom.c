/* $Id: multinom.c,v 1.2 2011/05/03 18:50:21 sfs Exp $ */
#include <cosi_rand/random.h>

void multinom(int nclass, int nitem, double prob[], int nbybin[]) {
  double x, probsum;
  int i, which;
  
  for (i = 0; i < nclass; i++) {
    nbybin[i]=0;
  }
  for(i = 0; i < nitem; i++) {
    probsum = prob[0];
    x = random_double();
    which = 0;
    /*    for (which = 0;  which <  */
    while( (x > probsum) && ( which<(nclass-1) ) )  {
      probsum += prob[++which];
    }
    nbybin[which]++;
  }
  return;
}
