/* $Id: rng.c,v 1.3 2011/05/06 23:13:09 sfs Exp $ */

#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <cosi_rand/mtwist.h>
#include <cosi_rand/random.h>

unsigned long 
seed_rng (void)
{
  return mt_seed();
}

void 
set_rng_seed(unsigned long rseed) {
  mt_seed32new(rseed);
}

double 
random_double (void) 
{
  return mt_drand();
}

